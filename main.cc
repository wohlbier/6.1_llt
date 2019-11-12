#include <assert.h>
#include <iostream>

extern "C" {
#include <emu_c_utils/hooks.h>
#include <io.h>
}

#ifndef THREAD_OVERSUBSCRIBE
#define THREAD_OVERSUBSCRIBE 1
#endif
#ifndef THREADS_PER_GC
#define THREADS_PER_GC 64
#endif
#ifndef GC_PER_NODELET
#define GC_PER_NODELET 3
#endif
#define THREADS_PER_NODELET \
    THREAD_OVERSUBSCRIBE * THREADS_PER_GC * GC_PER_NODELET

#include "local_arena_allocator.h"

#include <tuple>
#include <vector>

#include <cilk.h>
#include <memoryweb.h>

typedef long Index_t;
typedef long Scalar_t;
typedef std::vector<Index_t> IndexArray_t;
typedef std::tuple<Index_t, Scalar_t> Pair_t;
#if 1
typedef std::vector<Pair_t> Row_t;
#else
typedef std::vector<Pair_t, emu::local_arena_allocator<Pair_t>> Row_t;
#endif

typedef Row_t * pRow_t;
typedef pRow_t * ppRow_t;

static inline Index_t n_map(Index_t i) { return i % NODELETS(); }
static inline Index_t r_map(Index_t i) { return i / NODELETS(); }
// inverse mapping
static inline Index_t nr_inv(Index_t n, Index_t r)
{ return r * NODELETS() + n; }

/*
 * Overrides default new to always allocate replicated storage for instances
 * of this class. repl_new is intended to be used as a parent class for
 * distributed data structure types.
 */
class repl_new
{
public:
    // Overrides default new to always allocate replicated storage for
    // instances of this class
    static void *
    operator new(std::size_t sz)
    {
        return mw_mallocrepl(sz);
    }

    // Overrides default delete to safely free replicated storage
    static void
    operator delete(void * ptr)
    {
        mw_free(ptr);
    }
};

class rMatrix_t : public repl_new
{
public:
    static rMatrix_t * create(Index_t nrows)
    {
        return new rMatrix_t(nrows);
    }

    rMatrix_t() = delete;
    rMatrix_t(const rMatrix_t &) = delete;
    rMatrix_t & operator=(const rMatrix_t &) = delete;
    rMatrix_t(rMatrix_t &&) = delete;
    rMatrix_t & operator=(rMatrix_t &&) = delete;

    void build(Index_t nl_id,
               IndexArray_t::iterator i_it,
               IndexArray_t::iterator j_it,
               IndexArray_t::iterator v_it,
               Index_t nedges)
    {
        for (Index_t ix = 0; ix < nedges; ++ix)
        {
            setElement(*i_it, *j_it, *v_it);
            ++i_it; ++j_it; ++v_it; // increment iterators
        }

        // count max degree
        max_degree_ = 0;
        for (Index_t row_idx = 0; row_idx < nrows_per_nodelet_; ++row_idx)
        {
            Index_t irow = nr_inv(nl_id, row_idx); // absolute row index
            pRow_t r = getrow(irow);
            Index_t deg = r->size();
            max_degree_ = (deg > max_degree_) ? deg : max_degree_;
        }
    }

    void set_max_degree()
    {
        Index_t max = 0;
        for (Index_t i = 0; i < NODELETS(); ++i)
        {
            Index_t lmax = *(Index_t*)mw_get_nth(&this->max_degree_, i);
            max = (lmax > max) ? lmax : max;
        }
        mw_replicated_init(&this->max_degree_, max);
    }

    Index_t nrows() { return nrows_; }
    Index_t nrows() const { return nrows_; }

    Index_t nrows_nl() { return nrows_per_nodelet_; }
    Index_t nrows_nl() const { return nrows_per_nodelet_; }

    pRow_t getrow(Index_t i) { return rows_[n_map(i)] + r_map(i); }
    pRow_t getrow(Index_t i) const { return rows_[n_map(i)] + r_map(i); }

    Index_t * row_addr(Index_t i)
    {
        return (Index_t *)(rows_ + n_map(i));
    }

    Index_t max_degree() { return max_degree_; }

private:
    rMatrix_t(Index_t nrows) : nrows_(nrows)
    {
        nrows_per_nodelet_ = r_map(nrows_);
        if (n_map(nrows_) != 0) nrows_per_nodelet_++; // add empty rows

        rows_ = (ppRow_t)mw_malloc2d(NODELETS(),
                                     nrows_per_nodelet_ * sizeof(Row_t));

        // replicate the class across nodelets
        for (Index_t i = 1; i < NODELETS(); ++i)
        {
            memcpy(mw_get_nth(this, i), mw_get_nth(this, 0), sizeof(*this));
        }

        for (Index_t i = 0; i < NODELETS(); ++i)
        {
            cilk_migrate_hint(rows_ + i);
            cilk_spawn allocateRows(i);
        }
        cilk_sync;
    }

    void allocateRows(Index_t i)
    {
        for (Index_t row_idx = 0; row_idx < nrows_per_nodelet_; ++row_idx)
        {
            new(rows_[i] + row_idx) Row_t();
        }
    }

    void setElement(Index_t irow, Index_t icol, Index_t const &val)
    {
        pRow_t r = rows_[n_map(irow)] + r_map(irow);

        if (r->empty()) // empty row
        {
            r->push_back(std::make_tuple(icol, val));
        }
        else // insert into row
        {
            Row_t::iterator it = r->begin();
            while (it != r->end() and std::get<0>(*it) < icol)
            {
                ++it;
            }
            if (it == r->end())
            {
                r->push_back(std::make_tuple(icol, val));
            }
            else
            {
                it = r->insert(it, std::make_tuple(icol, val));
            }
        }
    }

    Index_t nrows_;
    Index_t nrows_per_nodelet_;
    ppRow_t rows_;
    Index_t max_degree_;
};
typedef rMatrix_t * prMatrix_t;

static inline
bool index_exists(pRow_t r, Index_t icol)
{
    bool result = false;
    Row_t::iterator rit = r->begin();
    while (rit != r->end())
    {
        if (icol == std::get<0>(*rit))
        {
            result = true;
            break;
        }
        ++rit;
    }
    return result;
}

static inline
bool dot(Scalar_t & ans, pRow_t a, pRow_t b) // no semiring
{
    bool result = false;
    Row_t::iterator ait = a->begin();
    Row_t::iterator bit = b->begin();

    ans = 0;
    while (ait != a->end() && bit != b->end())
    {
        Index_t a_idx = std::get<0>(*ait);
        Index_t b_idx = std::get<0>(*bit);

        if (a_idx == b_idx)
        {
            ans += std::get<1>(*ait) * std::get<1>(*bit);
            result = true;
            ++ait;
            ++bit;
        }
        else if (a_idx < b_idx)
        {
            ++ait;
        }
        else
        {
            ++bit;
        }
    }
    return result;
}

static inline
void row_kernel(Index_t irow,
                prMatrix_t C,
                prMatrix_t const M,
                prMatrix_t const A,
                prMatrix_t const B)
{
    // return for empty row of A
    if (!A->getrow(irow)) return;

    std::tuple<Index_t, Scalar_t> tmp;
    // loop over columns
    for (Index_t icol = 0; icol < A->nrows(); ++icol)
    {
        // continue for empty column of B
        if (!B->getrow(icol)) continue;
        // apply mask
        if (!index_exists(M->getrow(irow), icol)) continue;

        // compute the dot
        Scalar_t ans;
        if (dot(ans, A->getrow(irow), B->getrow(icol)))
        {
            std::get<0>(tmp) = icol;
            std::get<1>(tmp) = ans;
            C->getrow(irow)->push_back(tmp);
        }
    }
}

static inline
void multi_row_kernel(Index_t nl_id,
                      Index_t t,
                      Index_t nrow,
                      prMatrix_t C,
                      prMatrix_t const M,
                      prMatrix_t const A,
                      prMatrix_t const B)
{
    for (Index_t j = t*nrow; j < (t+1)*nrow; ++j)
    {
        // absolute row index
        Index_t irow = nr_inv(nl_id, j);
        row_kernel(irow, C, M, A, B);
    }

}

static inline
void ABT_Mask_NoAccum_kernel(
    Index_t nl_id,              // nodelet_id
    prMatrix_t C,               // output matrix
    prMatrix_t const M,         // mask matrix
    // SemiringT,               // semiring
    prMatrix_t const A,         // Input matrix 1
    prMatrix_t const B,         // Input matrix 2
    bool replace_flag = false)  // put the answer in place?
{
    // making use of the fact we know that B equals L^T

    // compute rows per thread
    Index_t threads_per_nodelet = THREADS_PER_NODELET;
    Index_t nrows_per_thread = A->nrows_nl() / threads_per_nodelet;
    // if nrows_nl < threads_per_nodelet, all rows are remainder rows
    Index_t nremainder_rows = A->nrows_nl() % threads_per_nodelet;

    // spawn threads_per_nodelet threads
    if (nrows_per_thread)
    {
        for (Index_t t = 0; t < threads_per_nodelet; ++t)
        {
            cilk_spawn multi_row_kernel(nl_id, t, nrows_per_thread,
                                        C, M, A, B);
        }
        cilk_sync;
    }

    // spawn nremainder_rows threads
    if (nremainder_rows)
    {
        Index_t offset = nrows_per_thread * threads_per_nodelet;

        for (Index_t t = 0; t < nremainder_rows; ++t)
        {
            // absolute row index
            Index_t irow = nr_inv(nl_id, t + offset);
            cilk_spawn row_kernel(irow, C, M, A, B);
        }
        cilk_sync;
    }
}

Scalar_t reduce(prMatrix_t A)
{
    Scalar_t sum = 0;
    for (Index_t irow = 0; irow < A->nrows(); ++irow)
    {
        pRow_t pArow = A->getrow(irow);
        Row_t::iterator ait = pArow->begin();
        while (ait != pArow->end())
        {
            sum += std::get<1>(*ait);
            ++ait;
        }
    }
    return sum;
}

void initialize(Index_t nl_id, std::string const & filename, prMatrix_t M,
                Index_t const nnodes, Index_t const nedges)
{
    Index_t tmp;
    FILE *infile = mw_fopen(filename.c_str(), "r", &tmp);
    mw_fread(&tmp, sizeof(Index_t), 1, infile);
    //assert(tmp == nnodes); // needed with 19.09
    mw_fread(&tmp, sizeof(Index_t), 1, infile);
    //assert(tmp == nedges); // needed with 19.09

    // thread local storage to read into
    IndexArray_t iL(nedges);
    IndexArray_t jL(nedges);
    mw_fread(reinterpret_cast<void *>(iL.data()),
             //sizeof(Index_t), iL.size(), infile);
             1, sizeof(Index_t)*iL.size(), infile); // bug work around
    mw_fread(reinterpret_cast<void *>(jL.data()),
             //sizeof(Index_t), jL.size(), infile);
             1, sizeof(Index_t) * jL.size(), infile); // bug work around
    mw_fclose(infile);

    // remove edges where i is a row not owned by this nodelet.
    IndexArray_t iL_nl;
    IndexArray_t jL_nl;
    Index_t nedges_nl = 0;

    for (Index_t e = 0; e < iL.size(); ++e)
    {
        Index_t i = iL[e];
        Index_t j = jL[e];
        if (n_map(i) == nl_id)
        {
            iL_nl.push_back(i);
            jL_nl.push_back(j);
            ++nedges_nl;
        }
    }

    // build matrix
    IndexArray_t v_nl(iL_nl.size(), 1);
    M->build(nl_id, iL_nl.begin(), jL_nl.begin(), v_nl.begin(), nedges_nl);
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Requires binary edge list." << std::endl;
        std::cerr << "Usage: ./llt input.bin" << std::endl;
        exit(1);
    }

    Index_t nnodes, nedges;
    std::string filename = std::string(argv[1]);

    // open file to get number of nodes and edges, then close
    FILE *infile = mw_fopen(filename.c_str(), "r", &nnodes);
    if (!infile)
    {
        fprintf(stderr, "Unable to open file: %s\n", filename.c_str());
        exit(1);
    }
    mw_fread(&nnodes, sizeof(Index_t), 1, infile);
    mw_fread(&nedges, sizeof(Index_t), 1, infile);
    mw_fclose(infile);

    std::cerr << "nnodes: " << nnodes << std::endl;
    std::cerr << "nedges: " << nedges << std::endl;

    prMatrix_t L = rMatrix_t::create(nnodes);

    // spawn threads on each nodelet to read and build
    for (Index_t i = 0; i < NODELETS(); ++i)
    {
        cilk_migrate_hint(L->row_addr(i));
        cilk_spawn initialize(i, filename, L, nnodes, nedges);
    }
    cilk_sync;
    L->set_max_degree();

    std::cerr << "Initialization complete." << std::endl;
    std::cerr << "Max degree: " << L->max_degree() << std::endl;

    // answer matrix
    prMatrix_t C = rMatrix_t::create(nnodes);

#ifdef __PROFILE__
    hooks_region_begin("6.1_llt");
#endif
    // solve L * L^T using ABT kernel
    for (Index_t i = 0; i < NODELETS(); ++i)
    {
        cilk_migrate_hint(L->row_addr(i));
        cilk_spawn ABT_Mask_NoAccum_kernel(i, C, L, L, L);
    }
    cilk_sync;

    // reduce
    std::cerr << "Start reduction." << std::endl;
    Scalar_t nTri = reduce(C);
    std::cerr << "nTri: " << nTri << std::endl;

    // clean up matrices
    delete L;
    delete C;

#ifdef __PROFILE__
    hooks_region_end();
#endif

    return 0;
}
