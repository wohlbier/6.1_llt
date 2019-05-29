#ifndef TYPES_HH
#define TYPES_HH

#ifndef THREADS_PER_NODELET
#define THREADS_PER_NODELET 64
#endif

#include <tuple>
#include <vector>

#include <cilk.h>

typedef long Index_t;
typedef long Scalar_t;
typedef std::vector<Index_t> IndexArray_t;
typedef std::vector<std::tuple<Index_t, Scalar_t>> Row_t;
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

    void build(IndexArray_t::iterator i_it,
               IndexArray_t::iterator j_it,
               IndexArray_t::iterator v_it,
               Index_t nedges)
    {
        for (Index_t ix = 0; ix < nedges; ++ix)
        {
            setElement(*i_it, *j_it, *v_it);
            ++i_it; ++j_it; ++v_it; // increment iterators
        }
    }

    Index_t nrows() { return nrows_; }
    Index_t nrows() const { return nrows_; }

    Index_t nrows_nl() { return nrows_per_nodelet_; }
    Index_t nrows_nl() const { return nrows_per_nodelet_; }

    pRow_t getrow(Index_t i) { return rows_[n_map(i)] + r_map(i); }
    pRow_t getrow(Index_t i) const { return rows_[n_map(i)] + r_map(i); }

    Index_t * row_addr(Index_t i)
    {
        return (Index_t *)(rows_ + i);
    }

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
        for (Index_t row_idx= 0; row_idx < nrows_per_nodelet_; ++row_idx)
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
            while (std::get<0>(*it) < icol and it != r->end())
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
};
typedef rMatrix_t * prMatrix_t;

#endif // TYPES_HH
