#ifndef ALGEBRA_HH
#define ALGEBRA_HH

#include "types.hh"

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
            C->getrow(irow)->push_back(std::make_tuple(icol, ans));
        }
    }
}

static inline
void ABT_Mask_NoAccum_kernel(
    prMatrix_t C,               // output matrix
    prMatrix_t const M,         // mask matrix
    // SemiringT,               // semiring
    prMatrix_t const A,         // Input matrix 1
    prMatrix_t const B,         // Input matrix 2
    bool replace_flag = false)  // put the answer in place?
{
    // making use of the fact we know that B equals L^T

    // spawn a thread for each row that the nodelet owns
    // spawn it locally, so no migrate hint
    for (Index_t i = 0; i < A->nrows_nl(); ++i)
    {
        // absolute row index
        Index_t irow = nr_inv(NODE_ID(), i);
        cilk_spawn row_kernel(irow, C, M, A, B);
    }
    cilk_sync;
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

#endif // ALGEBRA_HH
