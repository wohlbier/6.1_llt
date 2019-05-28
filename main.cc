#include <assert.h>
#include <iostream>
#include <tuple>
#include <vector>

#include <cilk.h>
#include <memoryweb.h>
#include <distributed.h>
extern "C" {
#include <emu_c_utils/layout.h>
#include <emu_c_utils/hooks.h>
}

#include "algebra.hh"
#include "types.hh"

void build(prMatrix_t L, prIndexArray_t riL, prIndexArray_t rjL)
{
    // create lists specific to this nodelet
    IndexArray_t iL, jL;
    Index_t nedgesL = 0;

    // loop over all edges and pick i nodes that will reside on this
    // nodelet.
    for (Index_t e = 0; e < riL->n(); e++)
    {
        Index_t i = riL->data(NODE_ID())[e];
        Index_t j = rjL->data(NODE_ID())[e];
        if (n_map(i) == NODE_ID())
        {
            iL.push_back(i);
            jL.push_back(j);
            ++nedgesL;
        }
    }
    IndexArray_t v(iL.size(), 1);
    L->build(iL.begin(), jL.begin(), v.begin(), nedgesL);
}

int main(int argc, char* argv[])
{
#ifdef __PROFILE__
    hooks_region_begin("6.1_llt");
#endif

    IndexArray_t iL;
    IndexArray_t jL;
    Index_t nedgesL = 0;
    Index_t max_id = 0;
    Index_t src, dst;

    FILE *infile = fopen(argv[1], "r");
    if (!infile)
    {
        fprintf(stderr, "Unable to open file: %s\n", argv[1]);
        exit(1);
    }

    // read edges in lower triangle of adjacency matrix
    while (!feof(infile))
    {
        fscanf(infile, "%ld %ld\n", &src, &dst);
        if (src > max_id) max_id = src;
        if (dst > max_id) max_id = dst;

        if (dst < src)
        {
            iL.push_back(src);
            jL.push_back(dst);
            ++nedgesL;
        }
    }
    fclose(infile);

    Index_t nnodes = max_id + 1;
    std::cerr << "num nodes: " << nnodes << std::endl;
    std::cerr << "num edges: " << nedgesL << std::endl;

    prIndexArray_t riL = rIndexArray_t::create(nedgesL);
    prIndexArray_t rjL = rIndexArray_t::create(nedgesL);

    // deep copy iL, jL into riL, rjL. each nodelet has entire edge list.
    for (Index_t i = 0; i < NODELETS(); ++i)
    {
        memcpy(riL->data(i), iL.data(), iL.size() * sizeof(Index_t));
        memcpy(rjL->data(i), jL.data(), jL.size() * sizeof(Index_t));
    }

    prMatrix_t L = rMatrix_t::create(nnodes);

    // spawn build functions on each nodelet
    for (Index_t i = 0; i < NODELETS(); ++i)
    {
        cilk_migrate_hint(riL->data(i));
        cilk_spawn build(L, riL, rjL);
    }
    cilk_sync;

    // clean up auxilliary arrays
    delete riL;
    delete rjL;

    rMatrix_t * C = rMatrix_t::create(nnodes);

    // solve L * L^T using ABT kernel
    for (Index_t i = 0; i < NODELETS(); ++i)
    {
        cilk_migrate_hint(L->row_addr(i));
        cilk_spawn ABT_Mask_NoAccum_kernel(C, L, L, L);
    }
    cilk_sync;

    // reduce
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
