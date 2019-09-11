#include <assert.h>
#include <iostream>

extern "C" {
#include <emu_c_utils/hooks.h>
#include <io.h>
}

#include "algebra.hh"
#include "types.hh"

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
             sizeof(Index_t)*iL.size(), 1, infile); // bug work around
    mw_fread(reinterpret_cast<void *>(jL.data()),
             //sizeof(Index_t), jL.size(), infile);
             sizeof(Index_t) * jL.size(), 1, infile); // bug work around
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
