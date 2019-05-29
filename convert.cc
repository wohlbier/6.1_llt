#include <fstream>
#include <iostream>
#include <vector>

typedef size_t Index_t;
typedef std::vector<Index_t> IndexArray_t;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Converts ASCII tsv edge list to binary file with only lower triangular edges." << std::endl;
        std::cerr << "Binary file has nnodes, nedges, iL(nedges), jL(nedges)." << std::endl;
        std::cerr << "Data type is size_t." << std::endl;
        std::cerr << "Usage: ./convert input.tsv" << std::endl;
        exit(1);
    }

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

    // replace .tsv with .bin for output file
    std::string ofname = std::string(argv[1]) + ".bin";
    size_t start_pos = ofname.find(".tsv");
    ofname.replace(start_pos, 4, "");
    std::cerr << "Writing file: " << ofname << std::endl;

    auto outfile = std::fstream(ofname, std::ios::out | std::ios::binary);
    outfile.write(reinterpret_cast<const char *>(&nnodes),  sizeof(Index_t));
    outfile.write(reinterpret_cast<const char *>(&nedgesL), sizeof(Index_t));
    outfile.write(reinterpret_cast<const char *>(iL.data()),
                  iL.size() * sizeof(Index_t));
    outfile.write(reinterpret_cast<const char *>(jL.data()),
                  jL.size() * sizeof(Index_t));
    outfile.close();

// Read and echo file
#if 0
    {
        Index_t nn, ne;

        auto infile = std::fstream(ofname, std::ios::in | std::ios::binary);
        infile.read(reinterpret_cast<char *>(&nn), sizeof(Index_t));
        infile.read(reinterpret_cast<char *>(&ne), sizeof(Index_t));
        IndexArray_t iL(ne), jL(ne);
        infile.read(reinterpret_cast<char *>(iL.data()),
                    iL.size() * sizeof(Index_t));
        infile.read(reinterpret_cast<char *>(jL.data()),
                    jL.size() * sizeof(Index_t));
        infile.close();

        std::cerr << "nn: " << nn << ", ne: " << ne << std::endl;

        for (Index_t i = 0; i < iL.size(); ++i)
        {
            std::cerr << "(iL, jL): (" << iL[i] << ", "
                      << jL[i] << ")" << std::endl;
        }
    }
#endif

    return 0;
}
