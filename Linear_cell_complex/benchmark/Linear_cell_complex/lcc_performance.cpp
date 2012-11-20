//== INCLUDES =================================================================


#include "lcc_performance.h"


//=============================================================================


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage:\n" << argv[0] << " <input-mesh>\n";
        exit(1);
    }

    LCC_performance().run(argv[1], "output_cgal.off");

    return 0;
}


//=============================================================================
