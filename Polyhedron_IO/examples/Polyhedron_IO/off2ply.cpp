#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

bool  verbose  = false;
bool  binary   = false;

int main(int argc, char *argv[])
{
    int n = 0; // number of filenames
    char *filename[2];
    bool help = false;
    for (int i = 1; i < argc; i++) { // check commandline options
        if ( strcmp( "-v", argv[i]) == 0)
            verbose = true;
        else if ( strcmp( "-b", argv[i]) == 0)
            binary = true;
        else if ( (strcmp( "-h", argv[i]) == 0) ||
                  (strcmp( "-help", argv[i]) == 0))
            help = true;
        else if ( n < 2 ) {
            filename[ n++] = argv[i];
        } else {
            ++n;
            break;
        }
    }

    if ((n > 2) || help) {
        if ( ! help)
            std::cerr << "Error: in parameter list" << std::endl;
        std::cerr << "Usage: " << argv[0] << " [<options>] [<infile> [<outfile>]]"
             << std::endl;
        std::cerr << "       convert a CGAL object (OFF) to (PLY) "
                "format." << std::endl;
        std::cerr << "       -v      verbose." << std::endl;
        std::cerr << "       -b      binary." << std::endl;
        exit( ! help);
    }

    CGAL::Verbose_ostream vout(verbose);
    vout << argv[0] << ": verbosity on." << std::endl;

    const char*  iname = "cin";
    std::istream*     p_in  = &(std::cin);
    std::ifstream     in;
    if ( n > 0) {
        in.open( filename[0]);
        p_in = &in;
        iname = filename[0];
    }
    if ( !*p_in) {
        std::cerr << argv[0] << ": error: cannot open file '"<< iname
             << "' for reading." <<std::endl;
        exit( 1);
    }

    const char*  oname = "cout";
    std::ostream*     p_out = &(std::cout);
    std::ofstream     out;
    if ( n > 1) {
        out.open( filename[1]);
        p_out = &out;
        oname = filename[1];
    }
    if ( !*p_out) {
        std::cerr << argv[0] << ": error: cannot open file '"<< oname
            << "' for writing." << std::endl;
        exit( 1);
    }

    Mesh sm1;
    *p_in >> sm1;

    CGAL::set_ascii_mode(*p_out);
    if(binary)
        CGAL::set_binary_mode(*p_out);

    CGAL::write_ply(*p_out, sm1);

    return 0;
}
