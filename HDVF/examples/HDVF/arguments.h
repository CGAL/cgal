#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <CGAL/HDVF/Hdvf.h>

using namespace std;
using namespace CGAL;
using namespace HDVF;

enum class Algorithm { HDVF, DualHDVF, PerHDVF };
enum class InputFormat { SIMP, OFF, CUB, PGM };
enum class StarFiltrStd { FiltrX, FiltrY, FiltrZ };

struct Options
{
    std::string outfile_root = "res", in_file ;
    Algorithm algorithm = Algorithm::HDVF ;
    InputFormat in_format = InputFormat::OFF ;
    int HDVF_opt = CGAL::HDVF::OPT_FULL ;
    int scalar = 0 ; // 0 : Z, n : Z/nZ
    bool with_export = true ; // Export reduction information to a file
    bool with_vtk_export = false ; // Export generators/cogenerators and PSC labels as vtk
    bool with_output = true ; // Consol output
    bool with_frame = true ; // Only for persistent homology (cub and pgm in_files)
    bool loop = false ; // Only for HDVF (enter an interaction loop for M, W, MW)
    double BB_ratio = 1.5 ; // Only for persistent homolgy (off in_files: size of the bounding sphere)
    bool primal = true ; // Only for cub_complexes (for dual construction, set primal to false)
    StarFiltrStd star_filtr = StarFiltrStd::FiltrZ ; // Only for per_hdvf (when filtrataion is a standard lower star filtration)
    bool random = false ; // Activate the computation of random perfect HDVF
    bool verbose = false ; // Activate verbose mode (matrix output after each A operation for HDVF computation)
    bool co_faces = false ; // Output co-faces of cohomology generators (or not)

    std::string comment() const
    {
        std::string str;
        if (algorithm == Algorithm::HDVF)
            str = "HDVF computation";
        else if (algorithm == Algorithm::DualHDVF)
            str = "Alexander duality HDVF computation";
        else if (algorithm == Algorithm::PerHDVF)
            str = "Persistent HDVF computation";
        else
            str = "undefined";
        return str;
    }
};

inline void usage() {
    cout << "usage:" << endl ;
    cout << "\talgorithm [options] in_file" << endl ;
    cout << "with:" << endl ;
    cout << "algorithm = hdvf (for hdvf (co)homology compuation)" << endl ;
    cout << "\t dual_hdvf (for hdvf Alexander duality)" << endl ;
    cout << "\t per_hdvf (for hdvf persistent (co)homology)" << endl ;
    cout << "and options:" << endl ;
    cout << "\t-outfile_root = xxx : sets the outfile root for outfiles names to xxx (default: res)" << endl ;
    cout << "\t-opt (bnd|f|g|full) : HDVF option (default: full)" << endl ;
    cout << "\t-s n : ring used to compute homology, 0: Z, n: Z/nZ (default: 0)" << endl ;
    cout << "\t-rand : compute a random perfect HDVF (default: false, computing a random HDVF is slower)" << endl ;
    cout << "\t-verbose : verbose HDVF computation, output reduction matrices after each A (default: false)" << endl ;
    cout << "\t-with_output / -no_output : output / do not output the HDVF and reduction to the console (default: with_output)" << endl ;
    cout << "\t-with_export / -no_export : export / do not export the HDVF and reduction to a file (default: with_export)" << endl ;
    cout << "\t-with_vtk_export / -no_vtk_export : export / do not export the HDVF and reduction to vtk (default: no_export)" << endl ;
    cout << "\t-cofaces_cohomology / -no_cofaces_cohomology : do not export co-faces of cohomology generators (default: co_faces_cohomology)" << endl ;
    cout << "For CubComplex only:" << endl ;
    cout << "\t-primal / -dual: build the primal / dual complex associated to data (default: primal)" << endl ;
    cout << "For Alexander duality only:" << endl ;
    cout << "\tFor SimpComplex:" << endl ;
    cout << "\t\t-BB_ratio x : ratio for the closing bounding sphere (default: 1.5)" << endl ;
    cout << "\tFor CubComplex:" << endl ;
    cout << "\t\t-with_frame / -no_frame : enclose / do not inclose the complex in a larger bounding box (default: with_frame, adds 1 pixel around)" << endl ;
    cout << "For persistent HDVF only:" << endl ;
    cout << "\t-lower (x|y|z) : compute lower star filtration with degree x, y or z" << endl ;
    cout << "\tFor other filtrations, define the filtration into the main and set #define OwnFiltration" << endl;
}

inline Options read_arguments_hdvf(int argc, char** argv) {
    Options opt_res ;

    // Minumum number of arguments: 2
    // HDVF_algorithm file

    if (argc <= 2)
    {
        cout << "HDVF wrong number of arguments. HDVF -h for help" << endl ;
        throw "HDVF wrong number of arguments" ;
    }
    // Read algorithm from argv[0] (last word of the path - separated by /)
    const string algo(argv[0]) ;
    string word1 ;
    {
        std::stringstream string_parse(algo);
        std::string tmp_word;

        std::getline(string_parse, tmp_word, '/') ;
        word1 = tmp_word ;
        while(std::getline(string_parse, tmp_word, '/'))
        {
            word1 = tmp_word ;
        }
    }

    if (word1.compare("hdvf") == 0)
        opt_res.algorithm = Algorithm::HDVF ;
    else if (word1.compare("dual_hdvf") == 0)
        opt_res.algorithm = Algorithm::DualHDVF ;
    else if (word1.compare("per_hdvf") == 0)
        opt_res.algorithm = Algorithm::PerHDVF ;
    else
    {
        cerr << "HDVF algorithm unknown. HDVF -h for help" << endl ;
        throw "HDVF algorithm unknown" ;
    }

    // Read in_file from argv[argc-1] ;
    opt_res.in_file = argv[argc-1] ;
    // Deduce the file type from the file name
    string word ;
    {
        // Get the file name (last word of the path - separated by /)
        std::stringstream string_parse(opt_res.in_file);
        std::string tmp_word;

        std::getline(string_parse, tmp_word, '/') ;
        word = tmp_word ;
        while(std::getline(string_parse, tmp_word, '/'))
        {
            word = tmp_word ;
        }
    }
    {
        // Get the extension (last word of the word - separated by .)
        string filename = word ;
        std::stringstream string_parse(filename);
        std::string tmp_word;

        std::getline(string_parse, tmp_word, '.') ;
        word = tmp_word ;
        while(std::getline(string_parse, tmp_word, '.'))
        {
            word = tmp_word ;
        }
    }
    // Assign in_format accordingly
    if (word.compare("simp")==0)
        opt_res.in_format = InputFormat::SIMP ;
    else if (word.compare("off")==0)
        opt_res.in_format = InputFormat::OFF ;
    else if (word.compare("cub")==0)
        opt_res.in_format = InputFormat::CUB ;
    else if (word.compare("pgm")==0)
        opt_res.in_format = InputFormat::PGM ;
    else
    {
        cerr << "HDVF input format unkown. HDVF -h for help" << endl ;
        throw "HDVF input format unknown" ;
    }


    // Read next options
    for(int i = 1; i<argc-1; ++i)
    {
        string tmp(argv[i]) ;
        if (tmp.compare("-h")==0)
        {
            usage() ;
        }
        else if (tmp.compare("-opt")==0) // HDVF option
        {
            if (i == argc-2)
            {
                cerr << "HDVF: -opt option requires an option. HDVF -h for help" << endl ;
                throw "HDVF: -opt option requires an option" ;
            }
            string tmp2(argv[++i]) ;
            if (tmp2.compare("bnd")==0)
                opt_res.HDVF_opt = OPT_BND ;
            else if (tmp2.compare("g")==0)
                opt_res.HDVF_opt = OPT_G ;
            else if (tmp2.compare("f")==0)
                opt_res.HDVF_opt = OPT_F ;
            else if (tmp2.compare("full")==0)
                opt_res.HDVF_opt = OPT_FULL ;
            else
            {
                cerr << "HDVF option unkown. HDVF -h for help" << endl ;
                throw "HDVF option unknown" ;
            }
        }
        else if (tmp.compare("-s")==0) // scalar
        {
            if (i == argc-2)
            {
                cerr << "HDVF: -s option requires an integer. HDVF -h for help" << endl ;
                throw "HDVF: -s option requires an integer" ;
            }
            string tmp2(argv[++i]) ;
            opt_res.scalar = stoi(tmp2) ;
        }
        else if (tmp.compare("-rand")==0) // random HDVF
        {
            opt_res.random = true ;
        }
        else if (tmp.compare("-verbose")==0) // verbose HDVF computation
        {
            opt_res.verbose = true ;
        }
        else if (tmp.compare("-BB_ratio")==0) // scalar
        {
            if (i == argc-2)
            {
                cerr << "HDVF: -BB_ratio option requires a float. HDVF -h for help" << endl ;
                throw "HDVF: -BB_ratio option requires a float" ;
            }
            else if ((opt_res.algorithm != Algorithm::PerHDVF) || (opt_res.in_format != InputFormat::OFF))
            {
                cerr << "HDVF: -BB_ratio only for per_HDVF with off input. HDVF -h for help" << endl ;
                throw "HDVF: -BB_ratio only for per_HDVF with off input" ;
            }
            string tmp2(argv[++i]) ;
            opt_res.BB_ratio = stod(tmp2) ;
        }
        else if (tmp.compare("-with_frame")==0) // scalar
        {
            if (! ((opt_res.in_format == InputFormat::CUB) || (opt_res.in_format == InputFormat::PGM)))
            {
                cerr << "HDVF: -with_frame only for cub or pgm input. HDVF -h for help" << endl ;
                throw "HDVF: -with_frame only for cub or pgm input" ;
            }
            opt_res.with_frame = true ;
        }
        else if (tmp.compare("-no_frame")==0) // scalar
        {
            if (! ((opt_res.in_format == InputFormat::CUB) || (opt_res.in_format == InputFormat::PGM)))
            {
                cerr << "HDVF: -no_frame only for cub or pgm input. HDVF -h for help" << endl ;
                throw "HDVF: -no_frame only for cub or pgm input" ;
            }
            opt_res.with_frame = false ;
        }
        else if (tmp.compare("-with_output")==0) // scalar
        {
            opt_res.with_output = true ;
        }
        else if (tmp.compare("-no_output")==0) // scalar
        {
            opt_res.with_output = false ;
        }
        else if (tmp.compare("-with_export")==0) // scalar
        {
            opt_res.with_export = true ;
        }
        else if (tmp.compare("-no_export")==0) // scalar
        {
            opt_res.with_export = false ;
        }
        else if (tmp.compare("-with_vtk_export")==0) // scalar
        {
            if (opt_res.in_format == InputFormat::SIMP)
            {
                cerr << "HDVF: -with_vtk_export requires a geometry. HDVF -h for help" << endl ;
                throw "HDVF: -with_vtk_export requires a geometry" ;
            }
            opt_res.with_vtk_export = true ;
        }
        else if (tmp.compare("-no_vtk_export")==0) // scalar
        {
            opt_res.with_vtk_export = false ;
        }
        else if (tmp.compare("-cofaces_cohomology")==0) // scalar
        {
            opt_res.co_faces = true ;
        }
        else if (tmp.compare("-no_cofaces_cohomology")==0) // scalar
        {
            opt_res.co_faces = false ;
        }
        else if (tmp.compare("-loop")==0) // scalar
        {
            if (!((opt_res.algorithm == Algorithm::HDVF) && (opt_res.HDVF_opt & OPT_FULL)))
            {
                cerr << "HDVF: -loop available only for HDVF computation and OPT_FULL required. HDVF -h for help" << endl ;
                throw "HDVF: -loop available only for HDVF computation and OPT_FULL required" ;
            }
            opt_res.loop = true ;
        }
        else if (tmp.compare("-primal")==0) // scalar
        {
            if ((opt_res.in_format != InputFormat::CUB) && (opt_res.in_format != InputFormat::PGM))
            {
                cerr << "HDVF: -primal available only for cubical data. HDVF -h for help" << endl ;
                throw "HDVF: -primal available only for cubical data" ;
            }
            opt_res.primal = true ;
        }
        else if (tmp.compare("-dual")==0) // scalar
        {
            if ((opt_res.in_format != InputFormat::CUB) && (opt_res.in_format != InputFormat::PGM))
            {
                cerr << "HDVF: -dual available only for cubical data. HDVF -h for help" << endl ;
                throw "HDVF: -dual available only for cubical data" ;
            }
            opt_res.primal = false ;
        }
        else if (tmp.compare("-lower")==0) // scalar
        {
            if (opt_res.algorithm != Algorithm::PerHDVF)
            {
                cerr << "HDVF: -lower available only for perHDVF for lower star filtration. HDVF -h for help" << endl ;
                throw "HDVF: -lower available only for perHDVF for lower star filtration" ;
            }
            if (i == argc-2)
            {
                cerr << "HDVF: -lower needs an argument x, y or z. HDVF -h for help" << endl ;
                throw "HDVF: -lower needs an argument x, y or z" ;
            }
            else
            {
                string tmp(argv[++i]) ;

                if (tmp.compare("x")==0)
                    opt_res.star_filtr = StarFiltrStd::FiltrX ;
                else if (tmp.compare("y")==0)
                    opt_res.star_filtr = StarFiltrStd::FiltrY ;
                else if (tmp.compare("z")==0)
                    opt_res.star_filtr = StarFiltrStd::FiltrZ ;
                else
                {
                    cerr << "HDVF: -lower needs an argument x, y or z. HDVF -h for help" << endl ;
                    throw "HDVF: -lower needs an argument x, y or z" ;
                }
            }
        }
        else // unknown option
        {
            string tmp("HDVF: ") ;
            tmp += argv[i] ;
            tmp += " unknown option. HDVF -f for help" ;
            cerr << tmp << ". HDVF -h for help" << endl ;
            throw tmp ;
        }
    }
    return opt_res ;
}
inline std::ostream& operator<< (std::ostream& out, const Options& opt)
{
    out << "outfile_root: " << opt.outfile_root << endl ;

    out << "in_file: " << opt.in_file << endl ;

    out << "algorithm: " << opt.comment() << endl ;

    out << "in_format: " ;
    if (opt.in_format == InputFormat::OFF)
        out << "OFF" ;
    else if (opt.in_format == InputFormat::SIMP)
        out << "SIMP" ;
    else if (opt.in_format == InputFormat::CUB)
        out << "CUB" ;
    else if (opt.in_format == InputFormat::PGM)
        out << "PGM" ;
    out << endl ;

    out << "HDVF_opt: " ;
    if (opt.HDVF_opt == OPT_FULL)
        out << "OPT_FULL" ;
    else if (opt.HDVF_opt == OPT_G)
        out << "OPT_G" ;
    else if (opt.HDVF_opt == OPT_F)
        out << "OPT_F" ;
    else if (opt.HDVF_opt == OPT_BND)
        out << "OPT_BND" ;
    out << endl ;

    out << "scalar: " << opt.scalar << endl ;
    out << "with_export: " << opt.with_export << endl ; // Export reduction information to a file
    out << "with_vtk_export: " << opt.with_vtk_export << endl ; // Export generators/cogenerators and PSC labels as vtk
    out << "co_faces: " << opt.co_faces << endl ; // Export cofaces of cogenerators
    out << "with_output: " << opt.with_output << endl ; // Consol output
    out << "with_frame: " << opt.with_frame << endl ; // Only for persistent homology (cub and pgm in_files)
    out << "BB_ratio: " << opt.BB_ratio << endl ; // Only for persistent homolgy (off in_files: size of the bounding sphere)
    out << "primal: " << opt.primal << endl ; // Only for cubical data (primal / dual = false)
    out << "start_filtr: " ;
    if (opt.star_filtr == StarFiltrStd::FiltrX)
        out << "x" << endl ;
    else if (opt.star_filtr == StarFiltrStd::FiltrY)
        out << "y" << endl ;
    else if (opt.star_filtr == StarFiltrStd::FiltrZ)
        out << "z" << endl ;
    return out ;
}

#endif // ARGUMENTS_H

