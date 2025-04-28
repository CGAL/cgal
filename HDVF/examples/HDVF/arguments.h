#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include "CGAL/HDVF/hdvf.hpp"

enum class Algorithm { HDVF, DualHDVF, PerHDVF };
enum class InputFormat { SIMP, OFF, CUB, PGM };
enum class StarFiltrStd { FiltrX, FiltrY, FiltrZ };

struct Options
{
    std::string outfile_root = "res", in_file ;
    Algorithm algorithm = Algorithm::HDVF ;
    InputFormat in_format = InputFormat::OFF ;
    int HDVF_opt = OPT_FULL ;
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

void usage() ;
Options read_arguments_hdvf(int argc, char** argv) ;
std::ostream& operator<< (std::ostream& out, const Options& opt) ;

#endif // ARGUMENTS_H

