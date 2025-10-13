// Persistent HDVF computation (command line version)
// --------
// Computes a "perfect" persistent HDVF from a given filtration
//      and provides a batch mode to specify arguments
// For help: per_hdvf -h
// --------
// A. Bac
// --------

#include <iostream>
#include <chrono>
#include <type_traits>
#include <typeinfo>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Simplex.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Cubical_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Hdvf_persistence.h>
#include <CGAL/OSM/OSM.h>
#include <CGAL/HDVF/hdvf_tools.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Cub_object_io.h>
#include "arguments.h"

namespace HDVF = CGAL::Homological_discrete_vector_field;

// ------- A ring
// For Z/nZ other than Z (ie. n=0) and Z/2Z, uncomment and set the following define properly

//#define SCALAR 5

using Kernel = CGAL::Simple_cartesian<double>;
using Traits = HDVF::Hdvf_traits_3<Kernel>;

// ------- Filtration

// For a "custom" filtration, please use the following filtration function and set the
//#define CUST_FILTRATION

// Degree type definition
#ifndef CUST_FILTRATION
using Degree = double ;
#else // Custom filtration
/// Define proper Degree below
using Degree = double ;
#endif

// Custom filtration
#ifdef CUST_FILTRATION

// Example returning x+z value for each vertex
template<typename Complex>
std::function<Degree(int)>  degree_function (const Complex& complex)
{
    std::function<Degree(int)> deg_fun_f = [&complex](int i)
    {
        const std::vector<double> Xi(complex.get_vertex_coords(i)) ;
        /// Define proper filtration value below
        return Xi.at(0)+Xi.at(2) ;
    } ;
    return deg_fun_f ;
}

#endif

template <typename MeshType, typename Complex>
void mesh_complex_output(const MeshType& mesh, const Complex& complex, const Options& options)
{
    if (options.with_output)
    {
        // Mesh
        std::cout << "----> mesh informations" << std::endl ;
        mesh.print_infos() ;

        // Complex
        std::cout << "----> complex informations" << std::endl ;
        std::cout << complex;
    }
}

template<typename Complex, typename Degree, typename FiltrationType >
HDVF::Hdvf_persistence<Complex, Degree, FiltrationType>& per_HDVF_comput (const Complex& complex, const FiltrationType& f, const Options &options)
{
    typedef typename Complex::Coefficient_ring CoefficientType;
    typedef HDVF::Hdvf_persistence<Complex, Degree, FiltrationType> HDVF_type ;
    HDVF_type& hdvf(*(new HDVF_type(complex, f, options.HDVF_opt, options.with_vtk_export)));

    std::cout << "----> START computing persistent homology" << std::endl ;
    hdvf.compute_perfect_hdvf() ;
    std::cout << "------> END computing persistent homology" << std::endl ;

    if (options.with_output)
    {
        std::cout << "----> perHDVF" << std::endl ;
        hdvf.print_hdvf_persistence_info(std::cout);
        std::cout << "----> reduction" << std::endl ;
        hdvf.insert_reduction() ;
        std::cout << "----> persistent diagram" << std::endl ;
        std::cout << hdvf ;
    }
    if (options.with_export)
    {
        std::cout << "----> exporting..." << std::endl ;
        std::string file(options.outfile_root+"_reduction.txt") ;
        std::ofstream out ( file, std::ios::out | std::ios::trunc);

        if ( not out . good () ) {
            std::cerr << "hdvf: with_export. Fatal Error:\n  " << file << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        out << "----> perHDVF" << std::endl ;
        hdvf.print_hdvf_persistence_info(out);
        out << "----> reduction" << std::endl ;
        hdvf.insert_reduction(out) ;

        out.close() ;

        std::string file_per(options.outfile_root+"_per.txt") ;
        std::ofstream out_per ( file_per, std::ios::out | std::ios::trunc);

        if ( not out_per . good () ) {
            std::cerr << "hdvf: with_export. Fatal Error:\n  " << file << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        out_per << "----> persistent diagram" << std::endl ;
        out_per << hdvf ;

        out_per.close() ;
    }
    return hdvf ;
}

template<typename Complex, typename Degree, typename FiltrationType>
FiltrationType& build_filtration(const Complex& complex, const Options& options)
{
    typedef FiltrationType Filtration;
    std::function<Degree(size_t)>  deg_function ;
#ifndef CUST_FILTRATION
    // Standard lower star filtration along x,y or z
    if (options.star_filtr == StarFiltrStd::FiltrX)
        deg_function = (degree_function<Complex,typename Traits::Point>(complex, Kernel::Compute_x_3())) ;
    else if (options.star_filtr == StarFiltrStd::FiltrY)
        deg_function = (degree_function<Complex,typename Traits::Point>(complex, Kernel::Compute_y_3())) ;
    else if (options.star_filtr == StarFiltrStd::FiltrZ)
        deg_function = (degree_function<Complex,typename Traits::Point>(complex, Kernel::Compute_z_3())) ;
    else
    {
        std::cout << "Unknown lower start filtration" << std::endl ;
        throw("Unknown lower start filtration") ;
    }
#else
    // Use degree_function defined above
    deg_function = (degree_function(complex)) ;
#endif

    std::cout << "----> START building filtration" << std::endl ;
    Filtration& f = *(new Filtration(complex, deg_function)) ;
    std::cout << "------> END building filtration" << std::endl ;
    return f ;
}


template <typename CoefficientType>
void main_code (const Options &options)
{
    /// SIMP format
    if (options.in_format == InputFormat::SIMP)
    {
//        using Complex = AbstractSimpComplex<CoefficientType>  ;
//        using HDVF_type = HDVF<CoefficientType, Complex> ;
//
//        // MeshObject
//        MeshObject mesh ;
//        mesh.read_simp(options.in_file) ;
//
//        // Complex
//        Complex complex(mesh);
//
//        mesh_complex_output<MeshObject, Complex>(mesh, complex, options) ;
//
//        // HDVF computation, export, output
//        HDVF_type hdvf(HDVF_comput<CoefficientType,Complex>(complex, options)) ;
//
//        // Export to vtk
//        // None for SIMP format
        std::cout << "not yet..." << std::endl ;
        throw("not yet") ;
    }
    /// OFF format
    else if (options.in_format == InputFormat::OFF)
    {
        using Complex = HDVF::Simplicial_chain_complex<CoefficientType, Traits> ;
        using FiltrationType = HDVF::Filtration_lower_star<Complex, Degree> ;
        using HDVF_type = HDVF::Hdvf_persistence<Complex, Degree, FiltrationType> ;

        // MeshObject
        HDVF::Mesh_object_io<Traits> mesh ;
        mesh.read_off(options.in_file) ;

        // Complex
        Complex complex(mesh);

        // Build filtration
        FiltrationType& f(build_filtration<Complex, Degree, FiltrationType>(complex, options)) ;

        mesh_complex_output<HDVF::Mesh_object_io<Traits>, Complex>(mesh, complex, options) ;

        // HDVF computation, export, output
        HDVF_type& hdvf(per_HDVF_comput<Complex,Degree, FiltrationType>(complex,f, options)) ;

        // Export to vtk
        if (options.with_vtk_export)
        {
            std::cout << "----> exporting to vtk" << std::endl ;
            CGAL::IO::write_VTK(hdvf, complex, options.outfile_root, options.co_faces) ;
        }
    }
    // CubComplex
    else if ((options.in_format == InputFormat::PGM) || (options.in_format == InputFormat::CUB))
    {
        using Complex = HDVF::Cubical_chain_complex<CoefficientType,Traits> ;
        using FiltrationType = HDVF::Filtration_lower_star<Complex, Degree> ;
        using HDVF_type = HDVF::Hdvf_persistence<Complex, Degree, FiltrationType> ;

        HDVF::Cub_object_io mesh ;
        typename Complex::Cubical_complex_primal_dual primal_dual(Complex::PRIMAL) ;
        if (options.primal)
        {
            if (options.in_format == InputFormat::PGM)
                mesh.read_pgm(options.in_file, true) ; // Read with Khalimsky coordinates (for primal)
            else
                mesh.read_cub(options.in_file, true) ; // Read with Khalimsky coordinates (for primal)
        }
        else // dual
        {
            if (options.in_format == InputFormat::PGM)
                mesh.read_pgm(options.in_file, false) ; // Read with pixel coordinates (for dual)
            else
                mesh.read_cub(options.in_file, false) ; // Read with pixel coordinates (for dual)
            primal_dual = Complex::DUAL ;
        }

        // Complex
        Complex complex(mesh, primal_dual);

        mesh_complex_output<HDVF::Cub_object_io, Complex>(mesh, complex, options) ;

        // Build filtration
        FiltrationType& f(build_filtration<Complex, Degree, FiltrationType>(complex, options)) ;

        // HDVF computation, export, output
        HDVF_type& hdvf(per_HDVF_comput<Complex,Degree,FiltrationType>(complex,f, options)) ;

        // Export to vtk
        if (options.with_vtk_export)
        {
            std::cout << "----> exporting to vtk" << std::endl ;
            CGAL::IO::write_VTK(hdvf, complex, options.outfile_root, options.co_faces) ;
        }
    }
}

// Main

int main(int argc, char **argv)
{
    if (argc <= 2)
        usage() ;
    else
    {
        for (int i=0;i<argc; ++i)
            std::cout << "arg " << i << " : " << argv[i] << std::endl ;

        Options options(read_arguments_hdvf(argc, argv)) ;
        std::cout << "options:" << std::endl << options ;

        // ----- Definition of the CoefficientType
#ifndef SCALAR
        if (options.scalar == 0)
        {
            using CoefficientType = int ;
            main_code<CoefficientType>(options) ;
        }
        else if (options.scalar == 2)
        {
            using CoefficientType = HDVF::Z2 ;
            main_code<CoefficientType>(options) ;
        }
        else
        {
            std::cerr << "Z" << options.scalar << " not instantiated, use the #define at line 27" << std::endl ;
        }
#else
        typedef HDVF::Zp<SCALAR,int,true> CoefficientType;
#endif
    }

    return 0 ;
}















