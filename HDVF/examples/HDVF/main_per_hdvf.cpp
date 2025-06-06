// HDVF off files
// --------
// Computes a "perfect" persistent HDVF from a given filtration
//      and provides a batch mode to specify arguments
// For help: hdvf -h
// --------
// A. Bac
// --------

#include <iostream>
#include <chrono>
#include <type_traits>
#include <typeinfo>
#include "CGAL/HDVF/Zp.hpp"
#include "CGAL/HDVF/Simplex.hpp"
#include "CGAL/HDVF/tools_io.hpp"
#include "CGAL/HDVF/Abstract_simplicial_chain_complex.hpp"
#include "CGAL/HDVF/SimpComplexTools.hpp"
#include "CGAL/HDVF/Cubical_chain_complex.hpp"
#include "CGAL/HDVF/CubComplexTools.hpp"
#include "CGAL/HDVF/Hdvf_persistence.h"
#include "CGAL/OSM/OSM.hpp"
#include "CGAL/HDVF/hdvf_tools.hpp"
#include "arguments.h"

using namespace CGAL ;
using namespace HDVF ;

// ------- A ring
// For Z/nZ other than Z (ie. n=0) and Z/2Z, uncomment and set the following define properly

//#define SCALAR 5

// ------- Filtration

// For a "custom" filtration, please use the following filtration function and set the
//#define CUST_FILTRATION

// Degree type definition
#ifndef CUST_FILTRATION
using DegType = double ;
#else // Custom filtration
/// Define proper DegType below
using DegType = double ;
#endif

// Custom filtration
#ifdef CUST_FILTRATION

// Example returning x+z value for each vertex
template<typename ComplexType>
std::function<DegType(int)>  deg_fun (const ComplexType& complex)
{
    std::function<DegType(int)> deg_fun_f = [&complex](int i)
    {
        const std::vector<double> Xi(complex.get_vertex_coords(i)) ;
        /// Define proper filtration value below
        return Xi.at(0)+Xi.at(2) ;
    } ;
    return deg_fun_f ;
}

#endif

template <typename MeshType, typename ComplexType>
void mesh_complex_output(const MeshType& mesh, const ComplexType& complex, const Options& options)
{
    if (options.with_output)
    {
        // Mesh
        cout << "----> mesh informations" << endl ;
        mesh.print_infos() ;
        
        // Complex
        cout << "----> complex informations" << endl ;
        complex.print_complex();
    }
}

template<typename CoefficientType, typename ComplexType, typename DegType, typename FiltrationType >
Hdvf_persistence<CoefficientType, ComplexType, DegType, FiltrationType>& per_HDVF_comput (const ComplexType& complex, const FiltrationType& f, const Options &options)
{
    typedef Hdvf_persistence<CoefficientType, ComplexType, DegType, FiltrationType> HDVFType ;
    HDVFType& hdvf(*(new HDVFType(complex, f, options.HDVF_opt, options.with_vtk_export, options.co_faces)));
    
    cout << "----> START computing persistent homology" << endl ;
    hdvf.compute_perfect_hdvf() ;
    cout << "------> END computing persistent homology" << endl ;
    
    if (options.with_output)
    {
        cout << "----> perHDVF" << endl ;
        hdvf.print_hdvf_persistence_info(cout);
        cout << "----> reduction" << endl ;
        hdvf.print_reduction() ;
        cout << "----> persistent diagram" << endl ;
        cout << hdvf ;
    }
    if (options.with_export)
    {
        cout << "----> exporting..." << endl ;
        string file(options.outfile_root+"_reduction.txt") ;
        std::ofstream out ( file, std::ios::out | std::ios::trunc);

        if ( not out . good () ) {
            std::cerr << "hdvf: with_export. Fatal Error:\n  " << file << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        out << "----> perHDVF" << endl ;
        hdvf.print_hdvf_persistence_info(out);
        out << "----> reduction" << endl ;
        hdvf.print_reduction(out) ;
        
        out.close() ;
        
        string file_per(options.outfile_root+"_per.txt") ;
        std::ofstream out_per ( file_per, std::ios::out | std::ios::trunc);

        if ( not out_per . good () ) {
            std::cerr << "hdvf: with_export. Fatal Error:\n  " << file << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        out_per << "----> persistent diagram" << endl ;
        out_per << hdvf ;
        
        out_per.close() ;
    }
    return hdvf ;
}

template<typename CoefficientType, typename ComplexType, typename DegType, typename FiltrationType>
FiltrationType& build_filtration(const ComplexType& complex, const Options& options)
{
    typedef FiltrationType Filtration;
    std::function<DegType(int)>  deg_function ;
#ifndef CUST_FILTRATION
    // Standard lower star filtration along x,y or z
    if (options.star_filtr == StarFiltrStd::FiltrX)
        deg_function = (deg_fun(complex, f_x)) ;
    else if (options.star_filtr == StarFiltrStd::FiltrY)
        deg_function = (deg_fun(complex, f_y)) ;
    else if (options.star_filtr == StarFiltrStd::FiltrZ)
        deg_function = (deg_fun(complex, f_z)) ;
    else
    {
        cout << "Unknown lower start filtration" << endl ;
        throw("Unknown lower start filtration") ;
    }
#else
    // Use deg_fun defined above
    deg_function = (deg_fun(complex)) ;
#endif
    
//    Filtration& f = *(new Filtration(complex)) ;
//    cout << "----> START building filtration" << endl ;
//    f.star_filtration(deg_function) ;
//    cout << "------> END building filtration" << endl ;
//    return f ;
    cout << "----> START building filtration" << endl ;
    Filtration& f = *(new Filtration(complex, deg_function)) ;
    cout << "------> END building filtration" << endl ;
    return f ;
}


template <typename CoefficientType>
void main_code (const Options &options)
{
    /// SIMP format
    if (options.in_format == InputFormat::SIMP)
    {
//        using ComplexType = AbstractSimpComplex<CoefficientType>  ;
//        using HDVFType = HDVF<CoefficientType, ComplexType> ;
//        
//        // MeshObject
//        MeshObject mesh ;
//        mesh.read_simp(options.in_file) ;
//        
//        // Complex
//        ComplexType complex(mesh);
//        
//        mesh_complex_output<MeshObject, ComplexType>(mesh, complex, options) ;
//        
//        // HDVF computation, export, output
//        HDVFType hdvf(HDVF_comput<CoefficientType,ComplexType>(complex, options)) ;
//        
//        // Export to vtk
//        // None for SIMP format
        cout << "not yet..." << endl ;
        throw("not yet") ;
    }
    /// OFF format
    else if (options.in_format == InputFormat::OFF)
    {
        using ComplexType = Simplicial_chain_complex<CoefficientType> ;
        using FiltrationType = Filtration_lower_star<CoefficientType, ComplexType, DegType> ;
        using HDVFType = Hdvf_persistence<CoefficientType, ComplexType,DegType, FiltrationType> ;
        
        // MeshObject
        Mesh_object mesh ;
        mesh.read_off(options.in_file) ;
        
        // Complex
        ComplexType complex(mesh, mesh.nodes);
        
        // Build filtration
        FiltrationType& f(build_filtration<CoefficientType, ComplexType, DegType, FiltrationType>(complex, options)) ;
        
        mesh_complex_output<Mesh_object, ComplexType>(mesh, complex, options) ;
        
        // HDVF computation, export, output
        HDVFType& hdvf(per_HDVF_comput<CoefficientType,ComplexType,DegType, FiltrationType>(complex,f, options)) ;
        
        // Export to vtk
        if (options.with_vtk_export)
        {
            cout << "----> exporting to vtk" << endl ;
            Per_Simp_output_vtk<CoefficientType,DegType,FiltrationType>(hdvf, complex, options.outfile_root) ;
        }
    }
    // CubComplex
    else if ((options.in_format == InputFormat::PGM) || (options.in_format == InputFormat::CUB))
    {
        using ComplexType = Cubical_chain_complex<CoefficientType> ;
        using FiltrationType = Filtration_lower_star<CoefficientType, ComplexType, DegType> ;
        using HDVFType = Hdvf_persistence<CoefficientType, ComplexType, DegType, FiltrationType> ;
        
        Cub_object mesh ;
        typename ComplexType::typeComplexCube primal_dual(ComplexType::PRIMAL) ;
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
            primal_dual = ComplexType::DUAL ;
        }
        
        // Complex
        ComplexType complex(mesh, primal_dual);
        
        mesh_complex_output<Cub_object, ComplexType>(mesh, complex, options) ;
        
        // Build filtration
        FiltrationType& f(build_filtration<CoefficientType, ComplexType, DegType, FiltrationType>(complex, options)) ;
        
        // HDVF computation, export, output
        HDVFType& hdvf(per_HDVF_comput<CoefficientType,ComplexType,DegType,FiltrationType>(complex,f, options)) ;
        
        // Export to vtk
        if (options.with_vtk_export)
        {
            cout << "----> exporting to vtk" << endl ;
            Per_Cub_output_vtk<CoefficientType, DegType,FiltrationType>(hdvf, complex, options.outfile_root) ;
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
            cout << "arg " << i << " : " << argv[i] << endl ;
        
        Options options(read_arguments_hdvf(argc, argv)) ;
        cout << "options:" << endl << options ;
        
        // ----- Definition of the CoefficientType
#ifndef SCALAR
        if (options.scalar == 0)
        {
            using CoefficientType = int ;
            main_code<CoefficientType>(options) ;
        }
        else if (options.scalar == 2)
        {
            using CoefficientType = Zp<2,int8_t> ;
            main_code<CoefficientType>(options) ;
        }
        else
        {
            cout << "Z" << options.scalar << " not instantiated, use the #define at line 27" << endl ;
        }
#else
        typedef Zp<SCALAR> CoefficientType;
#endif
    }
    
    return 0 ;
}









         
    



    
