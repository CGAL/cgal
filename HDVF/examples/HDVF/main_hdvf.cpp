// HDVF computation (command line version)
// --------
// Computes a "perfect" Hdvf and provides a batch mode to specify arguments
// For help: hdvf -h
// --------
// A. Bac
// --------

#include <iostream>
#include <chrono>
#include <type_traits>
#include <typeinfo>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Simplex.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Cubical_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>
#include <CGAL/OSM/Sparse_chain.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Cub_object_io.h>
#include <CGAL/HDVF/Hdvf_tools.h>
#include "arguments.h"


// ------- A ring
// For Z/nZ other than Z (ie. n=0) and Z/2Z, uncomment and set the following define properly

//#define SCALAR 5

// Hdvf computation, output and export

template <typename MeshType, typename ComplexType>
void mesh_complex_output(const MeshType& mesh, const ComplexType& complex, const Options& options)
{
    if (options.with_output)
    {
        // Mesh
        std::cout << "----> mesh informations" << std::endl ;
        mesh.print_infos() ;
        
        // Complex
        std::cout << "----> complex informations" << std::endl ;
        complex.print_complex();
    }
}

template <typename CoefficientType, typename ComplexType>
CGAL::HDVF::Hdvf<CoefficientType, ComplexType>& HDVF_comput (const ComplexType& complex, const Options &options)
{
    CGAL::HDVF::Hdvf<CoefficientType, ComplexType>& hdvf(*(new CGAL::HDVF::Hdvf<CoefficientType, ComplexType>(complex, options.HDVF_opt)));
    std::vector<CGAL::HDVF::PairCell> pairs ;
    if (!options.random)
        pairs = hdvf.compute_perfect_hdvf(options.verbose);
    else
        pairs = hdvf.compute_rand_perfect_hdvf(options.verbose);
    
    if (options.with_output)
    {
        std::cout << "----> pairs found by computePerfectHDVF" << std::endl ;
        std::cout << pairs ;
        std::cout << "----> reduction" << std::endl ;
        hdvf.print_reduction() ;
    }
    if (options.with_export)
    {
        std::string file(options.outfile_root+"_reduction.txt") ;
        std::ofstream out ( file, std::ios::out | std::ios::trunc);

        if ( not out . good () ) {
            std::cerr << "hdvf: with_export. Fatal Error:\n  " << file << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        out << "----> pairs found by computePerfectHDVF" << std::endl ;
        out << pairs ;
        out << "----> reduction" << std::endl ;
        hdvf.print_reduction(out) ;
        
        out.close() ;
    }
    return hdvf ;
}

template <typename CoefficientType>
void main_code (const Options &options)
{
    /// SIMP format
    if (options.in_format == InputFormat::SIMP)
    {
        using ComplexType = CGAL::HDVF::Abstract_simplicial_chain_complex<CoefficientType>  ;
        using HDVFType = CGAL::HDVF::Hdvf<CoefficientType, ComplexType> ;
        
        // MeshObject
        CGAL::HDVF::Mesh_object_io mesh ;
        mesh.read_simp(options.in_file) ;
        
        // Complex
        ComplexType complex(mesh);
        
        mesh_complex_output<CGAL::HDVF::Mesh_object_io, ComplexType>(mesh, complex, options) ;
        
        // Hdvf computation, export, output
        HDVFType hdvf(HDVF_comput<CoefficientType,ComplexType>(complex, options)) ;
        
        // Export to vtk
        // None for SIMP format
    }
    /// OFF format
    else if (options.in_format == InputFormat::OFF)
    {
        using ComplexType = CGAL::HDVF::Simplicial_chain_complex<CoefficientType> ;
        using HDVFType = CGAL::HDVF::Hdvf<CoefficientType, ComplexType> ;
        
        // MeshObject
        CGAL::HDVF::Mesh_object_io mesh ;
        mesh.read_off(options.in_file) ;
        
        // Complex
        ComplexType complex(mesh, mesh.get_nodes());
        
        mesh_complex_output<CGAL::HDVF::Mesh_object_io, ComplexType>(mesh, complex, options) ;
        
        // Hdvf computation, export, output
        HDVFType hdvf(HDVF_comput<CoefficientType,ComplexType>(complex, options)) ;
        
        // Loop on operations with vtk export
        if (options.loop)
        {
            auto output_vtk_simp = [options](HDVFType &hdvf, ComplexType& complex)
            {
                CGAL::HDVF::hdvf_geometric_chain_complex_output_vtk(hdvf, complex, options.outfile_root, options.co_faces) ;
            } ;
            
            CGAL::HDVF::interaction_loop<CoefficientType, ComplexType>(hdvf, complex, output_vtk_simp) ;
        }
        // Export to vtk
        else if (options.with_vtk_export)
        {
            std::cout << "----> exporting to vtk" << std::endl ;
            CGAL::HDVF::hdvf_geometric_chain_complex_output_vtk(hdvf, complex, options.outfile_root, options.co_faces) ;
        }
    }
    // CubComplex
    else if ((options.in_format == InputFormat::PGM) || (options.in_format == InputFormat::CUB))
    {
        using ComplexType = CGAL::HDVF::Cubical_chain_complex<CoefficientType> ;
        using HDVFType = CGAL::HDVF::Hdvf<CoefficientType, ComplexType> ;
        
        CGAL::HDVF::Cub_object_io mesh ;
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
        
        mesh_complex_output<CGAL::HDVF::Cub_object_io, ComplexType>(mesh, complex, options) ;
        
        // Hdvf computation, export, output
        HDVFType hdvf(HDVF_comput<CoefficientType,ComplexType>(complex, options)) ;
        
        // Loop on operations with vtk export
        if (options.loop)
        {
            auto output_vtk_cub = [options](HDVFType &hdvf, ComplexType& complex)
            {
                CGAL::HDVF::hdvf_geometric_chain_complex_output_vtk(hdvf, complex, options.outfile_root, options.co_faces) ;
            } ;
            
            CGAL::HDVF::interaction_loop<CoefficientType, ComplexType>(hdvf, complex, output_vtk_cub) ;
        }
        // Export to vtk
        else if (options.with_vtk_export)
        {
            std::cout << "----> exporting to vtk" << std::endl ;
            CGAL::HDVF::hdvf_geometric_chain_complex_output_vtk(hdvf, complex, options.outfile_root, options.co_faces) ;
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
        
        // TEST
        std::cout << "TEST" << std::endl ;
        CGAL::OSM::Sparse_chain<int, CGAL::OSM::COLUMN> gamma(4) ;
        std::cout << "before: " << gamma.is_null() << std::endl ;
        gamma.set_coef(2, 1) ;
        std::cout << "after: " << gamma.is_null() << std::endl ;
        gamma.set_coef(2, 0) ;
        std::cout << "after2: " << gamma.is_null() << std::endl ;
        std::cout << "END TEST" << std::endl ;
        
        // ----- Definition of the CoefficientType
#ifndef SCALAR
        if (options.scalar == 0)
        {
            using CoefficientType = int ;
            main_code<CoefficientType>(options) ;
        }
        else if (options.scalar == 2)
        {
//            using CoefficientType = CGAL::HDVF::Zp<2,int8_t> ;
            using CoefficientType = CGAL::HDVF::Z2 ;
            main_code<CoefficientType>(options) ;
        }
        else
        {
            std::cerr << "Z" << options.scalar << " not instantiated, use the #define at line 27" << std::endl ;
        }
#else
        typedef Zp<SCALAR> CoefficientType;
#endif
    }
    
    return 0 ;
}









         
    



    
