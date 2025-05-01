// Hdvf off files
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
#include "CGAL/Hdvf/Zp.hpp"
#include "CGAL/Hdvf/Simplex.hpp"
#include "CGAL/Hdvf/tools_io.hpp"
#include "CGAL/Hdvf/Abstract_simplicial_chain_complex.hpp"
#include "CGAL/Hdvf/SimpComplexTools.hpp"
#include "CGAL/Hdvf/Cubical_chain_complex.hpp"
#include "CGAL/Hdvf/CubComplexTools.hpp"
#include "CGAL/Hdvf/hdvf.hpp"
#include "CGAL/OSM/OSM.hpp"
#include "CGAL/Hdvf/hdvf_tools.hpp"
#include "arguments.h"

using namespace CGAL ;
using namespace HDVF ;

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
        cout << "----> mesh informations" << endl ;
        mesh.print_infos() ;
        
        // Complex
        cout << "----> complex informations" << endl ;
        complex.print_complex();
    }
}

template <typename CoefficientType, typename ComplexType>
Hdvf<CoefficientType, ComplexType>& HDVF_comput (const ComplexType& complex, const Options &options)
{
    Hdvf<CoefficientType, ComplexType>& hdvf(*(new Hdvf<CoefficientType, ComplexType>(complex, options.HDVF_opt)));
    std::vector<PairCell> pairs ;
    if (!options.random)
        pairs = hdvf.computePerfectHDVF(options.verbose);
    else
        pairs = hdvf.computeRandPerfectHDVF(options.verbose);
    
    if (options.with_output)
    {
        cout << "----> pairs found by computePerfectHDVF" << endl ;
        hdvf.print_pairs(pairs) ;
        cout << "----> reduction" << endl ;
        hdvf.print_reduction() ;
    }
    if (options.with_export)
    {
        string file(options.outfile_root+"_reduction.txt") ;
        std::ofstream out ( file, std::ios::out | std::ios::trunc);

        if ( not out . good () ) {
            std::cerr << "hdvf: with_export. Fatal Error:\n  " << file << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        out << "----> pairs found by computePerfectHDVF" << endl ;
        hdvf.print_pairs(pairs, out) ;
        out << "----> reduction" << endl ;
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
        using ComplexType = Abstract_simplicial_chain_complex<CoefficientType>  ;
        using HDVFType = Hdvf<CoefficientType, ComplexType> ;
        
        // MeshObject
        Mesh_object mesh ;
        mesh.read_simp(options.in_file) ;
        
        // Complex
        ComplexType complex(mesh);
        
        mesh_complex_output<Mesh_object, ComplexType>(mesh, complex, options) ;
        
        // Hdvf computation, export, output
        HDVFType hdvf(HDVF_comput<CoefficientType,ComplexType>(complex, options)) ;
        
        // Export to vtk
        // None for SIMP format
    }
    /// OFF format
    else if (options.in_format == InputFormat::OFF)
    {
        using ComplexType = Simplicial_chain_complex<CoefficientType> ;
        using HDVFType = Hdvf<CoefficientType, ComplexType> ;
        
        // MeshObject
        Mesh_object mesh ;
        mesh.read_off(options.in_file) ;
        
        // Complex
        ComplexType complex(mesh, mesh.nodes);
        
        mesh_complex_output<Mesh_object, ComplexType>(mesh, complex, options) ;
        
        // Hdvf computation, export, output
        HDVFType hdvf(HDVF_comput<CoefficientType,ComplexType>(complex, options)) ;
        
        // Loop on operations with vtk export
        if (options.loop)
        {
            auto output_vtk_simp = [options](HDVFType &hdvf, ComplexType& complex)
            {
                Simp_output_vtk<CoefficientType>(hdvf, complex, options.outfile_root) ;
            } ;
            
            interaction_loop<CoefficientType, ComplexType>(hdvf, complex, output_vtk_simp) ;
        }
        // Export to vtk
        else if (options.with_vtk_export)
        {
            cout << "----> exporting to vtk" << endl ;
            Simp_output_vtk<CoefficientType>(hdvf, complex, options.outfile_root) ;
        }
    }
    // CubComplex
    else if ((options.in_format == InputFormat::PGM) || (options.in_format == InputFormat::CUB))
    {
        using ComplexType = Cubical_chain_complex<CoefficientType> ;
        using HDVFType = Hdvf<CoefficientType, ComplexType> ;
        
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
        
        // Hdvf computation, export, output
        HDVFType hdvf(HDVF_comput<CoefficientType,ComplexType>(complex, options)) ;
        
        // Loop on operations with vtk export
        if (options.loop)
        {
            auto output_vtk_cub = [options](HDVFType &hdvf, ComplexType& complex)
            {
                Cub_output_vtk<CoefficientType>(hdvf, complex, options.outfile_root) ;
            } ;
            
            interaction_loop<CoefficientType, ComplexType>(hdvf, complex, output_vtk_cub) ;
        }
        // Export to vtk
        else if (options.with_vtk_export)
        {
            cout << "----> exporting to vtk" << endl ;
            Cub_output_vtk<CoefficientType>(hdvf, complex, options.outfile_root) ;
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









         
    



    
