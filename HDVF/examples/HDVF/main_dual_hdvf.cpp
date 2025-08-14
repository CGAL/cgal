// Dual HDVF computation (command line version)
// --------
// Computes a "perfect" dual HDVF from a object
//      (the object is embedded into a "ball")
//      and provides a batch mode to specify arguments
// For help: dual_hdvf -h
// --------
// A. Bac
// --------

#include <iostream>
#include <CGAL/HDVF/Simplex.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Cubical_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/hdvf_tools.h>
#include <CGAL/HDVF/Hdvf_duality.h>
#include <CGAL/OSM/OSM.h>
#include <CGAL/HDVF/Sub_sparse_matrix.h>
#include <CGAL/HDVF/hdvf_tools.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Cub_object_io.h>
#include <CGAL/HDVF/Icosphere_object_io.h>
#include <CGAL/HDVF/Tet_object_io.h>
#include "arguments.h"
#include <CGAL/HDVF/Zp.h>

// ------- A ring
// For Z/nZ other than Z (ie. n=0) and Z/2Z, uncomment and set the following define properly

//#define SCALAR 5


template <typename CoefficientType, typename MeshType, typename ComplexType>
void mesh_complex_output(const MeshType& mesh, const ComplexType& L, const CGAL::HDVF::Sub_chain_complex_mask<CoefficientType, ComplexType>& K, const Options& options)
{
    if (options.with_output)
    {
        // Mesh
        std::cout << "----> mesh informations" << std::endl ;
        mesh.print_infos() ;

        // Complex
        std::cout << "----> complex informations" << std::endl ;
        std::cout << "------> complex L" << std::endl ;
        std::cout << L;
        std::cout << "------> subcomplex K" << std::endl ;
        std::cout << K << std::endl ;
    }
}

inline std::ostream& dual_pairs_output(const std::vector<CGAL::HDVF::Pair_cells>& pairs, std::ostream& out=std::cout)
{
    out << "Pairs found by compute_perfect_hdvf:" << std::endl;
    for (const auto& pair : pairs) {
        out << "Sigma: " << pair.sigma << ", Tau: " << pair.tau << ", Dim: " << pair.dim << std::endl;
    }
    return out ;
}

template <typename CoefficientType, typename ComplexType>
void dual_HDVF_pair (CGAL::HDVF::Hdvf_duality<CoefficientType, ComplexType>& dual_hdvf, const Options &options)
{
    // Compute pairing
    std::vector<CGAL::HDVF::Pair_cells> pairs = dual_hdvf.compute_pairing_hdvf() ;

    if (options.with_output)
    {
        dual_pairs_output(pairs) ;
    }
    if (options.with_export)
    {
        std::string file(options.outfile_root+"_pairs.txt") ;
        std::ofstream out ( file, std::ios::out | std::ios::trunc);

        if ( not out . good () ) {
            std::cerr << "hdvf: with_export. Fatal Error:\n  " << file << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }

        dual_pairs_output(pairs, out) ;
        out.close() ;
    }
}

template <typename CoefficientType, typename ComplexType>
CGAL::HDVF::Hdvf_duality<CoefficientType, ComplexType>& dual_HDVF_comput (const ComplexType& L,  CGAL::HDVF::Sub_chain_complex_mask<CoefficientType, ComplexType>& K, const Options &options)
{
    using HDVFType = CGAL::HDVF::Hdvf_duality<CoefficientType, ComplexType> ;
    using SubCCType = CGAL::HDVF::Sub_chain_complex_mask<CoefficientType, ComplexType> ;

    HDVFType& hdvf(*(new HDVFType(L, K, options.HDVF_opt)));

    std::cout << "----> START computing dual HDVF" << std::endl ;
    if (options.random)
        hdvf.compute_rand_perfect_hdvf() ;
    else
        hdvf.compute_perfect_hdvf() ;
    std::cout << "------> END computing dual HDVF" << std::endl ;

    if (options.with_output)
    {
        std::cout << "----> reduction" << std::endl ;
        hdvf.insert_reduction() ;
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
        out << "----> reduction" << std::endl ;
        hdvf.insert_reduction(out) ;

        out.close() ;
    }
    return hdvf ;
}


template <typename CoefficientType>
void main_code (const Options &options)
{
#ifndef CUST_FILTRATION
    // Standard lower star filtration along x,y or z
    using DegreeType = double ;
#else
    // TODO
#endif

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
        std::cout << "not yet..." << std::endl ;
        throw("not yet") ;
    }
    /// OFF format
    else if (options.in_format == InputFormat::OFF)
    {
        using ComplexType = CGAL::HDVF::Simplicial_chain_complex<CoefficientType> ;
        using HDVFType = CGAL::HDVF::Hdvf_duality<CoefficientType, ComplexType> ;
        using ToolsType = CGAL::HDVF::Duality_simplicial_complex_tools<CoefficientType> ;
        using SubCCType = CGAL::HDVF::Sub_chain_complex_mask<CoefficientType, ComplexType> ;

        // MeshObject
        CGAL::HDVF::Mesh_object_io mesh ;
        mesh.read_off(options.in_file) ;

        // Complex
        ComplexType* complex = new ComplexType(mesh, mesh.get_nodes());

        // Build L (bounding sphere meshed with tetgen), K and L-K

        typename ToolsType::TripleRes t(ToolsType::simplicial_chain_complex_bb(*complex)) ;
        delete complex ;
        ComplexType& L(t.L) ;
        SubCCType& K(t.K) ;

        // Output/export mesh and complex

        mesh_complex_output<CoefficientType, CGAL::HDVF::Mesh_object_io, ComplexType>(mesh, L, K, options) ;

        // HDVF computation, export, output
        HDVFType& hdvf(dual_HDVF_comput<CoefficientType,ComplexType>(L, K, options)) ;

        // Export to vtk
        if (options.with_vtk_export)
        {
            std::cout << "----> exporting to vtk" << std::endl ;
            // K
            {
                hdvf.set_mask_K() ;
                hdvf_duality_geometric_chain_complex_output_vtk(hdvf, L, options.outfile_root+"_complex_K", options.co_faces) ;
            }
            // L-K
            {
                hdvf.set_mask_L_K() ;
                hdvf_duality_geometric_chain_complex_output_vtk(hdvf, L, options.outfile_root+"_cocomplex_L_K", options.co_faces) ;
            }
        }

        // Compute pairing
        dual_HDVF_pair<CoefficientType,ComplexType>(hdvf, options) ;
    }
    // CubComplex
    else if ((options.in_format == InputFormat::PGM) || (options.in_format == InputFormat::CUB))
    {
        using ComplexType = CGAL::HDVF::Cubical_chain_complex<CoefficientType> ;
        using HDVFType = CGAL::HDVF::Hdvf_duality<CoefficientType, ComplexType> ;
        using SubCCType = CGAL::HDVF::Sub_chain_complex_mask<CoefficientType, ComplexType> ;
        using ToolsType = CGAL::HDVF::Duality_cubical_complex_tools<CoefficientType> ;

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

        // Frame (add 1 voxel around data)
        if (options.with_frame)
        {
            mesh.frame() ;
        }

        // Complex
        ComplexType* complex = new ComplexType(mesh, primal_dual);

        // Build L, K and L-K

        std::pair<ComplexType&, SubCCType&> p(ToolsType::cubical_chain_complex_bb(*complex)) ;
        delete complex ;
        ComplexType &L(p.first) ;
        SubCCType &K(p.second) ;

        // Output/export mesh and complex

        mesh_complex_output<CoefficientType, CGAL::HDVF::Cub_object_io, ComplexType>(mesh, L, K, options) ;

        // HDVF computation, export, output
        HDVFType& hdvf(dual_HDVF_comput<CoefficientType,ComplexType>(L, K, options)) ;

        // Export to vtk
        if (options.with_vtk_export)
        {
            std::cout << "----> exporting to vtk" << std::endl ;
            // K
            {
                hdvf.set_mask_K() ;
                hdvf_duality_geometric_chain_complex_output_vtk(hdvf, L, options.outfile_root+"_complex_K", options.co_faces) ;
            }
            // L-K
            {
                hdvf.set_mask_L_K() ;
                hdvf_duality_geometric_chain_complex_output_vtk(hdvf, L, options.outfile_root+"_cocomplex_L_K", options.co_faces) ;
            }
        }

        // Compute pairing
        dual_HDVF_pair<CoefficientType,ComplexType>(hdvf, options) ;
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
            using CoefficientType = CGAL::HDVF::Zp<2,int8_t> ;
            main_code<CoefficientType>(options) ;
        }
        else
        {
            std::cerr << "Z" << options.scalar << " not instantiated, use the #define at line 27" << std::endl ;
        }
#else
        typedef CGAL::HDVF::Zp<SCALAR> CoefficientType;
#endif
    }

    return 0 ;
}

