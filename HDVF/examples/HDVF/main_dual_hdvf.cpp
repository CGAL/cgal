// HDVF off files
// --------
// Computes a "perfect" HDVF from a simple .off file
// Prints selected pairs for A and the resulting reduction
// --------
// A. Bac
// --------

#include <iostream>
#include "CGAL/HDVF/Simplex.hpp"
#include "CGAL/HDVF/tools_io.hpp"
#include "CGAL/HDVF/Abstract_Simplicial_chain_complex.hpp"
#include "CGAL/HDVF/SimpComplexTools.hpp"
#include "CGAL/HDVF/Cubical_chain_complex.hpp"
#include "CGAL/HDVF/CubComplexTools.hpp"
#include "CGAL/HDVF/hdvf_tools.hpp"
#include "CGAL/HDVF/hdvf_duality.hpp"
#include "CGAL/OSM/OSM.hpp"
#include "CGAL/HDVF/SubSparseMatrix.hpp"
#include "CGAL/HDVF/hdvf_tools.hpp"
#include "arguments.h"
#include "CGAL/HDVF/Zp.hpp"

// ------- A ring
// For Z/nZ other than Z (ie. n=0) and Z/2Z, uncomment and set the following define properly

//#define SCALAR 5

using namespace CGAL ;
using namespace HDVF ;

template <typename CoefficientType, typename MeshType, typename ComplexType>
void mesh_complex_output(const MeshType& mesh, const ComplexType& L, const Sub_chain_complex_mask<CoefficientType, ComplexType>& K, const Options& options)
{
    if (options.with_output)
    {
        // Mesh
        cout << "----> mesh informations" << endl ;
        mesh.print_infos() ;
        
        // Complex
        cout << "----> complex informations" << endl ;
        cout << "------> complex L" << endl ;
        L.print_complex();
        cout << "------> subcomplex K" << endl ;
        cout << K << endl ;
    }
}

inline ostream& dual_pairs_output(const std::vector<PairCell>& pairs, ostream& out=cout)
{
    out << "Pairs found by computeDualPerfectHDVF:" << std::endl;
    for (const auto& pair : pairs) {
        out << "Sigma: " << pair.sigma << ", Tau: " << pair.tau << ", Dim: " << pair.dim << std::endl;
    }
    return out ;
}

template <typename CoefficientType, typename ComplexType>
void dual_HDVF_pair (Hdvf_duality<CoefficientType, ComplexType>& dual_hdvf, const Options &options)
{
    // Compute pairing
    std::vector<PairCell> pairs = dual_hdvf.computePairingHDVF() ;
    
    if (options.with_output)
    {
        dual_pairs_output(pairs) ;
    }
    if (options.with_export)
    {
        string file(options.outfile_root+"_pairs.txt") ;
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
Hdvf_duality<CoefficientType, ComplexType>& dual_HDVF_comput (const ComplexType& L,  Sub_chain_complex_mask<CoefficientType, ComplexType>& K, const Options &options)
{
    using HDVFType = Hdvf_duality<CoefficientType, ComplexType> ;
    using SubCCType = Sub_chain_complex_mask<CoefficientType, ComplexType> ;
    
    HDVFType& hdvf(*(new HDVFType(L, K, options.HDVF_opt)));
    
    cout << "----> START computing dual HDVF" << endl ;
    hdvf.computeDualPerfectHDVF() ;
    cout << "------> END computing dual HDVF" << endl ;
    
    if (options.with_output)
    {
        cout << "----> reduction" << endl ;
        hdvf.print_reduction() ;
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
        out << "----> reduction" << endl ;
        hdvf.print_reduction(out) ;
        
        out.close() ;
    }
    return hdvf ;
}


template <typename CoefficientType>
void main_code (const Options &options)
{
#ifndef CUST_FILTRATION
    // Standard lower star filtration along x,y or z
    using DegType = double ;
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
        cout << "not yet..." << endl ;
        throw("not yet") ;
    }
    /// OFF format
    else if (options.in_format == InputFormat::OFF)
    {
        using ComplexType = Simplicial_chain_complex<CoefficientType> ;
        using HDVFType = Hdvf_duality<CoefficientType, ComplexType> ;
        using ToolsType = Duality_simplicial_complex_tools<CoefficientType> ;
        using SubCCType = Sub_chain_complex_mask<CoefficientType, ComplexType> ;

        // MeshObject
        Mesh_object mesh ;
        mesh.read_off(options.in_file) ;
        
        // Complex
        ComplexType* complex = new ComplexType(mesh, mesh.nodes);
        
        // Build L (bounding sphere meshed with tetgen), K and L-K
        
        typename ToolsType::TripleRes t(ToolsType::SimpComplexBB(*complex)) ;
        delete complex ;
        ComplexType& L(t.L) ;
        SubCCType& K(t.K) ;
        
        // Output/export mesh and complex
        
        mesh_complex_output<CoefficientType, Mesh_object, ComplexType>(mesh, L, K, options) ;
        
        // HDVF computation, export, output
        HDVFType& hdvf(dual_HDVF_comput<CoefficientType,ComplexType>(L, K, options)) ;
        
        // Export to vtk
        if (options.with_vtk_export)
        {
            cout << "----> exporting to vtk" << endl ;
            // K
            {
                hdvf.set_mask_K() ;
                Simp_output_vtk<CoefficientType, OSM::Sparse_chain, OSM::SubSparseMatrix>(hdvf, L, options.outfile_root+"_complex_K") ;
            }
            // L-K
            {
                hdvf.set_mask_L_K() ;
                Simp_output_vtk<CoefficientType, OSM::Sparse_chain, OSM::SubSparseMatrix>(hdvf, L, options.outfile_root+"_cocomplex_L_K") ;
            }
        }
        
        // Compute pairing
        dual_HDVF_pair<CoefficientType,ComplexType>(hdvf, options) ;
    }
    // CubComplex
    else if ((options.in_format == InputFormat::PGM) || (options.in_format == InputFormat::CUB))
    {
        using ComplexType = Cubical_chain_complex<CoefficientType> ;
        using HDVFType = Hdvf_duality<CoefficientType, ComplexType> ;
        using SubCCType = Sub_chain_complex_mask<CoefficientType, ComplexType> ;
        using ToolsType = Duality_cubical_complex_tools<CoefficientType> ;
        
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
        
        // Frame (add 1 voxel around data)
        if (options.with_frame)
        {
            mesh.frame() ;
        }
        
        // Complex
        ComplexType* complex = new ComplexType(mesh, primal_dual);
        
        // Build L, K and L-K
        
        std::pair<ComplexType&, SubCCType&> p(ToolsType::CubComplexBB(*complex)) ;
        delete complex ;
        ComplexType &L(p.first) ;
        SubCCType &K(p.second) ;
        
        // Output/export mesh and complex
        
        mesh_complex_output<CoefficientType, Cub_object, ComplexType>(mesh, L, K, options) ;
        
        // HDVF computation, export, output
        HDVFType& hdvf(dual_HDVF_comput<CoefficientType,ComplexType>(L, K, options)) ;
        
        // Export to vtk
        if (options.with_vtk_export)
        {
            cout << "----> exporting to vtk" << endl ;
            // K
            {
                hdvf.set_mask_K() ;
                Cub_output_vtk<CoefficientType, OSM::Sparse_chain, OSM::SubSparseMatrix>(hdvf, L, options.outfile_root+"_complex_K") ;
            }
            // L-K
            {
                hdvf.set_mask_L_K() ;
                Cub_output_vtk<CoefficientType, OSM::Sparse_chain, OSM::SubSparseMatrix>(hdvf, L, options.outfile_root+"_cocomplex_L_K") ;
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

