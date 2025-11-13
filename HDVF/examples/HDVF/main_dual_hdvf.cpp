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
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
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
#include <CGAL/HDVF/Z2.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

// ------- A ring
// For Z/nZ other than Z (ie. n=0) and Z/2Z, uncomment and set the following define properly

//#define SCALAR 5

using Kernel = CGAL::Simple_cartesian<double>;
using Traits = HDVF::Hdvf_traits_3<Kernel>;

template <typename MeshType, typename Complex>
void mesh_complex_output(const MeshType& mesh, const Complex& L, const HDVF::Sub_chain_complex_mask<Complex>& K, const Options& options)
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

inline std::ostream& dual_pairs_output(const std::vector<HDVF::Cell_pair>& pairs, std::ostream& out=std::cout)
{
    out << "Pairs found by compute_perfect_hdvf:" << std::endl;
    for (const auto& pair : pairs) {
        out << "Sigma: " << pair.sigma << ", Tau: " << pair.tau << ", Dim: " << pair.dim << std::endl;
    }
    return out ;
}

template <typename Complex>
void dual_HDVF_pair (HDVF::Hdvf_duality<Complex>& dual_hdvf, const Options &options)
{
    // Compute pairing
    std::vector<HDVF::Cell_pair> pairs = dual_hdvf.compute_pairing_hdvf() ;

    if (options.with_output)
    {
        dual_pairs_output(pairs) ;
    }
    if (options.with_export)
    {
        std::string file(options.outfile_root+"_pairs.txt") ;
        std::ofstream out ( file, std::ios::out | std::ios::trunc);

        if ( ! out . good () ) {
            std::cerr << "hdvf: with_export. Fatal Error:\n  " << file << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }

        dual_pairs_output(pairs, out) ;
        out.close() ;
    }
}

template <typename Complex>
HDVF::Hdvf_duality<Complex>& dual_HDVF_comput (const Complex& L,  HDVF::Sub_chain_complex_mask<Complex>& K, const Options &options)
{
    using Coefficient_ring = typename Complex::Coefficient_ring;
    using HDVF_type = HDVF::Hdvf_duality<Complex> ;
    using SubCCType = HDVF::Sub_chain_complex_mask<Complex> ;

    HDVF_type& hdvf(*(new HDVF_type(L, K, options.HDVF_opt)));

    std::cout << "----> START computing dual HDVF" << std::endl ;
    if (options.random)
        hdvf.compute_rand_perfect_hdvf() ;
    else
        hdvf.compute_perfect_hdvf() ;
    std::cout << "------> END computing dual HDVF" << std::endl ;

    if (options.with_output)
    {
        std::cout << "----> reduction" << std::endl ;
        hdvf.write_reduction() ;
    }
    if (options.with_export)
    {
        std::cout << "----> exporting..." << std::endl ;
        std::string file(options.outfile_root+"_reduction.txt") ;
        std::ofstream out ( file, std::ios::out | std::ios::trunc);

        if ( ! out . good () ) {
            std::cerr << "hdvf: with_export. Fatal Error:\n  " << file << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        out << "----> reduction" << std::endl ;
        hdvf.write_reduction(out) ;

        out.close() ;
    }
    return hdvf ;
}


template <typename Coefficient_ring>
void main_code (const Options &options)
{
#ifndef CUST_FILTRATION
    // Standard lower star filtration along x,y or z
    using Degree = double ;
#else
    // TODO
#endif

    /// SIMP format
    if (options.in_format == InputFormat::SIMP)
    {
//        using Complex = AbstractSimpComplex<Coefficient_ring>  ;
//        using HDVF_type = HDVF<Coefficient_ring, Complex> ;
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
//        HDVF_type hdvf(HDVF_comput<Coefficient_ring,Complex>(complex, options)) ;
//
//        // Export to vtk
//        // None for SIMP format
        std::cout << "not yet..." << std::endl ;
        throw("not yet") ;
    }
    /// OFF format
    else if (options.in_format == InputFormat::OFF)
    {
        using Complex = HDVF::Simplicial_chain_complex<Coefficient_ring,Traits> ;
        using HDVF_type = HDVF::Hdvf_duality<Complex> ;
        using ToolsType = HDVF::Duality_simplicial_complex_tools<Coefficient_ring, Traits> ;
        using SubCCType = HDVF::Sub_chain_complex_mask<Complex> ;

        // MeshObject
        HDVF::Mesh_object_io<Traits> mesh ;
        mesh.read_off(options.in_file) ;

        // Build L (bounding sphere meshed with tetgen), K and L-K

        typename ToolsType::Complex_duality_data t(ToolsType::dualize_complex(mesh)) ;
        Complex& L(t.L) ;
        SubCCType& K(t.K) ;

        // Output/export mesh and complex

        mesh_complex_output<HDVF::Mesh_object_io<Traits>, Complex>(mesh, L, K, options) ;

        // HDVF computation, export, output
        HDVF_type& hdvf(dual_HDVF_comput<Complex>(L, K, options)) ;

        // Export to vtk
        if (options.with_vtk_export)
        {
            std::cout << "----> exporting to vtk" << std::endl ;
            // K
            {
                hdvf.set_mask_K() ;
                CGAL::IO::write_VTK(hdvf, L, options.outfile_root+"_complex_K", options.co_faces) ;
            }
            // L-K
            {
                hdvf.set_mask_L_K() ;
                CGAL::IO::write_VTK(hdvf, L, options.outfile_root+"_cocomplex_L_K", options.co_faces) ;
            }
        }

        // Compute pairing
        dual_HDVF_pair<Complex>(hdvf, options) ;
    }
    // CubComplex
    else if ((options.in_format == InputFormat::PGM) || (options.in_format == InputFormat::CUB))
    {
        using Complex = HDVF::Cubical_chain_complex<Coefficient_ring, Traits> ;
        using HDVF_type = HDVF::Hdvf_duality<Complex> ;
        using SubCCType = HDVF::Sub_chain_complex_mask<Complex> ;
        using ToolsType = HDVF::Duality_cubical_complex_tools<Coefficient_ring, Traits> ;

        HDVF::Cub_object_io<Traits> mesh ;
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

        // Frame (add 1 voxel around data)
        if (options.with_frame)
        {
            mesh.frame() ;
        }

        // Complex
        Complex* complex = new Complex(mesh, primal_dual);

        // Build L, K and L-K

        std::pair<Complex&, SubCCType&> p(ToolsType::dualize_complex(*complex)) ;
        delete complex ;
        Complex &L(p.first) ;
        SubCCType &K(p.second) ;

        // Output/export mesh and complex

        mesh_complex_output<HDVF::Cub_object_io<Traits>, Complex>(mesh, L, K, options) ;

        // HDVF computation, export, output
        HDVF_type& hdvf(dual_HDVF_comput<Complex>(L, K, options)) ;

        // Export to vtk
        if (options.with_vtk_export)
        {
            std::cout << "----> exporting to vtk" << std::endl ;
            // K
            {
                hdvf.set_mask_K() ;
                CGAL::IO::write_VTK(hdvf, L, options.outfile_root+"_complex_K", options.co_faces) ;
            }
            // L-K
            {
                hdvf.set_mask_L_K() ;
                CGAL::IO::write_VTK(hdvf, L, options.outfile_root+"_cocomplex_L_K", options.co_faces) ;
            }
        }

        // Compute pairing
        dual_HDVF_pair<Complex>(hdvf, options) ;
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

        // ----- Definition of the Coefficient_ring
#ifndef SCALAR
        if (options.scalar == 0)
        {
            using Coefficient_ring = int ;
            main_code<Coefficient_ring>(options) ;
        }
        else if (options.scalar == 2)
        {
            using Coefficient_ring = HDVF::Z2 ;
            main_code<Coefficient_ring>(options) ;
        }
        else
        {
            std::cerr << "Z" << options.scalar << " not instantiated, use the #define at line 27" << std::endl ;
        }
#else
        typedef HDVF::Zp<SCALAR,int,true> Coefficient_ring;
#endif
    }

    return 0 ;
}

