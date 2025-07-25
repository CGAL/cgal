#include <iostream>
#include <fstream>
#include <chrono>
#include "CGAL/Hdvf/tools_io.hpp"
#include "CGAL/Hdvf/Cubical_chain_complex.hpp"
#include "CGAL/HDVF/Geometric_chain_complex_tools.h"
#include "CGAL/HDVF/Zp.hpp"
#include "CGAL/Hdvf/Hdvf.h"

#include "CGAL/OSM/OSM.hpp"

using namespace CGAL;
using namespace HDVF;

//typedef int CoefficientType;
typedef Zp<2> CoefficientType;

int main(int argc, char **argv)
{
    using ComplexType = Cubical_chain_complex<CoefficientType> ;
    using HDVFType = Hdvf<CoefficientType, ComplexType> ;
    
    if (argc != 2)
    {
        std::cout << "usage: test_CubObject cub_file" << std::endl;
    }
    else
    {
        // Load cub object
        Cub_object mesh ;
        mesh.read_cub(argv[1], true);
        
        mesh.print_infos();
        
        // Build complex (PRIMAL construction)
        ComplexType complex(mesh, ComplexType::PRIMAL);
        
        complex.print_complex();
        
        // Build empty HDVF
        HDVFType hdvf(complex, OPT_FULL) ;
        
        // Build HDVF step by step
        //        hdvf.A(7,0,1);
        //        hdvf.A(2,5,0);
        //        hdvf.A(3,8,0);
        hdvf.A(5,0,1);
        cout << "Step 1: is_perfect_hdvf " << hdvf.is_perfect_hdvf() << endl ;
        
        hdvf.compute_perfect_hdvf();
        //        hdvf.compute_rand_perfect_hdvf();
        cout << "Step 2: is_perfect_hdvf " << hdvf.is_perfect_hdvf() << endl ;
        
        // Output HDVF to console
        hdvf.print_matrices();
        hdvf.print_reduction();
        
        // Output HDVF to vtk
        hdvf_geometric_chain_complex_output_vtk<CoefficientType, ComplexType>(hdvf, complex, "res") ;
        
        // Test get_annotation
        
        // Compute the annotation of cycle1 in the homology basis
        std::vector<std::vector<size_t> > crit =  hdvf.get_flag(CRITICAL);
        std::vector<size_t> criticals = hdvf.get_flag_dim(CRITICAL, 1) ;
        HDVFType::CChain cycle1(hdvf.get_homology_chain(criticals.at(0), 1)) ;
        HDVFType::CChain annot1(hdvf.get_annotation(cycle1,1));
        cout << "Cycle1:" << cycle1 << endl ;
        cout << "Annotation of cycle 1: " << annot1 << endl ;
        
        // Compute the annotation of cycle2 (outer cycle) in the homology basis
        HDVFType::CChain cycle2(complex.nb_cells(1)) ;
        cycle2.set_coef(0, 1) ;
        cycle2.set_coef(1, 1) ;
        cycle2.set_coef(2, 1) ;
        cycle2.set_coef(3, 1) ;
        cycle2.set_coef(4, 1) ;
        cycle2.set_coef(9, 1) ;
        cycle2.set_coef(7, 1) ;
        cycle2.set_coef(10, 1) ;
        cycle2.set_coef(11, 1) ;
        cycle2.set_coef(12, 1) ;
        HDVFType::CChain annot2(hdvf.get_annotation(cycle2,1));
        cout << "Cycle2:" << cycle1 << endl ;
        ComplexType::chain_complex_chain_to_vtk(complex, "cycle2.vtk", cycle2, 1) ;
        cout << "Annotation of cycle 2: " << annot2 << endl ;
        
        // Test get_coannotation
        
        // Compute the co-annotation of cycle1 in the cohomology basis
        HDVFType::RChain cocycle1((hdvf.get_cohomology_chain(criticals.at(0), 1)).transpose()) ;
        HDVFType::RChain coannot1(hdvf.get_coannotation(cocycle1,1));
        cout << "Co-cycle1:" << cocycle1 << endl ;
        cout << "Co-annotation of co-cycle 1: " << coannot1 << endl ;
        
        
        // Compute the co-annotation of cycle2 (outer cycle) in the homology basis
        HDVFType::RChain cocycle2(complex.nb_cells(1)) ;
        cocycle2.set_coef(5, 1) ;
        cocycle2.set_coef(6, 1) ;
        HDVFType::RChain coannot2(hdvf.get_coannotation(cocycle2,1));
        cout << "Co-cycle2: " << cocycle2 << endl ;
        cout << "Co-annotation of co-cycle 2: " << coannot2 << endl ;
        
        // Test are_same_cycles
        HDVFType::CChain cycle3(cycle1) ;
        cycle3 += OSM::cget_column(complex.get_bnd_matrix(2), 0); // Add the boundary of the 2-cell
        cout << "Cycle3: " << cycle3 << endl ;
        ComplexType::chain_complex_chain_to_vtk(complex, "cycle3.vtk", cycle3, 1) ;
        cout << "are_same_cycles cycle1 and cycle3: " << hdvf.are_same_cycles(cycle1, cycle3, 1) << endl ;
        
        HDVFType::CChain cycle4(cycle3) ;
        cycle4 += hdvf.get_homology_chain(criticals.at(1), 1) ; // Cycle4: cycle3 + second hole
        cout << "Cycle4: " << cycle4 << endl ;
        ComplexType::chain_complex_chain_to_vtk(complex, "cycle4.vtk", cycle4, 1) ;
        cout << "are_same_cycles cycle1 and cycle4: " << hdvf.are_same_cycles(cycle1, cycle4, 1) << endl ;
        
        // Test are_same_cocycles
        HDVFType::RChain cocycle3(cocycle1) ;
        cocycle3 += OSM::get_row(complex.get_bnd_matrix(1), 0); // Add the coboundary of 0-cell
        cout << "Coycle3: " << cocycle3 << endl ;
        cout << "are_same_cocycles cocycle1 and cocycle3: " << hdvf.are_same_cocycles(cocycle1, cocycle3, 1) << endl ;
        
        HDVFType::RChain cocycle4(cocycle3) ;
        cocycle4 += (hdvf.get_cohomology_chain(criticals.at(1), 1).transpose()) ; // Cycle4: cycle3 + second cohomology generator
        cout << "Cocycle4: " << cocycle4 << endl ;
        cout << "are_same_cocycles cocycle1 and cocycle4: " << hdvf.are_same_cocycles(cocycle1, cocycle4, 1) << endl ;
        
        cout << "--------------" << endl ;
        
        typedef CGAL::OSM::Sparse_chain<int, OSM::COLUMN> CChain;
        typedef CGAL::OSM::Sparse_matrix<int, OSM::COLUMN> CMatrix;
        typedef CGAL::OSM::Sparse_matrix<int, OSM::ROW> RMatrix;
        // Create a column-major sparse matrix
        CMatrix M(5,4) ;

        // Fill coefficients
        OSM::set_coef(M, 0, 1, 1) ;
        OSM::set_coef(M, 0, 2, -1) ;
        OSM::set_coef(M, 2, 1, 2) ;

        // Iterate over non empty columns
         for(OSM::Bitboard::iterator it_col = M.begin(); it_col != M.end(); ++it_col)
            {
                cout << "col: " << *it_col << endl ;
                // Get a constant reference over the column (complexity O(1))
                const CChain& col(OSM::cget_column(M, *it_col));
                // Iterate over the column
                for (CChain::const_iterator it = col.begin(); it != col.end(); ++it)
                {
                    std::cout << "row: " << it->first << " - coef: " << it->second << std::endl ;
                }
            }
        // Direct output of the matrix with << operator
        std::cout << "M: " << M << std::endl;
        
        // Create a row-major sparse matrix
        RMatrix MM(5,4) ;

        // Fill coefficients
        OSM::set_coef(MM, 0, 1, 1) ;
        OSM::set_coef(MM, 0, 2, -1) ;
        OSM::set_coef(MM, 2, 1, 2) ;
        
        // Test write_matrix
        const string filename("test.osm") ;
        
        std::ofstream out ( filename, std::ios::out | std::ios::trunc);
        if ( not out . good () ) {
            std::cerr << "Out fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        
        write_matrix(M, out);
        write_matrix(MM, out);
        
        out.close();
        
        std::ifstream in ( filename );
        if ( not in . good () ) {
            std::cerr << "In fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        CMatrix M2 ;
        RMatrix MM2 ;
        
        read_matrix(M2, in) ;
        read_matrix(MM2, in) ;
        
        in.close();
        
        cout << "M:" << M << endl ;
        cout << "M2:" << M2 << endl ;
        cout << "MM:" << MM << endl ;
        cout << "MM2:" << MM2 << endl ;
        
        CMatrix M3 ;
        write_matrix(M3,cout);
        
        // Test save_hdvf
        
        // Test save_hdvf
        const string filename2("test.hdvf") ;
        
        std::ofstream out2 ( filename2, std::ios::out | std::ios::trunc);
        if ( not out2 . good () ) {
            std::cerr << "Out fatal Error:\n  " << filename2 << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        hdvf.save_hdvf_reduction(out2) ;
        
        out2.close() ;
        
        std::ifstream in2 ( filename2 );
        if ( not in2 . good () ) {
            std::cerr << "In fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        
        HDVFType hdvf2(complex);
        
        hdvf2.load_hdvf_reduction(in2) ;
        
        in2.close();
        
        cout << "=============== HDVF" << endl ;
        hdvf.print_reduction() ;
        cout << "=============== HDVF2" << endl ;
        hdvf2.print_reduction() ;
    }
    
    return 0;
}
