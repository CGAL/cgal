#include <iostream>
#include <fstream>
#include <chrono>
#include <CGAL/Hdvf/Cub_object_io.h>
#include <CGAL/Hdvf/Cubical_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/Hdvf/Hdvf.h>
#include <CGAL/OSM/OSM.h>

typedef int CoefficientType;
//typedef CGAL::HDVF::Zp<5> CoefficientType;
//typedef CGAL::HDVF::Z2 CoefficientType;

int main(int argc, char **argv)
{
    using ComplexType = CGAL::HDVF::Cubical_chain_complex<CoefficientType> ;
    using HDVFType = CGAL::HDVF::Hdvf<ComplexType> ;

    if (argc != 2)
    {
        std::cout << "usage: test_CubObject cub_file" << std::endl;
    }
    else
    {
        // Load cub object
        CGAL::HDVF::Cub_object_io mesh ;
        mesh.read_cub(argv[1], true);

        mesh.print_infos();

        // Build complex (PRIMAL construction)
        ComplexType complex(mesh, ComplexType::PRIMAL);

        std::cout << complex;

        // Build empty HDVF
        HDVFType hdvf(complex, CGAL::HDVF::OPT_FULL) ;

        // Build HDVF step by step
        //        hdvf.A(7,0,1);
        //        hdvf.A(2,5,0);
        //        hdvf.A(3,8,0);
        hdvf.A(5,0,1);
        std::cout << "Step 1: is_perfect_hdvf " << hdvf.is_perfect_hdvf() << std::endl ;

        hdvf.compute_perfect_hdvf();
        //        hdvf.compute_rand_perfect_hdvf();
        std::cout << "Step 2: is_perfect_hdvf " << hdvf.is_perfect_hdvf() << std::endl ;

        // Output HDVF to console
        hdvf.insert_matrices();
        hdvf.insert_reduction();

        // Output HDVF to vtk
        hdvf_geometric_chain_complex_output_vtk(hdvf, complex, "res") ;

        // Test get_annotation

        // Compute the annotation of cycle1 in the homology basis
        std::vector<std::vector<size_t> > crit =  hdvf.flag(CGAL::HDVF::CRITICAL);
        std::vector<size_t> criticals = hdvf.flag_dim(CGAL::HDVF::CRITICAL, 1) ;
        HDVFType::Column_chain cycle1(hdvf.homology_chain(criticals.at(0), 1)) ;
        HDVFType::Column_chain annot1(hdvf.get_annotation(cycle1,1));
        std::cout << "Cycle1:" << cycle1 << std::endl ;
        std::cout << "Annotation of cycle 1: " << annot1 << std::endl ;

        // Compute the annotation of cycle2 (outer cycle) in the homology basis
        HDVFType::Column_chain cycle2(complex.number_of_cells(1)) ;
        cycle2.set_coefficient(0, 1) ;
        cycle2.set_coefficient(1, 1) ;
        cycle2.set_coefficient(2, 1) ;
        cycle2.set_coefficient(3, 1) ;
        cycle2.set_coefficient(4, 1) ;
        cycle2.set_coefficient(9, 1) ;
        cycle2.set_coefficient(7, 1) ;
        cycle2.set_coefficient(10, 1) ;
        cycle2.set_coefficient(11, 1) ;
        cycle2.set_coefficient(12, 1) ;
        HDVFType::Column_chain annot2(hdvf.get_annotation(cycle2,1));
        std::cout << "Cycle2:" << cycle1 << std::endl ;
        ComplexType::chain_complex_chain_to_vtk(complex, "cycle2.vtk", cycle2, 1) ;
        std::cout << "Annotation of cycle 2: " << annot2 << std::endl ;

        // Test get_coannotation

        // Compute the co-annotation of cycle1 in the cohomology basis
        HDVFType::Row_chain cocycle1((hdvf.cohomology_chain(criticals.at(0), 1)).transpose()) ;
        HDVFType::Row_chain coannot1(hdvf.get_coannotation(cocycle1,1));
        std::cout << "Co-cycle1:" << cocycle1 << std::endl ;
        std::cout << "Co-annotation of co-cycle 1: " << coannot1 << std::endl ;


        // Compute the co-annotation of cycle2 (outer cycle) in the homology basis
        HDVFType::Row_chain cocycle2(complex.number_of_cells(1)) ;
        cocycle2.set_coefficient(5, 1) ;
        cocycle2.set_coefficient(6, 1) ;
        HDVFType::Row_chain coannot2(hdvf.get_coannotation(cocycle2,1));
        std::cout << "Co-cycle2: " << cocycle2 << std::endl ;
        std::cout << "Co-annotation of co-cycle 2: " << coannot2 << std::endl ;

        // Test are_same_cycles
        HDVFType::Column_chain cycle3(cycle1) ;
        cycle3 += CGAL::OSM::cget_column(complex.boundary_matrix(2), 0); // Add the boundary of the 2-cell
        std::cout << "Cycle3: " << cycle3 << std::endl ;
        ComplexType::chain_complex_chain_to_vtk(complex, "cycle3.vtk", cycle3, 1) ;
        std::cout << "are_same_cycles cycle1 and cycle3: " << hdvf.are_same_cycles(cycle1, cycle3, 1) << std::endl ;

        HDVFType::Column_chain cycle4(cycle3) ;
        cycle4 += hdvf.homology_chain(criticals.at(1), 1) ; // Cycle4: cycle3 + second hole
        std::cout << "Cycle4: " << cycle4 << std::endl ;
        ComplexType::chain_complex_chain_to_vtk(complex, "cycle4.vtk", cycle4, 1) ;
        std::cout << "are_same_cycles cycle1 and cycle4: " << hdvf.are_same_cycles(cycle1, cycle4, 1) << std::endl ;

        // Test are_same_cocycles
        HDVFType::Row_chain cocycle3(cocycle1) ;
        cocycle3 += CGAL::OSM::get_row(complex.boundary_matrix(1), 0); // Add the coboundary of 0-cell
        std::cout << "Coycle3: " << cocycle3 << std::endl ;
        std::cout << "are_same_cocycles cocycle1 and cocycle3: " << hdvf.are_same_cocycles(cocycle1, cocycle3, 1) << std::endl ;

        HDVFType::Row_chain cocycle4(cocycle3) ;
        cocycle4 += (hdvf.cohomology_chain(criticals.at(1), 1).transpose()) ; // Cycle4: cycle3 + second cohomology generator
        std::cout << "Cocycle4: " << cocycle4 << std::endl ;
        std::cout << "are_same_cocycles cocycle1 and cocycle4: " << hdvf.are_same_cocycles(cocycle1, cocycle4, 1) << std::endl ;

        std::cout << "--------------" << std::endl ;

        typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::COLUMN> Column_chain;
        typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::ROW> Row_chain;
        typedef CGAL::OSM::Sparse_matrix<int, CGAL::OSM::COLUMN> Column_matrix;
        typedef CGAL::OSM::Sparse_matrix<int, CGAL::OSM::ROW> Row_matrix;
        // Create a column-major sparse matrix
        Column_matrix M(5,4) ;

        // Fill coefficients
        CGAL::OSM::set_coefficient(M, 0, 1, 1) ;
        CGAL::OSM::set_coefficient(M, 0, 2, -1) ;
        CGAL::OSM::set_coefficient(M, 2, 1, 2) ;

        // Iterate over non empty columns
         for(CGAL::OSM::Bitboard::iterator it_col = M.begin(); it_col != M.end(); ++it_col)
            {
                std::cout << "col: " << *it_col << std::endl ;
                // Get a constant reference over the column (complexity O(1))
                const Column_chain& col(CGAL::OSM::cget_column(M, *it_col));
                // Iterate over the column
                for (Column_chain::const_iterator it = col.begin(); it != col.end(); ++it)
                {
                    std::cout << "row: " << it->first << " - coef: " << it->second << std::endl ;
                }
            }
        // Direct output of the matrix with << operator
        std::cout << "M: " << M << std::endl;

        // Create a row-major sparse matrix
        Row_matrix MM(5,4) ;

        // Fill coefficients
        CGAL::OSM::set_coefficient(MM, 0, 1, 1) ;
        CGAL::OSM::set_coefficient(MM, 0, 2, -1) ;
        CGAL::OSM::set_coefficient(MM, 2, 1, 2) ;

        // Test write_matrix
        const std::string filename("test.osm") ;

        std::ofstream out ( filename, std::ios::out | std::ios::trunc);
        if ( ! out . good () ) {
            std::cerr << "Out fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }

        write_matrix(M, out);
        write_matrix(MM, out);

        out.close();

        std::ifstream in ( filename );
        if ( ! in . good () ) {
            std::cerr << "In fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        Column_matrix M2 ;
        Row_matrix MM2 ;

        read_matrix(M2, in) ;
        read_matrix(MM2, in) ;

        in.close();

        out << "M:" << M << std::endl ;
        out << "M2:" << M2 << std::endl ;
        out << "MM:" << MM << std::endl ;
        out << "MM2:" << MM2 << std::endl ;

        Column_matrix M3 ;
        write_matrix(M3,out);

        // Test save_hdvf

        // Test save_hdvf
        const std::string filename2("test.hdvf") ;

        std::ofstream out2 ( filename2, std::ios::out | std::ios::trunc);
        if ( ! out2 . good () ) {
            std::cerr << "Out fatal Error:\n  " << filename2 << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        hdvf.insert_hdvf_reduction(out2) ;

        out2.close() ;

        std::ifstream in2 ( filename2 );
        if ( ! in2 . good () ) {
            std::cerr << "In fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }

        HDVFType hdvf2(complex);

        hdvf2.extract_hdvf_reduction(in2) ;

        in2.close();

        out << "=============== HDVF" << std::endl ;
        hdvf.insert_reduction() ;
        out << "=============== HDVF2" << std::endl ;
        hdvf2.insert_reduction() ;
    }

    // Test == operators

    // CChains
    HDVFType::Column_chain c1(4), c2(4), c3 ;
    c1.set_coefficient(1,1) ;
    c1.set_coefficient(3,-1) ;

    c2.set_coefficient(3,-1) ;
    c2.set_coefficient(1,1) ;

    c3 = c2 ;
    c3.set_coefficient(2, 1) ;

    std::cout << "c1: " << c1 << std::endl ;
    std::cout << "c2: " << c2 << std::endl ;
    std::cout << "c3: " << c3 << std::endl ;
    std::cout << "c1 == c2 : " << (c1 == c2) << std::endl ;
    std::cout << "c1 == c3 : " << (c1 == c3) << std::endl ;

    // CMatrices
    HDVFType::Column_matrix M1(3,4), M2(3,4), M3(3,4) ;
    CGAL::OSM::set_coefficient(M1, 0, 1, 1) ;
    CGAL::OSM::set_coefficient(M1, 0, 2, -1) ;
    CGAL::OSM::set_coefficient(M1, 1, 1, 2) ;
    CGAL::OSM::set_coefficient(M1, 2, 1, -2) ;

    CGAL::OSM::set_coefficient(M2, 0, 1, 1) ;
    CGAL::OSM::set_coefficient(M2, 2, 1, -2) ;
    CGAL::OSM::set_coefficient(M2, 0, 2, -1) ;
    CGAL::OSM::set_coefficient(M2, 1, 1, 2) ;

    M3 = M2 ;
    CGAL::OSM::set_coefficient(M3, 2, 2, 3) ;

    std::cout << "M1 == M2 : " << (M1 == M2) << std::endl ;
    std::cout << "M1 == M3 : " << (M1 == M3) << std::endl ;

    // RChains
    HDVFType::Row_chain cc1(4), cc2(4), cc3 ;
    cc1.set_coefficient(1,1) ;
    cc1.set_coefficient(3,-1) ;

    cc2.set_coefficient(3,-1) ;
    cc2.set_coefficient(1,1) ;

    cc3 = cc2 ;
    cc3.set_coefficient(2, 1) ;

    std::cout << "cc1: " << cc1 << std::endl ;
    std::cout << "cc2: " << cc2 << std::endl ;
    std::cout << "cc3: " << cc3 << std::endl ;
    std::cout << "cc1 == cc2 : " << (cc1 == cc2) << std::endl ;
    std::cout << "cc1 == cc3 : " << (cc1 == cc3) << std::endl ;

    // CMatrices
    HDVFType::Row_matrix MM1(3,4), MM2(3,4), MM3(3,4) ;
    CGAL::OSM::set_coefficient(MM1, 0, 1, 1) ;
    CGAL::OSM::set_coefficient(MM1, 0, 2, -1) ;
    CGAL::OSM::set_coefficient(MM1, 1, 1, 2) ;
    CGAL::OSM::set_coefficient(MM1, 2, 1, -2) ;

    CGAL::OSM::set_coefficient(MM2, 0, 1, 1) ;
    CGAL::OSM::set_coefficient(MM2, 2, 1, -2) ;
    CGAL::OSM::set_coefficient(MM2, 0, 2, -1) ;
    CGAL::OSM::set_coefficient(MM2, 1, 1, 2) ;

    MM3 = MM2 ;
    CGAL::OSM::set_coefficient(MM3, 2, 2, 3) ;

    std::cout << "MM1 == MM2 : " << (MM1 == MM2) << std::endl ;
    std::cout << "MM1 == MM3 : " << (MM1 == MM3) << std::endl ;

    return 0;
}
