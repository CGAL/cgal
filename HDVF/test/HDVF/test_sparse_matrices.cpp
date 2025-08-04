#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/OSM/OSM.h>

typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::COLUMN> CChain;
typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::ROW> RChain ;
typedef CGAL::OSM::Sparse_matrix<int, CGAL::OSM::COLUMN> CMatrix;
typedef CGAL::OSM::Sparse_matrix<int, CGAL::OSM::ROW> RMatrix;


int main(int argc, char **argv)
{
    std::cerr << "-- Test Sparse_matrices read/write" << std::endl;

    std::cerr << "----> Column matrices" << std::endl;

    CMatrix MC_rw(3,4);
    CGAL::OSM::set_coef(MC_rw, 0, 1, 1) ;
    CGAL::OSM::set_coef(MC_rw, 0, 2, -1) ;
    CGAL::OSM::set_coef(MC_rw, 1, 0, 2) ;
    CGAL::OSM::set_coef(MC_rw, 1, 2, -2) ;

    const std::string filename_col_save("test_col_mat.osm") ;

    std::ofstream out_col_save ( filename_col_save, std::ios::out | std::ios::trunc);
    if ( not out_col_save . good () ) {
        std::cerr << "Out fatal Error:\n  " << filename_col_save << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    write_matrix(MC_rw, out_col_save);

    out_col_save.close() ;

    CMatrix MC_rw2 ;

    std::ifstream in_col_save (filename_col_save);
    if ( not in_col_save . good () ) {
        std::cerr << "Out fatal Error:\n  " << filename_col_save << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }
    CGAL::OSM::read_matrix(MC_rw2, in_col_save);
    in_col_save.close() ;

    bool comp_col(MC_rw == MC_rw2) ;
    std::cerr << "saved and loaded matrices comparison: " << comp_col << std::endl ;
    assert(comp_col) ;

    std::cerr << "----> Row matrices" << std::endl;

    RMatrix MR_rw(3,4) ;
    CGAL::OSM::set_coef(MR_rw, 0, 1, 1) ;
    CGAL::OSM::set_coef(MR_rw, 0, 2, -1) ;
    CGAL::OSM::set_coef(MR_rw, 1, 0, 2) ;
    CGAL::OSM::set_coef(MR_rw, 1, 2, -2) ;

    const std::string filename_row_save("test_row_mat.osm") ;

    std::ofstream out_row_save ( filename_row_save, std::ios::out | std::ios::trunc);
    if ( not out_row_save . good () ) {
        std::cerr << "Out fatal Error:\n  " << filename_row_save << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    write_matrix(MR_rw, out_row_save);

    out_row_save.close() ;

    RMatrix MR_rw2 ;

    std::ifstream in_row_save (filename_row_save);
    if ( not in_row_save . good () ) {
        std::cerr << "Out fatal Error:\n  " << filename_row_save << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }
    CGAL::OSM::read_matrix(MR_rw2, in_row_save);
    in_row_save.close() ;

    bool comp_row(MR_rw == MR_rw2) ;
    std::cerr << "saved and loaded matrices comparison: " << comp_row << std::endl ;
    assert(comp_row) ;

    std::cerr << "-- Tests Sparse_matrices operations" << std::endl ;

    CChain a(4);
    RChain b(4);
    a.set_coef(1,1);
    a.set_coef(3,3);
    b.set_coef(0,2);
    b.set_coef(1,-3);

    CMatrix columnMajor = a * b, tmp = columnMajor;
    std::cout << columnMajor << std::endl ;

    CMatrix columnMajorRes(4,4) ;

    // How to read data in data/ directory ?

//    RMatrix rowMajor = a % b;
//
//    std::cout << "Création de matrice col-dominant:" << std::endl << columnMajor << std::endl << std::endl;
//    std::cout << "Suppr col 0:" << std::endl << tmp.delColumn(0) << std::endl << std::endl;
//    std::cout << "Suppr row 0:" << std::endl << tmp.delRow(1) << std::endl << std::endl;
//    std::cout << "Création de matrice row-dominant:" << std::endl << rowMajor << std::endl << std::endl;
//
//    std::cout << "Ajout de deux matrices CC: " << columnMajor + columnMajor << std::endl;
//    std::cout << "Soustraction de deux matrices CC: " << columnMajor - columnMajor << std::endl << std::endl;
//
//    std::cout << "Ajout de deux matrices RR: " << rowMajor + rowMajor << std::endl;
//    std::cout << "Soustraction de deux matrices RR: " << rowMajor - rowMajor << std::endl << std::endl;
//
//    std::cout << "Ajout de deux matrices CR: " << columnMajor + rowMajor << std::endl;
//    std::cout << "Soustraction de deux matrices CR: " << columnMajor - rowMajor << std::endl << std::endl;
//
//    std::cout << "Ajout de deux matrices RC: " << rowMajor + columnMajor << std::endl;
//    std::cout << "Soustraction de deux matrices RC: " << rowMajor - columnMajor << std::endl << std::endl;
//
//    std::cout << "Transposée de matrice:" << std::endl << "Col-dominant (devient row): " << columnMajor.transpose() << std::endl << "Row-dominant (devient col): " << rowMajor.transpose() << std::endl << std::endl;
//
//    std::cout << "Produits matriciels" << std::endl;
//
//    std::cout << "C*R = " << columnMajor * rowMajor << std::endl;
//    std::cout << "R*C = " << rowMajor * columnMajor << std::endl;
//    std::cout << "C*C = " << columnMajor * columnMajor << std::endl;
//    std::cout << "R*R = " << rowMajor * rowMajor << std::endl << std::endl;
//
//    std::cout << "C%R = " << columnMajor % rowMajor << std::endl;
//    std::cout << "R%C = " << rowMajor % columnMajor << std::endl;
//    std::cout << "C%C = " << columnMajor % columnMajor << std::endl;
//    std::cout << "R%R = " << rowMajor % rowMajor << std::endl << std::endl;
//
//    std::cout << "Récupération de colonne (col/row dom):" << OSM::getColumn(columnMajor, 0) << " " << OSM::getColumn(rowMajor, 0) << std::endl;
//    std::cout << "Récupération de ligne   (col/row dom):" << OSM::getRow(columnMajor, 1) << " " << OSM::getRow(rowMajor, 1) << std::endl << std::endl;
//
//    OSM::SparseMatrix colCol = columnMajor;
//    OSM::SparseMatrix colRow = rowMajor;
//    OSM::SparseMatrix rowCol = columnMajor;
//    OSM::SparseMatrix rowRow = rowMajor;
//
//    colCol *= columnMajor;
//    colRow *= rowMajor;
//    rowCol *= columnMajor;
//    rowRow *= rowMajor;
//
//    std::cout << "Multiplication de colonne (col/row dom):" << std::endl << colCol << std::endl << colRow << std::endl << std::endl;
//    std::cout << "Multiplication de ligne   (col/row dom):" << std::endl << rowCol << std::endl << rowRow << std::endl << std::endl;
//
//    cout << "-------- Test parcours (init par set_coef)" << endl ;
//    OSM::SparseMatrix<int, OSM::COLUMN> Mtest(10, 10) ;
//    Mtest.set_coef(1, 2, 1) ;
//    Mtest.set_coef(3, 2, 1) ;
//    Mtest.set_coef(1, 4, 1) ;
//    Mtest.set_coef(2, 5, 1) ;
//    for (OSM::Bitboard::iterator it  = Mtest.begin(); it != Mtest.end(); ++it)
//    {
//        cout << "col : " << *it << endl ;
//    }
//
    return 0;
}

