#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/OSM/OSM.h>

typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::COLUMN> Column_chain;
typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::ROW> Row_chain ;
typedef CGAL::OSM::Sparse_matrix<int, CGAL::OSM::COLUMN> Column_matrix;
typedef CGAL::OSM::Sparse_matrix<int, CGAL::OSM::ROW> Row_matrix;


int main(int argc, char **argv)
{

    std::cerr << "-- Tests Sparse_matrices operations" << std::endl ;

    std::cerr << "---- Test column chain * row chain -> columnMajor matrix" << std::endl ;
    
    Column_chain a(4);
    Row_chain b(4);
    a.set_coefficient(1,1);
    a.set_coefficient(3,3);
    b.set_coefficient(0,2);
    b.set_coefficient(1,-3);

    Column_matrix columnMajor = a * b;

    CGAL::OSM::write_matrix(columnMajor, "data/test_sparse_matrices/columnMajor.osm");
    
    std::cout << "a * b: " << columnMajor << std::endl ;
    
    std::cerr << "---- Test column chain % row chain -> rowMajor matrix" << std::endl ;
    
    Row_matrix rowMajor = a * b;

    std::cout << "a % b: " << rowMajor << std::endl ;
    
    CGAL::OSM::write_matrix(rowMajor, "data/test_sparse_matrices/rowMajor.osm");

    std::cerr << "---- Test column/row deletion" << std::endl ;
    
    std::cerr << "------ In Column_matrix" << std::endl ;
    
    {
        Column_matrix tmp(columnMajor);
        CGAL::OSM::del_column(tmp, 0) ;
        std::cout << "deletion column 0 in columnMajor: " << tmp << std::endl ;
        
        CGAL::OSM::write_matrix(tmp, "data/test_sparse_matrices/columnMajor_del_col.osm");
    }
    
    {
        Column_matrix tmp(columnMajor);
        CGAL::OSM::del_row(tmp, 1) ;
        std::cout << "deletion row 1 in columnMajor: " << tmp << std::endl ;
        
        CGAL::OSM::write_matrix(tmp, "data/test_sparse_matrices/columnMajor_del_row.osm");
    }
    
    std::cerr << "------ In Row_matrix" << std::endl ;
    
    {
        Row_matrix tmp(rowMajor);
        CGAL::OSM::del_row(tmp, 1) ;
        std::cout << "deletion row 1 in rowMajor: " << tmp << std::endl ;
        
        CGAL::OSM::write_matrix(tmp, "data/test_sparse_matrices/rowMajor_del_row.osm");
    }
    
    {
        Row_matrix tmp(rowMajor);
        CGAL::OSM::del_column(tmp, 0) ;
        std::cout << "deletion column 0 in rowMajor: " << tmp << std::endl ;
        
        CGAL::OSM::write_matrix(tmp, "data/test_sparse_matrices/rowMajor_del_column.osm");
    }
    
    std::cerr << "---- Test matrices sum" << std::endl ;
    
    std::cerr << "------ Test Column_matrix + Column_matrix" << std::endl ;
    
    {
        Column_matrix tmp(columnMajor), res ;
        CGAL::OSM::set_column(tmp, 0, a);
        
//        CGAL::OSM::del_coefficient(tmp, 0, 1);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
    }
    
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
//    cout << "-------- Test parcours (init par set_coefficient)" << endl ;
//    OSM::SparseMatrix<int, OSM::COLUMN> Mtest(10, 10) ;
//    Mtest.set_coefficient(1, 2, 1) ;
//    Mtest.set_coefficient(3, 2, 1) ;
//    Mtest.set_coefficient(1, 4, 1) ;
//    Mtest.set_coefficient(2, 5, 1) ;
//    for (OSM::Bitboard::iterator it  = Mtest.begin(); it != Mtest.end(); ++it)
//    {
//        cout << "col : " << *it << endl ;
//    }
//
    return 0;
}

