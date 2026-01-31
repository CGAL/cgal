#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/OSM/OSM.h>

// TODO: use also results to check += -= *= .....

typedef int Coefficient_ring;
typedef CGAL::OSM::Sparse_chain<Coefficient_ring, CGAL::OSM::COLUMN> Column_chain;
typedef CGAL::OSM::Sparse_chain<Coefficient_ring, CGAL::OSM::ROW> Row_chain ;
typedef CGAL::OSM::Sparse_matrix<Coefficient_ring, CGAL::OSM::COLUMN> Column_matrix;
typedef CGAL::OSM::Sparse_matrix<Coefficient_ring, CGAL::OSM::ROW> Row_matrix;


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

    std::cerr << "---- Test column chain * row chain -> columnMajor matrix" << std::endl ;

    Column_matrix columnMajor = a * b;

    CGAL::OSM::write_matrix(columnMajor, "data/test_sparse_matrices/columnMajor.osm");

    std::cout << "a * b: " << columnMajor << std::endl ;

    std::cerr << "---- Test column chain % row chain -> rowMajor matrix" << std::endl ;

    Row_matrix rowMajor = a % b;

    std::cout << "a % b: " << rowMajor << std::endl ;

    CGAL::OSM::write_matrix(rowMajor, "data/test_sparse_matrices/rowMajor.osm");

    std::cerr << "---- Test column/row deletion" << std::endl;

    std::cerr << "------ In Column_matrix" << std::endl ;

    {
        Column_matrix tmp(columnMajor);
        CGAL::OSM::remove_column(tmp, 0) ;
        std::cout << "deletion column 0 in columnMajor: " << tmp << std::endl ;

        CGAL::OSM::write_matrix(tmp, "data/test_sparse_matrices/columnMajor_del_col.osm");
    }

    {
        Column_matrix tmp(columnMajor);
        CGAL::OSM::remove_row(tmp, 1) ;
        std::cout << "deletion row 1 in columnMajor: " << tmp << std::endl ;

        CGAL::OSM::write_matrix(tmp, "data/test_sparse_matrices/columnMajor_del_row.osm");
    }

    std::cerr << "------ In Row_matrix" << std::endl ;

    {
        Row_matrix tmp(rowMajor);
        CGAL::OSM::remove_row(tmp, 1) ;
        std::cout << "deletion row 1 in rowMajor: " << tmp << std::endl ;

        CGAL::OSM::write_matrix(tmp, "data/test_sparse_matrices/rowMajor_del_row.osm");
    }

    {
        Row_matrix tmp(rowMajor);
        CGAL::OSM::remove_column(tmp, 0) ;
        std::cout << "deletion column 0 in rowMajor: " << tmp << std::endl ;

        CGAL::OSM::write_matrix(tmp, "data/test_sparse_matrices/rowMajor_del_column.osm");
    }

    // //////// Sum
    std::cerr << "---- Test matrices sum" << std::endl ;

    // COL + COL
    std::cerr << "------ Test Column_matrix + Column_matrix" << std::endl ;

    {
        Column_matrix tmp(columnMajor), res ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp + columnMajor ;
        std::cout << "sum columnMajor + otherColumnMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/ColCol_sum.osm");
    }


    // COL + ROW
    std::cerr << "------ Test Column_matrix + Row_matrix" << std::endl ;

    {
        Column_matrix tmp(columnMajor), res ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp + rowMajor ;
        std::cout << "sum otherColumnMajor + rowMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/ColRow_sum.osm");
    }

    // ROW + COL
    std::cerr << "------ Test Row_matrix + Column_matrix" << std::endl ;

    {
        Row_matrix tmp(rowMajor), res ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp + columnMajor ;
        std::cout << "sum otherRowMajor + columnMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/RowCol_sum.osm");
    }


    // ROW + ROW
    std::cerr << "------ Test Row_matrix + Row_matrix" << std::endl ;

    {
        Row_matrix tmp(rowMajor), res ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp + rowMajor ;
        std::cout << "sum otherRowMajor + rowMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/RowRow_sum.osm");
    }

    // //////// Subtract
    std::cerr << "---- Test matrices subtract" << std::endl ;

    // COL - COL
    std::cerr << "------ Test Column_matrix - Column_matrix" << std::endl ;

    {
        Column_matrix tmp(columnMajor), res ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp - columnMajor ;
        std::cout << "sub columnMajor - otherColumnMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/ColCol_sub.osm");
    }


    // COL - ROW
    std::cerr << "------ Test Column_matrix - Row_matrix" << std::endl ;

    {
        Column_matrix tmp(columnMajor), res ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp - rowMajor ;
        std::cout << "sub otherColumnMajor - rowMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/ColRow_sub.osm");
    }

    // ROW - COL
    std::cerr << "------ Test Row_matrix - Column_matrix" << std::endl ;

    {
        Row_matrix tmp(rowMajor), res ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp - columnMajor ;
        std::cout << "sub otherRowMajor - columnMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/RowCol_sub.osm");
    }


    // ROW - ROW
    std::cerr << "------ Test Row_matrix - Row_matrix" << std::endl ;

    {
        Row_matrix tmp(rowMajor), res ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp - rowMajor ;
        std::cout << "sub otherRowMajor - rowMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/RowRow_sub.osm");
    }

    // //////// Product * (column result)
    std::cerr << "---- Test matrices product *" << std::endl ;

    // COL * COL
    std::cerr << "------ Test Column_matrix * Column_matrix" << std::endl ;

    {
        Column_matrix res ;
        Column_matrix tmp(columnMajor) ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp * columnMajor ;
        std::cout << "prod columnMajor * otherColumnMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/ColCol_prod_col.osm");
    }


    // COL * ROW
    std::cerr << "------ Test Column_matrix * Row_matrix" << std::endl ;

    {
        Column_matrix res ;
        Column_matrix tmp(columnMajor) ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp * rowMajor ;
        std::cout << "prod otherColumnMajor * rowMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/ColRow_prod_col.osm");
    }

    // ROW * COL
    std::cerr << "------ Test Row_matrix * Column_matrix" << std::endl ;

    {
        Column_matrix res ;
        Row_matrix tmp(rowMajor) ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp * columnMajor ;
        std::cout << "prod otherRowMajor * columnMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/RowCol_prod_col.osm");
    }


    // ROW * ROW
    std::cerr << "------ Test Row_matrix * Row_matrix" << std::endl ;

    {
        Column_matrix res ;
        Row_matrix tmp(rowMajor) ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp * rowMajor ;
        std::cout << "prod otherRowMajor * rowMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/RowRow_prod_col.osm");
    }

    // //////// Product % (row result)
    std::cerr << "---- Test matrices product %" << std::endl ;

    // COL % COL
    std::cerr << "------ Test Column_matrix % Column_matrix" << std::endl ;

    {
        Row_matrix res ;
        Column_matrix tmp(columnMajor) ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp % columnMajor ;
        std::cout << "prod columnMajor % otherColumnMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/ColCol_prod_row.osm");
    }


    // COL % ROW
    std::cerr << "------ Test Column_matrix % Row_matrix" << std::endl ;

    {
        Row_matrix res ;
        Column_matrix tmp(columnMajor) ;

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp % rowMajor ;
        std::cout << "prod otherColumnMajor % rowMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/ColRow_prod_row.osm");
    }

    // ROW % COL
    std::cerr << "------ Test Row_matrix % Column_matrix" << std::endl ;

    {
        Row_matrix res;
        Row_matrix tmp(rowMajor);

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp % columnMajor ;
        std::cout << "prod otherRowMajor % columnMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/RowCol_prod_row.osm");
    }


    // ROW % ROW
    std::cerr << "------ Test Row_matrix % Row_matrix" << std::endl ;

    {
        Row_matrix res;
        Row_matrix tmp(rowMajor);

        CGAL::OSM::remove_coefficient(tmp, 1, 0);
        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;

        res = tmp % rowMajor ;
        std::cout << "prod otherRowMajor % rowMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/RowRow_prod_row.osm");
    }

    // //////// Transpose
    std::cerr << "---- Test matrices transpose" << std::endl ;

    // COL^t
    std::cerr << "------ Test Column_matrix transpose" << std::endl ;

    {
        Row_matrix res = columnMajor.transpose() ;

        std::cout << "transpose columnMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/Col_transpose.osm");
    }

    // ROW^t
    std::cerr << "------ Test Row_matrix transpose" << std::endl ;

    {
        Column_matrix res = rowMajor.transpose() ;

        std::cout << "transpose rowMajor: " << res << std::endl ;

        CGAL::OSM::write_matrix(res, "data/test_sparse_matrices/Row_transpose.osm");
    }

    return 0;
}

