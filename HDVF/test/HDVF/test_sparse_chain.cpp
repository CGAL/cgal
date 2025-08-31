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
    Column_chain c1(5), c2(5), c3(5);
    Row_chain r1(5), r2(5), r3(5) ;
    
    c1.set_coefficient(0, 1);
    c1.set_coefficient(2, -1);
    c2.set_coefficient(0, -1);
    c2.set_coefficient(1, 1);
    c3.set_coefficient(2, -1);
    c3.set_coefficient(0, 1);
    
    r1.set_coefficient(0, 1);
    r1.set_coefficient(2, -1);
    r2.set_coefficient(0, -1);
    r2.set_coefficient(1, 1);
    r3.set_coefficient(2, -1);
    r3.set_coefficient(0, 1);
    
    std::cerr << "-- Test Sparse_chain operator==" << std::endl;

    {
        std::cerr << "----> Column chains" << std::endl;
        bool comp_true(c1==c3) ;
        std::cerr << "Compare similar column chains: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(c1==c2) ;
        std::cerr << "Compare different column chains: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        bool comp_true(r1==r3) ;
        std::cerr << "Compare similar row chains: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(r1==r2) ;
        std::cerr << "Compare different row chains: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    
    std::cerr << "-- Test Sparse_chain operator+" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(5), res2(4);
        res.set_coefficient(2, -1);
        res.set_coefficient(1, 1);
        bool comp_true(c1+c2==res) ;
        std::cerr << "Test c1+c2: " << comp_true << std::endl ;
        assert(comp_true) ;
        res2.set_coefficient(2, -1);
        res2.set_coefficient(1, 1);
        bool comp_false(c1+c2==res2) ;
        std::cerr << "Compare c1+c2 with invalid size chain: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(5), res2(4);
        res.set_coefficient(2, -1);
        res.set_coefficient(1, 1);
        bool comp_true(r1+r2==res) ;
        std::cerr << "Test r1+r2: " << comp_true << std::endl ;
        assert(comp_true) ;
        res2.set_coefficient(2, -1);
        res2.set_coefficient(1, 1);
        bool comp_false(r1+r2==res2) ;
        std::cerr << "Compare r1+r2 with invalid size chain: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    
    std::cerr << "-- Test Sparse_chain operator-" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(5);
        res.set_coefficient(0, 2);
        res.set_coefficient(2, -1);
        res.set_coefficient(1, -1);
        bool comp_true(c1-c2==res) ;
        std::cerr << "Test c1-c2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(5);
        res.set_coefficient(0, 2);
        res.set_coefficient(2, -1);
        res.set_coefficient(1, -1);
        bool comp_true(r1-r2==res) ;
        std::cerr << "Test r1-r2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Sparse_chain operator* (coefficient*chain)" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(5);
        res.set_coefficient(0, 3);
        res.set_coefficient(2, -3);
        bool comp_true(3*c1==res) ;
        std::cerr << "Test 3*c1: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(5);
        res.set_coefficient(0, 3);
        res.set_coefficient(2, -3);
        bool comp_true(3*r1==res) ;
        std::cerr << "Test 3*r1: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Sparse_chain operator* (chain*coefficient)" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(5);
        res.set_coefficient(0, 3);
        res.set_coefficient(2, -3);
        bool comp_true(c1*3==res) ;
        std::cerr << "Test c1*3: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(5);
        res.set_coefficient(0, 3);
        res.set_coefficient(2, -3);
        bool comp_true(r1*3==res) ;
        std::cerr << "Test r1*3: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Sparse_chain operator/" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(c1), res2(5);
        res.set_coefficient(1, 3);
        bool comp_true(c1==(res/1)) ;
        std::cerr << "Test c1 == res/1: " << comp_true << std::endl ;
        assert(comp_true) ;
        res2.set_coefficient(0, 1);
        bool comp_true2(res2==(res/std::vector<size_t>({1,2}))) ;
        std::cerr << "Test res2 == res/{1,2}: " << comp_true2 << std::endl ;
        assert(comp_true2) ;
    }
    
    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(r1), res2(5);
        res.set_coefficient(1, 3);
        bool comp_true(r1==(res/1)) ;
        std::cerr << "Test r1 == res/1: " << comp_true << std::endl ;
        assert(comp_true) ;
        res2.set_coefficient(0, 1);
        bool comp_true2(res2==(res/std::vector<size_t>({1,2}))) ;
        std::cerr << "Test res2 == res/{1,2}: " << comp_true2 << std::endl ;
        assert(comp_true2) ;
    }
    
    
//        Row_matrix MR_rw(3,4);
//        CGAL::OSM::set_coefficient(MR_rw, 0, 1, 1) ;
//        CGAL::OSM::set_coefficient(MR_rw, 0, 2, -1) ;
//        CGAL::OSM::set_coefficient(MR_rw, 1, 0, 2) ;
//        CGAL::OSM::set_coefficient(MR_rw, 1, 2, -2) ;
//        
//        write_matrix(MR_rw, "test_row_mat.osm");
//        
//        Row_matrix MR_rw2 ;
//        
//        CGAL::OSM::read_matrix(MR_rw2, "test_row_mat.osm");
//        
//        bool comp_row(MR_rw == MR_rw2) ;
//        std::cerr << "saved and loaded matrices comparison: " << comp_row << std::endl ;
//        assert(comp_row) ;
//    }
//
//    std::cerr << "-- Tests Sparse_matrices operations" << std::endl ;
//
//    Column_chain a(4);
//    Row_chain b(4);
//    a.set_coefficient(1,1);
//    a.set_coefficient(3,3);
//    b.set_coefficient(0,2);
//    b.set_coefficient(1,-3);
//
//    std::cerr << "---- Test column chain * row chain -> columnMajor matrix" << std::endl ;
//    
//    Column_matrix columnMajor = a * b, tmp = columnMajor;
//    
//    {
//        Column_matrix tmp ;
//        CGAL::OSM::read_matrix(tmp, "data/test_sparse_matrices/columnMajor.osm") ;
//        bool bres(tmp == columnMajor);
//        std::cerr << "column chain * row chain: " << bres << std::endl;
//        assert(bres) ;
//    }
//    
//    std::cerr << "---- Test column chain % row chain -> rowMajor matrix" << std::endl ;
//    
//    Row_matrix rowMajor = a % b;
//
//    {
//        Row_matrix tmp ;
//        CGAL::OSM::read_matrix(tmp, "data/test_sparse_matrices/rowMajor.osm") ;
//        bool bres(tmp == rowMajor);
//        std::cerr << "column chain % row chain: " << bres << std::endl;
//        assert(bres) ;
//    }
//    
//    std::cerr << "---- Test column/row get" << std::endl ;
//    
//    std::cerr << "------ In Column_matrix" << std::endl ;
//    
//    {
//        Column_chain res(CGAL::OSM::get_column(columnMajor, 0)), tmp(4) ;
//        tmp.set_coefficient(1, 2) ;
//        tmp.set_coefficient(3, 6) ;
//        
//        bool bres (tmp == res) ;
//        std::cerr << "get_column 0 in columnMajor: " << bres << std::endl;
//        assert(bres) ;
//    }
//    
//    {
//        Row_chain res2(CGAL::OSM::get_row(columnMajor, 1)), tmp2(4) ;
//        tmp2.set_coefficient(0, 2) ;
//        tmp2.set_coefficient(1, -3) ;
//        
//        bool bres (tmp2 == res2) ;
//        std::cerr << "get_row 1 in columnMajor: " << bres << std::endl;
//        assert(bres);
//    }
//    
//    std::cerr << "------ In Row_matrix" << std::endl ;
//    
//    {
//        Column_chain res(CGAL::OSM::get_column(rowMajor, 0)), tmp(4) ;
//        tmp.set_coefficient(1, 2) ;
//        tmp.set_coefficient(3, 6) ;
//        
//        bool bres(tmp == res) ;
//        std::cerr << "get_column 0 in rowMajor: " << bres << std::endl;
//        assert(bres);
//    }
//    
//    {
//        Row_chain res2(CGAL::OSM::get_row(rowMajor, 1)), tmp2(4) ;
//        tmp2.set_coefficient(0, 2) ;
//        tmp2.set_coefficient(1, -3) ;
//        
//        bool bres(tmp2==res2);
//        std::cerr << "get_row 1 in columnMajor: " << bres << std::endl;
//        assert(bres);
//    }
//    
//    std::cerr << "---- Test column/row deletion" << std::endl;
//    
//    std::cerr << "------ In Column_matrix" << std::endl ;
//    
//    {
//        Column_matrix tmp(columnMajor), tmp2;
//        CGAL::OSM::del_column(tmp, 0) ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/columnMajor_del_col.osm");
//        
//        bool bres(tmp == tmp2);
//        std::cerr << "deletion column 0 in columnMajor: " << bres << std::endl ;
//        assert(bres);
//    }
//    
//    {
//        Column_matrix tmp(columnMajor), tmp2;
//        CGAL::OSM::del_row(tmp, 1) ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/columnMajor_del_row.osm");
//        bool bres(tmp == tmp2);
//        std::cerr << "deletion row 1 in columnMajor: " << bres << std::endl ;
//        assert(bres);
//    }
//    
//    std::cerr << "------ In Row_matrix" << std::endl ;
//    
//    {
//        Row_matrix tmp(rowMajor), tmp2;
//        CGAL::OSM::del_row(tmp, 1) ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/rowMajor_del_row.osm");
//        bool bres(tmp==tmp2);
//        std::cerr << "deletion row 1 in rowMajor: " << bres << std::endl ;
//        assert(bres);
//    }
//    
//    {
//        Row_matrix tmp(rowMajor), tmp2;
//        CGAL::OSM::del_column(tmp, 0) ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/rowMajor_del_column.osm");
//        bool bres(tmp == tmp2);
//        std::cerr << "deletion column 0 in rowMajor: " << bres << std::endl ;
//        assert(bres);
//    }
//    
//    // //////// Sum
//    std::cerr << "---- Test matrices sum" << std::endl ;
//    
//    // COL + COL
//    std::cerr << "------ Test Column_matrix + Column_matrix" << std::endl ;
//    
//    {
//        Column_matrix tmp(columnMajor), res, tmp2 ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp + columnMajor ;
//        Column_matrix res2(tmp);
//        res2 += columnMajor ;
//        
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/ColCol_sum.osm");
//        bool bres(res == tmp2) ;
//        std::cerr << "sum columnMajor + otherColumnMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(res2 == tmp2) ;
//        std::cerr << "sum columnMajor += otherColumnMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    // COL + ROW
//    std::cerr << "------ Test Column_matrix + Row_matrix" << std::endl ;
//    
//    {
//        Column_matrix tmp(columnMajor), res, tmp2 ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp + rowMajor ;
//        Column_matrix res2(tmp);
//        res2 += rowMajor ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/ColRow_sum.osm");
//        bool bres(tmp2 == res);
//        std::cerr << "sum otherColumnMajor + rowMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(tmp2 == res2);
//        std::cerr << "sum otherColumnMajor += rowMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    // ROW + COL
//    std::cerr << "------ Test Row_matrix + Column_matrix" << std::endl ;
//    
//    {
//        Row_matrix tmp(rowMajor), res, tmp2 ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp + columnMajor ;
//        Row_matrix res2(tmp);
//        res2 += columnMajor ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/RowCol_sum.osm");
//        bool bres(tmp2 == res);
//        std::cerr << "sum otherRowMajor + columnMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(tmp2 == res2);
//        std::cerr << "sum otherRowMajor += columnMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    
//    // ROW + ROW
//    std::cerr << "------ Test Row_matrix + Row_matrix" << std::endl ;
//    
//    {
//        Row_matrix tmp(rowMajor), res, tmp2 ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp + rowMajor ;
//        Row_matrix res2(tmp);
//        res2 += rowMajor;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/RowRow_sum.osm");
//        bool bres(tmp2 == res);
//        std::cerr << "sum otherRowMajor + rowMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(tmp2 == res2);
//        std::cerr << "sum otherRowMajor += rowMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    // //////// Subtract
//    std::cerr << "---- Test matrices subtract" << std::endl ;
//    
//    // COL - COL
//    std::cerr << "------ Test Column_matrix - Column_matrix" << std::endl ;
//    
//    {
//        Column_matrix tmp(columnMajor), res, tmp2 ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp - columnMajor ;
//        Column_matrix res2(tmp);
//        res2 -= columnMajor;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/ColCol_sub.osm");
//        bool bres(res == tmp2) ;
//        std::cerr << "sub columnMajor - otherColumnMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(res2 == tmp2) ;
//        std::cerr << "sub columnMajor -= otherColumnMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    
//    // COL - ROW
//    std::cerr << "------ Test Column_matrix - Row_matrix" << std::endl ;
//    
//    {
//        Column_matrix tmp(columnMajor), res, tmp2 ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp - rowMajor ;
//        Column_matrix res2(tmp);
//        res2 -= rowMajor;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/ColRow_sub.osm");
//        bool bres(res == tmp2) ;
//        std::cerr << "sub otherColumnMajor - rowMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(res2 == tmp2) ;
//        std::cerr << "sub otherColumnMajor -= rowMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    // ROW - COL
//    std::cerr << "------ Test Row_matrix - Column_matrix" << std::endl ;
//    
//    {
//        Row_matrix tmp(rowMajor), res, tmp2 ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp - columnMajor ;
//        Row_matrix res2(tmp);
//        res2 -= columnMajor;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/RowCol_sub.osm");
//        bool bres(res == tmp2);
//        std::cerr << "sub otherRowMajor - columnMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(res2 == tmp2);
//        std::cerr << "sub otherRowMajor -= columnMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    // ROW - ROW
//    std::cerr << "------ Test Row_matrix - Row_matrix" << std::endl ;
//    
//    {
//        Row_matrix tmp(rowMajor), res, tmp2 ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp - rowMajor ;
//        Row_matrix res2(tmp);
//        res2 -= rowMajor;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/RowRow_sub.osm");
//        bool bres(res == tmp2);
//        std::cerr << "sub otherRowMajor - rowMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(res2 == tmp2);
//        std::cerr << "sub otherRowMajor -= rowMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    // //////// Product * (column result)
//    std::cerr << "---- Test matrices product *" << std::endl ;
//    
//    // COL * COL
//    std::cerr << "------ Test Column_matrix * Column_matrix" << std::endl ;
//    
//    {
//        Column_matrix res, tmp2 ;
//        Column_matrix tmp(columnMajor) ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp * columnMajor ;
//        Column_matrix res2(tmp) ;
//        res2 *= columnMajor;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/ColCol_prod_col.osm");
//        bool bres(res == tmp2);
//        std::cerr << "prod columnMajor * otherColumnMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(res2 == tmp2);
//        std::cerr << "prod columnMajor *= otherColumnMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    // COL * ROW
//    std::cerr << "------ Test Column_matrix * Row_matrix" << std::endl ;
//    
//    {
//        Column_matrix res, tmp2 ;
//        Column_matrix tmp(columnMajor) ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp * rowMajor ;
//        Column_matrix res2(tmp);
//        res2 *= rowMajor;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/ColRow_prod_col.osm");
//        bool bres(res == tmp2);
//        std::cerr << "prod otherColumnMajor * rowMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(res2 == tmp2);
//        std::cerr << "prod otherColumnMajor *= rowMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    // ROW * COL
//    std::cerr << "------ Test Row_matrix * Column_matrix" << std::endl ;
//    
//    {
//        Column_matrix res, tmp2 ;
//        Row_matrix tmp(rowMajor) ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp * columnMajor ;
//        Column_matrix res2(tmp);
//        res2 *= columnMajor;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/RowCol_prod_col.osm");
//        bool bres(res == tmp2);
//        std::cerr << "prod otherRowMajor * columnMajor: " << bres << std::endl ;
//        assert(bres);
//        bool bres2(res2 == tmp2);
//        std::cerr << "prod otherRowMajor *= columnMajor: " << bres2 << std::endl ;
//        assert(bres2);
//    }
//    
//    
//    // ROW * ROW
//    std::cerr << "------ Test Row_matrix * Row_matrix" << std::endl ;
//    
//    {
//        Column_matrix res, tmp2 ;
//        Row_matrix tmp(rowMajor) ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp * rowMajor ;
//        Column_matrix res2(tmp) ;
//        res2 *= rowMajor;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/RowRow_prod_col.osm");
//        bool bres(res == tmp2);
//        std::cerr << "prod otherRowMajor * rowMajor: " << bres << std::endl;
//        assert(bres);
//        bool bres2(res2 == tmp2);
//        std::cerr << "prod otherRowMajor *= rowMajor: " << bres2 << std::endl;
//        assert(bres2);
//    }
//    
//    // //////// Product % (row result)
//    std::cerr << "---- Test matrices product %" << std::endl ;
//    
//    // COL % COL
//    std::cerr << "------ Test Column_matrix % Column_matrix" << std::endl ;
//    
//    {
//        Row_matrix res, tmp2 ;
//        Column_matrix tmp(columnMajor) ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp % columnMajor ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/ColCol_prod_row.osm");
//        bool bres(res == tmp2);
//        std::cerr << "prod columnMajor % otherColumnMajor: " << bres << std::endl ;
//        assert(bres);
//    }
//    
//    
//    // COL % ROW
//    std::cerr << "------ Test Column_matrix % Row_matrix" << std::endl ;
//    
//    {
//        Row_matrix res, tmp2 ;
//        Column_matrix tmp(columnMajor) ;
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp % rowMajor ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/ColRow_prod_row.osm");
//        bool bres(res == tmp2);
//        std::cerr << "prod otherColumnMajor % rowMajor: " << bres << std::endl ;
//        assert(bres);
//    }
//    
//    // ROW % COL
//    std::cerr << "------ Test Row_matrix % Column_matrix" << std::endl ;
//    
//    {
//        Row_matrix res, tmp2;
//        Row_matrix tmp(rowMajor);
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp % columnMajor ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/RowCol_prod_row.osm");
//        bool bres(res == tmp2);
//        std::cerr << "prod otherRowMajor % columnMajor: " << bres << std::endl ;
//        assert(bres);
//    }
//    
//    
//    // ROW % ROW
//    std::cerr << "------ Test Row_matrix % Row_matrix" << std::endl ;
//    
//    {
//        Row_matrix res, tmp2;
//        Row_matrix tmp(rowMajor);
//        
//        CGAL::OSM::del_coefficient(tmp, 1, 0);
//        CGAL::OSM::set_coefficient(tmp, 2, 2, -2) ;
//        
//        res = tmp % rowMajor ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/RowRow_prod_row.osm");
//        bool bres(res == tmp2);
//        std::cerr << "prod otherRowMajor % rowMajor: " << bres << std::endl ;
//        assert(bres);
//        
//    }
//    
//    // //////// Transpose
//    std::cerr << "---- Test matrices transpose" << std::endl ;
//    
//    // COL^t
//    std::cerr << "------ Test Column_matrix transpose" << std::endl ;
//    
//    {
//        Row_matrix res = columnMajor.transpose(), tmp2 ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/Col_transpose.osm");
//        bool bres(res == tmp2);
//        std::cerr << "transpose columnMajor: " << bres << std::endl ;
//        assert(bres);
//    }
//    
//    // ROW^t
//    std::cerr << "------ Test Row_matrix transpose" << std::endl ;
//    
//    {
//        Column_matrix res = rowMajor.transpose(), tmp2 ;
//        CGAL::OSM::read_matrix(tmp2, "data/test_sparse_matrices/Row_transpose.osm");
//        bool bres(res == tmp2);
//        std::cerr << "transpose rowMajor: " << bres << std::endl ;
//        assert(bres);
//    }
//    
//    // //////// Product scalar * matrix
//    std::cerr << "---- Test scalar - matrix and matrix - scalar product" << std::endl ;
//    
//    // COL
//    std::cerr << "------ Test 3 * Column_matrix and Column_matrix * 3" << std::endl ;
//    
//    {
//        Row_matrix tmp(rowMajor), tmp2(rowMajor), res1, res2, res3 ;
//        
//        // Build expected result
//        for (CGAL::OSM::Bitboard::iterator it = tmp2.begin(); it != tmp2.end(); ++it)
//        {
//            const Row_chain& tmp_row(CGAL::OSM::cget_row(tmp2, *it));
//            for (Row_chain::const_iterator it2=tmp_row.cbegin(); it2 != tmp_row.cend(); ++it2)
//            {
//                CGAL::OSM::set_coefficient(tmp2, *it, it2->first, 3*it2->second);
//            }
//        }
//        
//        res1 = 3*tmp ;
//        res2 = tmp*3 ;
//        res3 = tmp ;
//        res3 *= 3 ;
//        bool bres1(res1 == tmp2);
//        std::cerr << "prod 3*rowMajor: " << bres1 << std::endl ;
//        assert(bres1);
//        bool bres2(res2 == tmp2);
//        std::cerr << "prod rowMajor*3: " << bres2 << std::endl ;
//        assert(bres2);
//        bool bres3(res3 == tmp2);
//        std::cerr << "prod rowMajor *= 3: " << bres3 << std::endl ;
//        assert(bres3);
//    }
//    
//    // ROW
//    std::cerr << "------ Test 3 * Row_matrix and Row_matrix * 3" << std::endl ;
//    
//    {
//        Column_matrix tmp(columnMajor), tmp2(columnMajor), res1, res2, res3 ;
//        
//        // Build expected result
//        for (CGAL::OSM::Bitboard::iterator it = tmp2.begin(); it != tmp2.end(); ++it)
//        {
//            const Column_chain& tmp_col(CGAL::OSM::cget_column(tmp2, *it));
//            for (Column_chain::const_iterator it2=tmp_col.cbegin(); it2 != tmp_col.cend(); ++it2)
//            {
//                CGAL::OSM::set_coefficient(tmp2, it2->first, *it, 3*it2->second);
//            }
//        }
//        
//        res1 = 3*tmp ;
//        res2 = tmp*3 ;
//        res3 = tmp ;
//        res3 *= 3 ;
//        bool bres1(res1 == tmp2);
//        std::cerr << "prod 3*columnMajor: " << bres1 << std::endl ;
//        assert(bres1);
//        bool bres2(res2 == tmp2);
//        std::cerr << "prod columnMajor*3: " << bres2 << std::endl ;
//        assert(bres2);
//        bool bres3(res3 == tmp2);
//        std::cerr << "prod columnMajor *= 3: " << bres3 << std::endl ;
//        assert(bres3);
//    }
//
////    std::cout << "Multiplication de colonne (col/row dom):" << std::endl << colCol << std::endl << colRow << std::endl << std::endl;
////    std::cout << "Multiplication de ligne   (col/row dom):" << std::endl << rowCol << std::endl << rowRow << std::endl << std::endl;
////
////    cout << "-------- Test parcours (init par set_coefficient)" << endl ;
////    OSM::SparseMatrix<int, OSM::COLUMN> Mtest(10, 10) ;
////    Mtest.set_coefficient(1, 2, 1) ;
////    Mtest.set_coefficient(3, 2, 1) ;
////    Mtest.set_coefficient(1, 4, 1) ;
////    Mtest.set_coefficient(2, 5, 1) ;
////    for (OSM::Bitboard::iterator it  = Mtest.begin(); it != Mtest.end(); ++it)
////    {
////        cout << "col : " << *it << endl ;
////    }
////
//    return 0;
}

