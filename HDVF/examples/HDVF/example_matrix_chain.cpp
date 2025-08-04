#include <iostream>
#include <fstream>

#include <CGAL/OSM/OSM.h>

typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::COLUMN> CChain;
typedef CGAL::OSM::Sparse_matrix<int, CGAL::OSM::COLUMN> CMatrix;

int main ()
{
    // Create a column-major sparse matrix
    CMatrix M(5,4) ;

    // Fill coefficients
    CGAL::OSM::set_coef(M, 0, 1, 1) ;
    CGAL::OSM::set_coef(M, 0, 2, -1) ;
    CGAL::OSM::set_coef(M, 2, 1, 2) ;

    // Iterate over non empty columns
    for(CGAL::OSM::Bitboard::iterator it_col = M.begin(); it_col != M.end(); ++it_col)
    {
        std::cout << "col: " << *it_col << std::endl ;
        // Get a constant reference over the column (complexity O(1))
        const CChain& col(CGAL::OSM::cget_column(M, *it_col));
        // Iterate over the column
        for (CChain::const_iterator it = col.begin(); it != col.end(); ++it)
        {
            std::cout << "row: " << it->first << " - coef: " << it->second << std::endl ;
        }
    }
    // Direct output of the matrix with << operator
    std::cout << "M: " << M << std::endl;

    return 0 ;
}
