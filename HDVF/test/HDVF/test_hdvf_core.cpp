#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/HDVF/Simplex.h>

namespace HDVF = CGAL::Homological_discrete_vector_field

int main() {
    HDVF::Simplex s1({1,2}), s2({1,2,3});
    std::vector<HDVF::Simplex> bnd, bnd2 ;

    // Check s1
    std::cerr << "Check dimension s1" << std::endl ;
    assert(s1.dimension() == 1);
    std::cerr << "Check boundary s1" << std::endl ;
    bnd = s1.boundary() ;
    bnd2.push_back(HDVF::Simplex({2})) ;
    bnd2.push_back(HDVF::Simplex({1})) ;
    assert(bnd.size() == bnd2.size());
    for (int i=0; i<bnd.size(); ++i)
        assert(bnd.at(i) == bnd2.at(i)) ;

    // Check s2
    std::cerr << "Check dimension s2" << std::endl ;
    assert(s2.dimension() == 2);
    std::cerr << "Check boundary s2" << std::endl ;
    bnd = s2.boundary() ;
    bnd2.clear();
    bnd2.push_back(HDVF::Simplex({2,3})) ;
    bnd2.push_back(HDVF::Simplex({1,3})) ;
    bnd2.push_back(HDVF::Simplex({1,2})) ;
    assert(bnd.size() == bnd2.size());
    for (int i=0; i<bnd.size(); ++i)
        assert(bnd.at(i) == bnd2.at(i)) ;

    // Check Simplex() sort option
    std::cerr << "Check Simplex() sort option" << std::endl;
    HDVF::Simplex s3({3,2,1},true);
    assert(s2 == s3);

    // Check s1 < s2
    std::cerr << "Check s1 < s2" << std::endl;
    assert(s1 < s2);

    // Test iterator
    int i = 1 ;
    std::cerr << "Check iterator" << std::endl;
    for (HDVF::Simplex::const_iterator it = s2.cbegin(); it != s2.cend(); ++it)
        assert(*it == i++);

    return 0;
}


