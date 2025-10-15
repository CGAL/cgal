#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/HDVF/Zp.h>

namespace HDVF=CGAL::Homological_discrete_vector_field;

typedef HDVF::Zp<5, char, true> Coefficient_field;
typedef HDVF::Zp<12, char, false> Coefficient_ring;

typedef CGAL::Algebraic_structure_traits<Coefficient_field> ATCoefs_field;
typedef CGAL::Algebraic_structure_traits<Coefficient_ring> ATCoefs_ring;

int main() {
    Coefficient_field f1(1), f2(3), f0(0), fp(Coefficient_field::p_val());
    Coefficient_ring r1(1), r2(3), r0(0), rp(Coefficient_field::p_val());
    
    std::cerr << "-- Test Zp is_zero" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true(fp.is_zero()) ;
        std::cerr << "Test zero element: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(f1.is_zero()) ;
        std::cerr << "Test non zero element: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    
    std::cout << "n1==0 : " << ATCoefs_field::Is_zero()(f1) << std::endl;
    
    std::cout << "-1 : " << Coefficient_field(-1) << std::endl ;
    std::cout << "-1 == -1 : " << (Coefficient_field(-1) == -1) << std::endl ;
    std::cout << "4 == -1 : " << (Coefficient_field(4) == -1) << std::endl ;
    
    return 0;
}


