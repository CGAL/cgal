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
    Coefficient_field f1(1), f2(3), f0(0), fp((Coefficient_field()));
    Coefficient_ring r1(1), r2(3), r0(0), rp((Coefficient_ring()));
    
    std::cerr << "-- Test Zp ==" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true(Coefficient_field(2)==Coefficient_field(7)) ;
        std::cerr << "Test 2==7: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(Coefficient_field(2)==Coefficient_field(-2)) ;
        std::cerr << "Test 2==-2: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true(Coefficient_ring(2) == Coefficient_ring(14)) ;
        std::cerr << "Test 2==14: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(Coefficient_ring(2) == Coefficient_ring(-2)) ;
        std::cerr << "Test 2==-2: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    
    std::cerr << "-- Test Zp !=" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true(Coefficient_field(2)!=Coefficient_field(-2)) ;
        std::cerr << "Test 2 != -2: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(Coefficient_field(2)!=Coefficient_field(7)) ;
        std::cerr << "Test 2 != 7: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true(Coefficient_ring(2) != Coefficient_ring(-2)) ;
        std::cerr << "Test 2!=-2: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(Coefficient_ring(2) != Coefficient_ring(14)) ;
        std::cerr << "Test 2!=14: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    
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
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true(rp.is_zero()) ;
        std::cerr << "Test zero element: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(r1.is_zero()) ;
        std::cerr << "Test non zero element: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    
    std::cerr << "-- Test Zp unary operator+" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true((+Coefficient_field(2)) == Coefficient_field(2)) ;
        std::cerr << "Test +2 == 2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true((+Coefficient_ring(6)) == Coefficient_ring(6)) ;
        std::cerr << "Test +6 == 6: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp unary operator-" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true((-Coefficient_field(2)) == Coefficient_field(3)) ;
        std::cerr << "Test -2 == 3: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true((-Coefficient_ring(5)) == Coefficient_ring(7)) ;
        std::cerr << "Test -5 == 7: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp operator+" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true((Coefficient_field(2)+Coefficient_field(4)) == Coefficient_field(1)) ;
        std::cerr << "Test 2+4 == 1: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true((Coefficient_ring(6)+Coefficient_ring(8)) == Coefficient_ring(2)) ;
        std::cerr << "Test 6+8 == 2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp operator-" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true((Coefficient_field(2)-Coefficient_field(4)) == Coefficient_field(-2)) ;
        std::cerr << "Test 2-4 == -2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true((Coefficient_ring(6)-Coefficient_ring(8)) == Coefficient_ring(-2)) ;
        std::cerr << "Test 6-8 == -2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp operator*" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true((Coefficient_field(2)*Coefficient_field(4)) == Coefficient_field(3)) ;
        std::cerr << "Test 2*4 == 3: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true((Coefficient_ring(6)*Coefficient_ring(7)) == Coefficient_ring(6)) ;
        std::cerr << "Test 6*7 == 6: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp operator/" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true((Coefficient_field(4)/Coefficient_field(3)) == Coefficient_field(1)) ;
        std::cerr << "Test 4/3 == 1: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true((Coefficient_ring(8)/Coefficient_ring(3)) == Coefficient_ring(2)) ;
        std::cerr << "Test 8/3 == 2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp operator+=" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        Coefficient_field tmp(2);
        tmp += Coefficient_field(4);
        bool comp_true(tmp == Coefficient_field(1)) ;
        std::cerr << "Test 2+=4 == 1: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        Coefficient_ring tmp(6);
        tmp += Coefficient_ring(8);
        bool comp_true(tmp == Coefficient_ring(2)) ;
        std::cerr << "Test 6+=8 == 2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp operator-=" << std::endl;
    {
        Coefficient_field tmp(2);
        tmp -= Coefficient_field(4);
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true(tmp == Coefficient_field(-2)) ;
        std::cerr << "Test 2-=4 == -2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        Coefficient_ring tmp(6);
        tmp -= Coefficient_ring(8);
        bool comp_true(tmp == Coefficient_ring(-2)) ;
        std::cerr << "Test 6-=8 == -2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp operator*=" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        Coefficient_field tmp(2);
        tmp *= Coefficient_field(4);
        bool comp_true(tmp == Coefficient_field(3)) ;
        std::cerr << "Test 2*=4 == 3: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        Coefficient_ring tmp(6);
        tmp *= Coefficient_ring(7);
        bool comp_true(tmp == Coefficient_ring(6)) ;
        std::cerr << "Test 6*=7 == 6: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp operator/=" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        Coefficient_field tmp(4);
        tmp /= Coefficient_field(3);
        bool comp_true(tmp == Coefficient_field(1)) ;
        std::cerr << "Test 4/=3 == 1: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        Coefficient_ring tmp(8);
        tmp /= Coefficient_ring(3);
        bool comp_true(tmp == Coefficient_ring(2)) ;
        std::cerr << "Test 8/=3 == 2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    
    std::cerr << "-- Test Zp is_invertible" << std::endl;
    {
        std::cerr << "----> Field Z5" << std::endl;
        bool comp_true(Coefficient_field(2).is_invertible()) ;
        std::cerr << "Test is_invertible(2): " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(Coefficient_field(0).is_invertible()) ;
        std::cerr << "Test is_invertible(0): " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    {
        std::cerr << "----> Ring Z12" << std::endl;
        bool comp_true(Coefficient_ring(5).is_invertible()) ;
        std::cerr << "Test is_invertible(5): " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(Coefficient_ring(2).is_invertible()) ;
        std::cerr << "Test is_invertible(2): " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    
    std::cout << "n1==0 : " << ATCoefs_field::Is_zero()(f1) << std::endl;
    
    return 0;
}


