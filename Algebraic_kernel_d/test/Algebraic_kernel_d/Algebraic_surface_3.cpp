//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/config.h>

#include <sstream>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d/macros.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_surface_3.h>

template < class AK >
void test_Algebraic_surface_3() {

    CGAL_SNAP_AK_3_TYPEDEFS(AK);
    
    typedef Integer NT;

    typedef CGAL::Algebraic_surface_3< AK > Algebraic_surface_3;
    
    typename Algebraic_surface_3::Surface_cache cache;

    Poly_int3 t1, t2, t3;

    {
        Integer ra(10);
        Integer rb(5);
        
        Integer ra2 = ra*ra;
        Integer rb2 = rb*rb;
        
        Poly_int1 c0(Integer(0));
        Poly_int1 cp1(Integer(1));
        Poly_int1 cm4a2(Integer(-4)*ra2);
        Poly_int1 ca2(ra2);
        Poly_int1 cb2(rb2);
        
        Poly_int1 x2(Integer(0),Integer(0),Integer(1));
        
        Poly_int2 x2y2(x2,c0,cp1);
        
        Poly_int2 pt1_0 = x2y2 + Poly_int2(ca2 - cb2);
        Poly_int2 pt1_1(c0);
        Poly_int2 pt1_2(cp1);
        Poly_int2 pt2_0 = x2y2 * Poly_int2(cm4a2);
        
        Poly_int3 pt1(pt1_0, pt1_1, pt1_2);
        Poly_int3 pt2(pt2_0);
        
        t1 = pt1*pt1 + pt2;
    }
    
    //::CGAL::set_pretty_mode(std::cout);
    //std::cout << "t1: " << t1 << std::endl;
    
    Algebraic_surface_3 sf1(cache(t1));
    
    assert(sf1.f() == t1);
    
    {
        Integer ra(12);
        Integer rb(5);
        
        Integer ra2 = ra*ra;
        Integer rb2 = rb*rb;
        
        Poly_int1 c0(Integer(0));
        Poly_int1 cp1(Integer(1));
        Poly_int1 cm4a2(Integer(-4)*ra2);
        Poly_int1 ca2(ra2);
        Poly_int1 cb2(rb2);
        
        Poly_int1 x2(Integer(0),Integer(0),Integer(1));
        
        Poly_int2 x2y2(x2,c0,cp1);
        
        Poly_int2 pt1_0 = x2y2 + Poly_int2(ca2 - cb2);
        Poly_int2 pt1_1(c0);
        Poly_int2 pt1_2(cp1);
        Poly_int2 pt2_0 = x2y2 * Poly_int2(cm4a2);
        
        Poly_int3 pt1(pt1_0, pt1_1, pt1_2);
        Poly_int3 pt2(pt2_0);
        
        t2 = pt1*pt1 + pt2;
    }
    
    // cache
    //::CGAL::set_pretty_mode(std::cout);
    //std::cout << "t1: " << t1 << std::endl;
    
    Algebraic_surface_3 sf2(cache(t2));
        
    assert(sf2.f() == t2);


    Algebraic_surface_3 sf3(cache(t1));
    
    assert(sf3.id() == sf1.id());

    // operator<<

    std::stringstream strtest;
    std::stringstream strout;
    strtest << sf3.f() << std::endl;
    strout << sf3 << std::endl;
    assert(strout.str() == strtest.str());
    
    // TODO test for resultant/cofactors etc
}

int main() {
    #if CGAL_USE_LEDA
    {
        typedef CGAL::LEDA_arithmetic_kernel AK;;
        test_Algebraic_surface_3< AK >();
    }
#endif
#if LiS_HAVE_CORE
    {
        typedef CGAL::CORE_arithmetic_kernel AK;
        test_Algebraic_surface_3< AK >();
    }
#endif
    return EXIT_SUCCESS;
}
