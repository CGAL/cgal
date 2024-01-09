#include <iostream>
#include <CGAL/config.h>
#include <CGAL/Interval_nt.h>
#include <cassert>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Uncertain.h>
#include <limits>

#define CGAL_catch_error(expr, error)                            \
    {                                                            \
        bool b = false;                                          \
        try{(void) expr;}catch(error){ b = true;}                \
        if(!b) CGAL_error_msg( "Expr should throw exception");   \
    }

int main() {
{
    typedef CGAL::Interval_nt<true> NT;
    typedef CGAL::Field_with_sqrt_tag Tag;
    typedef CGAL::Tag_false Is_exact;

    CGAL::test_algebraic_structure<NT,Tag, Is_exact>();

    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(-15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6), NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(-15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(-15));

    CGAL::test_real_embeddable<NT>();
}
{
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);
    typedef CGAL::Interval_nt<false> NT;
    typedef CGAL::Field_with_sqrt_tag Tag;
    typedef CGAL::Tag_false Is_exact;
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(-15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6), NT(15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(-15));
    CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(-15));

    CGAL::test_real_embeddable<NT>();
}

    typedef CGAL::Interval_nt<true> Interval;

    { // DEFAULT CONSTRUCTOR I(0)
        Interval I; // TODO: Standard constructor does not initialize I with 0.
/*        assert(I==.0);
        assert(I.inf()==.0);
        assert(I.sup()==.0);*/
    }

    { // CONSTRUCTOR FROM INT
        Interval I(1);
        assert(I.inf()==1.0);
        assert(I.sup()==1.0);
        I = Interval(1,2);
        assert(I.inf()==1.0);
        assert(I.sup()==2.0);
    }
    { // CONSTRUCTOR FROM double
        Interval I(1.0);
        assert(I.inf()==1.0);
        assert(I.sup()==1.0);
        I= Interval(1.0,2.0);
        assert(I.inf()==1.0);
        assert(I.sup()==2.0);
    }
    { // assign
      // TODO: No assign available in Interval_nt
/*        Interval I;
        I.assign(2.0,3.0);
        assert(I.inf()==2.0);
        assert(I.sup()==3.0);
        I.assign(-2.0,-1.0);
        assert(I.inf()==-2.0);
        assert(I.sup()==-1.0); */
    }

    { //comparison
        Interval I,J;

        I=Interval(2);
        J=Interval(2);
        assert( (I==J));
        assert(!(I!=J));

        I=Interval(2);
        J=Interval(3);
        assert(!(I==J));
        assert( (I!=J));

        // I < J
        I=Interval(1,2);
        J=Interval(3,4);
        assert( (I<J));
        assert(!(J<I));
        assert( (I<=J));
        assert(!(J<=I));
        assert(!(I>J));
        assert( (J>I));
        assert(!(I>=J));
        assert( (J>=I));

        // OVERLAP
        I=Interval(1,3);
        J=Interval(2,4); // TODO: must explicitly convert to bool to get exception
        CGAL_catch_error((bool)(I==J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I!=J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I< J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I> J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I<=J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I>=J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(J> I),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(J> I),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(J>=I),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(J<=I),CGAL::Uncertain_conversion_exception&);

        // I<=J
        I=Interval(1,2);
        J=Interval(2,3);
        assert( (I<=J));
        assert( (J>=I));
        assert(!(I> J));
        assert(!(J< I));
        CGAL_catch_error((bool)(I==J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I!=J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I< J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(J> I),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I>=J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(J<=I),CGAL::Uncertain_conversion_exception&);

        // degenerated I
        I=Interval(1,1);
        J=Interval(1,1);
        assert( (I==J));
        assert(!(I!=J));
        assert(!(I> J));
        assert(!(I< J));
        assert( (I>=J));
        assert( (I<=J));

        // "I==J"
        I=Interval(1,2);
        J=Interval(1,2);
        CGAL_catch_error((bool)(I==J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I!=J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I< J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I> J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I>=J),CGAL::Uncertain_conversion_exception&);
        CGAL_catch_error((bool)(I<=J),CGAL::Uncertain_conversion_exception&);
    }
    { // Multiplication
        const double inf = std::numeric_limits<double>::infinity();
        Interval all = Interval::largest();
        Interval zero = 0;
        Interval a(0,1);
        Interval b(1,inf);
        Interval c(-1,inf);
        Interval d(0,inf);
        // For CGAL's purposes, [0,0]*anything is 0, but we tolerate larger intervals if they can be produced faster.
        for (Interval I : { all*zero, zero*all, all*0., 0.*all, b*zero, zero*b, -b*zero, zero*-b, c*zero, zero*c, -c*zero, zero*-c })
        {
          //std::cout << I << '\n';
          assert(I.inf()<=0 && I.sup()>=0);
        }
        // This should be [0,inf], but again we tolerate more.
        for (Interval I : { a*b, b*a, -a*-b, -b*-a, d*d, -d*-d })
        {
          //std::cout << I << '\n';
          assert(I.inf()<=0 && I.sup()==inf);
        }
        for (Interval I : { -a*b, -b*a, a*-b, b*-a, -d*d, d*-d })
        {
          //std::cout << I << '\n';
          assert(I.inf()==-inf && I.sup()>=0);
        }
        for (Interval I : { all*a, a*all, all*-a, -a*all })
        {
          //std::cout << I << '\n';
          assert(I.inf()==-inf && I.sup()==inf);
        }
        for (Interval I : { all*d, d*all, all*-d, -d*all })
        {
          //std::cout << I << '\n';
          assert(I.inf()==-inf && I.sup()==inf);
        }
    }
    {// external functions on Intervals
     // functions (abs, square, sqrt, pow)
        {//abs
            {
                Interval I=CGAL_NTS abs(Interval(2.0,3.0));
                assert(I.inf()==2.0);
                assert(I.sup()==3.0);
            }{
                Interval I=CGAL_NTS abs(Interval(-4.0,3.0));
                assert(I.inf()==0.0);
                assert(I.sup()==4.0);
            }{
                Interval I=CGAL_NTS abs(Interval(-4.0,-2.0));
                assert(I.inf()==2.0);
                assert(I.sup()==4.0);
            }
        }{//square
            {
                Interval I=CGAL_NTS square(Interval(2.0,3.0));
                assert(I.inf()==4.0);
                assert(I.sup()==9.0);
            }{
                Interval I=CGAL_NTS square(Interval(-4.0,3.0));
                assert(I.inf()==0.0);
                assert(I.sup()==16.0);
            }{
                Interval I=CGAL_NTS square(Interval(-4.0,-2.0));
                assert(I.inf()==4.0);
                assert(I.sup()==16.0);
            }
        }{//sqrt
            {
                Interval I=CGAL_NTS sqrt(Interval(2.0,3.0));
                assert(1.4<I.inf()&& I.inf()<1.5);
                assert(1.7<I.sup()&& I.sup()<1.8);
                assert(CGAL::in(CGAL_NTS sqrt(2.0),I));
                assert(CGAL::in(CGAL_NTS sqrt(3.0),I));
            }{
                Interval I=CGAL_NTS sqrt(Interval(-4.0,3.0));
                assert(I.inf()==0.0);
                assert(1.7<I.sup()&& I.sup()<1.8);
                assert(CGAL::in(sqrt(3.0),I));
            }{
                  // TODO: Throws no exception
//                CGAL_catch_error(sqrt(Interval(-4.0,-2.0)),
//                                ...);
            }
        }{//pow TODO: Not available for Interval_nt
/*            {
                Interval I=CGAL::pow(Interval(2.0,3.0),3);
                assert(I.inf()==8.0);
                assert(I.sup()==27.0);
            }{
                Interval I=CGAL::pow(Interval(-2.0,3.0),3);
                assert(I.inf()==-8.0);
                assert(I.sup()==27.0);
            }{
                Interval I=CGAL::pow(Interval(-4.0,-2.0),3);
                assert(I.inf()==-64.0);
                assert(I.sup()==-8.0);
            }*/
        }
    }{//  functions min max
        {//min
            {
                Interval I=(CGAL::min)(Interval(-2.0,-1.0),Interval(1.0,4.0));
                assert(I.inf()==-2.0);
                assert(I.sup()==-1.0);
            }{
                Interval I=(CGAL::min)(Interval(2.0,3.0),Interval(1.0,4.0));
                assert(I.inf()==1.0);
                assert(I.sup()==3.0);
            }
        }{//max
            {
                Interval I=(CGAL::max)(Interval(-2.0,-1.0),Interval(1.0,4.0));
                assert(I.inf()==1.0);
                assert(I.sup()==4.0);
            }{
                Interval I=(CGAL::max)(Interval(2.0,3.0),Interval(1.0,4.0));
                assert(I.inf()==2.0);
                assert(I.sup()==4.0);
            }
        }
    }{// functions width, median, singleton
        {//width TODO: Not available for Interval_nt
/*            assert(CGAL::width(Interval(2.0,2.0))==0.0);
            assert(CGAL::width(Interval(2.0,3.0))==1.0);
            assert(CGAL::width(Interval(-2.0,3.0))==5.0);*/
        }{//median TODO: Is simulated by to_double
            assert(CGAL_NTS to_double(Interval(2.0,2.0))==2.0);
            assert(CGAL_NTS to_double(Interval(2.0,3.0))==2.5);
            assert(CGAL_NTS to_double(Interval(-2.0,3.0))==0.5);
        }{ // TODO: Is simulated by is_point
            assert((Interval(3.0,3.0)).is_point() ==true);
            assert((Interval(2.0,3.0)).is_point() ==false);
        }
    }{// functions  in, in_zero
        {//in
            assert(CGAL::in(2.0,Interval( 2.0,2.0))==true);
            assert(CGAL::in(2.5,Interval( 2.0,3.0))==true);
            assert(CGAL::in(4.0,Interval(-2.0,3.0))==false);
        }{//in_zero TODO: Not available for Interval_nt
/*            assert(CGAL::in_zero(Interval( 2.0,2.0))==false);
            assert(CGAL::in_zero(Interval( 2.0,3.0))==false);
            assert(CGAL::in_zero(Interval(-2.0,3.0))==true);
            assert(CGAL::in_zero(Interval(-2.0,0.0))==true);
            assert(CGAL::in_zero(Interval( 0.0,3.0))==true);*/
        }
    }{ // functions equal, subset, proper_subset, overlap, intersect, hull
        {//equal TODO: Is simulated by is_same
            assert(Interval(2.0,2.0).is_same(Interval(2.0,2.0))==true);
            assert(Interval(1.0,2.0).is_same(Interval(1.0,2.0))==true);
            assert(Interval(0.0,2.0).is_same(Interval(1.0,2.0))==false);
        }{//subset TODO: Not available for Interval_nt
/*            assert(CGAL::subset(Interval(2.0,2.0),Interval(2.0,2.0))==true);
            assert(CGAL::subset(Interval(1.0,2.0),Interval(0.0,3.0))==true);
            assert(CGAL::subset(Interval(0.0,2.0),Interval(1.0,2.0))==false);*/
        }{//proper_subset TODO: Not available for Interval_nt
/*            assert(CGAL::proper_subset(Interval(2.0,2.0),Interval(2.0,2.0))
                     ==false);
            assert(CGAL::proper_subset(Interval(1.0,2.0),Interval(0.0,3.0))
                     ==true);
            assert(CGAL::proper_subset(Interval(0.0,2.0),Interval(1.0,2.0))
                     ==false);                          */
        }{//overlap TODO: Is simulated by do_overlap
            assert(Interval( 2,2).do_overlap(Interval(2,2))==true);
            assert(Interval(-1,3).do_overlap(Interval(1,2))==true);
            assert(Interval( 0,2).do_overlap(Interval(2,3))==true);
            assert(Interval( 2,3).do_overlap(Interval(3,4))==true);
            assert(Interval(-2,1).do_overlap(Interval(3,4))==false);
            assert(Interval( 2,3).do_overlap(Interval(5,6))==false);
        }{//intersect TODO: Not available for Interval_nt
/*            {
                Interval I=CGAL::intersect(Interval(2,2),Interval(2,2));
                assert(CGAL::equal(I,Interval(2.0,2.0)));
            }{
                Interval I=CGAL::intersect(Interval(0,2),Interval(1,3));
                assert(CGAL::equal(I,Interval(1,2)));
            }*/
        }{//hull TODO: Not available for Interval_nt
/*            {
                Interval I=CGAL::hull(Interval(2,2),Interval(2,2));
                assert(CGAL::equal(I,Interval(2,2)));
            }{
                Interval I=CGAL::hull(Interval(0,2),Interval(1,3));
                assert(CGAL::equal(I,Interval(0,3)));
            }{
                Interval I=CGAL::hull(Interval(-3,-1),Interval(1,3));
                assert(CGAL::equal(I,Interval(-3,3)));
            }*/
        }
    }

// IO ------------------------------------------------------------------------
/*    {
        // input/output

        Interval tmp1,tmp2;
        // tmp IS_GENERAL = sqrt 2
        tmp1 = Interval(1.5,2.1);
        std::ostringstream os;
        os << LiS::oformat(tmp1);
        std::istringstream is(os.str());
        is >> LiS::iformat(tmp2);
        assert(tmp1.inf()==tmp2.inf());
        assert(tmp1.sup()==tmp2.sup());
    }*/



  return 0;
}
