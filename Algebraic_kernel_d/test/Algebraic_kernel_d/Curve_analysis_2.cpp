// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_d/flags.h>

// Switches on/off tests for Sqrt-extension types
#if CGAL_ACK_WITH_ROTATIONS
#ifndef DO_SQRT_EXTENSION_TESTS
#define DO_SQRT_EXTENSION_TESTS 1
#endif
#endif

#if CGAL_ACK_USE_EXACUS
#include <AcX/Algebraic_curve_2.h>
#include <AcX/Algebraic_curve_pair_2.h>
#endif

#include <sstream>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h>

template<typename Poly_> Poly_ from_string(const char* s) {
    std::stringstream ss(s);
    Poly_ f;
    ss >> f;
    return f;
}

template<typename AlgebraicKernel_2>
int number_of_objects(typename AlgebraicKernel_2::Curve_analysis_2 c) {

    typedef typename AlgebraicKernel_2::Curve_analysis_2::Status_line_1
        Status_line_1;

    int vertical_arcs=0, non_vertical_arcs=0, isolated_vertices=0;

    int n = c.number_of_status_lines_with_event();

    for( int i = 0; i < n; i++ ) {
        const Status_line_1& status_line = c.status_line_at_event(i);
        if(status_line.covers_line()) {
            // vertical
            vertical_arcs += 1 + status_line.number_of_events();
        }
        for( int j = 0; j < status_line.number_of_events(); j++ ) {
            if(status_line.number_of_incident_branches(j).first == 0 &&
               status_line.number_of_incident_branches(j).second == 0) {
                isolated_vertices++;
            }
        }
    }
    for( int i = 0 ; i <= n; i++ ) {
        non_vertical_arcs += c.status_line_of_interval(i).number_of_events();
    }

    return vertical_arcs + non_vertical_arcs + isolated_vertices;
}



template<typename Arithmetic_kernel> void test_routine() {


    typedef typename Arithmetic_kernel::Rational Rational;
    typedef typename Arithmetic_kernel::Integer Integer;

    typedef Integer Coefficient;
    typedef typename
        CGAL::Polynomial_type_generator<Coefficient,1>::Type Poly_int1;
    typedef typename
        CGAL::Polynomial_type_generator<Coefficient,2>::Type Poly_int2;

    typedef CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi
        < Coefficient, Rational > Rep_class;
    typedef CGAL::internal::Bitstream_descartes
        < CGAL::internal::Bitstream_descartes_rndl_tree_traits
            < CGAL::internal::Bitstream_coefficient_kernel<Coefficient>
            >
        >
        Isolator;

    typedef CGAL::Algebraic_kernel_d_1<Coefficient,Rational,Rep_class, Isolator>
        Algebraic_kernel_d_1;

    typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real;



#if CGAL_ACK_USE_EXACUS
    typedef AcX::Algebraic_curve_2<Algebraic_kernel_d_1> Algebraic_curve_2;
    typedef AcX::Algebraic_curve_pair_2<Algebraic_curve_2>
        Algebraic_curve_pair_2;
    typedef CGAL::Algebraic_curve_kernel_2<Algebraic_curve_pair_2,
        Algebraic_kernel_d_1>
        Algebraic_kernel_d_2;
#else
    typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_d_1>
        Algebraic_kernel_d_2;
#endif

    Algebraic_kernel_d_2 kernel;

    typename Algebraic_kernel_d_2::Construct_curve_2 construct_curve_2
        = kernel.construct_curve_2_object();

    typedef typename Algebraic_kernel_d_2::Curve_analysis_2 Curve_analysis_2;

    typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

    Poly_int2 f;

    Curve_analysis_2 curve;

    Status_line_1 event;

    Rational eps(1,1000);

    {
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[1(0,P[1(0,2)(1,2)])(1,P[0(0,-3)])]"
                             << std::endl;
#endif
        f=from_string<Poly_int2>("P[1(0,P[1(0,2)(1,2)])(1,P[0(0,-3)])]");
        curve=construct_curve_2(f);
        assert(curve.number_of_status_lines_with_event()==0);
        assert(number_of_objects<Algebraic_kernel_d_2>(curve)==1);
        //assert(!curve.may_be_singular());

        CGAL::Box_parameter_space_2 loc;

        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,0)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (loc,
                curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,0)));
        assert( loc == CGAL::TOP_BOUNDARY);

    }
    {
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[3(0,P[2(0,-2)(2,2)])(1,P[1(1,-1)])(3,P[1(1,-6)])]" << std::endl;
#endif
        f=from_string<Poly_int2>("P[3(0,P[2(0,-2)(2,2)])(1,P[1(1,-1)])(3,P[1(1,-6)])]");
        ::CGAL::IO::set_pretty_mode(std::cout);
        curve=construct_curve_2(f);
        assert(curve.number_of_status_lines_with_event()==1);
        assert(number_of_objects<Algebraic_kernel_d_2>(curve)==2);
        event=curve.status_line_at_event(0);
        assert(event.number_of_events()==0);

        assert(event.number_of_branches_approaching_minus_infinity().first==0);
        assert(event.number_of_branches_approaching_minus_infinity().second
               ==1);
        assert(event.number_of_branches_approaching_plus_infinity().first==1);
        assert(event.number_of_branches_approaching_plus_infinity().second==0);
        CGAL::Box_parameter_space_2 loc;

        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,0)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (loc,
                curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,0)));
        assert( loc == CGAL::TOP_BOUNDARY);

    }
    {
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[10(1,P[9(4,18)(9,18)])(2,P[9(4,-9)(9,-9)])(3,P[7(2,-12)(7,-12)])(4,P[7(2,6)(3,-18)(7,6)])(5,P[5(3,9)(5,36)])(6,P[5(1,12)(5,-18)])(7,P[3(1,-6)(3,-42)])(8,P[3(3,21)])(9,P[1(1,12)])(10,P[1(1,-6)])]" << std::endl;
#endif
        f=from_string<Poly_int2>("P[10(1,P[9(4,18)(9,18)])(2,P[9(4,-9)(9,-9)])(3,P[7(2,-12)(7,-12)])(4,P[7(2,6)(3,-18)(7,6)])(5,P[5(3,9)(5,36)])(6,P[5(1,12)(5,-18)])(7,P[3(1,-6)(3,-42)])(8,P[3(3,21)])(9,P[1(1,12)])(10,P[1(1,-6)])]");
        curve=construct_curve_2(f);
        assert(curve.number_of_status_lines_with_event()==10);
        //    assert(curve.may_be_singular());

        event=curve.status_line_at_exact_x(Algebraic_real(1));
        assert(event.number_of_events()==6);
        event.refine_to(0,eps);
        assert(event.lower_bound(0) > Rational(-169,100));
        assert(event.upper_bound(0) < Rational(-168,100));
        assert(! event.is_event(0));
        event.refine_to(3,eps);
        assert(event.lower_bound(3) > Rational(122,100));
        assert(event.upper_bound(3) < Rational(123,100));
        assert(! event.is_event(0));

        event=curve.status_line_at_exact_x(Algebraic_real(0));
        assert(event.covers_line());
        assert(event.number_of_events()==3);
        event.refine_to(0,eps);
        assert(event.lower_bound(0) > Rational(-101,100));
        assert(event.upper_bound(0) < Rational(-99,100));
        assert(! event.is_event(0));
        event.refine_to(1,eps);
        assert(event.lower_bound(1) > Rational(-1,100));
        assert(event.upper_bound(1) < Rational(1,100));
        assert(event.is_event(1));
        assert(event.number_of_incident_branches(1).first==4);
        assert(event.number_of_incident_branches(1).second==4);
        event.refine_to(2,eps);
        assert(event.lower_bound(2) > Rational(199,100));
        assert(event.upper_bound(2) < Rational(201,100));
        assert(! event.is_event(0));

        event=curve.status_line_at_exact_x(Algebraic_real(Poly_int1(-3,0,1),0,2));
        assert(event.number_of_events()==6);
        event.refine_to(1,eps);
        assert(event.lower_bound(1) > Rational(-213,100));
        assert(event.upper_bound(1) < Rational(-212,100));
        assert(! event.is_event(0));
        event.refine_to(3,eps);
        assert(event.lower_bound(3) > Rational(199,100));
        assert(event.upper_bound(3) < Rational(201,100));
        assert(! event.is_event(0));

        CGAL::Box_parameter_space_2 loc;
        Algebraic_real y_coor;


        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,0)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,1)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (y_coor,
                curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,2)));
        assert( y_coor == Rational(0));

        assert(CGAL::assign
               (y_coor,
                curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,3)));
        assert( y_coor == Rational(2));

        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,4)));
        assert( loc == CGAL::TOP_BOUNDARY);

        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,5)));
        assert( loc == CGAL::TOP_BOUNDARY);


        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,0)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,1)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (y_coor,
                curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,2)));
        assert( y_coor == Rational(0));

        assert(CGAL::assign
               (y_coor,
                curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,3)));
        assert( y_coor == Rational(2));

        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,4)));
        assert( loc == CGAL::TOP_BOUNDARY);

        assert(CGAL::assign
               (loc,curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,5)));
        assert( loc == CGAL::TOP_BOUNDARY);


    }
    {
        curve=Curve_analysis_2();
#if !CGAL_ACK_USE_EXACUS
        assert(! curve.has_defining_polynomial());
#endif
    }
    {
        Poly_int2 f1
            = from_string<Poly_int2>("P[1(0,P[1(0,2)(1,2)])(1,P[0(0,-3)])]");
        Poly_int2 f2
            = from_string<Poly_int2>("P[1(0,P[1(0,3)(1,2)])(1,P[0(0,-3)])]");
        Curve_analysis_2 curve1=construct_curve_2(f1),
            curve2=construct_curve_2(f2);
        assert(!curve1.is_identical(curve2));
        assert(curve1.polynomial_2()!=curve2.polynomial_2());
        std::swap(curve1,curve2);
        assert(curve1.polynomial_2()!=curve2.polynomial_2());
    }
    {
        Poly_int2 f1
            = from_string<Poly_int2>("P[1(0,P[1(0,2)(1,2)])(1,P[0(0,-3)])]");
        Curve_analysis_2 curve1  = construct_curve_2(f1);
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[1(0,P[1(0,2)(1,2)])(1,P[0(0,-3)])]"
                             << std::endl;
#endif
        Curve_analysis_2 curve2(curve1);
        assert(curve1.is_identical(curve2));
        Poly_int2 new_f=from_string<Poly_int2>("P[1(0,P[1(0,3)(1,2)])(1,P[0(0,-3)])]");
    }
    {
        Poly_int2 f =from_string<Poly_int2>("P[2(0,P[2(1,3)(2,6)])(1,P[2(0,-3)(1,-11)(2,1)])(2,P[1(0,5)(1,-1)])]");
        Curve_analysis_2 curve = construct_curve_2(f);
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[2(0,P[2(1,3)(2,6)])(1,P[2(0,-3)(1,-11)(2,1)])(2,P[1(0,5)(1,-1)])]" << std::endl;
#endif
#if !CGAL_ACK_USE_EXACUS
        assert(curve.polynomial_2()==curve.primitive_polynomial_2());
        assert(number_of_objects<Algebraic_kernel_d_2>(curve)==4);
        Curve_analysis_2 sh_curve=curve.shear_primitive_part(2);
        assert(number_of_objects<Algebraic_kernel_d_2>(sh_curve)==7);
        assert(sh_curve.status_line_at_exact_x(Algebraic_real(-10)).number_of_events()==3);
        // Now, this should be cached
        sh_curve=curve.shear_primitive_part(2);
        assert(number_of_objects<Algebraic_kernel_d_2>(sh_curve)==7);
        assert(sh_curve.status_line_at_exact_x(Algebraic_real(-10)).number_of_events()==3);
#endif
    }

    { // More tests...just analyse some curves and compute their segments
        Poly_int2 f = from_string<Poly_int2>("P[8(0,P[8(0,24)(1,-8)(2,-162)(3,204)(4,106)(5,-340)(6,240)(7,-72)(8,8)])(1,P[6(0,-60)(1,8)(2,304)(3,-400)(4,148)(5,8)(6,-8)])(2,P[6(0,18)(1,80)(2,-165)(3,-132)(4,367)(5,-212)(6,38)])(3,P[4(0,-30)(1,-136)(2,264)(3,-72)(4,-26)])(4,P[4(0,-15)(1,36)(2,89)(3,-144)(4,49)])(5,P[2(0,30)(1,-24)(2,-6)])(6,P[2(0,-6)(1,-28)(2,22)])(8,P[0(0,3)])]");
        Curve_analysis_2 curve= construct_curve_2(f);
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[8(0,P[8(0,24)(1,-8)(2,-162)(3,204)(4,106)(5,-340)(6,240)(7,-72)(8,8)])(1,P[6(0,-60)(1,8)(2,304)(3,-400)(4,148)(5,8)(6,-8)])(2,P[6(0,18)(1,80)(2,-165)(3,-132)(4,367)(5,-212)(6,38)])(3,P[4(0,-30)(1,-136)(2,264)(3,-72)(4,-26)])(4,P[4(0,-15)(1,36)(2,89)(3,-144)(4,49)])(5,P[2(0,30)(1,-24)(2,-6)])(6,P[2(0,-6)(1,-28)(2,22)])(8,P[0(0,3)])]" << std::endl;
#endif
        assert(number_of_objects<Algebraic_kernel_d_2>(curve)==54);
        f = from_string<Poly_int2>("P[5(0,P[8(0,40)(1,-40)(2,-20)(3,20)(5,-8)(6,8)(7,4)(8,-4)])(1,P[8(0,-100)(1,100)(2,50)(3,-50)(5,20)(6,-20)(7,-10)(8,10)])(2,P[8(0,60)(1,-60)(2,-30)(3,30)(5,-12)(6,12)(7,6)(8,-6)])(3,P[8(0,-20)(1,20)(2,10)(3,-10)(5,4)(6,-4)(7,-2)(8,2)])(4,P[8(0,50)(1,-50)(2,-25)(3,25)(5,-10)(6,10)(7,5)(8,-5)])(5,P[8(0,-30)(1,30)(2,15)(3,-15)(5,6)(6,-6)(7,-3)(8,3)])]");
        curve=construct_curve_2(f);
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[5(0,P[8(0,40)(1,-40)(2,-20)(3,20)(5,-8)(6,8)(7,4)(8,-4)])(1,P[8(0,-100)(1,100)(2,50)(3,-50)(5,20)(6,-20)(7,-10)(8,10)])(2,P[8(0,60)(1,-60)(2,-30)(3,30)(5,-12)(6,12)(7,6)(8,-6)])(3,P[8(0,-20)(1,20)(2,10)(3,-10)(5,4)(6,-4)(7,-2)(8,2)])(4,P[8(0,50)(1,-50)(2,-25)(3,25)(5,-10)(6,10)(7,5)(8,-5)])(5,P[8(0,-30)(1,30)(2,15)(3,-15)(5,6)(6,-6)(7,-3)(8,3)])]" << std::endl;
#endif
        assert(number_of_objects<Algebraic_kernel_d_2>(curve)==31);
        f=from_string<Poly_int2>("P[8(0,P[10(5,80)(6,104)(7,44)(8,10)(9,4)(10,1)])(1,P[8(4,-80)(5,-36)(6,54)(7,32)(8,3)])(2,P[8(2,-80)(3,-104)(4,-72)(5,-66)(6,-8)(7,5)(8,1)])(3,P[7(1,80)(2,36)(3,-94)(4,-60)(5,-19)(6,-5)(7,-1)])(4,P[5(1,28)(2,56)(3,-10)(4,-30)(5,-6)])(5,P[5(0,40)(1,28)(2,16)(3,7)(4,-3)(5,-1)])(6,P[3(0,14)(1,24)(2,5)(3,1)])(7,P[2(0,-2)(1,4)(2,1)])(8,P[0(0,-1)])]");
        curve=construct_curve_2(f);
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[8(0,P[10(5,80)(6,104)(7,44)(8,10)(9,4)(10,1)])(1,P[8(4,-80)(5,-36)(6,54)(7,32)(8,3)])(2,P[8(2,-80)(3,-104)(4,-72)(5,-66)(6,-8)(7,5)(8,1)])(3,P[7(1,80)(2,36)(3,-94)(4,-60)(5,-19)(6,-5)(7,-1)])(4,P[5(1,28)(2,56)(3,-10)(4,-30)(5,-6)])(5,P[5(0,40)(1,28)(2,16)(3,7)(4,-3)(5,-1)])(6,P[3(0,14)(1,24)(2,5)(3,1)])(7,P[2(0,-2)(1,4)(2,1)])(8,P[0(0,-1)])]" << std::endl;
#endif
        assert(number_of_objects<Algebraic_kernel_d_2>(curve)==27);
        f=from_string<Poly_int2>("P[4(0,P[1(0,-10)(1,6)])(1,P[2(0,4)(1,8)(2,-6)])(2,P[2(0,-5)(1,-6)(2,5)])(3,P[3(0,2)(1,6)(2,1)(3,-3)])(4,P[3(1,-2)(2,-1)(3,1)])]");
        curve=construct_curve_2(f);
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[4(0,P[1(0,-10)(1,6)])(1,P[2(0,4)(1,8)(2,-6)])(2,P[2(0,-5)(1,-6)(2,5)])(3,P[3(0,2)(1,6)(2,1)(3,-3)])(4,P[3(1,-2)(2,-1)(3,1)])]" << std::endl;
#endif
        assert(number_of_objects<Algebraic_kernel_d_2>(curve)==18);
        f=from_string<Poly_int2>("P[7(0,P[7(3,-12)(4,-12)(5,7)(6,4)(7,-1)])(1,P[6(1,12)(2,24)(3,-31)(4,-11)(5,9)(6,1)])(2,P[5(0,-12)(1,24)(2,43)(3,-22)(4,-9)(5,5)])(3,P[5(0,-36)(1,14)(2,22)(4,-5)(5,1)])(4,P[4(0,-14)(1,-5)(2,-4)(3,5)(4,-1)])(5,P[3(0,9)(1,-6)(2,-5)(3,1)])(6,P[2(0,6)(1,-1)(2,-1)])(7,P[0(0,1)])]");
        curve=construct_curve_2(f);
    }
#if CGAL_ACK_WITH_ROTATIONS
#if DO_SQRT_EXTENSION_TESTS

    // Tests for Sqrt_extension
    {
        typedef CGAL::Sqrt_extension<Integer, Integer> Sqrt_extension;

        typedef Sqrt_extension Coefficient;
        typedef typename
            CGAL::Polynomial_type_generator<Coefficient,1>::Type Poly_sqrt1;
        typedef typename
            CGAL::Polynomial_type_generator<Coefficient,2>::Type Poly_sqrt2;

        typedef CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi
            < Coefficient, Rational > Rep_class;
        typedef CGAL::internal::Bitstream_descartes
            < CGAL::internal::Bitstream_descartes_rndl_tree_traits
                < CGAL::internal::Bitstream_coefficient_kernel<Coefficient>
                >
            >
        Isolator;

        typedef CGAL::Algebraic_kernel_d_1<Coefficient,Rational,Rep_class, Isolator>
            Algebraic_kernel_d_1_with_sqrt;

        typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_d_1_with_sqrt>
            Algebraic_kernel_d_2_with_sqrt;

        Algebraic_kernel_d_2_with_sqrt kernel;

        typename Algebraic_kernel_d_2_with_sqrt::Construct_curve_2
            sqrt_construct_curve_2
              = kernel.construct_curve_2_object();

        typedef typename Algebraic_kernel_d_2_with_sqrt::Curve_analysis_2
            Sqrt_curve_analysis_2;

        typedef typename Sqrt_curve_analysis_2::Algebraic_real_1
            Sqrt_algebraic_real;

        typedef typename Sqrt_curve_analysis_2::Status_line_1
            Sqrt_status_line_1;



        Sqrt_status_line_1 event;



        Poly_sqrt2 sqrt_f=from_string<Poly_sqrt2>("P[10(1,P[9(4,EXT[0,54,5])(9,EXT[0,54,5])])(2,P[9(4,EXT[0,-27,5])(9,EXT[0,-27,5])])(3,P[7(2,EXT[0,-36,5])(7,EXT[0,-36,5])])(4,P[7(2,EXT[0,18,5])(3,EXT[0,-54,5])(7,EXT[0,18,5])])(5,P[5(3,EXT[0,27,5])(5,EXT[0,108,5])])(6,P[5(1,EXT[0,36,5])(5,EXT[0,-54,5])])(7,P[3(1,EXT[0,-18,5])(3,EXT[0,-126,5])])(8,P[3(3,EXT[0,63,5])])(9,P[1(1,EXT[0,36,5])])(10,P[1(1,EXT[0,-18,5])])]");
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "P[10(1,P[9(4,EXT[0,54,5])(9,EXT[0,54,5])])(2,P[9(4,EXT[0,-27,5])(9,EXT[0,-27,5])])(3,P[7(2,EXT[0,-36,5])(7,EXT[0,-36,5])])(4,P[7(2,EXT[0,18,5])(3,EXT[0,-54,5])(7,EXT[0,18,5])])(5,P[5(3,EXT[0,27,5])(5,EXT[0,108,5])])(6,P[5(1,EXT[0,36,5])(5,EXT[0,-54,5])])(7,P[3(1,EXT[0,-18,5])(3,EXT[0,-126,5])])(8,P[3(3,EXT[0,63,5])])(9,P[1(1,EXT[0,36,5])])(10,P[1(1,EXT[0,-18,5])])]" << std::endl;
#endif
        Sqrt_curve_analysis_2 sqrt_curve= sqrt_construct_curve_2(sqrt_f);
        assert(sqrt_curve.number_of_status_lines_with_event()==10);
        //    assert(sqrt_curve.may_be_singular());

        event=sqrt_curve.status_line_at_exact_x(Sqrt_algebraic_real(1));
        assert(event.number_of_events()==6);
        event.refine_to(0,eps);
        assert(event.lower_bound(0) > Rational(-169,100));
        assert(event.upper_bound(0) < Rational(-168,100));
        assert(! event.is_event(0));
        event.refine_to(3,eps);
        assert(event.lower_bound(3) > Rational(122,100));
        assert(event.upper_bound(3) < Rational(123,100));
        assert(! event.is_event(0));

        event=sqrt_curve.status_line_at_exact_x(Sqrt_algebraic_real(0));
        assert(event.covers_line());
        assert(event.number_of_events()==3);
        event.refine_to(0,eps);
        assert(event.lower_bound(0) > Rational(-101,100));
        assert(event.upper_bound(0) < Rational(-99,100));
        assert(! event.is_event(0));
        event.refine_to(1,eps);
        assert(event.lower_bound(1) > Rational(-1,100));
        assert(event.upper_bound(1) < Rational(1,100));
        assert(event.is_event(1));
        assert(event.number_of_incident_branches(1).first==4);
        assert(event.number_of_incident_branches(1).second==4);
        event.refine_to(2,eps);
        assert(event.lower_bound(2) > Rational(199,100));
        assert(event.upper_bound(2) < Rational(201,100));
        assert(! event.is_event(0));

        event=sqrt_curve.status_line_at_exact_x(Sqrt_algebraic_real(Poly_sqrt1(-3,0,1),0,2));
        assert(event.number_of_events()==6);
        event.refine_to(1,eps);
        assert(event.lower_bound(1) > Rational(-213,100));
        assert(event.upper_bound(1) < Rational(-212,100));
        assert(! event.is_event(0));
        event.refine_to(3,eps);
        assert(event.lower_bound(3) > Rational(199,100));
        assert(event.upper_bound(3) < Rational(201,100));
        assert(! event.is_event(0));


        CGAL::Box_parameter_space_2 loc;
        Sqrt_algebraic_real y_coor;


        assert(CGAL::assign
               (loc,sqrt_curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,0)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (loc,sqrt_curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,1)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (y_coor,
                sqrt_curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,2)));
        assert( y_coor == Rational(0));

        assert(CGAL::assign
               (y_coor,
                sqrt_curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,3)));
        assert( y_coor == Rational(2));

        assert(CGAL::assign
               (loc,sqrt_curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,4)));
        assert( loc == CGAL::TOP_BOUNDARY);

        assert(CGAL::assign
               (loc,sqrt_curve.asymptotic_value_of_arc(CGAL::LEFT_BOUNDARY,5)));
        assert( loc == CGAL::TOP_BOUNDARY);


        assert(CGAL::assign
               (loc,sqrt_curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,0)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (loc,sqrt_curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,1)));
        assert( loc == CGAL::BOTTOM_BOUNDARY);

        assert(CGAL::assign
               (y_coor,
                sqrt_curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,2)));
        assert( y_coor == Rational(0));

        assert(CGAL::assign
               (y_coor,
                sqrt_curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,3)));
        assert( y_coor == Rational(2));

        assert(CGAL::assign
               (loc,sqrt_curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,4)));
        assert( loc == CGAL::TOP_BOUNDARY);

        assert(CGAL::assign
               (loc,sqrt_curve.asymptotic_value_of_arc(CGAL::RIGHT_BOUNDARY,5)));
        assert( loc == CGAL::TOP_BOUNDARY);

    }
#endif
#endif

}


int main() {

#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
    test_routine<CGAL::LEDA_arithmetic_kernel>();
#else
    std::cerr << "LEDA tests skipped" << std::endl;
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
    test_routine<CGAL::CORE_arithmetic_kernel>();
#else
    std::cerr << "CORE tests skipped" << std::endl;
#endif
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
    test_routine<CGAL::GMP_arithmetic_kernel>();
#else
    std::cerr << "GMP tests skipped" << std::endl;
#endif
    return 0;
}
