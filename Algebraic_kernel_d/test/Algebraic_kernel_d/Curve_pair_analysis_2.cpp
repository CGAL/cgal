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


#include <sstream>

// In the testsuite we do not want that pairs of Curve_analysis_2 are ordered
#define CGAL_ALGEBRAIC_KERNEL_DONT_SWAP 1

#if CGAL_ACK_USE_EXACUS
#include <AcX/Algebraic_curve_2.h>
#include <AcX/Algebraic_curve_pair_2.h>
#endif

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

template <typename Arithmetic_kernel>
void test_routine() {


    typedef typename Arithmetic_kernel::Rational Rational;
    typedef typename Arithmetic_kernel::Integer Integer;

    typedef Integer Coefficient;
    typedef typename
        CGAL::Polynomial_type_generator<Coefficient,1>::Type Poly_1;
    CGAL_USE_TYPE(Poly_1);
    typedef typename
        CGAL::Polynomial_type_generator<Coefficient,2>::Type Poly_2;

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
    CGAL_USE_TYPE(Algebraic_real);

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

    typedef typename Algebraic_kernel_d_2::Curve_analysis_2 Curve_analysis_2;

    typedef typename Algebraic_kernel_d_2::Curve_pair_analysis_2
        Curve_pair_analysis_2;

    typename Algebraic_kernel_d_2::Construct_curve_2
        construct_curve_2
        = kernel.construct_curve_2_object();

    typename Algebraic_kernel_d_2::Construct_curve_pair_2
        construct_curve_pair_2
        = kernel.construct_curve_pair_2_object();


    {
        Poly_2 f=from_string<Poly_2>("P[4(0,P[4(3,-1)(4,2)])(2,P[1(1,1)])(4,P[0(0,1)])]");
        Poly_2 g=from_string<Poly_2>("P[4(0,P[4(4,1)])(1,P[2(2,1)])(3,P[0(0,-1)])(4,P[0(0,2)])]");
        Curve_analysis_2 c1=construct_curve_2(f);
        Curve_analysis_2 c2=construct_curve_2(g);
        Curve_pair_analysis_2 curve_pair=construct_curve_pair_2(c1,c2);
        assert(curve_pair.number_of_status_lines_with_event()==10);
        typedef typename Curve_pair_analysis_2::Status_line_1 Status_line_1;
#if CGAL_ACK_USE_EXACUS
        typedef SoX::Index_triple Triple;
#else
        typedef CGAL::internal::Event_indices<int> Triple;
#endif
        int i;
        {
            i=0;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==0);
        }
        {
            i=0;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==0);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==0);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==1);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
        }
        {
            i=1;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
        }
        {
            i=1;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==0);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==3);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==2);
        }
        {
            i=2;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==4);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==2);
            assert(slice.curves_at_event(3).first==-1);
            assert(slice.curves_at_event(3).second==3);
        }
        {
            i=2;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==0);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==0);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==6);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
            assert(slice.curves_at_event(4).first==-1);
            assert(slice.curves_at_event(4).second==2);
            assert(slice.curves_at_event(5).first==-1);
            assert(slice.curves_at_event(5).second==3);
        }
        {
            i=3;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==8);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==0);
            assert(slice.curves_at_event(3).first==-1);
            assert(slice.curves_at_event(3).second==1);
            assert(slice.curves_at_event(4).first==2);
            assert(slice.curves_at_event(4).second==-1);
            assert(slice.curves_at_event(5).first==3);
            assert(slice.curves_at_event(5).second==-1);
            assert(slice.curves_at_event(6).first==-1);
            assert(slice.curves_at_event(6).second==2);
            assert(slice.curves_at_event(7).first==-1);
            assert(slice.curves_at_event(7).second==3);
        }
        {
            i=3;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==0);
            assert(triple.ffy==-1);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(slice.is_intersection());
            assert(slice.number_of_events()==7);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==0);
            assert(slice.curves_at_event(3).first==-1);
            assert(slice.curves_at_event(3).second==1);
            assert(slice.curves_at_event(4).first==2);
            assert(slice.curves_at_event(4).second==-1);
            assert(slice.curves_at_event(5).first==3);
            assert(slice.curves_at_event(5).second==2);
            assert(slice.curves_at_event(6).first==-1);
            assert(slice.curves_at_event(6).second==3);

            assert(slice.multiplicity_of_intersection(5)==1);
        }
        {
            i=4;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==8);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==0);
            assert(slice.curves_at_event(3).first==-1);
            assert(slice.curves_at_event(3).second==1);
            assert(slice.curves_at_event(4).first==2);
            assert(slice.curves_at_event(4).second==-1);
            assert(slice.curves_at_event(5).first==-1);
            assert(slice.curves_at_event(5).second==2);
            assert(slice.curves_at_event(6).first==3);
            assert(slice.curves_at_event(6).second==-1);
            assert(slice.curves_at_event(7).first==-1);
            assert(slice.curves_at_event(7).second==3);
        }
        {
            i=4;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==1);
            assert(triple.ffy==1);
            assert(triple.ggy==2);
            //assert(slice.event_of_curve(i,0)==1);
            //assert(slice.event_of_curve(i,1)==2);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(slice.is_intersection());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
        }
        {
            i=5;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==6);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==0);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
            assert(slice.curves_at_event(4).first==-1);
            assert(slice.curves_at_event(4).second==2);
            assert(slice.curves_at_event(5).first==-1);
            assert(slice.curves_at_event(5).second==3);
        }
        {
            i=5;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==3);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==5);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==0);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
            assert(slice.curves_at_event(4).first==-1);
            assert(slice.curves_at_event(4).second==2);
        }
        {
            i=6;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==4);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==0);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
        }
        {
            i=6;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==2);
            assert(triple.ffy==-1);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(slice.is_intersection());
            assert(slice.number_of_events()==3);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
            assert(slice.curves_at_event(2).first==1);
            assert(slice.curves_at_event(2).second==-1);
            assert(slice.multiplicity_of_intersection(0)==1);
        }
        {
            i=7;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==4);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
        }
        {
            i=7;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==4);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==4);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==3);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==1);
            assert(slice.curves_at_event(2).second==-1);
        }
        {
            i=8;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==2);
             assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
        }
        {
            i=8;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==2);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==2);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==1);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
        }
        {
            i=9;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==0);
        }
        {
            i=9;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==3);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==3);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==0);
        }
        {
            i=10;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==0);
        }

    }
    {
        Poly_2 f=from_string<Poly_2>("P[4(0,P[1(0,1)(1,1)])(4,P[0(0,-1)])]");
        Poly_2 g=from_string<Poly_2>("P[2(0,P[2(0,-5)(2,6)])(2,P[0(0,4)])]");
        Curve_analysis_2 ca1=construct_curve_2(f),
            ca2=construct_curve_2(g);

        Curve_pair_analysis_2 curve_pair=construct_curve_pair_2(ca1,ca2);
        assert(curve_pair.number_of_status_lines_with_event()==7);
        typedef typename Curve_pair_analysis_2::Status_line_1 Status_line_1;
#if CGAL_ACK_USE_EXACUS
        typedef SoX::Index_triple Triple;
#else
        typedef CGAL::internal::Event_indices<int> Triple;
#endif
        int i;
        {
            i=0;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==0);
        }
        {
            i=0;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==0);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==0);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==1);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
        }
        {
            i=1;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
        }
        {
            i=1;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==0);
            assert(triple.ffy==-1);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
        }
        {
            i=2;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
        }
        {
            i=2;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==0);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==0);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==3);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==1);
            assert(slice.curves_at_event(2).second==-1);

        }
        {
            i=3;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==4);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
        }
        {
            i=3;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==1);
            assert(triple.ffy==-1);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(slice.is_intersection());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==1);
            assert(slice.multiplicity_of_intersection(0)==1);
            assert(slice.multiplicity_of_intersection(1)==1);
        }
        {
            i=4;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==4);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==0);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==1);
            assert(slice.curves_at_event(2).second==-1);
            assert(slice.curves_at_event(3).first==-1);
            assert(slice.curves_at_event(3).second==1);
        }
        {
            i=4;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==2);
            assert(triple.ffy==-1);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(slice.is_intersection());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==1);
            assert(slice.multiplicity_of_intersection(0)==1);
            assert(slice.multiplicity_of_intersection(1)==1);
        }
        {
            i=5;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==4);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
        }
        {
            i=5;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==3);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==1);
            assert(slice.curves_at_event(2).second==-1);
        }
        {
            i=6;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
        }
        {
            i=6;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==3);
            assert(triple.ffy==-1);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
        }
        {
            i=7;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
        }
    }
    {
        Poly_2 f=from_string<Poly_2>("P[4(0,P[4(3,-1)(4,2)])(2,P[1(1,1)])(4,P[0(0,1)])]");
        Poly_2 g=from_string<Poly_2>("P[4(0,P[4(0,12)(1,-380)(2,4200)(3,-18000)(4,20000)])(2,P[1(0,-4000)(1,40000)])(4,P[0(0,160000)])]");
        Curve_analysis_2 ca1=construct_curve_2(f),
            ca2=construct_curve_2(g);
        Curve_pair_analysis_2 curve_pair=construct_curve_pair_2(ca1,ca2);
        //assert(curve_pair.number_of_status_lines_with_event()==11);


    }
#if CGAL_ACK_WITH_ROTATIONS
#if DO_SQRT_EXTENSION_TESTS
    // Sqrt extension:
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

        typedef CGAL::Algebraic_kernel_d_1< Coefficient,Rational,
                                          Rep_class, Isolator >
            Algebraic_kernel_d_1_with_sqrt;

        Poly_sqrt2 f,g;
        f=from_string<Poly_sqrt2>("P[4(0,P[4(3,EXT[0,-2,7])(4,EXT[0,4,7])])(2,P[1(1,EXT[0,2,7])])(4,P[0(0,EXT[0,2,7])])]");
        g=from_string<Poly_sqrt2>("P[4(0,P[4(4,EXT[0,2,7])])(1,P[2(2,EXT[0,2,7])])(3,P[0(0,EXT[0,-2,7])])(4,P[0(0,EXT[0,4,7])])]");

        typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_d_1_with_sqrt>
            Algebraic_kernel_d_2;

        typedef typename Algebraic_kernel_d_2::Curve_analysis_2 Curve_analysis_2;
        typedef typename Algebraic_kernel_d_2::Curve_pair_analysis_2
            Curve_pair_analysis_2;

        Algebraic_kernel_d_2 kernel;

        typename Algebraic_kernel_d_2::Construct_curve_2
            construct_curve_2
            = kernel.construct_curve_2_object();

        typename Algebraic_kernel_d_2::Construct_curve_pair_2
            construct_curve_pair_2
            = kernel.construct_curve_pair_2_object();


        Curve_analysis_2 c1=construct_curve_2(f),
            c2=construct_curve_2(g);

        Curve_pair_analysis_2 curve_pair=construct_curve_pair_2(c1,c2);
        assert(curve_pair.number_of_status_lines_with_event()==10);
        typedef typename Curve_pair_analysis_2::Status_line_1 Status_line_1;
        typedef CGAL::internal::Event_indices<int> Triple;
        int i;

        {
            i=0;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==0);
        }
        {
            i=0;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==0);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==0);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==1);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
        }
        {
            i=1;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
        }
        {
            i=1;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==0);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==3);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==2);
        }
        {
            i=2;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==4);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==2);
            assert(slice.curves_at_event(3).first==-1);
            assert(slice.curves_at_event(3).second==3);
        }
        {
            i=2;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==0);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==0);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==6);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
            assert(slice.curves_at_event(4).first==-1);
            assert(slice.curves_at_event(4).second==2);
            assert(slice.curves_at_event(5).first==-1);
            assert(slice.curves_at_event(5).second==3);
        }
        {
            i=3;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==8);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==0);
            assert(slice.curves_at_event(3).first==-1);
            assert(slice.curves_at_event(3).second==1);
            assert(slice.curves_at_event(4).first==2);
            assert(slice.curves_at_event(4).second==-1);
            assert(slice.curves_at_event(5).first==3);
            assert(slice.curves_at_event(5).second==-1);
            assert(slice.curves_at_event(6).first==-1);
            assert(slice.curves_at_event(6).second==2);
            assert(slice.curves_at_event(7).first==-1);
            assert(slice.curves_at_event(7).second==3);
        }
        {
            i=3;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==0);
            assert(triple.ffy==-1);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(slice.is_intersection());
            assert(slice.number_of_events()==7);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==0);
            assert(slice.curves_at_event(3).first==-1);
            assert(slice.curves_at_event(3).second==1);
            assert(slice.curves_at_event(4).first==2);
            assert(slice.curves_at_event(4).second==-1);
            assert(slice.curves_at_event(5).first==3);
            assert(slice.curves_at_event(5).second==2);
            assert(slice.curves_at_event(6).first==-1);
            assert(slice.curves_at_event(6).second==3);

            assert(slice.multiplicity_of_intersection(5)==1);
        }
        {
            i=4;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==8);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==0);
            assert(slice.curves_at_event(3).first==-1);
            assert(slice.curves_at_event(3).second==1);
            assert(slice.curves_at_event(4).first==2);
            assert(slice.curves_at_event(4).second==-1);
            assert(slice.curves_at_event(5).first==-1);
            assert(slice.curves_at_event(5).second==2);
            assert(slice.curves_at_event(6).first==3);
            assert(slice.curves_at_event(6).second==-1);
            assert(slice.curves_at_event(7).first==-1);
            assert(slice.curves_at_event(7).second==3);
        }
        {
            i=4;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==1);
            assert(triple.ffy==1);
            assert(triple.ggy==2);
            //assert(slice.event_of_curve(i,0)==1);
            //assert(slice.event_of_curve(i,1)==2);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(slice.is_intersection());
            assert(slice.number_of_events()==2);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
        }
        {
            i=5;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==6);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==0);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
            assert(slice.curves_at_event(4).first==-1);
            assert(slice.curves_at_event(4).second==2);
            assert(slice.curves_at_event(5).first==-1);
            assert(slice.curves_at_event(5).second==3);
        }
        {
            i=5;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==3);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==5);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==0);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
            assert(slice.curves_at_event(4).first==-1);
            assert(slice.curves_at_event(4).second==2);
        }
        {
            i=6;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==4);
            assert(slice.curves_at_event(0).first==-1);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==0);
            assert(slice.curves_at_event(1).second==-1);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
        }
        {
            i=6;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==2);
            assert(triple.ffy==-1);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(slice.is_intersection());
            assert(slice.number_of_events()==3);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==0);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==1);
            assert(slice.curves_at_event(2).first==1);
            assert(slice.curves_at_event(2).second==-1);
            assert(slice.multiplicity_of_intersection(0)==1);
        }
        {
            i=7;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==4);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==-1);
            assert(slice.curves_at_event(2).second==1);
            assert(slice.curves_at_event(3).first==1);
            assert(slice.curves_at_event(3).second==-1);
        }
        {
            i=7;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==-1);
            assert(triple.ggy==4);
            //assert(slice.event_of_curve(i,0)==-1);
            //assert(slice.event_of_curve(i,1)==4);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==3);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==-1);
            assert(slice.curves_at_event(1).second==0);
            assert(slice.curves_at_event(2).first==1);
            assert(slice.curves_at_event(2).second==-1);
        }
        {
            i=8;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==2);
             assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
            assert(slice.curves_at_event(1).first==1);
            assert(slice.curves_at_event(1).second==-1);
        }
        {
            i=8;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==2);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==2);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==1);
            assert(slice.curves_at_event(0).first==0);
            assert(slice.curves_at_event(0).second==-1);
        }
        {
            i=9;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==0);
        }
        {
            i=9;
            const Status_line_1& slice=curve_pair.status_line_at_event(i);
            Triple triple = curve_pair.event_indices(i);
            assert(triple.fg==-1);
            assert(triple.ffy==3);
            assert(triple.ggy==-1);
            //assert(slice.event_of_curve(i,0)==3);
            //assert(slice.event_of_curve(i,1)==-1);
            assert(slice.index()==i);
            assert(slice.is_event());
            assert(! slice.is_intersection());
            assert(slice.number_of_events()==0);
        }
        {
            i=10;
            const Status_line_1& slice=curve_pair.status_line_of_interval(i);
            assert(! slice.is_event());
            assert(slice.number_of_events()==0);
        }

    }
#endif
#endif
}


int main() {
#ifdef NDEBUG
    std::cout << "Assertions switched off!" << std::endl;
    return 0;
#endif
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
    test_routine<CGAL::LEDA_arithmetic_kernel>();
#else
    std::cout << "LEDA tests skipped!" << std::endl;
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
    test_routine<CGAL::CORE_arithmetic_kernel>();
#else
    std::cout << "CORE tests skipped!" << std::endl;
#endif
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
    test_routine<CGAL::GMP_arithmetic_kernel>();
#else
    std::cout << "GMP tests skipped!" << std::endl;
#endif
    return 0;
}
