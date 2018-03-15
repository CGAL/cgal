#define  CGAL_ACK_DEBUG_FLAG 1

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



#include <CGAL/basic.h>

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

template <typename Arithmetic_kernel> 
void test_routine() {

        
    typedef typename Arithmetic_kernel::Rational Rational;
    typedef typename Arithmetic_kernel::Integer Integer;
    
    typedef Integer Coefficient;
  
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


    typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_d_1> 
        Algebraic_kernel_d_2;

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
        Poly_2 f=from_string<Poly_2>("P[4(0,P[1(0,1)(1,1)])(4,P[0(0,-1)])]");
        Poly_2 g=from_string<Poly_2>("P[2(0,P[2(0,-5)(2,6)])(2,P[0(0,4)])]");
        Curve_analysis_2 ca1=construct_curve_2(f), 
            ca2=construct_curve_2(g);

        Curve_pair_analysis_2 curve_pair=construct_curve_pair_2(ca1,ca2);
        assert(curve_pair.number_of_status_lines_with_event()==7);
        typedef typename Curve_pair_analysis_2::Status_line_1 Status_line_1;

        typedef CGAL::internal::Event_indices<int> Triple;

        int i=0;

        Triple triple = curve_pair.event_indices(i);
        assert(triple.fg==-1);
        assert(triple.ffy==0);
        assert(triple.ggy==-1);

    }



}


int main() {


  test_routine<CGAL::CORE_arithmetic_kernel>();

    return 0;
}
