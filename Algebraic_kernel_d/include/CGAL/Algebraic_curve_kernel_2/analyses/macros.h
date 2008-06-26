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


#ifndef CGAL_ACK_MACROS_H
#define CGAL_ACK_MACROS_H 1

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

#define CGAL_ACK_SNAP_ALGEBRAIC_CURVE_TYPEDEFS \
    typedef typename Algebraic_kernel_2::Algebraic_kernel_1 \
        Algebraic_kernel_1;\
    typedef typename Algebraic_kernel_1::Coefficient Coefficient;\
    typedef typename Algebraic_kernel_1::Boundary Boundary; \
    typedef typename CGAL::Get_arithmetic_kernel<Boundary> \
        ::Arithmetic_kernel Arithmetic_kernel; \
    typedef typename Arithmetic_kernel::Integer Integer; \
    typedef typename Algebraic_kernel_1::Algebraic_real_1 X_coordinate_1;   \
    typedef X_coordinate_1  Y_coordinate_1; \
    typedef CGAL::Polynomial< Coefficient > Polynomial_1; \
    typedef CGAL::Polynomial< Polynomial_1 > Polynomial_2; \
    typedef typename Algebraic_kernel_1::Solve_1 Solve_1;       \
    typedef \
    CGAL::CGALi::Bitstream_descartes_traits_on_vert_line<Polynomial_1,   \
                                                    X_coordinate_1,      \
                                                    Integer >          \
    Bitstream_traits; \
    typedef CGAL::CGALi::Bitstream_descartes_bfs<Bitstream_traits> \
        Bitstream_descartes;                                           \
    typedef CGAL::CGALi::Status_line_CA_1< Curve_analysis_2 > Status_line_1; \


#define AcX_TIMERS \
    CGAL::Timer no_timer;

  // This is not really a macro, but it would be if I knew how to put #ifdef in a macro...
  template <typename OutputStream>
    void compiler_flag_info(OutputStream& out) {
    out << ">>>>>>>>>>>>>>> Compilation info >>>>>>>>>>>>>>>" << std::endl; 
    std::cout << "AcX_USE_CORE " << std::flush;
#if AcX_USE_CORE 
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 
    std::cout << "AcX_USE_LEDA " << std::flush;
#if AcX_USE_LEDA 
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 
    std::cout << "LiS_HAVE_NTL " << std::flush;
#ifdef LiS_HAVE_NTL
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 
    std::cout << "AcX_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS " << std::flush;
#if AcX_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 
    out << "NDEBUG " << std::flush;
#ifdef NDEBUG
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 
    out << "AcX_USE_MAPLE_FOR_MODULAR_RESULTANT " << std::flush;
#if AcX_USE_MAPLE_FOR_MODULAR_RESULTANT
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 
    out << "AcX_SPEED_UP_FOR_REGULAR_CURVES  " << std::flush;
#if AcX_SPEED_UP_FOR_REGULAR_CURVES
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 
#ifdef AcX_SPEED_UP_FOR_DEGREE_GREATER_EQUAL
    out << "AcX_SPEED_UP_FOR_DEGREE_GREATER_EQUAL " << std::flush;
    out << AcX_SPEED_UP_FOR_DEGREE_GREATER_EQUAL << std::endl;
#else
    out << "AcX_SPEED_UP_FOR_DEGREE_GREATER_EQUAL  " << std::flush;
    out << "undefined" << std::endl;
#endif
    out << "AcX_STATIC_SEED  " << std::flush;
#if AcX_STATIC_SEED
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 
    out << "AcX_CHECK_POLYNOMIALS_FOR_COPRIMABILITY  " << std::flush;
#if AcX_CHECK_POLYNOMIALS_FOR_COPRIMABILITY
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 
    out << "AcX_NO_ARC_FLIP " << std::flush;
#if AcX_NO_ARC_FLIP
    out << "yes" << std::endl;
#else 
    out << "no" << std::endl;
#endif 


    out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl; 
  }
  
CGAL_END_NAMESPACE

#endif //CGAL_ACK_MACROS_H
