// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ACK_MACROS_H
#define CGAL_ACK_MACROS_H 1

/*!\file include/CGAL/Algebraic_kernel_d/macros.d
 * \brief Macro definitions wrt algebraic kernels
 */

#include <CGAL/config.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Sqrt_extension.h>

CGAL_BEGIN_NAMESPACE

#define CGAL_ACK_SNAP_ALGEBRAIC_CURVE_KERNEL_2_TYPEDEFS \
    typedef typename Algebraic_kernel_2::Algebraic_kernel_1 \
        Algebraic_kernel_1;\
    typedef typename Algebraic_kernel_1::Coefficient Coefficient;\
    typedef typename Algebraic_kernel_1::Boundary Boundary; \
    typedef typename CGAL::Get_arithmetic_kernel<Boundary> \
        ::Arithmetic_kernel Arithmetic_kernel; \
    typedef typename Arithmetic_kernel::Integer Integer; \
    typedef typename Algebraic_kernel_1::Algebraic_real_1 Algebraic_real_1;   \
    typedef typename Algebraic_kernel_2::Polynomial_1 Polynomial_1; \
    typedef typename Algebraic_kernel_2::Polynomial_2 Polynomial_2;      \
    typedef typename Algebraic_kernel_1::Solve_1 Solve_1;       \
    typedef CGAL::CGALi::Bitstream_descartes_rndl_tree_traits \
        < CGAL::CGALi::Bitstream_coefficient_kernel_at_alpha \
              < typename Polynomial_traits_d<Polynomial_2>::Coefficient_type, \
                Algebraic_real_1 \
              > \
        > \
        Bitstream_traits; \
    typedef CGAL::CGALi::Bitstream_descartes<Bitstream_traits>     \
        Bitstream_descartes;                                           \
    typedef CGAL::CGALi::Status_line_CA_1< Curve_analysis_2 > Status_line_1; \



#define CGAL_SNAP_AK_3_TYPEDEFS(Arithmetic_kernel)    \
    CGAL_SNAP_ARITHMETIC_KERNEL_TYPEDEFS(Arithmetic_kernel) \
    typedef CGAL::Polynomial< Integer > Poly_int1; \
    typedef CGAL::Polynomial< Poly_int1 > Poly_int2; \
    typedef CGAL::Polynomial< Poly_int2 > Poly_int3; \
    typedef CGAL::Polynomial< Rational > Poly_rat1; \
    typedef CGAL::Polynomial< Poly_rat1 > Poly_rat2; \
    typedef CGAL::Polynomial< Poly_rat2 > Poly_rat3; \
    typedef CGAL::Sqrt_extension< Rational, Integer > Extn; \
    typedef CGAL::Sqrt_extension< Extn, Extn > Nested_extn; \
    typedef CGAL::Polynomial< Extn > Poly_extn1; \
    typedef CGAL::Polynomial< Poly_extn1 > Poly_extn2; \
    typedef CGAL::Polynomial< Poly_extn2 > Poly_extn3; \
    typedef CGAL::Polynomial< Nested_extn > Poly_nested_extn1; \
    typedef CGAL::Polynomial< Poly_nested_extn1 > Poly_nested_extn2; \
    typedef CGAL::Polynomial< Poly_nested_extn2 > Poly_nested_extn3; \
// end #define CGAL_SNAP_AK_3_TYPEDEFS(AT)


CGAL_END_NAMESPACE

#endif //CGAL_ACK_MACROS_H
// EOF
