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

#define CGAL_ACK_SNAP_ALGEBRAIC_CURVE_KERNEL_2_TYPEDEFS \
    typedef typename Algebraic_kernel_2::Algebraic_kernel_1 \
        Algebraic_kernel_1;\
    typedef typename Algebraic_kernel_1::Coefficient Coefficient;\
    typedef typename Algebraic_kernel_1::Boundary Boundary; \
    typedef typename CGAL::Get_arithmetic_kernel<Boundary> \
        ::Arithmetic_kernel Arithmetic_kernel; \
    typedef typename Arithmetic_kernel::Integer Integer; \
    typedef typename Algebraic_kernel_1::Algebraic_real_1 Algebraic_real_1;   \
        typedef CGAL::Polynomial< Coefficient > Polynomial_1; \
    typedef CGAL::Polynomial< Polynomial_1 > Polynomial_2; \
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


CGAL_END_NAMESPACE

#endif //CGAL_ACK_MACROS_H
