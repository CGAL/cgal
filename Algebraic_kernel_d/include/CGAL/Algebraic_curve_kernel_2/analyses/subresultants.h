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
#if 0

#ifndef CGAL_ACK_SUBRESULTANTS
#define CGAL_ACK_SUBRESULTANTS 1

#include<vector>
#include<list>

#include <CGAL/Polynomial.h>
#include <CGAL/Algebraic_curve_kernel_2/analyses/subresultants.h>

#if AcX_USE_MAPLE_FOR_MODULAR_RESULTANT
#if !AcX_MAPLE_USE_WORKSHEET

#include <NiX/Maple_session.h>

#endif
#endif

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  /*!
   * \brief Computes the resultant of \c A and \c B
   *
   * The default method is to use Polynomial remainder sequence method. With 
   * the compiler flag AcX_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS, one can 
   * compute the resultant with the Bezout matrix instead.
   *
   * A modular algorithm for resultant computation would certainly improve the
   * performance.
   */
  template<typename NT>
    NT resultant(CGAL::Polynomial<NT> A,
		 CGAL::Polynomial<NT> B) {
#if AcX_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS
    return NiX::hybrid_bezout_subresultant(A,B);
#else
#if AcX_USE_MAPLE_FOR_MODULAR_RESULTANT

    return NiX::maple_res2(A,B);

#else
    if(B.degree()==0) {
      if(B.is_zero()) {
	return NT(0);
      }
      else {
	return CGAL::ipower(B[0],A.degree());
      }
    }
    std::vector<NT> sres;
    NiX::prs_principal_subresultants(A,B,std::back_inserter(sres));
    return sres[0];
#endif
#endif     
  }

}
      
CGAL_END_NAMESPACE

#endif

#endif // if 0
