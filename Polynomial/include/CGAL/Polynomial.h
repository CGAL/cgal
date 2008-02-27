// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Seel <seel@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file NiX/Polynomial.h
 *  \brief Defines class NiX::Polynomial.
 *  
 *  Polynomials in one variable (or more, by recursion)
 */

#ifndef CGAL_POLYNOMIAL_H
#define CGAL_POLYNOMIAL_H

#include <cstdarg>
#include <cctype>
#include <vector>
#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/mpl/if.hpp>
#include <CGAL/Flattening_iterator.h>
//#include <NiX/Modular.h>

#include <CGAL/Exponent_vector.h>

#include <boost/static_assert.hpp>

#ifdef CGAL_USE_LEDA
#if CGAL_LEDA_VERSION >= 500
#include <LEDA/core/array.h>
#else
#include <LEDA/array.h>
#endif
#endif // CGAL_USE_LEDA

#include <CGAL/Polynomial/Polynomial_type.h>
#include <CGAL/Polynomial/Algebraic_structure_traits.h>
#include <CGAL/Polynomial/Real_embeddable_traits.h>
#include <CGAL/Polynomial/Fraction_traits.h>
#include <CGAL/Polynomial/Scalar_factor_traits.h>
#include <CGAL/Polynomial/Modular_traits.h>
#include <CGAL/Polynomial/Coercion_traits.h>

// TODO: Are these still includes necessary?
#include <CGAL/Polynomial/polynomial_gcd.h> // used above for NT_traits<Poly...>::Gcd
#include <CGAL/Polynomial/prs_resultant.h>  // for compatibility

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial/polynomial_utils.h>



CGAL_BEGIN_NAMESPACE


// Internally, Polynomials also need this:

//! used internally for data exchanged between nesting levels of polynomials
// this traits-class provides
//   a) poly_nesting_depth: a counter for identifying a polynomials nesting level
//   b) Innermost_coefficient: a type which defines the numbertype of the
//        innermost nesting level, which is not a polynomial itself
//   c) Innermost_lcoeff: returns the leading coefficient of the polynomial in its
//       common sense
//   d) Innermost_coefficient_to_polynomial: transforms an innermost
//       numbertype into a polynomial like a typecast

 
// fwd Polynomial_traits_d
//template <typename Polynomial_d> class Polynomial_traits_d;




//
// Non-Member Functions
//



// The following former parts of Polynomial.h have been moved to the new files
// <NiX/polynomial_gcd.h> (items 2,3,5) and <NiX/prs_resultant.h> (item 4):
//
// 2) gcd (basic form without cofactors)
// 3) extended gcd computation (with cofactors)
// 4) resultant computation from polynomial remainder sequences (PRS)
// 5) square-free factorization


// } // namespace NiX

// Cofraction_traits added by Michael Hemmer
/*namespace NiX{
namespace Intern{
    template <class NT, class TAG> class Cofraction_traits_base;

    template <class NT_, class TAG>
    class Cofraction_traits_base<NiX::Polynomial<NT_>, TAG > {
        typedef NT_ NT;
    public:
        typedef Polynomial<NT_>                        Numerator_type;
        typedef ::LiS::False_tag                       Is_composable;
        typedef ::LiS::Null_tag                       Denominator_type;
        typedef ::LiS::Null_tag                       Type;
        typedef ::LiS::Null_tag                       Compose;
    };
    
    template <class NT_>
    class Cofraction_traits_base<NiX::Polynomial<NT_>, LiS::True_tag > {
        typedef NT_ NT;
        typedef NiX::Cofraction_traits<NT_> CFT_NT;
    public:
        typedef Polynomial<NT>                             Numerator_type;
        typedef ::LiS::True_tag                            Is_composable;
        typedef typename CFT_NT::Denominator_type          Denominator_type;
        typedef Polynomial<typename CFT_NT::Type> Type;
        
        class Compose {
        public:
            //! first argument type
            typedef Numerator_type   first_argument_type;
            //! second argument type
            typedef Denominator_type second_argument_type;
            //! result type
            typedef Type    result_type;
            //! Compose fraction
            Type operator() (Numerator_type num, 
                                      Denominator_type den 
                                      = Denominator_type(1)){
                Type tmp1; NiX::convert_to(num,tmp1);
                Type tmp2; NiX::convert_to(den,tmp2);
                return tmp1/tmp2;
            }
        };
    };
} //namespace Intern
*/
/*! \ingroup NiX_Cofraction_traits_specs
 *  \brief Specialization of Cofraction_traits for NiX::Polynomial<NT>.
 */
//template<class NT>
/*class Cofraction_traits<Polynomial<NT> > :
    public Intern::Cofraction_traits_base<
      Polynomial<NT>,
      typename Cofraction_traits<NT>::Is_composable>{
    //nothing new
};*/ 





CGAL_END_NAMESPACE


//
// trailing documentation
//

// Literature reference
//
// [Akritas, 1989]
//      Alkiviadis G. Akritas
//      Elements of Computer Algebra With Applications
//      Wiley, New York, 1989.
//
// [Cohen, 1993]
//      Cohen, Henri
//      A Course in Computational Algebraic Number Theory
//      Springer GTM 138, 1993
//
// [Cox et al, 1997]
//      David Cox; John Little; Donal O'Shea
//      Ideals, Varieties, and Algorithms
//      2nd ed., Springer UTM, 1997
//
// [Geddes et al, 1992]
//      Geddes, Keith O. and Czapor, Stephen R. and Labahn, George
//      Algorithms for Computer Algebra
//      Kluwer, 1992
//
// [Mignotte, 1992]
//      Mignotte, Maurice
//      Mathematics for Computer Algebra
//      Springer, 1992
//
// [PARI]
//      (a computer algebra system by Henri Cohen and collaborators)
//      http://www.parigp-home.de/

#endif  // NiX_POLYNOMIAL_H

// EOF
