// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@informatik.uni-mainz.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


#ifndef CGAL_POLYNOMIAL_REAL_EMBEDDABLE_TRAITS_H
#define CGAL_POLYNOMIAL_REAL_EMBEDDABLE_TRAITS_H

#include <CGAL/basic.h>

namespace CGAL {

namespace internal {
template< class Polynomial , class TAG> class Real_embeddable_traits_poly_base;

template< class NT , class TAG> class Real_embeddable_traits_poly_base< Polynomial<NT>, TAG >
  : public INTERN_RET::Real_embeddable_traits_base< Polynomial<NT> , CGAL::Tag_false > {};

// Real embeddable traits
// TODO: Polynomials aren't Real_embeddable! But for debugging and testing
//       reasons, the real embeddable functors are provided.
template< class NT > class Real_embeddable_traits_poly_base< Polynomial<NT>, CGAL::Tag_true >
  : public INTERN_RET::Real_embeddable_traits_base< Polynomial<NT> , CGAL::Tag_false > {
public:

    typedef Tag_false Is_real_embeddable;

    class Abs {
    public:
        typedef Polynomial<NT> argument_type;
        typedef Polynomial<NT> result_type;
        Polynomial<NT> operator()( const Polynomial<NT>& x ) const {
            return x.abs();
        }
    };

    class Sgn {
    public:
        typedef Polynomial<NT>              argument_type;
        typedef CGAL::Sign        result_type;
        CGAL::Sign operator()( const Polynomial<NT>& x ) const {
            return x.sign();
        }
    };

    class Compare {
    public:
        typedef Polynomial<NT>                    first_argument_type;
        typedef Polynomial<NT>                    second_argument_type;
        typedef CGAL::Comparison_result result_type;

        CGAL::Comparison_result operator()(
                const Polynomial<NT>& x,
                const Polynomial<NT>& y ) const {
            return x.compare(y);
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Polynomial<NT>,
                CGAL::Comparison_result )
    };

    class To_double {
    public:
        typedef typename Real_embeddable_traits<NT>::To_double NT_to_double;
        typedef Polynomial<typename NT_to_double::result_type> result_type;
        typedef Polynomial<NT> argument_type;
        result_type operator()( const Polynomial<NT>& x ) const {
            CGAL_precondition(x.degree() >= 0);
            NT_to_double to_double;
            return result_type(
                    ::boost::make_transform_iterator(x.begin(),to_double),
                    ::boost::make_transform_iterator(x.end()  ,to_double));
        }
    };

    class To_interval {
    public:
        typedef typename Real_embeddable_traits<NT>::To_interval NT_to_interval;
        typedef Polynomial<typename NT_to_interval::result_type> result_type;
        typedef Polynomial<NT> argument_type;
        result_type operator()( const Polynomial<NT>& x ) const {
            CGAL_precondition( x.degree() >= 0 );
            NT_to_interval to_interval;
            return result_type(
                    ::boost::make_transform_iterator(x.begin(),to_interval),
                    ::boost::make_transform_iterator(x.end()  ,to_interval));
        }
    };
};
} // namespace internal

template <typename NT>
class Real_embeddable_traits<Polynomial<NT> >
  :public internal::Real_embeddable_traits_poly_base<
  Polynomial<NT>,
  typename Real_embeddable_traits<typename internal::Innermost_coefficient_type<NT>::Type>::Is_real_embeddable>
{};

} //namespace CGAL
#endif // CGAL_POLYNOMIAL_REAL_EMBEDDABLE_TRAITS_H
