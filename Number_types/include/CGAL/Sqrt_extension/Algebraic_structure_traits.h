// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_ALGEBRAIC_STRUCTURE_TRAITS_H
#define CGAL_SQRT_EXTENSION_ALGEBRAIC_STRUCTURE_TRAITS_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// Algebraic structure traits
template< class Type, class Algebraic_type >
class Sqrt_extension_algebraic_structure_traits_base;

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                        CGAL::Integral_domain_without_division_tag >
  : public Algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_without_division_tag > {
  public:
    typedef CGAL::Integral_domain_without_division_tag Algebraic_category;

    class Simplify
      : public Unary_function< Type&, void > {
      public:
        typedef void result_type;
        typedef Type& argument_type;

        void operator()( Type& x ) const {
          x.simplify();
        }
    };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Integral_domain_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_without_division_tag > {
  public:
    typedef CGAL::Integral_domain_tag Algebraic_category;

    class Integral_division
      : public Binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                                        const Type& y ) const {
          return x/y;
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Unique_factorization_domain_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Euclidean_ring_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Field_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  public:
    typedef Field_tag Algebraic_category;

    class Unit_part
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return( x == Type(0) ? Type(1) : x );
        }
    };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Field_with_sqrt_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Field_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                CGAL::Field_with_kth_root_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      // TODO: Why not Fiel_tag?
                                      CGAL::Field_with_sqrt_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                CGAL::Field_with_root_of_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      // TODO: Why not Fiel_tag?
                                      CGAL::Field_with_sqrt_tag > {
  // Nothing new
};

} // namespace CGALi


template< class COEFF, class ROOT>
class Algebraic_structure_traits< Sqrt_extension< COEFF, ROOT > >
    : public CGALi::Sqrt_extension_algebraic_structure_traits_base<
      Sqrt_extension< COEFF, ROOT >,
      typename Algebraic_structure_traits< COEFF >::Algebraic_category > {
public:
    typedef Sqrt_extension< COEFF, ROOT > Type;

    // Tag_true if COEFF and ROOT are exact
    typedef typename ::boost::mpl::if_c<
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<ROOT >::Is_exact,::CGAL::Tag_true>::value )&&
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<COEFF>::Is_exact,::CGAL::Tag_true>::value )
           ,::CGAL::Tag_true,::CGAL::Tag_false>::type Is_exact;

    typedef typename Algebraic_structure_traits<COEFF>::Is_numerical_sensitive
    Is_numerical_sensitive;
};

CGAL_END_NAMESPACE

#endif
