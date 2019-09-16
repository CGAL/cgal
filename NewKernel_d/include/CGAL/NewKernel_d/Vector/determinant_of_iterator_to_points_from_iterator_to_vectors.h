// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Marc Glisse

#ifndef CGAL_VECTOR_DET_ITER_PTS_ITER_VEC_H
#define CGAL_VECTOR_DET_ITER_PTS_ITER_VEC_H
#include <functional>
#include <CGAL/transforming_iterator.h>
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class LA, class Dim_=typename LA::Dimension,
	 class Max_dim_=typename LA::Max_dimension,
	 bool = LA::template Property<Has_determinant_of_iterator_to_points_tag>::value,
	 bool = LA::template Property<Has_determinant_of_iterator_to_vectors_tag>::value>
struct Add_determinant_of_iterator_to_points_from_iterator_to_vectors : LA {
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_iterator_to_points_from_iterator_to_vectors<LA2> Other;
    };
};

template <class LA, class Dim_,class Max_dim_>
struct Add_determinant_of_iterator_to_points_from_iterator_to_vectors
<LA, Dim_, Max_dim_, false, true> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_iterator_to_points_from_iterator_to_vectors<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_iterator_to_points_tag, D> :
    boost::true_type {};

  // TODO: use std::minus, boost::bind, etc
  template<class T> struct Minus_fixed {
    T const& a;
    Minus_fixed(T const&a_):a(a_){}
    T operator()(T const&b)const{return b-a;}
  };
  template<class Iter>
  static NT determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Minus_fixed<Vector> f(a);
    return LA::determinant_of_iterator_to_vectors(make_transforming_iterator(first,f),make_transforming_iterator(end,f));
  }
  template<class Iter>
  static Sign sign_of_determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Minus_fixed<Vector> f(a);
    return LA::sign_of_determinant_of_iterator_to_vectors(make_transforming_iterator(first,f),make_transforming_iterator(end,f));
  }
};

}
#endif
