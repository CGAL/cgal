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
//
// Author(s)     : Marc Glisse

#ifndef CGAL_VECTOR_DET_ITER_PTS_PTS_H
#define CGAL_VECTOR_DET_ITER_PTS_PTS_H
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class LA, class Dim_=typename LA::Dimension,
	 class Max_dim_=typename LA::Max_dimension,
	 bool = LA::template Property<Has_determinant_of_iterator_to_points_tag>::value,
	 bool = LA::template Property<Has_determinant_of_points_tag>::value>
struct Add_determinant_of_iterator_to_points_from_points : LA {
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_iterator_to_points_from_points<LA2> Other;
    };
};

//FIXME: Use variadics and boost so it works in any dimension.
template <class LA, class Max_dim_>
struct Add_determinant_of_iterator_to_points_from_points
<LA, Dimension_tag<2>, Max_dim_, false, true> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_iterator_to_points_from_points<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_iterator_to_points_tag, D> :
    boost::true_type {};

  template<class Iter>
  static NT determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; CGAL_assertion(++first==end);
    return LA::determinant_of_points(a,b,c);
  }
  template<class Iter>
  static Sign sign_of_determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; CGAL_assertion(++first==end);
    return LA::sign_of_determinant_of_points(a,b,c);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_iterator_to_points_from_points
<LA, Dimension_tag<3>, Max_dim_, false, true> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_iterator_to_points_from_points<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_iterator_to_points_tag, D> :
    boost::true_type {};

  template<class Iter>
  static NT determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; ++first;
    Vector const&d=*first; CGAL_assertion(++first==end);
    return LA::determinant_of_points(a,b,c,d);
  }
  template<class Iter>
  static Sign sign_of_determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; ++first;
    Vector const&d=*first; CGAL_assertion(++first==end);
    return LA::sign_of_determinant_of_points(a,b,c,d);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_iterator_to_points_from_points
<LA, Dimension_tag<4>, Max_dim_, false, true> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_iterator_to_points_from_points<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_iterator_to_points_tag, D> :
    boost::true_type {};

  template<class Iter>
  static NT determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; ++first;
    Vector const&d=*first; ++first;
    Vector const&e=*first; CGAL_assertion(++first==end);
    return LA::determinant_of_points(a,b,c,d,e);
  }
  template<class Iter>
  static Sign sign_of_determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; ++first;
    Vector const&d=*first; ++first;
    Vector const&e=*first; CGAL_assertion(++first==end);
    return LA::sign_of_determinant_of_points(a,b,c,d,e);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_iterator_to_points_from_points
<LA, Dimension_tag<5>, Max_dim_, false, true> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_iterator_to_points_from_points<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_iterator_to_points_tag, D> :
    boost::true_type {};

  template<class Iter>
  static NT determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; ++first;
    Vector const&d=*first; ++first;
    Vector const&e=*first; ++first;
    Vector const&f=*first; CGAL_assertion(++first==end);
    return LA::determinant_of_points(a,b,c,d,e,f);
  }
  template<class Iter>
  static Sign sign_of_determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; ++first;
    Vector const&d=*first; ++first;
    Vector const&e=*first; ++first;
    Vector const&f=*first; CGAL_assertion(++first==end);
    return LA::sign_of_determinant_of_points(a,b,c,d,e,f);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_iterator_to_points_from_points
<LA, Dimension_tag<6>, Max_dim_, false, true> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_iterator_to_points_from_points<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_iterator_to_points_tag, D> :
    boost::true_type {};

  template<class Iter>
  static NT determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; ++first;
    Vector const&d=*first; ++first;
    Vector const&e=*first; ++first;
    Vector const&f=*first; ++first;
    Vector const&g=*first; CGAL_assertion(++first==end);
    return LA::determinant_of_points(a,b,c,d,e,f,g);
  }
  template<class Iter>
  static Sign sign_of_determinant_of_iterator_to_points(Iter const&first, Iter const&end){
    Vector const&a=*first; ++first;
    Vector const&b=*first; ++first;
    Vector const&c=*first; ++first;
    Vector const&d=*first; ++first;
    Vector const&e=*first; ++first;
    Vector const&f=*first; ++first;
    Vector const&g=*first; CGAL_assertion(++first==end);
    return LA::sign_of_determinant_of_points(a,b,c,d,e,f,g);
  }
};

}
#endif
