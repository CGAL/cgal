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

#ifndef CGAL_ALLOWED_INCLUSION
#error Must not include this header directly
#endif
#if !defined(CGAL_TAG) \
  || ! defined(CGAL_CLASS) \
  || ! defined(CGAL_FUNC) \
  || ! defined(CGAL_SIGN_FUNC) \
  || ! defined(CGAL_SHIFT)

#error Forgot one macro
#endif

namespace CGAL {

template <class LA, class Dim_=typename LA::Dimension,
	 class Max_dim_=typename LA::Max_dimension,
	 bool=LA::template Property<CGAL_TAG>::value>
struct CGAL_CLASS : LA {
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef CGAL_CLASS<LA2> Other;
    };
};

template <class LA, class Max_dim_>
struct CGAL_CLASS
<LA, Dimension_tag<2+CGAL_SHIFT>, Max_dim_, false> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef CGAL_CLASS<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<CGAL_TAG, D> :
    boost::true_type {};

  static NT CGAL_FUNC(Vector const&a, Vector const&b){
    return CGAL::determinant_of_vectors<NT>(a,b);
  }
  template <class V1, class V2>
  static Sign CGAL_SIGN_FUNC(V1 const&a, V2 const&b){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b);
  }
};

template <class LA, class Max_dim_>
struct CGAL_CLASS
<LA, Dimension_tag<3+CGAL_SHIFT>, Max_dim_, false> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef CGAL_CLASS<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<CGAL_TAG, D> :
    boost::true_type {};

  static NT CGAL_FUNC(Vector const&a, Vector const&b,
      Vector const&c){
    return CGAL::determinant_of_vectors<NT>(a,b,c);
  }
  static Sign CGAL_SIGN_FUNC(Vector const&a, Vector const&b,
      Vector const&c){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b,c);
  }
};

template <class LA, class Max_dim_>
struct CGAL_CLASS
<LA, Dimension_tag<4+CGAL_SHIFT>, Max_dim_, false> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef CGAL_CLASS<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<CGAL_TAG, D> :
    boost::true_type {};

  static NT CGAL_FUNC(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d){
    return CGAL::determinant_of_vectors<NT>(a,b,c,d);
  }
  static Sign CGAL_SIGN_FUNC(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b,c,d);
  }
};

template <class LA, class Max_dim_>
struct CGAL_CLASS
<LA, Dimension_tag<5+CGAL_SHIFT>, Max_dim_, false> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef CGAL_CLASS<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<CGAL_TAG, D> :
    boost::true_type {};

  static NT CGAL_FUNC(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e){
    return CGAL::determinant_of_vectors<NT>(a,b,c,d,e);
  }
  static Sign CGAL_SIGN_FUNC(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b,c,d,e);
  }
};

template <class LA, class Max_dim_>
struct CGAL_CLASS
<LA, Dimension_tag<6+CGAL_SHIFT>, Max_dim_, false> : LA {
  typedef typename LA::NT NT;
  typedef typename LA::Vector Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef CGAL_CLASS<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<CGAL_TAG, D> :
    boost::true_type {};

  static NT CGAL_FUNC(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e, Vector const&f){
    return CGAL::determinant_of_vectors<NT>(a,b,c,d,e,f);
  }
  static Sign CGAL_SIGN_FUNC(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e, Vector const&f){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b,c,d,e,f);
  }
};

}
