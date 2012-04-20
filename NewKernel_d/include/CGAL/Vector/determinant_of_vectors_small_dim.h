#ifndef CGAL_VECTOR_DETVEC_SMALL_H
#define CGAL_VECTOR_DETVEC_SMALL_H
#include <CGAL/functor_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/determinant_of_vectors.h>

namespace CGAL {

template <class LA, class Dim_=typename LA::Dimension,
	 class Max_dim_=typename LA::Max_dimension,
	 bool=LA::template Property<Has_determinant_of_vectors_tag>::value>
struct Add_determinant_of_vectors_small_dim : LA {
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_vectors_small_dim<LA2> Other;
    };
};

template <class LA, class Max_dim_>
struct Add_determinant_of_vectors_small_dim
<LA, Dimension_tag<2>, Max_dim_, false> : LA {
  using typename LA::NT;
  using typename LA::Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_vectors_small_dim<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_vectors_tag, D> :
    boost::true_type {};

  static NT determinant_of_vectors(Vector const&a, Vector const&b){
    return CGAL::determinant_of_vectors<NT>(a,b);
  }
  template <class V1, class V2>
  static Sign sign_of_determinant_of_vectors(V1 const&a, V2 const&b){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_vectors_small_dim
<LA, Dimension_tag<3>, Max_dim_, false> : LA {
  using typename LA::NT;
  using typename LA::Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_vectors_small_dim<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_vectors_tag, D> :
    boost::true_type {};

  static NT determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c){
    return CGAL::determinant_of_vectors<NT>(a,b,c);
  }
  static Sign sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b,c);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_vectors_small_dim
<LA, Dimension_tag<4>, Max_dim_, false> : LA {
  using typename LA::NT;
  using typename LA::Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_vectors_small_dim<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_vectors_tag, D> :
    boost::true_type {};

  static NT determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d){
    return CGAL::determinant_of_vectors<NT>(a,b,c,d);
  }
  static Sign sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b,c,d);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_vectors_small_dim
<LA, Dimension_tag<5>, Max_dim_, false> : LA {
  using typename LA::NT;
  using typename LA::Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_vectors_small_dim<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_vectors_tag, D> :
    boost::true_type {};

  static NT determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e){
    return CGAL::determinant_of_vectors<NT>(a,b,c,d,e);
  }
  static Sign sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b,c,d,e);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_vectors_small_dim
<LA, Dimension_tag<6>, Max_dim_, false> : LA {
  using typename LA::NT;
  using typename LA::Vector;
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_vectors_small_dim<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_vectors_tag, D> :
    boost::true_type {};

  static NT determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e, Vector const&f){
    return CGAL::determinant_of_vectors<NT>(a,b,c,d,e,f);
  }
  static Sign sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d, Vector const&e, Vector const&f){
    return CGAL::sign_of_determinant_of_vectors<NT>(a,b,c,d,e,f);
  }
};

}
#endif
