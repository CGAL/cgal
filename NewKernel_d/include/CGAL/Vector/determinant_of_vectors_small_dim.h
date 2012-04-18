#ifndef CGAL_VECTOR_DETVEC_SMALL_H
#define CGAL_VECTOR_DETVEC_SMALL_H
#include <CGAL/functor_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/determinant.h>
#include <CGAL/predicates/sign_of_determinant.h>

namespace CGAL {

template <class LA, class Dim_=LA::Dimension, class Max_dim_=LA::Max_dimension,
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
  template< class D2, class D3=D2 >
    struct Rebind_dimension {
      typedef typename LA::template Rebind_dimension<D2,D3> LA2;
      typedef Add_determinant_of_vectors_small_dim<LA2> Other;
    };
  template<class P,class=void> struct Property : LA::template Property<P> {};
  template<class D> struct Property<Has_determinant_of_vectors_tag, D> :
    boost::true_type {};

  template <class V1, class V2>
  static NT determinant_of_vectors(V1 const&a, V2 const&b){
    return determinant(a[0],a[1],b[0],b[1]);
  }
  template <class V1, class V2>
  static Sign sign_of_determinant_of_vectors(V1 const&a, V2 const&b){
    return sign_of_determinant(a[0],a[1],b[0],b[1]);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_vectors_small_dim
<LA, Dimension_tag<3>, Max_dim_, false> : LA {
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
    return determinant(a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2]);
  }
  static Sign sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c){
    return sign_of_determinant(a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2]);
  }
};

template <class LA, class Max_dim_>
struct Add_determinant_of_vectors_small_dim
<LA, Dimension_tag<4>, Max_dim_, false> : LA {
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
    return determinant(
	a[0],a[1],a[2],a[3],
	b[0],b[1],b[2],b[3],
	c[0],c[1],c[2],c[3],
	d[0],d[1],d[2],d[3]);
  }
  static Sign sign_of_determinant_of_vectors(Vector const&a, Vector const&b,
      Vector const&c, Vector const&d){
    return sign_of_determinant(
	a[0],a[1],a[2],a[3],
	b[0],b[1],b[2],b[3],
	c[0],c[1],c[2],c[3],
	d[0],d[1],d[2],d[3]);
  }
};

//TODO: Go up to 6. First check that it won't be done differently eventually.

}
#endif
