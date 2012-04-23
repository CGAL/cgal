#ifndef CGAL_VECTOR_2INT_H
#define CGAL_VECTOR_2INT_H

#include <stdint.h>
#include <CGAL/array.h>
#include <CGAL/Dimension.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/NT_converter.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/determinant_of_vectors.h>



namespace CGAL {

  struct Vector_2_int {
#if 1
    typedef double NT; // try lying a bit
    typedef int32_t NT1; // what is really stored
    // (sign_of_)determinant_of_vectors needs adapting for unsigned NT1
    //typedef unsigned int NTu1;
    typedef int_least64_t NT2; // longer type for computations
    //typedef uint_least64_t NTu2;
#else
    typedef long double NT;
    typedef int64_t NT1;
    //typedef uint64_t NTu1;
    typedef __int128 NT2;
    //typedef unsigned __int128 NTu2;
#endif

    typedef Dimension_tag<2> Dimension;
    typedef Dimension_tag<2> Max_dimension;
    // No Rebind_dimension, this is a building block
    template<class,bool=true> struct Property : boost::false_type {};
    //template<bool b> struct Property<Has_vector_plus_minus_tag,b>
    //  : boost::true_type {};
    template<bool b> struct Property<Has_determinant_of_vectors_tag,b>
      : boost::true_type {};
    //template<bool b> struct Property<Has_determinant_of_points_tag,b>
    //  : boost::true_type {};
    // Advertise somehow that the sign_of_determinant* are exact?

    typedef cpp0x::array<NT1,2> Vector;
    struct Construct_vector {
      struct Dimension {
	Vector operator()(unsigned d) const {
	  CGAL_assertion(d==2);
	  return Vector();
	}
      };

      // TODO (for all constructors): check that input fits in NT1...
      struct Iterator {
	template<typename Iter>
	  Vector operator()(unsigned d,Iter const& f,Iter const& e) const {
	    CGAL_assertion(d==2);
	    NT1 x0 = *f;
	    NT1 x1 = *++f;
	    CGAL_assertion(++f==e);
	    Vector a = { x0, x1 };
	    return a;
	  }
      };

      struct Iterator_and_last {
	template<typename Iter,typename T>
	  Vector operator()(unsigned d,Iter const& f,Iter const& e,double t) const {
	    CGAL_assertion(d==2);
	    Vector a = { static_cast<NT1>(*f), t };
	    CGAL_assertion(++f==e);
	    return a;
	  }
      };

      struct Values {
	  Vector operator()(NT1 a,NT1 b) const {
	    Vector r = { a, b };
	    return r;
	  }
      };

      /*
	 // Maybe safer not to provide it
      struct Values_divide {
	Vector operator()(double h,double a,double b) const {
	  Vector r = { a/h, b/h };
	  return r;
	}
      };
      */
    };

    // Since we lie about NT, be consistent about it
    typedef transforming_iterator<NT_converter<NT1,NT>,NT1 const*> Vector_const_iterator;
    static inline Vector_const_iterator vector_begin(Vector const&a){
      return Vector_const_iterator(a.begin());
    }
    static inline Vector_const_iterator vector_end(Vector const&a){
      return Vector_const_iterator(a.end());
    }
    static inline unsigned size_of_vector(Vector){
      return 2;
    }

    // for unsigned NT1, check what changes to do.
    // return NT instead?
    static NT2 determinant_of_vectors(Vector a, Vector b) {
      return CGAL::determinant_of_vectors<NT2>(a,b);
    }
    static CGAL::Sign sign_of_determinant_of_vectors(Vector a, Vector b) {
      return CGAL::sign_of_determinant_of_vectors<NT2>(a,b);
    }

#if 0
    // WARNING: FIXME: Completely broken as is
    // TODO: put an assertion at construction that abs(coord)<INTMAX/sqrt(2) or something
    static NT determinant_of_points(Vector a, Vector b, Vector c) {
      NT2 a0=a[0];    NT2 a1=a[1];
      NT2 x0=b[0]-a0; NT2 x1=b[1]-a1;
      NT2 y0=c[0]-a0; NT2 y1=c[1]-a1;
      return CGAL::determinant<NT>(x0,x1,y0,y1);
    }
    static CGAL::Sign sign_of_determinant_of_points(Vector a, Vector b, Vector c) {
      NT2 a0=a[0];    NT2 a1=a[1];
      NT2 x0=b[0]-a0; NT2 x1=b[1]-a1;
      NT2 y0=c[0]-a0; NT2 y1=c[1]-a1;
      return CGAL::compare(x0*y1,x1*y0);
    }
#endif
  };

}
#endif
