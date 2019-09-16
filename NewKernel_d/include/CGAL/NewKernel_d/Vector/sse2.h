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

#ifndef CGAL_VECTOR_SSE2_H
#define CGAL_VECTOR_SSE2_H

// Check what needs adapting for clang, intel and microsoft
#if !defined __SSE2__ || (__GNUC__ * 100 + __GNUC_MINOR__ < 408)
#error Requires SSE2 and gcc 4.8+
#endif
#include <x86intrin.h> // FIXME: other platforms call it differently

#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/enum.h> // CGAL::Sign
#include <CGAL/number_utils.h> // CGAL::sign



namespace CGAL {

  struct Sse_vector_2 {
    typedef double NT;
    typedef Dimension_tag<2> Dimension;
    typedef Dimension_tag<2> Max_dimension;
    // No Rebind_dimension, this is a building block
    template<class,bool=true> struct Property : boost::false_type {};
    template<bool b> struct Property<Has_vector_plus_minus_tag,b>
      : boost::true_type {};
    /* MAYBE?
       template<bool b> struct Property<Has_vector_scalar_ops_tag,b>
       : boost::true_type {};
       */
    template<bool b> struct Property<Has_determinant_of_vectors_tag,b>
      : boost::true_type {};
    template<bool b> struct Property<Has_dot_product_tag,b>
      : boost::true_type {};

    typedef __m128d Vector;
    struct Construct_vector {
      struct Dimension {
	// Initialize with NaN?
	Vector operator()(unsigned d) const {
	  CGAL_assertion(d==2);
	  return Vector();
	}
      };

      struct Iterator {
	template<typename Iter>
	  Vector operator()(unsigned d,Iter const& f,Iter const& e) const {
	    CGAL_assertion(d==2);
	    double x0 = *f;
	    double x1 = *++f;
	    CGAL_assertion(++f==e);
	    Vector a = { x0, x1 };
	    return a;
	  }
      };

      struct Iterator_and_last {
	template<typename Iter,typename T>
	  Vector operator()(unsigned d,Iter const& f,Iter const& e,double t) const {
	    CGAL_assertion(d==2);
	    Vector a = { *f, t };
	    CGAL_assertion(++f==e);
	    return a;
	  }
      };

      struct Values {
	  Vector operator()(double a,double b) const {
	    Vector r = { a, b };
	    return r;
	  }
      };

      struct Values_divide {
	Vector operator()(double h,double a,double b) const {
	  // {a,b}/{h,h} is probably slower
	  Vector r = { a/h, b/h };
	  return r;
	}
      };
    };

    typedef double const* Vector_const_iterator;
    static inline Vector_const_iterator vector_begin(Vector const&a){
      return (Vector_const_iterator)(&a);
    }
    static inline Vector_const_iterator vector_end(Vector const&a){
      return (Vector_const_iterator)(&a)+2;
    }
    static inline unsigned size_of_vector(Vector){
      return 2;
    }
    public:

    static double determinant_of_vectors(Vector a, Vector b) {
      __m128d c = _mm_shuffle_pd (b, b, 1); // b1, b0
      __m128d d = a * c; // a0*b1, a1*b0
#ifdef __SSE3__
      __m128d e = _mm_hsub_pd (d, d);
      return e[0];
#else
      return d[0]-d[1];
#endif
    }
    static CGAL::Sign sign_of_determinant_of_vectors(Vector a, Vector b) {
      return CGAL::sign(determinant_of_vectors(a,b));
    }

    static double dot_product(Vector a,Vector b){
#ifdef __SSE4_1__
      return _mm_dp_pd (a, b, 1+16+32)[0];
#else
      __m128d p = a * b;
#if defined __SSE3__
      __m128d s = _mm_hadd_pd (p, p);
      return s[0];
#else
      return p[0]+p[1];
#endif
#endif
    };
  };

}
#endif
