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

#ifndef CGAL_VECTOR_AVX4_H
#define CGAL_VECTOR_AVX4_H

#if !defined __AVX__ || (__GNUC__ * 100 + __GNUC_MINOR__ < 408)
#error Requires AVX and gcc 4.8+
#endif
#include <x86intrin.h>

#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/enum.h> // CGAL::Sign
#include <CGAL/number_utils.h> // CGAL::sign



namespace CGAL {

  struct Avx_vector_4 {
    typedef double NT;
    typedef Dimension_tag<4> Dimension;
    typedef Dimension_tag<4> Max_dimension;
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
    template<bool b> struct Property<Has_determinant_of_vectors_omit_last_tag,b>
      : boost::true_type {};

    typedef __m256d Vector;
    struct Construct_vector {
      struct Dimension {
	// Initialize with NaN?
	Vector operator()(unsigned d) const {
	  CGAL_assertion(d==4);
	  return Vector();
	}
      };

      struct Iterator {
	template<typename Iter>
	  Vector operator()(unsigned d,Iter const& f,Iter const& e) const {
	    CGAL_assertion(d==4);
	    double x0 = *f;
	    double x1 = *++f;
	    double x2 = *++f;
	    double x3 = *++f;
	    CGAL_assertion(++f==e);
	    Vector a = { x0, x1, x2, x3 };
	    return a;
	  }
      };

      struct Iterator_and_last {
	template<typename Iter,typename T>
	  Vector operator()(unsigned d,Iter const& f,Iter const& e,double t) const {
	    CGAL_assertion(d==4);
	    double x0 = *f;
	    double x1 = *++f;
	    double x2 = *++f;
	    CGAL_assertion(++f==e);
	    Vector a = { x0, x1, x2, t };
	    return a;
	  }
      };

      struct Values {
	  Vector operator()(double a,double b,double c,double d) const {
	    Vector r = { a, b, c, d };
	    return r;
	  }
      };

      struct Values_divide {
	Vector operator()(double h,double a,double b,double c,double d) const {
	  // {a,b,c,d}/{h,h,h,h} should be roughly the same
	  Vector r = { a/h, b/h, c/h, d/h };
	  return r;
	}
      };
    };

    public:
    typedef double const* Vector_const_iterator;
    static inline Vector_const_iterator vector_begin(Vector const&a){
      return (Vector_const_iterator)(&a);
    }
    static inline Vector_const_iterator vector_end(Vector const&a){
      return (Vector_const_iterator)(&a)+4;
    }
    static inline unsigned size_of_vector(Vector){
      return 4;
    }
    static inline double dot_product(__m256d x, __m256d y){
      __m256d p=x*y;
      __m256d z=_mm256_hadd_pd(p,p);
      return z[0]+z[2];
    }
    private:
    static inline __m256d avx_sym(__m256d x){
#if 0
      return __builtin_shuffle(x,(__m256i){2,3,0,1});
#else
      return _mm256_permute2f128_pd(x,x,1);
#endif
    }
    static inline __m256d avx_left(__m256d x){
#if 0
      return __builtin_shuffle(x,(__m256i){1,2,3,0});
#else
#ifdef __AVX2__
      return _mm256_permute4x64_pd(x,1+2*4+3*16+0*64);
#else
      __m256d s = _mm256_permute2f128_pd(x,x,1);
      return _mm256_shuffle_pd(x,s,5);
#endif
#endif
    }
    static inline __m256d avx_right(__m256d x){
#if 0
      return __builtin_shuffle(x,(__m256i){3,0,1,2});
#else
#ifdef __AVX2__
      return _mm256_permute4x64_pd(x,3+0*4+1*16+2*64);
#else
      __m256d s = _mm256_permute2f128_pd(x,x,1);
      return _mm256_shuffle_pd(s,x,5);
#endif
#endif
    }
    static inline double avx_altprod(__m256d x, __m256d y){
      __m256d p=x*y;
      __m256d z=_mm256_hsub_pd(p,p);
      return z[0]+z[2];
    }
    public:
    static double
      determinant_of_vectors(Vector a, Vector b, Vector c, Vector d) {
      __m256d x=a*avx_left(b)-avx_left(a)*b;
      __m256d yy=a*avx_sym(b);
      __m256d y=yy-avx_sym(yy);
      __m256d z0=x*avx_sym(c);
      __m256d z1=avx_left(x)*c;
      __m256d z2=y*avx_left(c);
      __m256d z=z0+z1-z2;
      return avx_altprod(z,avx_right(d));
    }
    static CGAL::Sign
      sign_of_determinant_of_vectors(Vector a, Vector b, Vector c, Vector d) {
      return CGAL::sign(determinant_of_vectors(a,b,c,d));
    }

    private:
    static inline __m256d avx3_right(__m256d x){
#if 0
      return __builtin_shuffle(x,(__m256i){2,0,1,3}); // can replace 3 with anything
#else
#ifdef __AVX2__
      return _mm256_permute4x64_pd(x,2+0*4+1*16+3*64);
#else
      __m256d s = _mm256_permute2f128_pd(x,x,1);
      return _mm256_shuffle_pd(s,x,12);
#endif
#endif
    }
    public:
    static inline double dot_product_omit_last(__m256d x, __m256d y){
      __m256d p=x*y;
      __m128d q=_mm256_extractf128_pd(p,0);
      double z=_mm_hadd_pd(q,q)[0];
      return z+p[2];
    }
    // Note: without AVX2, is it faster than the scalar computation?
    static double
      determinant_of_vectors_omit_last(Vector a, Vector b, Vector c) {
      __m256d x=a*avx3_right(b)-avx3_right(a)*b;
      return dot_product_omit_last(c,avx3_right(x));
    }
    static CGAL::Sign
      sign_of_determinant_of_vectors_omit_last(Vector a, Vector b, Vector c) {
      return CGAL::sign(determinant_of_vectors_omit_last(a,b,c));
    }

  };

}
#endif
