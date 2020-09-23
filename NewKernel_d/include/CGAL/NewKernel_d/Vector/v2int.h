// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_VECTOR_2INT_H
#define CGAL_VECTOR_2INT_H

#include <stdint.h>
#include <cmath>
#include <CGAL/array.h>
#include <CGAL/Dimension.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/NT_converter.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/determinant_of_vectors.h>
#include <CGAL/NewKernel_d/functor_tags.h>


// What are the pros and cons of having NT be int vs double?

namespace CGAL {
  struct Vector_2_int_prop1 {
    typedef double NT; // try lying a bit
    typedef int32_t NT1; // what is really stored
    typedef int32_t NT1b; // slightly longer
    typedef int_fast64_t NT2; // longer type for computations
    typedef int_fast64_t NT2b; // slightly longer
    bool check_limits(int32_t x){return std::abs(x)<(1<<30);}
    // TODO: find nice bounds
  };
#ifdef __SIZEOF_INT128__
  struct Vector_2_int_prop2 {
    typedef double NT;
    typedef int32_t NT1;
    typedef int_fast64_t NT1b;
    typedef int_fast64_t NT2;
    typedef __int128 NT2b;
    bool check_limits(int32_t){return true;}
    // take a template/int64_t input and still check the limits?
  };
  struct Vector_2_int_prop3 {
    typedef long double NT;
    typedef int64_t NT1;
    typedef int64_t NT1b;
    typedef __int128 NT2;
    typedef __int128 NT2b;
    enum { has_limit=true };
    bool check_limits(int32_t x){return std::abs(x)<(1L<<62);}
    // TODO: find nice bounds
  };
#endif

  template<class Prop=Vector_2_int_prop1>
  struct Vector_2_int : Prop {
    using typename Prop::NT;
    using typename Prop::NT1;
    using typename Prop::NT1b;
    using typename Prop::NT2;
    using typename Prop::NT2b;
    using Prop::check_limits;

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

    typedef std::array<NT1,2> Vector;
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
            CGAL_assertion (++f == e);
            CGAL_assertion (check_limits(x0) && check_limits(x1));
            Vector a = { x0, x1 };
            return a;
          }
      };

      struct Iterator_and_last {
        template<typename Iter,typename T>
          Vector operator()(unsigned d,Iter const& f,Iter const& e,double t) const {
            CGAL_assertion(d==2);
            NT1 x = *f;
            CGAL_assertion (++f == e);
            CGAL_assertion (check_limits(x) && check_limits(t));
            Vector a = { x, t };
            return a;
          }
      };

      struct Values {
          Vector operator()(NT1 a,NT1 b) const {
            CGAL_assertion (check_limits(a) && check_limits(b));
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
    // return NT or NT2?
    static NT determinant_of_vectors(Vector a, Vector b) {
      return CGAL::determinant_of_vectors<NT2>(a,b);
    }
    static CGAL::Sign sign_of_determinant_of_vectors(Vector a, Vector b) {
      return CGAL::sign_of_determinant_of_vectors<NT2>(a,b);
    }

    static NT determinant_of_points(Vector a, Vector b, Vector c) {
      // could be faster to convert to NT directly
      NT1b a0=a[0];    NT1b a1=a[1];
      NT1b x0=b[0]-a0; NT1b x1=b[1]-a1;
      NT1b y0=c[0]-a0; NT1b y1=c[1]-a1;
      return CGAL::determinant<NT>(x0,x1,y0,y1);
    }
    static CGAL::Sign sign_of_determinant_of_points(Vector a, Vector b, Vector c) {
      NT1b a0=a[0];    NT1b a1=a[1];
      NT1b x0=b[0]-a0; NT1b x1=b[1]-a1;
      NT2b y0=c[0]-a0; NT2b y1=c[1]-a1;
      return CGAL::compare(x0*y1,x1*y0);
    }
  };

}
#endif
