// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SKIN_SURFACE_QUADRATIC_SURFACE_H
#define CGAL_SKIN_SURFACE_QUADRATIC_SURFACE_H

#include <CGAL/Weighted_point.h>

CGAL_BEGIN_NAMESPACE

template < class Kernel >
class Skin_surface_quadratic_surface_3 {
public:
  typedef Kernel                          K;
  typedef typename K::Point_3             Point;
  typedef typename K::Vector_3            Vector;
  typedef typename K::Segment_3           Segment;
  typedef typename K::RT                  RT;
  typedef Weighted_point<Point, RT> Weighted_point;

  Skin_surface_quadratic_surface_3(RT Qinput[], Point p, RT c)
    : p(p), c(c) 
  {
    for (int i=0; i<6; i++) Q[i] = Qinput[i];
  }

  template <class Input_point>
  RT value(Input_point const &x) const {
    typedef Cartesian_converter<typename Input_point::R, K> Converter;
    
    Vector v = Converter()(x) - p;

    return 
      v.x()*(Q[0]*v.x()) +
      v.y()*(Q[1]*v.x()+Q[2]*v.y()) +
      v.z()*(Q[3]*v.x()+Q[4]*v.y()+Q[5]*v.z()) +
      c;
  }


  Vector gradient(Point const &x) {
    std::cout << "NGHK: NOT YET IMPLEMENTED" << std::endl;
    // NGHK: TODO:
    return (x-p);
  }
	
  /// Construct the intersection point with the segment (p0,p1)
  template <class Input_point>
  Input_point to_surface(const Input_point &p0, const Input_point &p1) {
    // NGHK: We use the kernel trick again: DOCUMENT!!!!
    typedef Cartesian_converter<typename Input_point::R, K> Converter;
    Converter conv;

    Input_point pp0 = p0, pp1 = p1, mid;
    double sq_d = to_double(squared_distance(pp0,pp1));
    
    if (value(conv(pp1)) < value(conv(pp0))) {
      std::swap(pp0, pp1);
    }
    
    while (sq_d > 1.e-10) {
      mid = midpoint(pp0,pp1);
      if (value(conv(mid)) > 0) {
	pp1 = mid;
      } else {
	pp0 = mid;
      }
      sq_d /= 4;
    }
    return mid;
  };
private:
  RT Q[6]; 
  Point p; 
  RT c;
};

// template < class Kernel >
// class Skin_surface_quadratic_surface_3 {
// public:
//   typedef Kernel                          K;
//   typedef typename K::Point_3             Point;
//   typedef typename K::Vector_3            Vector;
//   typedef typename K::Segment_3           Segment;
//   typedef typename K::RT                  RT;
//   typedef Weighted_point<Point, RT> Weighted_point;

//   Skin_surface_quadratic_surface_3() {
//   }
//   virtual ~Skin_surface_quadratic_surface_3() {};
//   /// Construct the intersection point with the segment (p0,p1)
//   template <class Input_point>
//   Input_point to_surface(const Input_point &p0, const Input_point &p1) {
//     // NGHK: We use the kernel trick again: DOCUMENT!!!!
//     typedef Cartesian_converter<typename Input_point::R, K> Converter;
//     Converter conv;

//     Input_point pp0 = p0, pp1 = p1, mid;
//     double sq_d = to_double(squared_distance(pp0,pp1));
    
//     if (value(conv(pp1)) < value(conv(pp0))) {
//       std::swap(pp0, pp1);
//     }
    
//     while (sq_d > 1.e-10) {
//       mid = midpoint(pp0,pp1);
//       if (value(conv(mid)) > 0) {
// 	pp1 = mid;
//       } else {
// 	pp0 = mid;
//       }
//       sq_d /= 4;
//     }
//     return mid;
//   };
// //   virtual Point to_surface(Point const &p, Vector const &v) = 0;
//   inline Point to_surface(const Segment &s) {
//     return to_surface(s.source(), s.target());
//   }

// //   // Gradient descent
// //   virtual Point to_surface(Point const &p0) = 0;

//   // compute the function value of p
//   template <class Input_point>
//   RT value(const Input_point &p) const {
//     typedef Cartesian_converter<typename Input_point::R, K> Converter;
//     return _value(Converter()(p));
//   }
//   virtual RT _value(const Point &p) const = 0;
//   // Compute the gradient in p
//   virtual Vector gradient(const Point &p) = 0; 

// //   // NGHK: REMOVE: 
//   virtual int dimension() const = 0;

// //   // return a continuous density function on the skin surface:
// //   virtual RT sq_density(const Point &p) = 0;

// };


// /* *************************************************
//  *                  Hyperboloid
//  *
//  * Point p.
//  * Let y = (p-wp)*t/|t|
//  *     x = |p-wp|^2 - y^2
//  *
//  * orient*((1-s) x^2 - s y^2 + s(1-s) W = 0
//  * orient*((1-s)|p-wp|^2 - y^2 + s(1-s)W = 0
//  *
//  ***************************************************/

// template < class Kernel >
// class Skin_surface_hyperboloid_3 : public Skin_surface_quadratic_surface_3<Kernel> {
// public:
//   typedef Skin_surface_quadratic_surface_3<Kernel> Parent; 
//   typedef typename Parent::Point              Point;
//   typedef typename Parent::Vector             Vector;
//   typedef typename Parent::RT                 RT;
//   typedef typename Parent::Weighted_point     Weighted_point;

//   Skin_surface_hyperboloid_3(Weighted_point wp, Vector t, RT s, int orient)
//     : Skin_surface_quadratic_surface_3<Kernel>(), wp(wp), t(t), s(s), orient(orient) {
//     assert((orient == -1) || (orient == 1));
//     sq_t = t*t;
//   }
	
//   RT _value(Point const &x) const {
//     Vector dir = x-wp;
//     RT tmp = dir*t;
//     tmp = tmp*tmp/sq_t;

//     return orient * (dir*dir/s - tmp/(s*(1-s))) + wp.weight();
//   }
// //   Point to_surface(Point const &p0, Point const &p1) {
// //     assert(value(p0) * value(p1) <= 0);

// //     Vector p0c = p0-wp;
// //     Vector p1p0 = p1-p0;
// //     RT sq_d = p1p0*p1p0;
		
// //     RT p0_sym = p0c*t;
// //     RT p1_sym = (p1-wp)*t;
// //     RT d_sym = p0_sym-p1_sym;

// //     RT den = ((1-s)*sq_d - d_sym*d_sym/sq_t);
// //     RT top = -((1-s)*(p0c*p1p0) + p0_sym*d_sym/sq_t) / den;
// //     RT extr_val =
// //       orient*((1-s)*p0c*p0c-p0_sym*p0_sym/sq_t-s*(1-s)*wp.weight());
			
// //     {
// //       Point extr = p0+top*p1p0;
// //       Vector dir = extr-wp;
// //       RT tmp = dir*t; tmp *= tmp/sq_t;
// //       extr_val = orient*(-tmp + (1-s)*dir*dir + s*(1-s)*wp.weight());
// //     }
// //     RT d = sqrt(-orient*extr_val/den);

// //     RT t, t1;
// //     t = top + d; t1 = top - d;
// //     if ((2*t1-1)*(2*t1-1) < (2*t-1)*(2*t-1)) t = t1;

// //     if (t < 0) {
// //       std::cerr << "Hl[" << dimension() <<"] " <<  t << "\n";
// //       return p0;
// //     } else if (t > 1) {
// //       std::cerr << "Hh[" << dimension() <<"] " <<  t << "\n";
// //       return  p1;
// //     } else {
// //       assert (std::abs(value(p0 + t*p1p0)) < 0.001);
// //       return p0 + t*p1p0;
// //     }
// //   }
//   Point to_surface(Point const &p, Vector const &v){
//     Vector pc = p-wp;
//     Point p1 = p + v;
//     RT sq_d = v*v;
		
//     RT p_sym = pc*t;
//     RT p1_sym = (p1-wp)*t;
//     RT d_sym = p_sym-p1_sym;

//     RT den = ((1-s)*sq_d - d_sym*d_sym/sq_t);
//     RT top = -((1-s)*(pc*v) + p_sym*d_sym/sq_t) / den;
//     RT extr_val =
//       orient*((1-s)*pc*pc-p_sym*p_sym/sq_t-s*(1-s)*wp.weight());
			
//     {
//       Point extr = p+top*v;
//       Vector dir = extr-wp;
//       RT tmp = dir*t; tmp *= tmp/sq_t;
//       extr_val = orient*(-tmp + (1-s)*dir*dir + s*(1-s)*wp.weight());
//     }
//     RT d = -orient*extr_val/den;
//     if (d<0) return p;
//     d = sqrt(d);

//     RT t, t1;
//     t = top + d; t1 = top - d;
//     if (t1*t1 < t*t) t = t1;

//     assert (std::abs(value(p + t*v)) < 0.001);
//     return p + t*v;
//   }
//   Point to_surface(Point const &p0) {
//     return to_surface(p0, gradient(p0));
//   }
//   Vector gradient(Point const &p) {
//     // -s x + (1-s) y 
//     Vector v = p - wp;
//     Vector vt = (v*t)/(t*t)*t;
//     return orient*((1-s)*v - vt);
//   }

//   int dimension() const {
//     if (orient == 1) return 1; else return 2;
//   }

//   RT sq_density(Point const &p) {
//     Vector v = p - wp;
//     RT vt = v*t;
//     RT scale = 1-s;
//     if (s < RT(.5)) {
//       scale = (scale-s)/(scale*scale);
//     } else {
//       scale = (s-scale)/(s*s);
//     }
//     assert(scale >= 0);
		
//     return squared_distance(wp, p) - scale*vt*vt/(t*t);
//   }
// private:
//   Weighted_point wp;
//   Vector t;
//   RT s, sq_t;
//   int orient;
	
// };

CGAL_END_NAMESPACE

#endif // CGAL_SKIN_SURFACE_QUADRATIC_SURFACE_H
