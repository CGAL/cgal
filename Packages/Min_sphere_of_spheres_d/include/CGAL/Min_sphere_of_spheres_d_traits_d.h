// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind.
//
// Every use of CGAL requires a license.
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation.
//
// Commercial licenses
// - Please check the CGAL web site http://www.cgal.org/index2.html for
//   availability.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// chapter       : $CGAL_Chapter: Optimisation $
// file          : include/CGAL/Min_sphere_of_spheres_d_traits_d.h
// package       : Min_sphere_of_spheres_d (1.10)
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Kaspar Fischer
// maintainer    : Kaspar Fischer <fischerk@inf.ethz.ch>
// coordinator   : ETH Zurich (Kaspar Fischer)
//
// implementation: dD Smallest Enclosing Sphere of Spheres
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

#ifndef CGAL_MIN_SPHERE_OF_SPHERES_D_TRAITS_D_H
#define CGAL_MIN_SPHERE_OF_SPHERES_D_TRAITS_D_H

#include <CGAL/Cartesian_d.h>
#include <CGAL/Weighted_point.h>

namespace CGAL {

  template<typename K_,                      // kernel
    typename FT_,                            // number type
    int Dim_,                                // dimension
    typename UseSqrt_ = Tag_false,           // whether to use square-roots
    typename Algorithm_ = Default_algorithm> // algorithm to use
  class Min_sphere_of_spheres_d_traits_d {
  public: // types:
    typedef FT_ FT;
    typedef FT_ Weight;
    typedef CGAL::Point_d<K_> Point;
    typedef CGAL::Weighted_point<Point,Weight> Sphere;
    typedef typename Point::Cartesian_const_iterator Cartesian_const_iterator;
    typedef UseSqrt_ Use_square_roots;
    typedef Algorithm_ Algorithm;

  public: // constants:
    static const int D = Dim_;               // dimension

  public: // accessors:
    static inline const FT radius(const Sphere& s) {
      return s.weight();
    }

    static inline Cartesian_const_iterator
      center_cartesian_begin(const Sphere& s) {
      return s.cartesian_begin();
    }
  };

} // namespace CGAL

#endif // CGAL_MIN_SPHERE_OF_SPHERES_D_TRAITS_D_H
