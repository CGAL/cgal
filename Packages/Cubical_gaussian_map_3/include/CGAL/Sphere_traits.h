// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
//
// Author(s)     : Kapelushnik Lior <liorkape@post.tau.ac.il>

/*! \file
 * spherical arrangements of none intersecting arcs of great circles on a sphere
 */

#ifndef CGAL_SPHERE_TRAITS_H
#define CGAL_SPHERE_TRAITS_H

#include <CGAL/Sphere_arc.h>

#include <ostream>

CGAL_BEGIN_NAMESPACE

/*
  traits class for a sphere,
  holds curves that represents arcs of great circles on a sphere

  Kernel_ - the kernel with number types for spherical arcs endpoints directions
 */
template <class Kernel_>
class Sphere_traits : public Kernel_ {
public:
  // spherical traits types decleration
  typedef Kernel_                Kernel;
  typedef typename Kernel::Direction_3    Direction_3;
  typedef typename Kernel::Vector_3      Vector_3;

  typedef typename Kernel::Point_3      Point_3;
  // the curve is a spherical arc on a great circle
  typedef Sphere_arc<Kernel_>          X_monotone_curve_2;
  typedef X_monotone_curve_2           Curve_2;

private:
};

CGAL_END_NAMESPACE

#endif
