// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_ARR_PLANE_3_h
#define CGAL_ARR_PLANE_3_h

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Construct and maintain a plane in 3D that contains the origin. A plane is
 * defined by 4 coordinates, namely a, b, c, and d. When the plane contains
 * the origin, the last coordinate d vanishes. This simplifies arithmetic
 * expressions involved with the construction of the plane and other
 * operations applied to planes.
 */

#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/Kernel/solve.h>

namespace CGAL {

/*! A plane that contains the origin extended with a few methods */
template <class Kernel_>
class Arr_plane_3 {
public:
  typedef Kernel_                       Kernel;
  typedef Arr_plane_3<Kernel>           Self;
  typedef typename Kernel::Point_2      Point_2;
  typedef typename Kernel::Point_3      Point_3;
  typedef typename Kernel::Vector_3     Vector_3;
  typedef typename Kernel::Direction_3  Direction_3;
  typedef typename Kernel::Line_3       Line_3;
  typedef typename Kernel::FT           FT;

private:
  /*! The x coefficient */
  FT m_a;

  /*! The y coefficient */
  FT m_b;

  /*! The z coefficient */
  FT m_c;
  
public:
  /*! Default Constructor */
  Arr_plane_3() : m_a(0), m_b(0), m_c(0) {}

  /*! Constructor */
  Arr_plane_3(int a, int b, int c) : m_a(a), m_b(b), m_c(c) {}

  /*! Constructor */
  Arr_plane_3(typename Kernel::Plane_3 p)
  {
    CGAL_precondition_code(Kernel kernel;);
    CGAL_precondition_code(typename Kernel::Point_3 orig = kernel.construct_point_3_object()(ORIGIN););
    CGAL_precondition(kernel.has_on_3_object()(p, orig));

    m_a = p.a(); m_b = p.b(); m_c = p.c() ;
  }

  /*! Constructor */
  Arr_plane_3(const Point_3 & p, const Point_3 & r)
  {
    FT prx = r.x() - p.x();
    FT pry = r.y() - p.y();
    FT prz = r.z() - p.z();
    m_a = r.y() * prz - pry * r.z();
    m_b = r.z() * prx - prz * r.x();
    m_c = r.x() * pry - prx * r.y();
  }

  /*! Obtain the x coefficient */
  const FT & a() const { return m_a; }

  /*! Obtain the y coefficient */
  const FT & b() const { return m_b; }

  /*! Obtain the z coefficient */
  const FT & c() const { return m_c; }
  
  /*! Obtain the i-th coefficient of the plane
   * \param i the index of the coefficient
   * \return the i-th coefficient
   */
  FT operator[](unsigned int i) const
  {
    CGAL_assertion(i < 3);
    return (i == 0) ? m_a : ((i == 1) ? m_b : m_c);
  }

  bool equal(Self plane) const
  {
    return ((a() == plane.a()) &&
            (b() == plane.b()) &&
            (c() == plane.c()));
  }

  /*! Convert to kernel's plane */
  operator typename Kernel::Plane_3 () const
  {
    Kernel kernel;
    return kernel.construct_plane_3_object() (m_a, m_b, m_c, 0);
  }

  /*! Compute the image point of the projection of p under an affine
   * transformation, which maps the plane onto the xy-plane, with the
   * z-coordinate removed.
   * \param p the point
   * \return the image point
   */
  Point_2 to_2d(const Point_3 & p) const
  {
    Kernel kernel;
    typename Kernel::Plane_3 base_plane(m_a, m_b, m_c, 0);
    Vector_3 v1 = kernel.construct_base_vector_3_object()(base_plane, 1);    
    Vector_3 v2 = kernel.construct_base_vector_3_object()(base_plane, 2);    
    Vector_3 v3 = kernel.construct_orthogonal_vector_3_object()(base_plane);

    FT denom = v2[1]*v3[0]*v1[2] - v2[0]*v3[1]*v1[2] + v3[2]*v2[0]*v1[1] +
               v2[2]*v3[1]*v1[0] - v3[0]*v2[2]*v1[1] - v2[1]*v3[2]*v1[0];
    FT x = - (v2[1]*v3[2]*p[0] - v2[1]*v3[0]*p[2] + v3[0]*v2[2]*p[1] +
              v2[0]*v3[1]*p[2] - v3[2]*v2[0]*p[1] - v2[2]*v3[1]*p[0]) / denom;
    FT y = (v1[1]*v3[2]*p[0] - v1[1]*v3[0]*p[2] - v3[1]*p[0]*v1[2] +
            v3[1]*v1[0]*p[2] + p[1]*v3[0]*v1[2] - p[1]*v3[2]*v1[0]) / denom;
    
    return Point_2(x, y);
  }

  /*! Compute a 3d point p_3 coincident to the plane, such that the image point
   * of the projection of p_3 under an affine transformation, which maps the
   * plane onto the a given axis-parallel plane is a given 2d point.
   * \param p_2 the image point
   * \param i the index of the axis-parallel plane. 0, 1, or 2 indicate the
   * yz-, zx-, and xy-plane respectively
   * \return the coincident point
   */
  Point_3 to_3d(const Point_2 & p_2, unsigned int i) const
  {
    CGAL_assertion(i < 3);

#if 0
    std::cout << "(a, b, c): " << a() << "," << b() << "," << c()
              << std::endl;
#endif
    
    if (i == 0) {
      CGAL_assertion(m_a != 0);
      FT y = p_2.x();
      FT z = p_2.y();
      FT x = -(m_b * y + m_c * z) / m_a;
      Point_3 p_3(x, y, z);
      return p_3;
    }

    if (i == 1) {
      CGAL_assertion(m_b != 0);
      FT z = p_2.x();
      FT x = p_2.y();
      FT y = -(m_a * x + m_c * z) / m_b;
      Point_3 p_3(x, y, z);
      return p_3;
    }

    // if (i == 2) return Base::to_3d(p_2);
    CGAL_assertion(m_c != 0);
    FT x = p_2.x();
    FT y = p_2.y();
    FT z = -(m_a * x + m_b * y) / m_c;
    Point_3 p_3(x, y, z);
    return p_3;
  }

  /*! Determine the relative position of a point and the plane
   * \param p the point
   * \return ON_ORIENTED_BOUNDARY, ON_POSITIVE_SIDE, or ON_NEGATIVE_SIDE,
   * determined by the position of p relative to the oriented plane.
   */
  Oriented_side oriented_side(const Point_3 & p) const
  {
    return CGAL_NTS sign(m_a*p.x() + m_b*p.y() + m_c*p.z());
  }
};

/*! Intersect 2 planes
 * \param plane1 the first plane
 * \param plane2 the second plane
 * \return a geometric object that represents the intersection. Could be
 * the line of intersection, or a plane in case plane1 and plane2 coincide.
 */
template <class Kernel>
CGAL::Object intersect(const Arr_plane_3<Kernel> & plane1,
                       const Arr_plane_3<Kernel> & plane2)
{
  typedef typename Kernel::Point_3      Point_3;
  typedef typename Kernel::Direction_3  Direction_3;
  typedef typename Kernel::Line_3       Line_3;
  typedef typename Kernel::FT           FT;

  // We know that the plane goes throgh the origin
  const FT & a1 = plane1.a();
  const FT & b1 = plane1.b();
  const FT & c1 = plane1.c();

  const FT & a2 = plane2.a();
  const FT & b2 = plane2.b();
  const FT & c2 = plane2.c();

  FT det = a1*b2 - a2*b1;
  if (det != 0) {
    Point_3 is_pt = Point_3(0, 0, 0, det);
    Direction_3 is_dir = Direction_3(b1*c2 - c1*b2, a2*c1 - a1*c2, det);
    return make_object(Line_3(is_pt, is_dir));
  }

  det = a1*c2 - a2*c1;
  if (det != 0) {
    Point_3 is_pt = Point_3(0, 0, 0, det);
    Direction_3 is_dir = Direction_3(c1*b2 - b1*c2, det, a2*b1 - a1*b2);
    return make_object(Line_3(is_pt, is_dir));
  }
  det = b1*c2 - c1*b2;
  if (det != 0) {
    Point_3 is_pt = Point_3(0, 0, 0, det);
    Direction_3 is_dir = Direction_3(det, c1*a2 - a1*c2, a1*b2 - b1*a2);
    return make_object(Line_3(is_pt, is_dir));
  }

  // degenerate case
  return make_object(plane1);
}

/*! Compute the image point of the projection of p under an affine
 * transformation, which maps the plane onto the xy-plane, with the
 * z-coordinate removed.
 * \param plane the plane
 * \param p the point
 * \return the image point
 */
template <class Kernel>
typename Kernel::Point_2
construct_projected_xy_point(const Arr_plane_3<Kernel> & plane,
                             const typename Kernel::Point_3 & p)
{
  return plane.to_2d(p);
}

/*! Export a plane to an output stream
 * \param os the output stream
 * \param plane the plane
 * \return the output stream
 */
template <class Kernel>
inline std::ostream & operator<<(std::ostream & os,
                                 const Arr_plane_3<Kernel> & plane)
{
  os << plane[0] << ", " << plane[1] << ", " << plane[2];
  return os;
}  

} //namespace CGAL

#endif
