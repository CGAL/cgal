// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_3_TORUS_SHAPE_H
#define CGAL_SHAPE_DETECTION_3_TORUS_SHAPE_H

#include "Shape_base.h"
#include <cmath>
#include <CGAL/Circle_2.h>

/*!
 \file Torus_shape.h
 */

namespace CGAL {
    /*!
     \ingroup PkgPointSetShapeDetection3
     \brief Torus_shape implements Shape_base. The torus is represented by the symmetry axis, its center on the axis and the major and minor radii.     */
  template <class Sd_traits>
  class Torus_shape : public Shape_base<Sd_traits> {
  public:
    /// \cond SKIP_IN_MANUAL
    typedef typename Sd_traits::Input_iterator Input_iterator;
     ///< random access iterator for input data.
    typedef typename Sd_traits::Point_pmap Point_pmap;
     ///< property map to access the location of an input point.
    typedef typename Sd_traits::Normal_pmap Normal_pmap;
     ///< property map to access the unoriented normal of an input point.
    typedef typename Sd_traits::Geom_traits::FT FT; ///< number type.
    typedef typename Sd_traits::Geom_traits::Point_3 Point; ///< point type.
    typedef typename Sd_traits::Geom_traits::Vector_3 Vector; ///< vector type.
    typedef typename Sd_traits::Geom_traits::Point_2 Point_2;
     ///< 2D point type used during construction.
    typedef typename Sd_traits::Geom_traits::Circle_2 Circle;
     ///< cricle type used during construction.
    /// \endcond

    Torus_shape() : Shape_base<Sd_traits>() {}
      
    /*!
      Direction of symmetry axis.
     */
    Vector axis() const {
      return m_axis;
    }

    /*!
      Center point on symmetry axis.
     */
    Point center() const {
      return m_center;
    }
      
    /*!
      Major radius of the torus.
     */
    FT major_radius() const {
      return m_majorRad;
    }
      
    /*!
      Minor radius of the torus.
      */
    FT minor_radius() const {
      return m_minorRad;
    }

    /// \cond SKIP_IN_MANUAL

    /*!
      Helper function to write center point, symmetry axis
      and the two radii into a string.
     */
    std::string info() const {
      std::stringstream sstr;
      sstr << "Type: torus center(" << m_center.x() << ", " << m_center.y();
      sstr << ", " << m_center.z() << ") axis(" << m_axis.x() << ", ";
      sstr << m_axis.y() << ", " << m_axis.z() << ") major radius = ";
      sstr << m_majorRad << " minor radius = " << m_minorRad << " #Pts: ";
      sstr << this->m_indices.size();

      return sstr.str();
    }

    /*!
      Computes squared Euclidean distance from query point to the shape.
      */
    FT squared_distance(const Point &p) const {
      const Vector d = p - m_center;
      
	    // height over symmetry plane
      const FT height = d * m_axis;

      // distance from axis in plane
      const FT l = sqrt(d * d - height * height);

      // inPlane distance from circle
      const FT l2 = m_majorRad - l;

      // distance from torus
      const FT squared_dist = sqrt(height * height + l2 * l2) - m_minorRad;

      return squared_dist * squared_dist;
    }
    /// \endcond

  protected:
    /// \cond SKIP_IN_MANUAL
    void create_shape(const std::vector<std::size_t> &indices) {
      std::vector<Point> p;
      std::vector<Vector> n;

      p.resize(indices.size());
      n.resize(indices.size());

      for (std::size_t i = 0;i<indices.size();i++) {
        p[i] = this->point(indices[i]);
        n[i] = this->normal(indices[i]);
      }

      // Implemented method from 'Geometric least-squares fitting of spheres, cylinders, cones and tori' by G. Lukacs,A.D. Marshall, R. R. Martin
      double a01 = CGAL::cross_product(n[0], n[1]) * n[2];
      double b01 = CGAL::cross_product(n[0], n[1]) * n[3];
      double a0 = CGAL::cross_product(p[2] - p[1], n[0]) * n[2];
      double b0 = CGAL::cross_product(p[3] - p[1], n[0]) * n[3];
      double a1 = CGAL::cross_product(p[0] - p[2], n[1]) * n[2];
      double b1 = CGAL::cross_product(p[0] - p[3], n[1]) * n[3];
      double a = CGAL::cross_product(p[0] - p[2], p[1] - p[0]) * n[2];
      double b = CGAL::cross_product(p[0] - p[3], p[1] - p[0]) * n[3];

      double div = 1.0 / (b1 * a01 - b01 * a1);
      double r = ((a01 * b + b1 * a0 - b0 * a1 - b01 * a)) * div * 0.5;
      double q = (b * a0 - b0 * a) * div;

      FT root = r * r - q;
      if (r * r - q < 0)
        root = 0;

      double y1 = -r - sqrt(root);
      double y2 = -r + sqrt(root);
      double x1 = -(a1 * y1 + a) / (a01 * y1 + a0);
      double x2 = -(a1 * y2 + a) / (a01 * y2 + a0);

      // 1. center + axis
      FT majorRad1, minorRad1, dist1 = FLT_MAX;
      Point c1;
      Vector axis1;
      if (is_finite(x1) && is_finite(y1)) {
        c1 = p[0] + n[0] * x1;
        axis1 = c1 - (p[1] + n[1] * y1);

        FT l = axis1.squared_length();
        if (l > 0.00001 && l == l) {
          axis1 = axis1 / sqrt(l);
          dist1 = getCircle(c1, axis1, p, majorRad1, minorRad1);
        }
      }

      // 2. center + axis
      FT majorRad2 = 0, minorRad2 = 0, dist2 = FLT_MAX;
      Point c2;
      Vector axis2;
      if (is_finite(x2) && is_finite(y2)) {
        c2 = p[0] + n[0] * x2;
        axis2 = c2 - (p[1] + n[1] * y2);

        FT l = axis2.squared_length();
        if (l > 0.00001 && l == l) {
          axis2 = axis2 / sqrt(l);
          dist2 = getCircle(c2, axis2, p, majorRad2, minorRad2);
        }
      }

      if (dist1 < dist2) {
        m_center = c1;
        m_axis = axis1;
        m_majorRad = majorRad1;
        m_minorRad = sqrt(minorRad1);
      }
      else {
        m_center = c2;
        m_axis = axis2;
        m_majorRad = majorRad2;
        m_minorRad = sqrt(minorRad2);
      }

      //validate points and normals
      for (std::size_t i = 0;i<indices.size();i++) {
        // check distance
        if (squared_distance(p[i]) > this->m_epsilon) {
          this->m_isValid = false;
          return;
        }

        // check normal deviation
        Vector d = p[i] - m_center;
        // height over symmetry plane
        //FT p = d * m_axis;
        // distance from axis in plane
        //FT l = sqrt(d * d - p * p);

        Vector in_plane = CGAL::cross_product(m_axis,
          CGAL::cross_product(m_axis, d));
        if (in_plane * d < 0)
          in_plane = -in_plane;
        in_plane = in_plane / sqrt(in_plane.squared_length());

        d = p[i]
          - (m_center + in_plane * m_majorRad);
        d = d / sqrt(d.squared_length());
        if (abs(d * n[i]) < this->m_normal_threshold) {
          this->m_isValid = false;
          return;
        }
      }

      this->m_isValid = true;
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) {
      for (std::size_t i = 0;i<indices.size();i++) {
        Point po = this->point(indices[i]);
        Vector d = po - m_center;
        // height over symmetry plane
        const FT p = d * m_axis;
        // distance from axis in plane
        FT l = sqrt(d * d - p * p);

        // inPlane distance from circle
        const FT l2 = m_majorRad - l;

        // distance from torus
        l = sqrt(p * p + l2 * l2) - m_minorRad;
        dists[i] = l * l;
      }
    }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                               std::vector<FT> &angles) const {
      for (std::size_t i = 0;i<indices.size();i++) {
        Vector d = this->point(indices[i]) - m_center;

        Vector in_plane = CGAL::cross_product(m_axis,
                                              CGAL::cross_product(m_axis, d));
        if (in_plane * d < 0)
          in_plane = -in_plane;
        in_plane = in_plane / sqrt(in_plane.squared_length());

        d = this->point(indices[i]) - (m_center + in_plane * m_majorRad);
        d = d / sqrt(d.squared_length());
        angles[i] = abs(d * this->normal(indices[i]));
      }
    }

    FT cos_to_normal(const Point &p, const Vector &n) const {
      Vector d = p - m_center;

      Vector in_plane = CGAL::cross_product(m_axis,
                                           CGAL::cross_product(m_axis, d));
      if (in_plane * d < 0)
        in_plane = -in_plane;
      in_plane = in_plane / sqrt(in_plane.squared_length());

      d = p - (m_center + in_plane * m_majorRad);
      d = d / sqrt(d.squared_length());
      return abs(d * n);
    }
      
    virtual std::size_t required_samples() const {
        return 4;
    }

    virtual bool supports_connected_component() const {
      return false;
    }

    virtual bool wraps_u() const {
      return false;
    }

    virtual bool wraps_v() const {
      return false;
    }

  private:
    FT getCircle(Point &center, const Vector &axis, std::vector<Point> p, FT &majorRad, FT &minorRad) const {
      // create spin image
      std::vector<Point_2> pts;
      pts.resize(p.size());
      for (unsigned int i = 0;i<p.size();i++) {
        Vector d = p[i] - center;
        FT e = d * axis;
        FT f = d * d - e * e;
        if (f <= 0)
          pts[i] = Point_2(e, 0);
        else
          pts[i] = Point_2(e, sqrt(d * d - e * e));
      }

      if (CGAL::collinear(pts[0], pts[1], pts[2])) {
        return (std::numeric_limits<FT>::max)();
      }

      Circle c(pts[0], pts[1], pts[2]);
      minorRad = c.squared_radius();
      majorRad = c.center().y();
      center = center + c.center().x() * axis;

      return abs((pts[3] - c.center()).squared_length() - c.squared_radius());
    }

    Point m_center;
    Vector m_axis;
    FT m_majorRad;
    FT m_minorRad;
    /// \endcond
  };
}
#endif