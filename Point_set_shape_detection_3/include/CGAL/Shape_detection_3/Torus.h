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

#ifndef CGAL_SHAPE_DETECTION_3_TORUS_H
#define CGAL_SHAPE_DETECTION_3_TORUS_H

#include <CGAL/Shape_detection_3/Shape_base.h>
#include <CGAL/Circle_2.h>
#include <CGAL/number_utils.h>
#include <cmath>

/*!
 \file Torus.h
 */

namespace CGAL {
  namespace Shape_detection_3 {
    /*!
      \ingroup PkgPointSetShapeDetection3Shapes
      \brief Torus implements Shape_base. The torus is represented by the
      symmetry axis, its center on the axis and the major and minor radii.
      \tparam Traits a model of `EfficientRANSACTraits` with the additional
              requirement that the types `Traits::Point_2` and `Traits::Circle_2`
              are provided.
     */
  template <class Traits>
  class Torus : public Shape_base<Traits> {
  public:
    /// \cond SKIP_IN_MANUAL
    typedef typename Traits::Point_map Point_map;
     ///< property map to access the location of an input point.
    typedef typename Traits::Normal_map Normal_map;
     ///< property map to access the unoriented normal of an input point.
    typedef typename Traits::FT FT; ///< number type.
    typedef typename Traits::Point_3 Point_3; ///< point type.
    typedef typename Traits::Vector_3 Vector_3; ///< vector type.
    typedef typename Traits::Point_2 Point_2;
     ///< 2D point type used during construction.
    typedef typename Traits::Circle_2 Circle;
     ///< cricle type used during construction.
    /// \endcond

    Torus() : Shape_base<Traits>() {}
      
    /*!
      Direction of symmetry axis.
     */
    Vector_3 axis() const {
      return m_axis;
    }

    /*!
      Center point on symmetry axis.
     */
    Point_3 center() const {
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
    FT squared_distance(const Point_3 &p) const {
      const Vector_3 d = p - m_center;
      
	    // height over symmetry plane
      const FT height = d * m_axis;

      // distance from axis in plane
      const FT l = CGAL::sqrt(d * d - height * height);

      // inPlane distance from circle
      const FT l2 = m_majorRad - l;

      // distance from torus
      const FT squared_dist = CGAL::sqrt(height * height + l2 * l2) - m_minorRad;

      return squared_dist * squared_dist;
    }
    /// \endcond

  protected:
    /// \cond SKIP_IN_MANUAL
    void create_shape(const std::vector<std::size_t> &indices) {
      std::vector<Point_3> p;
      std::vector<Vector_3> n;

      p.resize(indices.size());
      n.resize(indices.size());

      for (std::size_t i = 0;i<indices.size();i++) {
        p[i] = this->point(indices[i]);
        n[i] = this->normal(indices[i]);
      }

      // Implemented method from 'Geometric least-squares fitting of spheres, cylinders, cones and tori' by G. Lukacs,A.D. Marshall, R. R. Martin
      const FT a01 = CGAL::cross_product(n[0], n[1]) * n[2];
      const FT b01 = CGAL::cross_product(n[0], n[1]) * n[3];
      const FT a0 = CGAL::cross_product(p[2] - p[1], n[0]) * n[2];
      const FT b0 = CGAL::cross_product(p[3] - p[1], n[0]) * n[3];
      const FT a1 = CGAL::cross_product(p[0] - p[2], n[1]) * n[2];
      const FT b1 = CGAL::cross_product(p[0] - p[3], n[1]) * n[3];
      const FT a = CGAL::cross_product(p[0] - p[2], p[1] - p[0]) * n[2];
      const FT b = CGAL::cross_product(p[0] - p[3], p[1] - p[0]) * n[3];

      FT div = (b1 * a01 - b01 * a1);
      if (div == 0)
        return;

      div = (FT)1.0 / div;
      const FT r = ((a01 * b + b1 * a0 - b0 * a1 - b01 * a)) * div * (FT)0.5;
      const FT q = (b * a0 - b0 * a) * div;

      FT root = r * r - q;
      if (r * r - q < 0)
        root = 0;

      const FT y1 = -r - CGAL::sqrt(root);
      const FT y2 = -r + CGAL::sqrt(root);
      FT x1 = (a01 * y1 + a0);
      FT x2 =  (a01 * y2 + a0);

      if (x1 == 0 || x2 == 0)
        return;

      x1 = -(a1 * y1 + a) / x1;
      x2 = -(a1 * y2 + a) / x2;

      // 1. center + axis
      FT majorRad1 = FLT_MAX, minorRad1 = FLT_MAX, dist1 = FLT_MAX;
      Point_3 c1;
      Vector_3 axis1;
      if (is_finite(x1) && is_finite(y1)) {
        c1 = p[0] + n[0] * x1;
        axis1 = c1 - (p[1] + n[1] * y1);

        FT l = axis1.squared_length();
        if (l > (FT)0.00001 && l == l) {
          axis1 = axis1 / CGAL::sqrt(l);
          dist1 = getCircle(c1, axis1, p, majorRad1, minorRad1);
        }
      }

      // 2. center + axis
      FT majorRad2 = 0, minorRad2 = 0, dist2 = FLT_MAX;
      Point_3 c2;
      Vector_3 axis2;
      if (is_finite(x2) && is_finite(y2)) {
        c2 = p[0] + n[0] * x2;
        axis2 = c2 - (p[1] + n[1] * y2);

        FT l = axis2.squared_length();
        if (l > (FT)0.00001 && l == l) {
          axis2 = axis2 / CGAL::sqrt(l);
          dist2 = getCircle(c2, axis2, p, majorRad2, minorRad2);
        }
      }

      if (dist1 < dist2) {
        m_center = c1;
        m_axis = axis1;
        m_majorRad = majorRad1;
        m_minorRad = CGAL::sqrt(minorRad1);
      }
      else {
        m_center = c2;
        m_axis = axis2;
        m_majorRad = majorRad2;
        m_minorRad = CGAL::sqrt(minorRad2);
      }

      // Drop if shape is probably sphere
      if (m_majorRad < this->m_epsilon) {
        this->m_is_valid = false;

        return;
      }

      //validate points and normals
      for (std::size_t i = 0;i<indices.size();i++) {
        // check distance
        if (squared_distance(p[i]) > this->m_epsilon) {
          this->m_is_valid = false;
          return;
        }

        // check normal deviation
        Vector_3 d = p[i] - m_center;

        Vector_3 in_plane = CGAL::cross_product(m_axis,
          CGAL::cross_product(m_axis, d));
        if (in_plane * d < 0)
          in_plane = -in_plane;

        FT length = CGAL::sqrt(in_plane.squared_length());
        if (length == 0)
          return;

        in_plane = in_plane / length;

        d = p[i] - (m_center + in_plane * m_majorRad);

        length = CGAL::sqrt(d.squared_length());
        if (length == 0)
          return;

        d = d / length;
        if (CGAL::abs(d * n[i]) < this->m_normal_threshold) {
          this->m_is_valid = false;
          return;
        }
      }

      this->m_is_valid = true;
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) const {
      for (std::size_t i = 0;i<indices.size();i++) {
        Point_3 po = this->point(indices[i]);
        Vector_3 d = po - m_center;
        // height over symmetry plane
        const FT p = d * m_axis;
        // distance from axis in plane
        FT l = CGAL::sqrt(d * d - p * p);

        // inPlane distance from circle
        const FT l2 = m_majorRad - l;

        // distance from torus
        l = CGAL::sqrt(p * p + l2 * l2) - m_minorRad;
        dists[i] = l * l;
      }
    }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                               std::vector<FT> &angles) const {
      for (std::size_t i = 0;i<indices.size();i++) {
        Vector_3 d = this->point(indices[i]) - m_center;

        Vector_3 in_plane = CGAL::cross_product(m_axis,
                                              CGAL::cross_product(m_axis, d));
        if (in_plane * d < 0)
          in_plane = -in_plane;

        FT length = (FT) CGAL::sqrt(in_plane.squared_length());

        // If length is 0 the point is on the axis, maybe in the apex. We
        // accept any normal for that position.
        if (length == 0) {
          angles[i] = (FT)1.0;
          continue;
        }

        in_plane = in_plane / CGAL::sqrt(in_plane.squared_length());

        d = this->point(indices[i]) - (m_center + in_plane * m_majorRad);
        d = d / CGAL::sqrt(d.squared_length());
        angles[i] = CGAL::abs(d * this->normal(indices[i]));
      }
    }

    FT cos_to_normal(const Point_3 &p, const Vector_3 &n) const {
      Vector_3 d = p - m_center;

      Vector_3 in_plane = CGAL::cross_product(m_axis,
                                           CGAL::cross_product(m_axis, d));
      if (in_plane * d < 0)
        in_plane = -in_plane;
      
      float length = CGAL::sqrt(in_plane.squared_length());

      // If length is 0 the point is on the axis, maybe in the apex. We
      // accept any normal for that position.
      if (length == 0) {
        return (FT)1.0;
      }

      in_plane = in_plane / CGAL::sqrt(in_plane.squared_length());

      d = p - (m_center + in_plane * m_majorRad);
      d = d / CGAL::sqrt(d.squared_length());

      return CGAL::abs(d * n);
    }
      
    virtual std::size_t minimum_sample_size() const {
        return 4;
    }

    virtual bool supports_connected_component() const {
      return false;
    }

  private:
    FT getCircle(Point_3 &center, const Vector_3 &axis, std::vector<Point_3> p, FT &majorRad, FT &minorRad) const {
      // create spin image
      std::vector<Point_2> pts;
      pts.resize(p.size());
      for (unsigned int i = 0;i<p.size();i++) {
        Vector_3 d = p[i] - center;
        FT e = d * axis;
        FT f = d * d - e * e;
        if (f <= 0)
          pts[i] = Point_2(e, (FT) 0);
        else
          pts[i] = Point_2(e, CGAL::sqrt(d * d - e * e));
      }

      if (CGAL::collinear(pts[0], pts[1], pts[2])) {
        return (std::numeric_limits<FT>::max)();
      }

      Circle c(pts[0], pts[1], pts[2]);
      minorRad = c.squared_radius();
      majorRad = c.center().y();
      center = center + c.center().x() * axis;

      return CGAL::abs((pts[3] - c.center()).squared_length() - c.squared_radius());
    }

    Point_3 m_center;
    Vector_3 m_axis;
    FT m_majorRad;
    FT m_minorRad;
    /// \endcond
  };
}
}
#endif
