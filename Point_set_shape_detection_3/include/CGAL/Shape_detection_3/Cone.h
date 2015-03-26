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

#ifndef CGAL_SHAPE_DETECTION_3_CONE_H
#define CGAL_SHAPE_DETECTION_3_CONE_H

#include <CGAL/Shape_detection_3/Shape_base.h>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/*!
 \file Cone.h
 */


namespace CGAL {
  namespace Shape_detection_3 {
  /*!
   \brief Cone implements Shape_base.
    The cone is represented by its apex, the axis and the opening angle.
    This representation models an open infinite single-cone.
   \ingroup PkgPointSetShapeDetection3
   */

  template <class ERTraits>
  class Cone : public Shape_base<ERTraits> {
  public:
    /// \cond SKIP_IN_MANUAL
    typedef typename ERTraits::Input_iterator Input_iterator;
     ///< random access iterator for input data.
    typedef typename ERTraits::Point_pmap Point_pmap;
     ///< property map to access the location of an input point.
    typedef typename ERTraits::Normal_pmap Normal_pmap;
     ///< property map to access the unoriented normal of an input point.
    typedef typename ERTraits::Geom_traits::FT FT; ///< number type.
    typedef typename ERTraits::Geom_traits::Point_3 Point;///< point type.
    typedef typename ERTraits::Geom_traits::Vector_3 Vector;///< vector type.
    /// \endcond

	
    Cone() : Shape_base<ERTraits>() {}
      
    /*!
      The opening angle between the axis and the surface of the cone.
     */
    FT angle() const {
        return m_angle;
    }
    
    /*!
      The apex of the cone.
     */
    Point apex() const {
        return m_apex;
    }
    
    /*!
      The axis points from the apex into the cone.
     */
    Vector axis() const {
        return m_axis;
    }
    
    /*!
      Helper function to write apex, axis and angle of the cone and
      number of assigned points into a string.
     */
    /// \cond SKIP_IN_MANUAL
    std::string info() const {
        std::stringstream sstr;
        
        sstr << "Type: cone apex: (" << m_apex.x() << ", " << m_apex.y();
        sstr << ", " << m_apex.z() << ") axis: (" << m_axis.x() << ", ";
        sstr << m_axis.y() << ", " << m_axis.z() << ") angle:" << m_angle;
        sstr << " #Pts: " << this->m_indices.size();
        
        return sstr.str();
    }

    /*!
    Computes squared Euclidean distance from query point to the shape.
    */ 
    FT squared_distance(const Point &p) const {
      Vector toApex = p - m_apex;
      FT a = toApex.squared_length();

      // projection on axis
      FT b = toApex * m_axis;

      // distance to axis
      if (a - b * b <= 0)
        return 0;

      FT l = CGAL::sqrt(a - b * b);
      FT c = m_cos_ang * l;
      FT d = m_neg_sin_ang * b;

      // far on other side?
      return (b < 0 && c - d < 0) ? a : CGAL::abs(c + d) * CGAL::abs(c + d);
    }
    /// \endcond

  protected:
      /// \cond SKIP_IN_MANUAL
    virtual void create_shape(const std::vector<std::size_t> &indices) {
      Point p1 = this->point(indices[0]);
      Point p2 = this->point(indices[1]);
      Point p3 = this->point(indices[2]);

      Vector n1 = this->normal(indices[0]);
      Vector n2 = this->normal(indices[1]);
      Vector n3 = this->normal(indices[2]);

      // first calculate intersection of three planes -> apex

      Vector lineDir = CGAL::cross_product(n1, n2);
      FT length = sqrt(lineDir.squared_length());
      if (length == 0)
        return;

      lineDir = lineDir * (FT)1.0 / length;

      // lineDir not normalized direction of intersection lines
      //  of two planes (p1, n1) and (p2, n2)
      // get point on line by moving point p1 onto line
      Vector orthLineInPlane = CGAL::cross_product(n1, lineDir);
      length = sqrt(orthLineInPlane.squared_length());
      if (length == 0)
        return;

      orthLineInPlane = orthLineInPlane * (FT)1.0 / length;

      // distance of p1 to (p2, n2)
      FT d = (p1 - CGAL::ORIGIN) * n2 - (p2 - CGAL::ORIGIN) * n2;
      // projection of orthLineInPlane onto p2
      FT l = orthLineInPlane * n2;
      if (l == 0)
        return;

      Point pointOnLine = p1 - (d/l) * orthLineInPlane;

      // distance of pLineDir to (p3, n3)
      d = (pointOnLine - CGAL::ORIGIN) * n3 - (p3 - CGAL::ORIGIN) * n3;
      l = lineDir * n3;
      if (l == 0)
        return;

      m_apex = pointOnLine - (d/l) * lineDir;

      // 2. find axis
      Vector v1 = p1 - m_apex;
      length = sqrt(v1.squared_length());
      if (length == 0)
        return;
      v1 = v1 * (FT)1.0 / length;
      Point c1 = m_apex + v1;

      Vector v2 = p2 - m_apex;
      length = sqrt(v2.squared_length());
      if (length == 0)
        return;
      v2 = v2 * (FT)1.0 / length;
      Point c2 = m_apex + v2;

      Vector v3 = p3 - m_apex;
      length = sqrt(v3.squared_length());
      if (length == 0)
        return;
      v3 = v3 * (FT)1.0 / length;
      Point c3 = m_apex + v3;

      m_axis = CGAL::cross_product(c1 - c2, c1 - c3);
      m_axis = (orthLineInPlane * m_axis < 0) ? -m_axis : m_axis;
      length = CGAL::sqrt(m_axis.squared_length());
      if (length == 0)
        return;
      m_axis = m_axis * (FT)1.0 / length;

      m_angle = acos(v1 * m_axis) + acos(v2 * m_axis) + acos(v3 * m_axis);
      m_angle /= (FT)3.0;
      if (m_angle < 0 || m_angle > M_PI / (FT)2.12)
        return;

      m_neg_sin_ang = -sin(m_angle);
      m_cos_ang = cos(m_angle);

      this->m_is_valid = true;
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) {
      for (std::size_t i = 0;i<indices.size();i++) {
          Vector to_apex = this->point(indices[i]) - m_apex;

          FT a = to_apex.squared_length();

          // projection on axis
          FT b = to_apex * m_axis;

          // distance to axis
          FT l = CGAL::sqrt(a - b * b);
          FT c = m_cos_ang * l;
          FT d = m_neg_sin_ang * b;

          // far on other side?
          dists[i] = (b < 0 && c - d < 0) ? a : abs(c + d) * abs(c + d);
        }
      }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                               std::vector<FT> &angles) const {
      for (std::size_t i = 0;i<indices.size();i++) {
          // construct vector orthogonal to axis in direction of the point
        Vector a = this->point(indices[i]) - m_apex;

          Vector b = CGAL::cross_product(m_axis, 
                                         CGAL::cross_product(m_axis, a));
          b = (a * b < 0) ? -b : b;
          FT length = CGAL::sqrt(b.squared_length());

          if (length == 0) {
            angles[i] = (FT)1.0;
            continue;
          }

          b = b * (FT)1.0 / length;
          b = m_cos_ang * b + m_neg_sin_ang * m_axis;

          angles[i] = CGAL::abs(this->normal(indices[i]) * b);
        }
      }

    virtual FT cos_to_normal(const Point &p, const Vector &n) const {
      // construct vector orthogonal to axis in direction of the point
      Vector a = p - m_apex;
      Vector b = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, a));
      b = (a * b < 0) ? -b : b;
      FT length = sqrt(b.squared_length());
      if (length == 0) {
        return (FT)1.0;
      }

      b = b * (FT)1.0 / length;
      b = m_cos_ang * b + m_neg_sin_ang * m_axis;

      return abs(n * b);
    }
      
    virtual std::size_t minimum_sample_size() const {
          return 3;
      }

    virtual bool supports_connected_component() const {
      return false;
    }

    // U is longitude
    virtual bool wraps_u() const {
      return true;
    }

    // V is between caps
    virtual bool wraps_v() const {
      return false;
    }

  private:
    FT m_angle;
    Point m_apex;
    Vector m_axis;
    FT m_neg_sin_ang, m_cos_ang;
      /// \endcond
  };
}
}
#endif
