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
#include <CGAL/number_utils.h>
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
    \tparam Traits a model of `EfficientRANSACTraits`
   \ingroup PkgPointSetShapeDetection3Shapes
   */

  template <class Traits>
  class Cone : public Shape_base<Traits> {
  public:
    /// \cond SKIP_IN_MANUAL
    typedef typename Traits::Point_map Point_map;
     ///< property map to access the location of an input point.
    typedef typename Traits::Normal_map Normal_map;
     ///< property map to access the unoriented normal of an input point.
    typedef typename Traits::FT FT; ///< number type.
    typedef typename Traits::Point_3 Point_3;///< point type.
    typedef typename Traits::Vector_3 Vector_3;///< vector type.
    /// \endcond

	
    Cone() : Shape_base<Traits>() {}
      
    /*!
      The opening angle between the axis and the surface of the cone.
     */
    FT angle() const {
        return m_angle;
    }
    
    /*!
      The apex of the cone.
     */
    Point_3 apex() const {
        return m_apex;
    }
    
    /*!
      The axis points from the apex into the cone.
     */
    Vector_3 axis() const {
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
    FT squared_distance(const Point_3 &p) const {
      Vector_3 toApex = p - m_apex;
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
      Point_3 p1 = this->point(indices[0]);
      Point_3 p2 = this->point(indices[1]);
      Point_3 p3 = this->point(indices[2]);

      Vector_3 n1 = this->normal(indices[0]);
      Vector_3 n2 = this->normal(indices[1]);
      Vector_3 n3 = this->normal(indices[2]);

      // first calculate intersection of three planes -> apex

      Vector_3 lineDir = CGAL::cross_product(n1, n2);
      FT length = CGAL::sqrt(lineDir.squared_length());
      if (length == 0)
        return;

      lineDir = lineDir * (FT)1.0 / length;

      // lineDir not normalized direction of intersection lines
      //  of two planes (p1, n1) and (p2, n2)
      // get point on line by moving point p1 onto line
      Vector_3 orthLineInPlane = CGAL::cross_product(n1, lineDir);
      length = CGAL::sqrt(orthLineInPlane.squared_length());
      if (length == 0)
        return;

      orthLineInPlane = orthLineInPlane * (FT)1.0 / length;

      // distance of p1 to (p2, n2)
      FT d = (p1 - CGAL::ORIGIN) * n2 - (p2 - CGAL::ORIGIN) * n2;
      // projection of orthLineInPlane onto p2
      FT l = orthLineInPlane * n2;
      if (l == 0)
        return;

      Point_3 pointOnLine = p1 - (d/l) * orthLineInPlane;

      // distance of pLineDir to (p3, n3)
      d = (pointOnLine - CGAL::ORIGIN) * n3 - (p3 - CGAL::ORIGIN) * n3;
      l = lineDir * n3;
      if (l == 0)
        return;

      m_apex = pointOnLine - (d/l) * lineDir;

      // 2. find axis
      Vector_3 v1 = p1 - m_apex;
      length = CGAL::sqrt(v1.squared_length());
      if (length == 0)
        return;
      v1 = v1 * (FT)1.0 / length;
      Point_3 c1 = m_apex + v1;

      Vector_3 v2 = p2 - m_apex;
      length = CGAL::sqrt(v2.squared_length());
      if (length == 0)
        return;
      v2 = v2 * (FT)1.0 / length;
      Point_3 c2 = m_apex + v2;

      Vector_3 v3 = p3 - m_apex;
      length = CGAL::sqrt(v3.squared_length());
      if (length == 0)
        return;
      v3 = v3 * (FT)1.0 / length;
      Point_3 c3 = m_apex + v3;

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
          Vector_3 to_apex = this->point(indices[i]) - m_apex;

          FT a = to_apex.squared_length();

          // projection on axis
          FT b = to_apex * m_axis;

          // distance to axis
          FT l = CGAL::sqrt(a - b * b);
          FT c = m_cos_ang * l;
          FT d = m_neg_sin_ang * b;

          // far on other side?
          dists[i] = 
            (b < 0 && c - d < 0) ? a : CGAL::abs(c + d) * CGAL::abs(c + d);
        }
      }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                               std::vector<FT> &angles) const {
      for (std::size_t i = 0;i<indices.size();i++) {
          // construct vector orthogonal to axis in direction of the point
        Vector_3 a = this->point(indices[i]) - m_apex;

          Vector_3 b = CGAL::cross_product(m_axis, 
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

    virtual FT cos_to_normal(const Point_3 &p, const Vector_3 &n) const {
      // construct vector orthogonal to axis in direction of the point
      Vector_3 a = p - m_apex;
      Vector_3 b = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, a));
      b = (a * b < 0) ? -b : b;
      FT length = CGAL::sqrt(b.squared_length());
      if (length == 0) {
        return (FT)1.0;
      }

      b = b * (FT)1.0 / length;
      b = m_cos_ang * b + m_neg_sin_ang * m_axis;

      return CGAL::abs(n * b);
    }
      
    virtual std::size_t minimum_sample_size() const {
          return 3;
      }

    virtual void parameters(const std::vector<std::size_t> &indices,
      std::vector<std::pair<FT, FT> > &parameterSpace,
      FT &cluster_epsilon,
      FT min[2],
      FT max[2]) const {

        // gap suchen? und dann shiften, um wrap zu vermeiden? vielleicht einfach später mit rad multiplizieren

    }

    virtual bool supports_connected_component() const {
      return false;
    }

  private:
    FT m_angle;
    Point_3 m_apex;
    Vector_3 m_axis;
    FT m_neg_sin_ang, m_cos_ang;
      /// \endcond
  };
}
}
#endif
