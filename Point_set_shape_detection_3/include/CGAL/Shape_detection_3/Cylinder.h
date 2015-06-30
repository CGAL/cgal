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

#ifndef CGAL_SHAPE_DETECTION_3_CYLINDER_H
#define CGAL_SHAPE_DETECTION_3_CYLINDER_H

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
 \file Cylinder.h
 */

namespace CGAL {
  namespace Shape_detection_3 {
    /*!
     \brief Cylinder implements Shape_base. The cylinder is represented
     by the axis, i.e. a point and direction, and the radius. The cylinder is
     unbounded, thus caps are not modelled.
     \tparam Traits a model of `EfficientRANSACTraits` with the additional 
             requirement that the type `Traits::Line_3` is provided.
     \ingroup PkgPointSetShapeDetection3Shapes
     */
  template <class Traits>
  class Cylinder : public Shape_base<Traits> {
  public:
    /// \cond SKIP_IN_MANUAL
    typedef typename Traits::Point_map Point_map;
     ///< property map to access the location of an input point.
    typedef typename Traits::Normal_map Normal_map;
     ///< property map to access the unoriented normal of an input point.
    typedef typename Traits::Vector_3 Vector_3; ///< vector type.
    typedef typename Traits::Point_3 Point_3; ///< point type.
    typedef typename Traits::FT FT; ///< number type.
    /// \endcond

    typedef typename Traits::Line_3 Line_3; ///< line type.

  public:
    Cylinder() : Shape_base<Traits>() {}

    /*!
      Axis of the cylinder.
     */
    Line_3 axis() const {
      return m_axis;
    }

    /*!
      Radius of the cylinder.
     */
    FT radius() const {
      return m_radius;
    }

    /// \cond SKIP_IN_MANUAL

    /*!
      Helper function to write axis and radius of the cylinder and number of assigned points into a string.
     */
    std::string info() const {
      std::stringstream sstr;
      Point_3 c = m_axis.point();
      Vector_3 a = m_axis.to_vector();

      sstr << "Type: cylinder center: (" << c.x() << ", " << c.y() << ", " << c.z() << ") axis: (" << a.x() << ", " << a.y() << ", " << a.z() << ") radius:" << m_radius
        << " #Pts: " <<  this->m_indices.size();

      return sstr.str();
    }

    /*!
      Computes squared Euclidean distance from query point to the shape.
      */ 
    FT squared_distance(const Point_3 &p) const {
      Vector_3 a = m_axis.to_vector();
      a = a * ((FT)1.0 / CGAL::sqrt(a.squared_length()));

      Vector_3 v = p - m_point_on_axis;
      v = v - ((v * a) * a);
      FT d = CGAL::sqrt(v.squared_length()) - m_radius;

      return d * d;
    }
    /// \endcond

  protected:
      /// \cond SKIP_IN_MANUAL
    virtual void create_shape(const std::vector<std::size_t> &indices) {
      Point_3 p1 = this->point(indices[0]);
      Point_3 p2 = this->point(indices[1]);

      Vector_3 n1 = this->normal(indices[0]);
      Vector_3 n2 = this->normal(indices[1]);

      Vector_3 axis = CGAL::cross_product(n1, n2);
      FT axisL = CGAL::sqrt(axis.squared_length());
      if (axisL < (FT)0.0001) {
        return;
      }
      axis = axis * (FT(1.0) / axisL);

      // establish two directions in the plane axis * x = 0, 
      // whereas xDir is the projected n1
      Vector_3 xDir = n1 - (n1 * axis) * axis;
      xDir = xDir * ((FT)1.0 / CGAL::sqrt(xDir.squared_length()));
      Vector_3 yDir = CGAL::cross_product(axis, xDir);
      yDir = yDir * ((FT)1.0 / CGAL::sqrt(yDir.squared_length()));

      FT n2x = n2 * yDir;
      FT n2y = -n2 * xDir;

      Vector_3 dist = p2 - p1;

      FT Ox = xDir * dist;
      FT Oy = yDir * dist;

      FT lineDist = n2x * Ox + n2y * Oy;

      m_radius = lineDist / n2x;
      m_point_on_axis = p1 + m_radius * xDir;
      m_radius = CGAL::abs(m_radius);

      m_axis = Line_3(m_point_on_axis, axis);

      if (squared_distance(p1) > this->m_epsilon ||
          (cos_to_normal(p1, n1) < this->m_normal_threshold))
        return;

      this->m_is_valid = true;
    }

    void parameters(const std::vector<std::size_t> &indices,
                    std::vector<std::pair<FT, FT> > &parameterSpace,                    
                    FT min[2],
                    FT max[2]) const {
      Vector_3 d1 = Vector_3((FT) 0, (FT) 0, (FT) 1);
      Vector_3 a = m_axis.to_vector();
      a = a * ((FT)1.0 / CGAL::sqrt(a.squared_length()));

      Vector_3 d2 = CGAL::cross_product(a, d1);
      FT l = d2.squared_length();
      if (l < (FT)0.0001) {
        d1 = Vector_3((FT) 1, (FT) 0, (FT) 0);
        d2 = CGAL::cross_product(m_axis.to_vector(), d1);
        l = d2.squared_length();
      }
      d2 = d2 / CGAL::sqrt(l);

      d1 = CGAL::cross_product(m_axis.to_vector(), d2);
      FT length = CGAL::sqrt(d1.squared_length());
      if (length == 0)
        return;

      d1 = d1 * (FT)1.0 / length;

      // first one separate for initializing min/max
      Vector_3 vec = this->point(indices[0]) - m_point_on_axis;
      FT v = vec * a;
      vec = vec - ((vec * a) * a);
      length = CGAL::sqrt(vec.squared_length());
      vec = vec * (FT)1.0 / length;

      FT a1 = acos(vec * d1);
      FT a2 = acos(vec * d2);

      FT u = FT((a2 < M_PI_2) ? 2 * M_PI - a1 : a1) * m_radius;

      parameterSpace[0] = std::pair<FT, FT>(u, v);

      min[1] = max[1] = v;

      for (std::size_t i = 0;i<indices.size();i++) {
        Vector_3 vec = this->point(indices[i]) - m_point_on_axis;
        FT v = vec * a;
        vec = vec - ((vec * a) * a);
        length = CGAL::sqrt(vec.squared_length());
        vec = vec * (FT)1.0 / length;

        FT a1 = acos(vec * d1);
        FT a2 = acos(vec * d2);

        FT u = FT((a2 < M_PI_2) ? 2 * M_PI - a1 : a1) * m_radius;
        min[1] = (std::min<FT>)(min[1], v);
        max[1] = (std::max<FT>)(max[1], v);

        parameterSpace[i] = std::pair<FT, FT>(u, v);
      }

      // Due to wrapping, the u parameter 'rotation around axis' always needs
      // to be the full extend.
      min[0] = 0;
      max[0] = FT(M_PI * 2.0 * m_radius);
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) {
      Vector_3 a = m_axis.to_vector();
      a = a * ((FT)1.0 / CGAL::sqrt(a.squared_length()));
      for (std::size_t i = 0;i<indices.size();i++) {
          Vector_3 v = this->point(indices[i]) - m_point_on_axis;

          v = v - ((v * a) * a);
          dists[i] = CGAL::sqrt(v.squared_length()) - m_radius;
          dists[i] = dists[i] * dists[i];
        }
      }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                              std::vector<FT> &angles) const {
      Vector_3 a = m_axis.to_vector();
      a = a * ((FT)1.0 / CGAL::sqrt(a.squared_length()));
      for (std::size_t i = 0;i<indices.size();i++) {
          Vector_3 v = this->point(indices[i]) - m_point_on_axis;

          v = v - ((v * a) * a);
          FT length = CGAL::sqrt(v.squared_length());
          if (length == 0) {
            angles[i] = (FT)1.0;
            continue;
          }
          v = v * (FT)1.0 / length;
          angles[i] = CGAL::abs(v * this->normal(indices[i]));
      }
    }

    FT cos_to_normal(const Point_3 &p, const Vector_3 &n) const {
      Vector_3 a = m_axis.to_vector();
      a = a * (FT)1.0 / CGAL::sqrt(a.squared_length());

      Vector_3 v = p - m_point_on_axis;
      v = v - ((v * a) * a);

      FT length = CGAL::sqrt(v.squared_length());
      if (length == 0)
       return (FT)1.0;

      v = v * (FT)1.0 / length;
      return CGAL::abs(v * n);
    }

    virtual std::size_t minimum_sample_size() const {
      return 2;
    }

    virtual bool supports_connected_component() const {
      return true;
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
    FT m_radius;
    Line_3 m_axis;
    Point_3 m_point_on_axis;
      
    /// \endcond
  };
}
}
#endif
