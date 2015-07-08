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
      this->m_wrap_u = true;
    }

    virtual void parameters(const std::vector<std::size_t> &indices,
                            std::vector<std::pair<FT, FT> > &parameterSpace,
                            FT &cluster_epsilon,
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

      FT a1 = vec * d1;
      a1 = (a1 < (FT) -1.0) ? (FT) -1.0 : ((a1 > (FT) 1.0) ? (FT) 1.0 : a1);
      a1 = acos(a1);
      FT a2 = vec * d2;
      a2 = (a2 < (FT) -1.0) ? (FT) -1.0 : ((a2 > (FT) 1.0) ? (FT) 1.0 : a2);
      a2 = acos(a2);

      FT u = FT((a2 < CGAL_M_PI_2) ? 2 * CGAL_PI - a1 : a1) * m_radius;

      parameterSpace[0] = std::pair<FT, FT>(u, v);

      min[0] = max[0] = u;
      min[1] = max[1] = v;

      for (std::size_t i = 0;i<indices.size();i++) {
        vec = this->point(indices[i]) - m_point_on_axis;
        v = vec * a;
        vec = vec - ((vec * a) * a);
        length = CGAL::sqrt(vec.squared_length());
        vec = vec * (FT)1.0 / length;

        a1 = vec * d1;
        a1 = (a1 < (FT) -1.0) ? (FT) -1.0 : ((a1 > (FT) 1.0) ? (FT) 1.0 : a1);
        a1 = acos(a1);
        a2 = vec * d2;
        a2 = (a2 < (FT) -1.0) ? (FT) -1.0 : ((a2 > (FT) 1.0) ? (FT) 1.0 : a2);
        a2 = acos(a2);

        u = FT((a2 < CGAL_M_PI_2) ? 2 * CGAL_PI - a1 : a1) * m_radius;
        min[0] = (std::min<FT>)(min[0], u);
        max[0] = (std::max<FT>)(max[0], u);

        min[1] = (std::min<FT>)(min[1], v);
        max[1] = (std::max<FT>)(max[1], v);

        parameterSpace[i] = std::pair<FT, FT>(u, v);
      }

      // Is close to wrapping around?
      FT diff_to_full_range = min[0] + FT(CGAL_PI * 2.0 * m_radius) - max[0];
      if (diff_to_full_range < cluster_epsilon) {
        m_wrap_u = true;
        FT frac = (max[0] - min[0]) / cluster_epsilon;
        FT trunc = floor(frac);
        frac = frac - trunc;

        if (frac < (FT) 0.5) {
          cluster_epsilon = (max[0] - min[0]) / (trunc - (FT) 0.01);
        }
      }
      else m_wrap_u = false;
    }
    
    // The u coordinate corresponds to the rotation around the axis and
    // therefore needs to be wrapped around.
    virtual void post_wrap(const std::vector<unsigned int> &bitmap,
      const std::size_t &u_extent,
      const std::size_t &v_extent,
      std::vector<unsigned int> &labels) const {
        if (!m_wrap_u)
          return;

        // handle top index separately
        unsigned int nw = bitmap[u_extent - 1];
        unsigned int l = bitmap[0];

        // Special case v_extent is just 1
        if (v_extent == 1) {
          if (nw && nw != l)
            update_label(labels, (std::max<unsigned int>)(nw, l), l = (std::min<unsigned int>)(nw, l));

          return;
        }

        unsigned int w = bitmap[2 * u_extent - 1];
        unsigned int sw;

        if (l) {
          if (nw && nw != l)
            update_label(labels, (std::max<unsigned int>)(nw, l), l = (std::min<unsigned int>)(nw, l));
          else if (w && w != l)
            update_label(labels, (std::max<unsigned int>)(w, l), l = (std::min<unsigned int>)(w, l));
        }
        
        // handle mid indices
        for (std::size_t y = 1;y<v_extent - 1;y++) {
          l = bitmap[y * u_extent];
          if (!l)
            continue;

          nw = bitmap[y * u_extent - 1];
          w = bitmap[(y + 1) * u_extent - 1];
          sw = bitmap[(y + 2) * u_extent - 1];

          if (nw && nw != l)
            update_label(labels, (std::max<unsigned int>)(nw, l), l = (std::min<unsigned int>)(nw, l));
          if (w && w != l)
            update_label(labels, (std::max<unsigned int>)(w, l), l = (std::min<unsigned int>)(w, l));
          else if (sw && sw != l)
            update_label(labels, (std::max<unsigned int>)(sw, l), l = (std::min<unsigned int>)(sw, l));
        }

        // handle last index
        l = bitmap[(v_extent - 1) * u_extent];
        if (!l)
          return;

        nw = bitmap[(v_extent - 1) * u_extent - 1];
        w = bitmap[u_extent * v_extent - 1];

        if (nw && nw != l)
          update_label(labels, (std::max<unsigned int>)(nw, l), l = (std::min<unsigned int>)(nw, l));
        else if (w && w != l)
          update_label(labels, (std::max<unsigned int>)(w, l), l = (std::min<unsigned int>)(w, l));
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) const {
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

  private:
    FT m_radius;
    Line_3 m_axis;
    Point_3 m_point_on_axis;
    mutable bool m_wrap_u;
      
    /// \endcond
  };
}
}
#endif
