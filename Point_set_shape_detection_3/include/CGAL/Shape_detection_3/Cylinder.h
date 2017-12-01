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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_3_CYLINDER_H
#define CGAL_SHAPE_DETECTION_3_CYLINDER_H

#include <CGAL/license/Point_set_shape_detection_3.h>


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
             requirement for cylinders (see `EfficientRANSACTraits` documentation).
     \ingroup PkgPointSetShapeDetection3Shapes
     */
  template <class Traits>
  class Cylinder : public Shape_base<Traits> {
    using Shape_base<Traits>::update_label;

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
      Point_3 c = this->constr_point_on(m_axis);
      Vector_3 a = this->constr_vec(m_axis);

      sstr << "Type: cylinder center: (" << this->get_x(c) << ", " << this->get_y(c) << ", " << this->get_z(c) << ") axis: (" << this->get_x(a) << ", " << this->get_y(a) << ", " << this->get_z(a) << ") radius:" << m_radius
        << " #Pts: " <<  this->m_indices.size();

      return sstr.str();
    }

    /*!
      Computes squared Euclidean distance from query point to the shape.
      */ 
    FT squared_distance(const Point_3 &p) const {
      Vector_3 a = this->constr_vec(m_axis);
      a = this->scale(a, (FT)1.0 / CGAL::sqrt(this->sqlen(a)));

      Vector_3 v = this->constr_vec(m_point_on_axis, p);
      v = this->sum_vectors(v, this->scale(a, -this->scalar_pdct(v, a)));
      FT d = CGAL::sqrt(this->sqlen(v)) - m_radius;

      return d * d;
    }
    /// \endcond

  protected:
      /// \cond SKIP_IN_MANUAL
    
    // ------------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------------
    using Shape_base<Traits>::constr_vec;
    Vector_3 constr_vec(const Line_3& l) const
    { return this->m_traits.construct_vector_3_object()(l); }
    Point_3 constr_point_on(const Line_3& l) const
    { return this->m_traits.construct_point_on_3_object()(l, 0); }

    virtual void create_shape(const std::vector<std::size_t> &indices) {
      Point_3 p1 = this->point(indices[0]);
      Point_3 p2 = this->point(indices[1]);

      Vector_3 n1 = this->normal(indices[0]);
      Vector_3 n2 = this->normal(indices[1]);

      Vector_3 axis = this->cross_pdct(n1, n2);
      FT axisL = CGAL::sqrt(this->sqlen(axis));
      if (axisL < (FT)0.0001) {
        return;
      }
      axis = this->scale(axis, FT(1.0) / axisL);

      // establish two directions in the plane axis * x = 0, 
      // whereas xDir is the projected n1
      Vector_3 xDir = this->sum_vectors(
        n1, this->scale(axis, -this->scalar_pdct(n1, axis)));
      xDir = this->scale(xDir, (FT)1.0 / CGAL::sqrt(this->sqlen(xDir)));
      Vector_3 yDir = this->cross_pdct(axis, xDir);
      yDir = this->scale(yDir, (FT)1.0 / CGAL::sqrt(this->sqlen(yDir)));

      FT n2x = this->scalar_pdct(n2, yDir);
      FT n2y = -this->scalar_pdct(n2, xDir);

      Vector_3 dist = this->constr_vec(p1, p2);

      FT Ox = this->scalar_pdct(xDir, dist);
      FT Oy = this->scalar_pdct(yDir, dist);

      FT lineDist = n2x * Ox + n2y * Oy;

      m_radius = lineDist / n2x;
      m_point_on_axis = this->transl(p1, this->scale(xDir, m_radius));
      m_radius = CGAL::abs(m_radius);

      m_axis = this->m_traits.construct_line_3_object()(m_point_on_axis, axis);

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
      Vector_3 d1 = this->constr_vec(
        ORIGIN, this->constr_pt(FT(0), FT(0), FT(1)));
      Vector_3 a = this->constr_vec(m_axis);
      a = this->scale(a, (FT)1.0 / CGAL::sqrt(this->sqlen(a)));

      Vector_3 d2 = this->cross_pdct(a, d1);
      FT l = this->sqlen(d2);
      if (l < (FT)0.0001) {
        d1 = this->constr_vec(ORIGIN, this->constr_pt(FT(1), FT(0), FT(0)));
        d2 = this->cross_pdct(this->constr_vec(m_axis), d1);
        l = this->sqlen(d2);
      }
      d2 = this->scale(d2, FT(1) / CGAL::sqrt(l));

      d1 = this->cross_pdct(this->constr_vec(m_axis), d2);
      FT length = CGAL::sqrt(this->sqlen(d1));
      if (length == 0)
        return;

      d1 = this->scale(d1, (FT)1.0 / length);

      // first one separate for initializing min/max
      Vector_3 vec = this->constr_vec(m_point_on_axis, this->point(indices[0]));
      FT v = this->scalar_pdct(vec, a);
      vec = this->sum_vectors(vec, this->scale(a, -this->scalar_pdct(vec, a)));
      length = CGAL::sqrt(this->sqlen(vec));
      vec = this->scale(vec, (FT)1.0 / length);

      FT a1 = this->scalar_pdct(vec, d1);
      a1 = (a1 < (FT) -1.0) ? (FT) -1.0 : ((a1 > (FT) 1.0) ? (FT) 1.0 : a1);
      a1 = acos(a1);
      FT a2 = this->scalar_pdct(vec, d2);
      a2 = (a2 < (FT) -1.0) ? (FT) -1.0 : ((a2 > (FT) 1.0) ? (FT) 1.0 : a2);
      a2 = acos(a2);

      FT u = FT((a2 < CGAL_M_PI_2) ? 2 * CGAL_PI - a1 : a1) * m_radius;

      parameterSpace[0] = std::pair<FT, FT>(u, v);

      min[0] = max[0] = u;
      min[1] = max[1] = v;

      for (std::size_t i = 0;i<indices.size();i++) {
        vec = this->constr_vec(m_point_on_axis, this->point(indices[i]));
        v = this->scalar_pdct(vec, a);
        vec = this->sum_vectors(vec, this->scale(a, -this->scalar_pdct(vec, a)));
        length = CGAL::sqrt(this->sqlen(vec));
        vec = this->scale(vec, (FT)1.0 / length);

        a1 = this->scalar_pdct(vec, d1);
        a1 = (a1 < (FT) -1.0) ? (FT) -1.0 : ((a1 > (FT) 1.0) ? (FT) 1.0 : a1);
        a1 = acos(a1);
        a2 = this->scalar_pdct(vec, d2);
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

        if (frac < 1)
          return;

        FT trunc = floor(frac);
        frac = frac - trunc;

        if (frac < (FT) 0.5) {
          cluster_epsilon = (max[0] - min[0]) / (trunc * FT(0.99999));
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
          if (nw && nw != l) {
            l = (std::min<unsigned int>)(nw, l);
            update_label(labels, (std::max<unsigned int>)(nw, l), l);
          }

          return;
        }

        unsigned int w = bitmap[2 * u_extent - 1];
        unsigned int sw;

        if (l) {
          if (nw && nw != l) {
            l = (std::min<unsigned int>)(nw, l);
            update_label(labels, (std::max<unsigned int>)(nw, l), l);
          }
          else if (w && w != l) {
            l = (std::min<unsigned int>)(w, l);
            update_label(labels, (std::max<unsigned int>)(w, l), l);
          }
        }
        
        // handle mid indices
        for (std::size_t y = 1;y<v_extent - 1;y++) {
          l = bitmap[y * u_extent];
          if (!l)
            continue;

          nw = bitmap[y * u_extent - 1];
          w = bitmap[(y + 1) * u_extent - 1];
          sw = bitmap[(y + 2) * u_extent - 1];

          if (nw && nw != l) {
            l = (std::min<unsigned int>)(nw, l);
            update_label(labels, (std::max<unsigned int>)(nw, l), l);
          }
          if (w && w != l) {
            l = (std::min<unsigned int>)(w, l);
            update_label(labels, (std::max<unsigned int>)(w, l), l);
          }
          else if (sw && sw != l) {
            l = (std::min<unsigned int>)(sw, l);
            update_label(labels, (std::max<unsigned int>)(sw, l), l);
          }
        }

        // handle last index
        l = bitmap[(v_extent - 1) * u_extent];
        if (!l)
          return;

        nw = bitmap[(v_extent - 1) * u_extent - 1];
        w = bitmap[u_extent * v_extent - 1];

        if (nw && nw != l) {
          l = (std::min<unsigned int>)(nw, l);
          update_label(labels, (std::max<unsigned int>)(nw, l), l);
        }
        else if (w && w != l) {
          l = (std::min<unsigned int>)(w, l);
          update_label(labels, (std::max<unsigned int>)(w, l), l);
        }
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) const {
      Vector_3 a = this->constr_vec(m_axis);
      a = this->scale(a, (FT)1.0 / CGAL::sqrt(this->sqlen(a)));
      for (std::size_t i = 0;i<indices.size();i++) {
        Vector_3 v = this->constr_vec(m_point_on_axis, this->point(indices[i]));

        v = this->sum_vectors(v, this->scale(a, -this->scalar_pdct(v, a)));
        dists[i] = CGAL::sqrt(this->sqlen(v)) - m_radius;
          dists[i] = dists[i] * dists[i];
        }
      }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                              std::vector<FT> &angles) const {
      Vector_3 a = this->constr_vec(m_axis);
      a = this->scale(a, (FT)1.0 / CGAL::sqrt(this->sqlen(a)));
      for (std::size_t i = 0;i<indices.size();i++) {
        Vector_3 v = this->constr_vec(m_point_on_axis, this->point(indices[i]));

        v = this->sum_vectors(v, this->scale(a, -this->scalar_pdct(v, a)));
        FT length = CGAL::sqrt(this->sqlen(v));
          if (length == 0) {
            angles[i] = (FT)1.0;
            continue;
          }
        v = this->scale(v, (FT)1.0 / length);
        angles[i] = CGAL::abs(this->scalar_pdct(v, this->normal(indices[i])));
      }
    }

    FT cos_to_normal(const Point_3 &p, const Vector_3 &n) const {
      Vector_3 a = this->constr_vec(m_axis);
      a = this->scale(a, (FT)1.0 / CGAL::sqrt(this->sqlen(a)));

      Vector_3 v = this->constr_vec(m_point_on_axis, p);
      v = this->sum_vectors(v, this->scale(a, -this->scalar_pdct(v, a)));

      FT length = CGAL::sqrt(this->sqlen(v));
      if (length == 0)
       return (FT)1.0;

      v = this->scale(v, (FT)1.0 / length);
      return CGAL::abs(this->scalar_pdct(v, n));
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
