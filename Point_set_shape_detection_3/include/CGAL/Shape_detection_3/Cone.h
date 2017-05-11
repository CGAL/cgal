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

#include <CGAL/license/Point_set_shape_detection_3.h>


#include <CGAL/Shape_detection_3/Shape_base.h>
#include <CGAL/number_utils.h>
#include <cmath>

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
    using Shape_base<Traits>::update_label;

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


    Cone() : Shape_base<Traits>(), m_wrap(false) {}
      
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
        
        sstr << "Type: cone apex: (" << this->get_x(m_apex) << ", " << this->get_y(m_apex);
        sstr << ", " << this->get_z(m_apex) << ") axis: (" << this->get_x(m_axis) << ", ";
        sstr << this->get_y(m_axis) << ", " << this->get_z(m_axis) << ") angle:" << m_angle;
        sstr << " #Pts: " << this->m_indices.size();
        
        return sstr.str();
    }

    /*!
    Computes squared Euclidean distance from query point to the shape.
    */ 
    FT squared_distance(const Point_3 &p) const {
      Vector_3 toApex = this->constr_vec(m_apex, p);
      FT a = this->sqlen(toApex);

      // projection on axis
      FT b = this->scalar_pdct(toApex, m_axis);

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

      Vector_3 lineDir = this->cross_pdct(n1, n2);
      FT length = CGAL::sqrt(this->sqlen(lineDir));
      if (length == 0)
        return;

      lineDir = this->scale(lineDir, (FT)1.0 / length);

      // lineDir not normalized direction of intersection lines
      //  of two planes (p1, n1) and (p2, n2)
      // get point on line by moving point p1 onto line
      Vector_3 orthLineInPlane = this->cross_pdct(n1, lineDir);
      length = CGAL::sqrt(this->sqlen(orthLineInPlane));
      if (length == 0)
        return;

      orthLineInPlane = this->scale(orthLineInPlane, (FT)1.0 / length);

      // distance of p1 to (p2, n2)
      FT d = this->scalar_pdct(this->constr_vec(CGAL::ORIGIN, p1), n2)
        - this->scalar_pdct(this->constr_vec(CGAL::ORIGIN, p2), n2);
      // projection of orthLineInPlane onto p2
      FT l = this->scalar_pdct(orthLineInPlane, n2);
      if (l == 0)
        return;

      Point_3 pointOnLine = this->transl(
        p1, this->scale(orthLineInPlane, -d/l));

      // distance of pLineDir to (p3, n3)
      d = this->scalar_pdct(this->constr_vec(CGAL::ORIGIN, pointOnLine), n3)
        - this->scalar_pdct(this->constr_vec(CGAL::ORIGIN, p3), n3);
      l = this->scalar_pdct(lineDir, n3);
      if (l == 0)
        return;

      m_apex = this->transl(pointOnLine, this->scale(lineDir, -d/l));

      // 2. find axis
      Vector_3 v1 = this->constr_vec(m_apex, p1);
      length = CGAL::sqrt(this->sqlen(v1));
      if (length == 0)
        return;
      v1 = this->scale(v1, (FT)1.0 / length);
      Point_3 c1 = this->transl(m_apex, v1);

      Vector_3 v2 = this->constr_vec(m_apex, p2);
      length = CGAL::sqrt(this->sqlen(v2));
      if (length == 0)
        return;
      v2 = this->scale(v2, (FT)1.0 / length);
      Point_3 c2 = this->transl(m_apex, v2);

      Vector_3 v3 = this->constr_vec(m_apex, p3);
      length = CGAL::sqrt(this->sqlen(v3));
      if (length == 0)
        return;
      v3 = this->scale(v3, (FT)1.0 / length);
      Point_3 c3 = this->transl(m_apex, v3);

      m_axis = this->cross_pdct(this->constr_vec(c2, c1), this->constr_vec(c3, c1));
      m_axis = (this->scalar_pdct(orthLineInPlane, m_axis) < 0) ? 
        this->scale(m_axis, FT(-1)) : m_axis;
      length = CGAL::sqrt(this->sqlen(m_axis));
      if (length == 0)
        return;
      m_axis = this->scale(m_axis, (FT)1.0 / length);

      m_angle = acos(this->scalar_pdct(v1, m_axis)) + acos(this->scalar_pdct(v2, m_axis)) + acos(this->scalar_pdct(v3, m_axis));
      m_angle /= (FT)3.0;
      if (m_angle < 0 || m_angle > CGAL_PI / (FT)2.12)
        return;

      m_neg_sin_ang = -sin(m_angle);
      m_cos_ang = cos(m_angle);

      this->m_is_valid = true;
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) const {
      for (std::size_t i = 0;i<indices.size();i++) {
          Vector_3 to_apex = this->constr_vec(m_apex, this->point(indices[i]));

          FT a = this->sqlen(to_apex);

          // projection on axis
          FT b = this->scalar_pdct(to_apex, m_axis);

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
        Vector_3 a = this->constr_vec(m_apex, this->point(indices[i]));

          Vector_3 b = this->cross_pdct(m_axis, 
                                         this->cross_pdct(m_axis, a));
          b = (this->scalar_pdct(a, b) < 0) ? this->scale(b, FT(-1)) : b;
          FT length = CGAL::sqrt(this->sqlen(b));

          if (length == 0) {
            angles[i] = (FT)1.0;
            continue;
          }

          b = this->scale(b, (FT)1.0 / length);
          b = this->sum_vectors(
            this->scale(b, m_cos_ang), 
            this->scale(m_axis, m_neg_sin_ang));

          angles[i] = CGAL::abs(this->scalar_pdct(this->normal(indices[i]), b));
        }
      }

    virtual FT cos_to_normal(const Point_3 &p, const Vector_3 &n) const {
      // construct vector orthogonal to axis in direction of the point
      Vector_3 a = this->constr_vec(m_apex, p);
      Vector_3 b = this->cross_pdct(m_axis, this->cross_pdct(m_axis, a));
      b = (this->scalar_pdct(a, b) < 0) ? this->scale(b, FT(-1)) : b;
      FT length = CGAL::sqrt(this->sqlen(b));
      if (length == 0) {
        return (FT)1.0;
      }

      b = this->scale(b, (FT)1.0 / length);
      b = this->sum_vectors(
        this->scale(b, m_cos_ang), 
        this->scale(m_axis, m_neg_sin_ang));

      return CGAL::abs(this->scalar_pdct(n, b));
    }
      
    virtual std::size_t minimum_sample_size() const {
          return 3;
      }

    virtual void post_wrap(const std::vector<unsigned int> &bitmap,
      const std::size_t &u_extent,
      const std::size_t &v_extent,
      std::vector<unsigned int> &labels) const {
      if (!m_wrap)
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

    virtual void parameters(const std::vector<std::size_t> &indices,
      std::vector<std::pair<FT, FT> > &parameterSpace,
      FT &cluster_epsilon,
      FT min[2],
      FT max[2]) const {

      // Create basis d1, d2
      Vector_3 d1 = this->constr_vec(
        ORIGIN, this->constr_pt(FT(0), FT(0), FT(1)));
      Vector_3 d2 = this->cross_pdct(m_axis, d1);
      FT l = this->sqlen(d2);
      if (l < (FT)0.0001) {
        d1 = this->constr_vec(ORIGIN, this->constr_pt(FT(1), FT(0), FT(0)));
        d2 = this->cross_pdct(m_axis, d1);
        l = this->sqlen(d2);
      }
      d2 = this->scale(d2, FT(1) / CGAL::sqrt(l));

      d1 = this->cross_pdct(m_axis, d2);
      l = CGAL::sqrt(this->sqlen(d1));
      if (l == 0)
        return;

      d1 = this->scale(d1, (FT)1.0 / l);

      if (m_angle > CGAL_M_PI_4) {
        // Projection onto a disk preserving distance to apex

        m_wrap = false;

        // First index separately to initialize min/max
        Vector_3 d = this->constr_vec(m_apex, this->point(indices[0]));
        FT l = this->scalar_pdct(d, m_axis) / m_cos_ang;
        FT u = this->scalar_pdct(d, d1);
        FT v = this->scalar_pdct(d, d2);
        FT l2 = CGAL::sqrt(u * u + v * v);
        u = u * l/l2;
        v = v * l/l2;
        min[0] = max[0] = u;
        min[1] = max[1] = v;
        parameterSpace[0] = std::pair<FT, FT>(u, v);

        for (std::size_t i = 1;i<indices.size();i++) {
          d = this->constr_vec(m_apex, this->point(indices[i]));
          l = this->scalar_pdct(d, m_axis) / m_cos_ang;
          u = this->scalar_pdct(d, d1);
          v = this->scalar_pdct(d, d2);
          l2 = CGAL::sqrt(u * u + v * v);
          u = u * l/l2;
          v = v * l/l2;

          min[0] = (std::min<FT>)(min[0], u);
          max[0] = (std::max<FT>)(max[0], u);

          min[1] = (std::min<FT>)(min[1], v);
          max[1] = (std::max<FT>)(max[1], v);

          parameterSpace[i] = std::pair<FT, FT>(u, v);
        }
      }
      else {
        // Map onto triangle.
        // u coordinate is arclength
        // v coordinate is distance to apex

        Vector_3 d = this->constr_vec(m_apex, this->point(indices[0]));
        FT v = this->scalar_pdct(d, m_axis) / m_cos_ang;
        FT phi = atan2(this->scalar_pdct(d, d2), this->scalar_pdct(d, d1));
        FT u = FT(phi + CGAL_PI);
        FT avg_v = v;

        min[0] = max[0] = u;
        min[1] = max[1] = v;
        parameterSpace[0] = std::pair<FT, FT>(u, v);
        for (std::size_t i = 1;i<indices.size();i++) {
          d = this->constr_vec(m_apex, this->point(indices[i]));
          v = this->scalar_pdct(d, m_axis) / m_cos_ang;
          phi = atan2(this->scalar_pdct(d, d2), this->scalar_pdct(d, d1));
          u = FT(phi + CGAL_PI);

          min[0] = (std::min<FT>)(min[0], u);
          max[0] = (std::max<FT>)(max[0], u);

          min[1] = (std::min<FT>)(min[1], v);
          max[1] = (std::max<FT>)(max[1], v);

          avg_v += v;

          parameterSpace[i] = std::pair<FT, FT>(u, v);
        }

        // Scale u parameter by average circumference to arc length
        avg_v /= indices.size();
        const FT scale = -m_neg_sin_ang * avg_v;

        m_wrap = (min[0] + 2 * CGAL_PI - max[0]) * scale < cluster_epsilon;

        for (std::size_t i = 0;i<parameterSpace.size();i++) {
          std::pair<FT, FT> p = parameterSpace[i];
          parameterSpace[i] = std::pair<FT, FT>(p.first * scale, p.second);
        }

        min[0] *= scale;
        max[0] *= scale;
      }
    }

    virtual bool supports_connected_component() const {
      return true;
    }

  private:
    FT m_angle;
    Point_3 m_apex;
    Vector_3 m_axis;
    FT m_neg_sin_ang, m_cos_ang;
    mutable bool m_wrap;
      /// \endcond
  };
}
}
#endif
