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

#ifndef CGAL_SHAPE_DETECTION_3_SPHERE_H
#define CGAL_SHAPE_DETECTION_3_SPHERE_H

#include <CGAL/Shape_detection_3/Shape_base.h>
#include <CGAL/number_utils.h>
#include <cmath>


/*!
 \file Sphere.h
 */

namespace CGAL {
  namespace Shape_detection_3 {
    /*!
     \ingroup PkgPointSetShapeDetection3Shapes
     \brief Sphere implements Shape_base. The sphere is represented by its center and the radius.
     \tparam Traits a model of `EfficientRANSACTraits` with the additional 
             requirement for spheres (see `EfficientRANSACTraits` documentation).
     */
  template <class Traits>
  class Sphere : public Shape_base<Traits> {
    using Shape_base<Traits>::update_label;

  public:    
    /// \cond SKIP_IN_MANUAL
    typedef typename Traits::Point_map Point_map;
      ///< property map to access the location of an input point.
    typedef typename Traits::Normal_map Normal_map;
      ///< property map to access the unoriented normal of an input point.
    typedef typename Traits::Vector_3 Vector_3;
      ///< vector type.
    typedef typename Traits::Sphere_3 Sphere_3;
      ///< sphere type.
    typedef typename Traits::FT FT;
      ///< number type.
    typedef typename Traits::Point_3 Point_3;
      ///< point type.
    /// \endcond

  public:
    Sphere() : Shape_base<Traits>() {}

    /*!
      Conversion operator to convert to `Sphere_3` type.
     */
    operator Sphere_3() const {
      return m_sphere;
    }

    /*!
      Access to the center.
     */
    Point_3 center() const {
      return this->sph_center(m_sphere);
    }
      
    /*!
      Access to the radius of the sphere.
     */
    FT radius() const {
      return CGAL::sqrt(this->sqradius(m_sphere));
    }

    /// \cond SKIP_IN_MANUAL
    /*!
      Computes the squared Euclidean distance from query point to the shape.
      */
    FT squared_distance(const Point_3 &p) const {
      const FT d = CGAL::sqrt(
        this->sqlen(this->constr_vec(
          p, this->sph_center(m_sphere)))) 
        - CGAL::sqrt(this->sqradius(m_sphere));
      return d * d;
    }
    
    /*!
      Helper function to write center, 
      radius of the sphere and number of assigned points into a string.
     */
    std::string info() const {
      std::stringstream sstr;
      Point_3 c = this->sph_center(m_sphere);
      FT r = CGAL::sqrt(this->sqradius(m_sphere));

      sstr << "Type: sphere center: (" << this->get_x(c) << ", " << this->get_y(c);
      sstr << ", " << this->get_z(c) << ") radius:" << r;
      sstr << " #Pts: " <<  this->m_indices.size();

      return sstr.str();
    }
    /// \endcond
  protected:
      /// \cond SKIP_IN_MANUAL

    // ------------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------------
    Sphere_3 constr_sphere(const Point_3& c, FT r) const
    { return this->m_traits.construct_sphere_3_object()(c, r); }
    Point_3 sph_center(const Sphere_3& s) const
    { return this->m_traits.construct_center_3_object()(s); }
    FT sqradius(const Sphere_3& s) const
    { return this->m_traits.compute_squared_radius_3_object()(s); }

    void create_shape(const std::vector<std::size_t> &indices) {
      Point_3 p1 = this->point(indices[0]);
      Point_3 p2 = this->point(indices[1]);
      Point_3 p3 = this->point(indices[2]);

      Vector_3 n1 = this->normal(indices[0]);
      Vector_3 n2 = this->normal(indices[1]);
      Vector_3 n3 = this->normal(indices[2]);


      // Determine center: select midpoint of shortest line segment
      // between p1 and p2. Implemented from "3D game engine design" by Eberly 2001
      Vector_3 diff = this->constr_vec(p2, p1);
      FT a = this->scalar_pdct(n1, n1);
      FT b = -this->scalar_pdct(n1, n2);
      FT c = this->scalar_pdct(n2, n2);
      FT d = this->scalar_pdct(n1, diff);

      FT det = CGAL::abs(a * c - b * b);

      // degenerated when nearly parallel
      if (det < (FT)0.00001) {
        this->m_is_valid = false;
        return;
      }

      FT e = -this->scalar_pdct(n2, diff);
      FT invDet = (FT) 1.0 / det;
      FT s = (b * e - c * d) * invDet;
      FT t = (d * b - a * e) * invDet;

      Vector_3 v_transl = this->sum_vectors(
        this->constr_vec(CGAL::ORIGIN, this->transl(p1, this->scale(n1, s))),
        this->constr_vec(CGAL::ORIGIN, this->transl(p2, this->scale(n2, t))));
      Point_3 center = this->transl(
        CGAL::ORIGIN, this->scale(v_transl, (FT)0.5));

      Vector_3 v1 = (this->constr_vec(center, p1));
      Vector_3 v2 = (this->constr_vec(center, p2));
      FT d1 = CGAL::sqrt(this->sqlen(v1));
      FT d2 = CGAL::sqrt(this->sqlen(v2));

      if (CGAL::abs(d1 - d2) > (FT)2.0 * this->m_epsilon) {
        this->m_is_valid = false;
        return;
      }

      v1 = this->scale(v1, (FT)1.0 / d1);
      v2 = this->scale(v2, (FT)1.0 / d2);

      if (this->scalar_pdct(n1, v1) < this->m_normal_threshold ||
          this->scalar_pdct(n2, v2) < this->m_normal_threshold) {
        this->m_is_valid = false;
        return;
      }

      Vector_3 v3 = this->constr_vec(center, p3);
      FT d3 = CGAL::sqrt(this->sqlen(v3));
      v3 = this->scale(v3, (FT)1.0 / d3);

      FT radius = (d1 + d2) * (FT)0.5;

      if (CGAL::abs(d3 - radius) > this->m_epsilon ||
          this->scalar_pdct(n3, v3) < this->m_normal_threshold) {
        this->m_is_valid = false;
        return;
      }

      this->m_is_valid = true;

      m_sphere = this->constr_sphere(center, radius * radius);
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) const {

      FT radius = CGAL::sqrt(this->sqradius(m_sphere));

      for (std::size_t i = 0;i<indices.size();i++) {
        dists[i] = CGAL::sqrt(this->sqlen(this->constr_vec(
          this->sph_center(m_sphere), this->point(indices[i]))))
          - radius;

        dists[i] = dists[i] * dists[i];
      }
    }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                               std::vector<FT> &angles) const {
      for (std::size_t i = 0;i<indices.size();i++) {
        Vector_3 n = this->constr_vec(
          this->point(indices[i]),
          this->sph_center(m_sphere));

        FT length = CGAL::sqrt(this->sqlen(n));
        if (length == 0) {
          angles[i] = (FT)1.0;
          continue;
        }

        n = this->scale(n, (FT)1.0 / length);
        angles[i] = CGAL::abs(this->scalar_pdct(this->normal(indices[i]), n));
      }
    }

    virtual FT cos_to_normal(const Point_3 &p, const Vector_3 &n) const {
      Vector_3 sphere_normal = this->constr_vec(p, this->sph_center(m_sphere));
      FT length = (FT)(CGAL::sqrt(this->sqlen(n)));
      if (length == 0)
        return 1;

      sphere_normal = this->scale(sphere_normal, (FT)1.0 / length);
      return CGAL::abs(this->scalar_pdct(sphere_normal, n));
    }
      
    virtual std::size_t minimum_sample_size() const {
      return 3;
    }

    // Maps to the range [-1,1]^2
    static void concentric_mapping(FT phi, FT proj, FT rad, FT &x, FT &y) {
      phi = (phi < FT(-CGAL_M_PI_4)) ? phi + FT(2 * CGAL_PI) : phi;
      proj = (proj < FT(-1.0)) ? FT(-1.0) : ((proj > FT(1.0)) ? FT(1.0) : proj);
      FT r = FT(acos(double(CGAL::abs(proj)))) / FT(CGAL_M_PI_2);

      FT a = 0, b = 0;
      if (phi < FT(CGAL_M_PI_4)) {
        a = r;
        b = phi * r / FT(CGAL_M_PI_4);
      }
      else if (phi < FT(3.0 * CGAL_M_PI_4)) {
        a = -FT(phi - CGAL_M_PI_2) * r / FT(CGAL_M_PI_4);
        b = r;
      }
      else if (phi < FT(5.0 * CGAL_M_PI_4)) {
        a = -r;
        b = (phi - FT(CGAL_PI)) * (-r) / FT(CGAL_M_PI_4);
      }
      else {
        a = (phi - 3 * FT(CGAL_M_PI_2)) * r / FT(CGAL_M_PI_4);
        b = -r;
      }

      x = a;
      y = b;

      // Map into hemisphere
      if (proj >= 0)
        y += 1;
      else
        y = -1 - y;

      // Scale to surface distance
      x = FT(x * CGAL_M_PI_2 * rad);
      y = FT(y * CGAL_M_PI_2 * rad);
    }

    virtual void parameters(const std::vector<std::size_t> &indices,
                            std::vector<std::pair<FT, FT> > &parameterSpace,
                            FT &cluster_epsilon,
                            FT min[2],
                            FT max[2]) const {
      Vector_3 axis = this->constr_vec();
      FT rad = radius();
      // Take average normal as axis
      for (std::size_t i = 0;i<indices.size();i++)
        axis = this->sum_vectors(axis, this->normal(indices[i]));
      axis = this->scale(axis, FT(1) / CGAL::sqrt(this->sqlen(axis)));

      // create basis d1, d2
      Vector_3 d1 = this->constr_vec(
        ORIGIN, this->constr_pt(FT(0), FT(0), FT(1)));
      Vector_3 d2 = this->cross_pdct(axis, d1);
      FT l = this->sqlen(d2);
      if (l < (FT)0.0001) {
        d1 = this->constr_vec(ORIGIN, this->constr_pt(FT(1), FT(0), FT(0)));
        d2 = this->cross_pdct(axis, d1);
        l = this->sqlen(d2);
      }
      d2 = this->scale(d2, FT(1) / CGAL::sqrt(l));

      d1 = this->cross_pdct(axis, d2);
      l = CGAL::sqrt(this->sqlen(d1));
      if (l == 0)
        return;

      d1 = this->scale(d1, (FT)1.0 / l);

      // Process first point separately to initialize min/max
      Vector_3 vec = this->constr_vec(
        this->sph_center(m_sphere), this->point(indices[0]));

      // sign indicates northern or southern hemisphere
      FT proj = (this->scalar_pdct(axis, vec)) / rad;
      FT phi = atan2(this->scalar_pdct(vec, d2), this->scalar_pdct(vec, d1));
      FT x = FT(0), y = FT(0);
      concentric_mapping(phi, proj, rad, x, y);
      CGAL_assertion( x==x && y==y); // check not nan's

      min[0] = max[0] = x;
      min[1] = max[1] = y;

      parameterSpace[0] = std::pair<FT, FT>(x, y);

      for (std::size_t i = 1;i<indices.size();i++) {
        Vector_3 vec = this->constr_vec(
          this->sph_center(m_sphere), this->point(indices[i]));
        // sign indicates northern or southern hemisphere
        proj = (this->scalar_pdct(axis, vec)) / rad;
        phi = atan2(this->scalar_pdct(vec, d2), this->scalar_pdct(vec, d1));

        concentric_mapping(phi, proj, rad, x, y);
        CGAL_assertion( x==x && y==y); // check not nan's

        min[0] = (std::min<FT>)(min[0], x);
        max[0] = (std::max<FT>)(max[0], x);

        min[1] = (std::min<FT>)(min[1], y);
        max[1] = (std::max<FT>)(max[1], y);

        parameterSpace[i] = std::pair<FT, FT>(x, y);
      }

      // Is close to wrapping around? Check all three directions separately
      m_wrap_right = abs(max[0] - CGAL_M_PI_2 * rad) < (cluster_epsilon * 0.5);
      m_wrap_left = abs(min[0] + CGAL_M_PI_2 * rad) < (cluster_epsilon * 0.5);

      FT diff_top = CGAL::abs(-FT(CGAL_PI * rad) - min[1]) 
        + FT(CGAL_PI * rad) - max[1];
      m_wrap_top = diff_top < cluster_epsilon;

      if (m_wrap_top || m_wrap_left || m_wrap_right) {
        FT fl = FT(floor((CGAL_PI * rad) / cluster_epsilon));

        if (fl > 0.9) {
          FT adjusted_cf = FT(CGAL_PI * rad) / fl;

          if ( (adjusted_cf < (2 * cluster_epsilon)))
            cluster_epsilon = adjusted_cf;
        }

        // center bitmap at equator
        FT required_space = ceil(
          (std::max<FT>)(CGAL::abs(min[1]), max[1]) / cluster_epsilon)
          * cluster_epsilon;
        min[1] = -required_space;
        max[1] = required_space;
      }

      m_equator = std::size_t((abs(min[1])) / cluster_epsilon - 0.5);
    }

    virtual void post_wrap(const std::vector<unsigned int> &bitmap,
      const std::size_t &u_extent,
      const std::size_t &v_extent,
      std::vector<unsigned int> &labels) const {
        unsigned int l;
        unsigned int nw, n, ne;
        if (m_wrap_top && v_extent > 2) {
          // Handle first index separately.
          l = bitmap[0];
          if (l) {
            n = bitmap[(v_extent - 1) * u_extent];

            if (u_extent == 1) {
              if (n && l != n) {
                l = (std::min<unsigned int>)(n, l);
                update_label(labels, (std::max<unsigned int>)(n, l), l);
                return;
              }
            }

            ne = bitmap[(v_extent - 1) * u_extent + 1];

            if (n && n != l) {
              l = (std::min<unsigned int>)(n, l);
              update_label(labels, (std::max<unsigned int>)(n, l), l);
            }
            else if (ne && ne != l) {
              l = (std::min<unsigned int>)(ne, l);
              update_label(labels, (std::max<unsigned int>)(ne, l), l);
            }
          }

          for (std::size_t i = 1;i<u_extent - 1;i++) {
            l = bitmap[i];
            if (!l)
              continue;

            nw = bitmap[(v_extent - 1) * u_extent + i - 1];
            n = bitmap[(v_extent - 1) * u_extent + i];
            ne = bitmap[(v_extent - 1) * u_extent + i + 1];

            if (nw && nw != l) {
              l = (std::min<unsigned int>)(nw, l);
              update_label(labels, (std::max<unsigned int>)(nw, l), l);
            }
            if (n && n != l) {
              l = (std::min<unsigned int>)(n, l);
              update_label(labels, (std::max<unsigned int>)(n, l), l);
            }
            else if (ne && ne != l) {
              l = (std::min<unsigned int>)(ne, l);
              update_label(labels, (std::max<unsigned int>)(ne, l), l);
            }
          }

          // Handle last index separately
          l = bitmap[u_extent - 1];
          if (l) {
            n = bitmap[u_extent * v_extent - 1];
            nw = bitmap[u_extent * v_extent - 2];

            if (n && n != l) {
              l = (std::min<unsigned int>)(n, l);
              update_label(labels, (std::max<unsigned int>)(n, l), l);
            }
            else if (nw && nw != l) {
              l = (std::min<unsigned int>)(nw, l);
              update_label(labels, (std::max<unsigned int>)(nw, l), l);
            }
          }
        }

        // Walk upwards on the right side in the northern hemisphere
        if (m_wrap_right && v_extent > 2) {
          // First index
          l = bitmap[(m_equator + 1) * u_extent - 1];
          unsigned int ws = bitmap[(m_equator + 3) * u_extent - 1];
          
          if (l && ws && l != ws) {
            l = (std::min<unsigned int>)(ws, l);
            update_label(labels, (std::max<unsigned int>)(ws, l), l);
          }

          for (std::size_t i = 1;i<(v_extent>>1) - 1;i++) {
            l = bitmap[(m_equator - i + 1) * u_extent - 1];
            if (!l)
              continue;

            unsigned int wn = bitmap[(m_equator + i) * u_extent - 1];
            unsigned int w = bitmap[(m_equator + i + 1) * u_extent - 1];
            ws = bitmap[(m_equator + i + 2) * u_extent - 1];

            if (wn && wn != l) {
              l = (std::min<unsigned int>)(wn, l);
              update_label(labels, (std::max<unsigned int>)(wn, l), l);
            }
            if (w && w != l) {
              l = (std::min<unsigned int>)(w, l);
              update_label(labels, (std::max<unsigned int>)(w, l), l);
            }
            else if (ws && ws != l) {
              l = (std::min<unsigned int>)(ws, l);
              update_label(labels, (std::max<unsigned int>)(ws, l), l);
            }
          }

          // Last index
          l = bitmap[u_extent - 1];
          if (l) {
            unsigned int w = bitmap[u_extent * v_extent - 1];
            unsigned int wn = bitmap[(v_extent - 1) * u_extent - 1];

            if (w && w != l) {
              l = (std::min<unsigned int>)(w, l);
              update_label(labels, (std::max<unsigned int>)(w, l), l);
            }
            else if (wn && wn != l) {
              l = (std::min<unsigned int>)(wn, l);
              update_label(labels, (std::max<unsigned int>)(wn, l), l);
            }
          }
        }

        if (m_wrap_left && v_extent > 2) {
          // First index
          l = bitmap[(m_equator) * u_extent];
          unsigned int es = bitmap[(m_equator + 2) * u_extent];

          if (l && l != es) {
            l = (std::min<unsigned int>)(es, l);
            update_label(labels, (std::max<unsigned int>)(es, l), l);
          }

          for (std::size_t i = 1;i<(v_extent>>1) - 1;i++) {
            l = bitmap[(m_equator - i) * u_extent];
            if (!l)
              continue;

            unsigned int en = bitmap[(m_equator + i) * u_extent];
            unsigned int e = bitmap[(m_equator + i + 1) * u_extent];
            es = bitmap[(m_equator + i + 2) * u_extent];

            if (en && en != l) {
              l = (std::min<unsigned int>)(en, l);
              update_label(labels, (std::max<unsigned int>)(en, l), l);
            }
            if (e && e != l) {
              l = (std::min<unsigned int>)(e, l);
              update_label(labels, (std::max<unsigned int>)(e, l), l);
            }
            else if (es && es != l) {
              l = (std::min<unsigned int>)(es, l);
              update_label(labels, (std::max<unsigned int>)(es, l), l);
            }
          }

          // Last index
          l = bitmap[0];
          if (l) {
            unsigned int w = bitmap[(v_extent - 1) * u_extent];
            unsigned int wn = bitmap[(v_extent - 2) * u_extent];

            if (w && w != l) {
              l = (std::min<unsigned int>)(w, l);
              update_label(labels, (std::max<unsigned int>)(w, l), l);
            }
            else if (wn && wn != l) {
              l = (std::min<unsigned int>)(wn, l);
              update_label(labels, (std::max<unsigned int>)(wn, l), l);
            }
          }
        }
    }

    virtual bool supports_connected_component() const {
      return true;
    }
    
  private:
    Sphere_3 m_sphere;
    mutable bool m_wrap_right, m_wrap_top, m_wrap_left;
    mutable std::size_t m_equator;
/// \endcond
  };
}
}
#endif
