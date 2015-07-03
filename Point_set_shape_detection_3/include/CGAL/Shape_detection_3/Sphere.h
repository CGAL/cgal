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
     \tparam Traits a model of `EfficientRANSACTraits`  with the additional 
             requirement that the type `Traits::Sphere_3` is provided.
     */
  template <class Traits>
  class Sphere : public Shape_base<Traits> {
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
      return m_sphere.center();
    }
      
    /*!
      Access to the radius of the sphere.
     */
    FT radius() const {
      return CGAL::sqrt(m_sphere.squared_radius());
    }

    /// \cond SKIP_IN_MANUAL
    /*!
      Computes the squared Euclidean distance from query point to the shape.
      */
    FT squared_distance(const Point_3 &p) const {
      const FT d = CGAL::sqrt((m_sphere.center() - p).squared_length()) - CGAL::sqrt(m_sphere.squared_radius());
      return d * d;
    }
    
    /*!
      Helper function to write center, 
      radius of the sphere and number of assigned points into a string.
     */
    std::string info() const {
      std::stringstream sstr;
      Point_3 c = m_sphere.center();
      FT r = CGAL::sqrt(m_sphere.squared_radius());

      sstr << "Type: sphere center: (" << c.x() << ", " << c.y();
      sstr << ", " << c.z() << ") radius:" << r;
      sstr << " #Pts: " <<  this->m_indices.size();

      return sstr.str();
    }
    /// \endcond
  protected:
      /// \cond SKIP_IN_MANUAL
    void create_shape(const std::vector<std::size_t> &indices) {
      Point_3 p1 = this->point(indices[0]);
      Point_3 p2 = this->point(indices[1]);
      Point_3 p3 = this->point(indices[2]);

      Vector_3 n1 = this->normal(indices[0]);
      Vector_3 n2 = this->normal(indices[1]);
      Vector_3 n3 = this->normal(indices[2]);


      // Determine center: select midpoint of shortest line segment
      // between p1 and p2. Implemented from "3D game engine design" by Eberly 2001
      Vector_3 diff = p1 - p2;
      FT a = n1 * n1;
      FT b = -(n1 * n2);
      FT c = n2 * n2;
      FT d = n1 * diff;

      FT det = CGAL::abs(a * c - b * b);

      // degenerated when nearly parallel
      if (det < (FT)0.00001) {
        this->m_is_valid = false;
        return;
      }

      FT e = -n2 * diff;
      FT invDet = (FT) 1.0 / det;
      FT s = (b * e - c * d) * invDet;
      FT t = (d * b - a * e) * invDet;

      Point_3 center = CGAL::ORIGIN + (FT)0.5 * (((p1 + s * n1) - CGAL::ORIGIN)
                     + ((p2 + t * n2) - CGAL::ORIGIN));

      Vector_3 v1 = (p1 - center);
      Vector_3 v2 = (p2 - center);
      FT d1 = CGAL::sqrt(v1.squared_length());
      FT d2 = CGAL::sqrt(v2.squared_length());

      if (CGAL::abs(d1 - d2) > (FT)2.0 * this->m_epsilon) {
        this->m_is_valid = false;
        return;
      }

      v1 = v1 * ((FT)1.0 / d1);
      v2 = v2 * ((FT)1.0 / d2);

      if (n1 * v1 < this->m_normal_threshold ||
          n2 * v2 < this->m_normal_threshold) {
        this->m_is_valid = false;
        return;
      }

      Vector_3 v3 = (p3 - center);
      FT d3 = CGAL::sqrt(v3.squared_length());
      v3 = v3 * ((FT)1.0 / d3);

      FT radius = (d1 + d2) * (FT)0.5;

      if (CGAL::abs(d3 - radius) > this->m_epsilon ||
          n3 * v3 < this->m_normal_threshold) {
        this->m_is_valid = false;
        return;
      }

      this->m_is_valid = true;

      m_sphere = Sphere_3(center, radius * radius);
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) {

      FT radius = CGAL::sqrt(m_sphere.squared_radius());

      for (std::size_t i = 0;i<indices.size();i++) {
        dists[i] = CGAL::sqrt((m_sphere.center()
          - this->point(indices[i])).squared_length())
          - radius;

        dists[i] = dists[i] * dists[i];
      }
    }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                               std::vector<FT> &angles) const {
      for (std::size_t i = 0;i<indices.size();i++) {
        Vector_3 n = m_sphere.center() - this->point(indices[i]);

        FT length = CGAL::sqrt(n.squared_length());
        if (length == 0) {
          angles[i] = (FT)1.0;
          continue;
        }

        n = n * (FT)1.0 / length;
        angles[i] = CGAL::abs(this->normal(indices[i]) * n);
      }
    }

    virtual FT cos_to_normal(const Point_3 &p, const Vector_3 &n) const {
      Vector_3 sphere_normal = m_sphere.center() - p;
      FT length = (FT)(CGAL::sqrt(n.squared_length()));
      if (length == 0)
        return 1;

      sphere_normal = sphere_normal * ((FT)1.0 / length);
      return CGAL::abs(sphere_normal * n);
    }
      
    virtual std::size_t minimum_sample_size() const {
      return 3;
    }

    // Maps to the range [-1,1]^2
    static void concentric_mapping(FT phi, FT proj, FT rad, FT &x, FT &y) {
      phi = (phi < FT(-M_PI_4)) ? phi + FT(2 * M_PI) : phi;
      proj = (proj < FT(-1.0)) ? FT(-1.0) : ((proj > FT(1.0)) ? FT(1.0) : proj);
      FT r = FT(acos(double(CGAL::abs(proj)))) / FT(M_PI_2);

      FT a = 0, b = 0;
      if (phi < FT(M_PI_4)) {
        a = r;
        b = phi * r / FT(M_PI_4);
      }
      else if (phi < FT(3.0 * M_PI_4)) {
        a = -FT(phi - M_PI_2) * r / FT(M_PI_4);
        b = r;
      }
      else if (phi < FT(5.0 * M_PI_4)) {
        a = -r;
        b = (phi - FT(M_PI)) * (-r) / FT(M_PI_4);
      }
      else {
        a = (phi - 3 * FT(M_PI_2)) * r / FT(M_PI_4);
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
      x = FT(x * M_PI_2 * rad);
      y = FT(y * M_PI_2 * rad);
    }

    virtual void parameters(const std::vector<std::size_t> &indices,
                            std::vector<std::pair<FT, FT> > &parameterSpace,
                            FT &cluster_epsilon,
                            FT min[2],
                            FT max[2]) const {
      Vector_3 axis;
      FT rad = radius();
      // Take average normal as axis
      for (std::size_t i = 0;i<indices.size();i++)
        axis = axis + this->normal(indices[i]);
      axis = axis / (CGAL::sqrt(axis.squared_length()));

      // create basis d1, d2
      Vector_3 d1 = Vector_3((FT) 0, (FT) 0, (FT) 1);
      Vector_3 d2 = CGAL::cross_product(axis, d1);
      FT l = d2.squared_length();
      if (l < (FT)0.0001) {
        d1 = Vector_3((FT) 1, (FT) 0, (FT) 0);
        d2 = CGAL::cross_product(axis, d1);
        l = d2.squared_length();
      }
      d2 = d2 / CGAL::sqrt(l);

      d1 = CGAL::cross_product(axis, d2);
      l = CGAL::sqrt(d1.squared_length());
      if (l == 0)
        return;

      d1 = d1 * (FT)1.0 / l;

      // Process first point separately to initialize min/max
      Vector_3 vec = this->point(indices[0]) - m_sphere.center();

      // sign indicates northern or southern hemisphere
      FT proj = (axis * vec) / rad;
      FT phi = atan2(vec * d2, vec * d1);
      FT x = FT(0), y = FT(0);
      concentric_mapping(phi, proj, rad, x, y);

      min[0] = max[0] = x;
      min[1] = max[1] = y;

      parameterSpace[0] = std::pair<FT, FT>(x, y);

      for (std::size_t i = 1;i<indices.size();i++) {
        Vector_3 vec = this->point(indices[i]) - m_sphere.center();
        // sign indicates northern or southern hemisphere
        proj = (axis * vec) / rad;
        phi = atan2(vec * d2, vec * d1);

        concentric_mapping(phi, proj, rad, x, y);

        min[0] = (std::min<FT>)(min[0], x);
        max[0] = (std::max<FT>)(max[0], x);

        min[1] = (std::min<FT>)(min[1], y);
        max[1] = (std::max<FT>)(max[1], y);

        parameterSpace[i] = std::pair<FT, FT>(x, y);
      }

      // Is close to wrapping around? Check all three directions separately
      m_wrap_right = abs(max[0] - M_PI_2 * rad) < (cluster_epsilon * 0.5);
      m_wrap_left = abs(min[0] + M_PI_2 * rad) < (cluster_epsilon * 0.5);

      FT diff_top = CGAL::abs(-FT(M_PI * rad) - min[1]) 
        + FT(M_PI * rad) - max[1];
      m_wrap_top = diff_top < cluster_epsilon;

      if (m_wrap_top || m_wrap_left || m_wrap_right) {
        cluster_epsilon = FT(M_PI * rad)
          / FT(floor((M_PI * rad) / cluster_epsilon));

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
                update_label(labels, (std::max<unsigned int>)(n, l),
                             l = (std::min<unsigned int>)(n, l));
                return;
              }
            }

            ne = bitmap[(v_extent - 1) * u_extent + 1];

            if (n && n != l)
              update_label(labels, (std::max<unsigned int>)(n, l),
                           l = (std::min<unsigned int>)(n, l));
            else if (ne && ne != l)
              update_label(labels, (std::max<unsigned int>)(ne, l),
                           l = (std::min<unsigned int>)(ne, l));
          }

          for (std::size_t i = 1;i<u_extent - 1;i++) {
            l = bitmap[i];
            if (!l)
              continue;

            nw = bitmap[(v_extent - 1) * u_extent + i - 1];
            n = bitmap[(v_extent - 1) * u_extent + i];
            ne = bitmap[(v_extent - 1) * u_extent + i + 1];

            if (nw && nw != l)
              update_label(labels, (std::max<unsigned int>)(nw, l),
                           l = (std::min<unsigned int>)(nw, l));
            if (n && n != l)
              update_label(labels, (std::max<unsigned int>)(n, l),
                           l = (std::min<unsigned int>)(n, l));
            else if (ne && ne != l)
              update_label(labels, (std::max<unsigned int>)(ne, l),
                           l = (std::min<unsigned int>)(ne, l));
          }

          // Handle last index separately
          l = bitmap[u_extent - 1];
          if (l) {
            n = bitmap[u_extent * v_extent - 1];
            nw = bitmap[u_extent * v_extent - 2];

            if (n && n != l)
              update_label(labels, (std::max<unsigned int>)(n, l), 
                           l = (std::min<unsigned int>)(n, l));
            else if (nw && nw != l)
              update_label(labels, (std::max<unsigned int>)(nw, l),
                           l = (std::min<unsigned int>)(nw, l));
          }
        }

        // Walk upwards on the right side in the northern hemisphere
        if (m_wrap_right && v_extent > 2) {
          // First index
          l = bitmap[(m_equator + 1) * u_extent - 1];
          unsigned int ws = bitmap[(m_equator + 3) * u_extent - 1];
          
          if (l && ws && l != ws)
            update_label(labels, (std::max<unsigned int>)(ws, l),
                         l = (std::min<unsigned int>)(ws, l));

          for (std::size_t i = 1;i<(v_extent>>1) - 1;i++) {
            l = bitmap[(m_equator - i + 1) * u_extent - 1];
            if (!l)
              continue;

            unsigned int wn = bitmap[(m_equator + i) * u_extent - 1];
            unsigned int w = bitmap[(m_equator + i + 1) * u_extent - 1];
            ws = bitmap[(m_equator + i + 2) * u_extent - 1];

            if (wn && wn != l)
              update_label(labels, (std::max<unsigned int>)(wn, l),
                           l = (std::min<unsigned int>)(wn, l));
            if (w && w != l)
              update_label(labels, (std::max<unsigned int>)(w, l),
                           l = (std::min<unsigned int>)(w, l));
            else if (ws && ws != l)
              update_label(labels, (std::max<unsigned int>)(ws, l),
                           l = (std::min<unsigned int>)(ws, l));
          }

          // Last index
          l = bitmap[u_extent - 1];
          if (l) {
            unsigned int w = bitmap[u_extent * v_extent - 1];
            unsigned int wn = bitmap[(v_extent - 1) * u_extent - 1];

            if (w && w != l)
              update_label(labels, (std::max<unsigned int>)(w, l),
                           l = (std::min<unsigned int>)(w, l));
            else if (wn && wn != l)
              update_label(labels, (std::max<unsigned int>)(wn, l),
                                    l = (std::min<unsigned int>)(wn, l));
          }
        }

        if (m_wrap_left && v_extent > 2) {
          // First index
          l = bitmap[(m_equator) * u_extent];
          unsigned int es = bitmap[(m_equator + 2) * u_extent];

          if (l && l != es)
            update_label(labels, (std::max<unsigned int>)(es, l),
                         l = (std::min<unsigned int>)(es, l));

          for (std::size_t i = 1;i<(v_extent>>1) - 1;i++) {
            l = bitmap[(m_equator - i) * u_extent];
            if (!l)
              continue;

            unsigned int en = bitmap[(m_equator + i) * u_extent];
            unsigned int e = bitmap[(m_equator + i + 1) * u_extent];
            es = bitmap[(m_equator + i + 2) * u_extent];

            if (en && en != l)
              update_label(labels, (std::max<unsigned int>)(en, l),
                           l = (std::min<unsigned int>)(en, l));
            if (e && e != l)
              update_label(labels, (std::max<unsigned int>)(e, l),
                           l = (std::min<unsigned int>)(e, l));
            else if (es && es != l)
              update_label(labels, (std::max<unsigned int>)(es, l),
                           l = (std::min<unsigned int>)(es, l));
          }

          // Last index
          l = bitmap[0];
          if (l) {
            unsigned int w = bitmap[(v_extent - 1) * u_extent];
            unsigned int wn = bitmap[(v_extent - 2) * u_extent];

            if (w && w != l)
              update_label(labels, (std::max<unsigned int>)(w, l),
                           l = (std::min<unsigned int>)(w, l));
            else if (wn && wn != l)
              update_label(labels, (std::max<unsigned int>)(wn, l),
                           l = (std::min<unsigned int>)(wn, l));
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
