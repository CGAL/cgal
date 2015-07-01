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

    virtual bool supports_connected_component() const {
      return false;
    }

  private:
    Sphere_3 m_sphere;
/// \endcond
  };
}
}
#endif
