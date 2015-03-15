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

#ifndef CGAL_SHAPE_DETECTION_3_SPHERE_SHAPE_H
#define CGAL_SHAPE_DETECTION_3_SPHERE_SHAPE_H

#include "Shape_base.h"

/*!
 \file Sphere_shape.h
 */

namespace CGAL {
    /*!
     \ingroup PkgPointSetShapeDetection3
     \brief Sphere_shape implements Shape_base. The sphere is represented by its center and the radius.
     */
  template <class Sd_traits>
  class Sphere_shape : public Shape_base<Sd_traits> {
  public:    
    /// \cond SKIP_IN_MANUAL
    typedef typename Sd_traits::Input_iterator Input_iterator;
      ///< random access iterator for input data.
    typedef typename Sd_traits::Point_pmap Point_pmap;
      ///< property map to access the location of an input point.
    typedef typename Sd_traits::Normal_pmap Normal_pmap;
      ///< property map to access the unoriented normal of an input point.
    typedef typename Sd_traits::Geom_traits::Vector_3 Vector;
      ///< vector type.
    typedef typename Sd_traits::Geom_traits::Sphere_3 Sphere;
      ///< sphere type.
    typedef typename Sd_traits::Geom_traits::FT FT;
      ///< number type.
    typedef typename Sd_traits::Geom_traits::Point_3 Point;
      ///< point type.
    /// \endcond

  public:
    Sphere_shape() :  Shape_base<Sd_traits>() {}

    /*!
      Conversion operator to convert to common Sphere_3 type.
     */
    operator Sphere() const {
      return m_sphere;
    }

    /*!
      Access to the center.
     */
    Point center() const {
      return m_sphere.center();
    }
      
    /*!
      Access to the radius of the sphere.
     */
    FT radius() const {
      return sqrt(m_sphere.squared_radius());
    }

    /// \cond SKIP_IN_MANUAL
    /*!
      Computes the squared Euclidean distance from query point to the shape.
      */
    FT squared_distance(const Point &p) const {
      const FT d = sqrt((m_sphere.center() - p).squared_length()) - sqrt(m_sphere.squared_radius());
      return d * d;
    }
    
    /*!
      Helper function to write center, 
      radius of the sphere and number of assigned points into a string.
     */
    std::string info() const {
      std::stringstream sstr;
      Point c = m_sphere.center();
      FT r = sqrt(m_sphere.squared_radius());

      sstr << "Type: sphere center: (" << c.x() << ", " << c.y();
      sstr << ", " << c.z() << ") radius:" << r;
      sstr << " #Pts: " <<  this->m_indices.size();

      return sstr.str();
    }
    /// \endcond
  protected:
      /// \cond SKIP_IN_MANUAL
    void create_shape(const std::vector<std::size_t> &indices) {
      Point p1 = this->point(indices[0]);
      Point p2 = this->point(indices[1]);
      Point p3 = this->point(indices[2]);

      Vector n1 = this->normal(indices[0]);
      Vector n2 = this->normal(indices[1]);
      Vector n3 = this->normal(indices[2]);


      // Determine center: select midpoint of shortest line segment
      // between p1 and p2. Implemented from "3D game engine design" by Eberly 2001
      Vector diff = p1 - p2;
      FT a = n1 * n1;
      FT b = -(n1 * n2);
      FT c = n2 * n2;
      FT d = n1 * diff;

      FT det = abs(a * c - b * b);

      // degenerated when nearly parallel
      if (det < (FT)0.00001) {
        this->m_isValid = false;
        return;
      }

      FT e = -n2 * diff;
      FT invDet = 1.0 / det;
      FT s = (b * e - c * d) * invDet;
      FT t = (d * b - a * e) * invDet;

      Point center = CGAL::ORIGIN + 0.5 * (((p1 + s * n1) - CGAL::ORIGIN)
                     + ((p2 + t * n2) - CGAL::ORIGIN));

      Vector v1 = (p1 - center);
      Vector v2 = (p2 - center);
      FT d1 = sqrt(v1.squared_length());
      FT d2 = sqrt(v2.squared_length());

      if (abs(d1 - d2) > (FT)2.0 * this->m_epsilon) {
        this->m_isValid = false;
        return;
      }

      v1 = v1 * (1.0 / d1);
      v2 = v2 * (1.0 / d2);

      if (n1 * v1 < this->m_normal_threshold ||
          n2 * v2 < this->m_normal_threshold) {
        this->m_isValid = false;
        return;
      }

      Vector v3 = (p3 - center);
      FT d3 = sqrt(v3.squared_length());
      v3 = v3 * (1.0 / d3);

      FT radius = (d1 + d2) * 0.5;

      if (abs(d3 - radius) > this->m_epsilon ||
          n3 * v3 < this->m_normal_threshold) {
        this->m_isValid = false;
        return;
      }

      m_sphere = Sphere(center, radius * radius);
    }

    virtual void squared_distance(const std::vector<std::size_t> &indices,
                                  std::vector<FT> &dists) {

      FT radius = sqrt(m_sphere.squared_radius());

      for (std::size_t i = 0;i<indices.size();i++) {
        dists[i] = sqrt((m_sphere.center()
          - this->point(indices[i])).squared_length())
          - radius;

        dists[i] = dists[i] * dists[i];
      }
    }

    virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                               std::vector<FT> &angles) const {
      for (std::size_t i = 0;i<indices.size();i++) {
        Vector n = m_sphere.center() - this->point(indices[i]);

        n = n * (1.0 / (sqrt(n.squared_length())));
        angles[i] = abs(this->normal(indices[i]) * n);
      }
    }

    virtual FT cos_to_normal(const Point &p, const Vector &n) const {
      Vector sphere_normal = m_sphere.center() - p;
      sphere_normal = sphere_normal * ((FT)1.0 / (CGAL::sqrt(n.squared_length())));
      return abs(sphere_normal * n);
    }
      
    virtual std::size_t required_samples() const {
      return 3;
    }

    virtual bool supports_connected_component() const {
      return false;
    }

    // U is longitude
    virtual bool wraps_u() const {
      return true;
    }

    // V is latitude
    virtual bool wraps_v() const {
      return false;
    }

  private:
    Sphere m_sphere;
/// \endcond
  };
}
#endif
