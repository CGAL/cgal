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

#ifndef CGAL_SHAPE_DETECTION_3_CYLINDER_SHAPE_H
#define CGAL_SHAPE_DETECTION_3_CYLINDER_SHAPE_H

#include "Shape_base.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

// CODE REVIEW
// fix naming, eg _p
// check all degenerate cases
// use (FT) in front of const such as (FT)1.0
// remove commented lines

/*!
 \file Cylinder_shape.h
 */

namespace CGAL {
    /*!
     \brief Cylinder_shape implements Shape_base. 
            The cylinder is parameterized by the axis,
            i.e. point and direction, and the radius.
     \ingroup PkgPointSetShapeDetection3
     */
  template <class Sd_traits>
  class Cylinder_shape : public Shape_base<Sd_traits> {
  public:
    /// \cond SKIP_IN_MANUAL
    typedef typename Sd_traits::Input_iterator Input_iterator; ///< random access iterator for input data.
    typedef typename Sd_traits::Point_pmap Point_pmap;  ///< property map to access the location of an input point.
    typedef typename Sd_traits::Normal_pmap Normal_pmap;  ///< property map to access the unoriented normal of an input point.
    typedef typename Sd_traits::Geom_traits::Vector_3 Vector;///< vector type.
    /// \endcond

    typedef typename Sd_traits::Geom_traits::FT FT; ///< number type.
    typedef typename Sd_traits::Geom_traits::Line_3 Line; ///< line type.
    typedef typename Sd_traits::Geom_traits::Point_3 Point;///< point type.

  public:
    Cylinder_shape() : Shape_base<Sd_traits>() {}
      /*!
       Helper function to write axis and radius of the cylinder and number of assigned points into a string.
       */
    std::string info() const {
      std::stringstream sstr;
      Point c = m_axis.point();
      Vector a = m_axis.to_vector();

      sstr << "Type: cylinder center: (" << c.x() << ", " << c.y() << ", " << c.z() << ") axis: (" << a.x() << ", " << a.y() << ", " << a.z() << ") radius:" << m_radius
        << " #Pts: " <<  this->m_indices.size();

      return sstr.str();
    }
      
      /*!
        Axis of the cylinder.
       */
    Line axis() const {
      return m_axis;
    }
      /*!
        Radius of the cylinder.
       */
    FT radius() const {
      return m_radius;
    }


	

    /*!
      Computes squared Euclidean distance from query point to the shape.
      */ 
    FT squared_distance(const Point &_p) const {
      Vector a = m_axis.to_vector();
      a = a * (1.0 / sqrt(a.squared_length()));
      Vector v = _p - m_point_on_axis;
      v = v - ((v * a) * a);
      FT d = sqrt(v.squared_length()) - m_radius;
      return d * d;
    }

  protected:
      /// \cond SKIP_IN_MANUAL
    virtual void create_shape(const std::vector<size_t> &indices) {
      Point p1 = get(this->m_point_pmap, *(this->m_first + indices[0]));
      Point p2 = get(this->m_point_pmap, *(this->m_first + indices[1]));
      //Point p3 = get(this->m_point_pmap, (this->m_first + indices[2]));

      Vector n1 = get(this->m_normal_pmap, *(this->m_first + indices[0]));
      Vector n2 = get(this->m_normal_pmap, *(this->m_first + indices[1]));
      //Vector n3 = get(this->m_normal_pmap, (this->m_first + indices[2]));

      Vector axis = CGAL::cross_product(n1, n2);
      FT axisL = sqrt(axis.squared_length());
      if (axisL < 0.001) {
        this->m_isValid = false;
        return;
      }
      axis = axis * (1.0 / axisL);

      // establish two directions in the plane axis * x = 0, 
      // whereas xDir is the projected n1
      Vector xDir = n1 - (n1 * axis) * axis;
      xDir = xDir * (1.0 / sqrt(xDir.squared_length()));
      Vector yDir = CGAL::cross_product(axis, xDir);
      yDir = yDir * (1.0 / sqrt(yDir.squared_length()));

      FT n2x = n2 * yDir;
      FT n2y = -n2 * xDir;

      Vector dist = p2 - p1;

      FT Ox = xDir * dist;
      FT Oy = yDir * dist;

      FT lineDist = n2x * Ox + n2y * Oy;

      m_radius = lineDist / n2x;
      m_point_on_axis = p1 + m_radius * xDir;
      m_radius = abs(m_radius);

      m_axis = Line(m_point_on_axis, axis);

      if (squared_distance(p1) > this->m_epsilon ||
          (cos_to_normal(p1, n1) < this->m_normal_threshold)) {
        this->m_isValid = false;
        return;
      }
    }

    void parameters(std::vector<std::pair<FT, FT> > &parameterSpace,
                    const std::vector<size_t> &indices,
                    FT min[2],
                    FT max[2]) const {

      Vector d1 = Vector(0, 0, 1);
      Vector a = m_axis.to_vector();
      a = a * (1.0 / sqrt(a.squared_length()));

      Vector d2 = CGAL::cross_product(a, d1);
      FT l = d2.squared_length();
      if (l < 0.0001) {
        d1 = Vector(1, 0, 0);
        d2 = CGAL::cross_product(m_axis.to_vector(), d1);
        l = d2.squared_length();
        if (l < 0.0001) {
          std::cout << "Cylinder::pointOnPrimitive() \
                       construction failed!" << std::endl;
        }
      }
      d2 = d2 / sqrt(l);

      d1 = CGAL::cross_product(m_axis.to_vector(), d2);
      d1 = d1 * (1.0 / sqrt(d1.squared_length()));

      // 1.0 / circumfence
      FT c = 1.0 / 2 * M_PI * m_radius;

      // first one separate for initializing min/max
      Vector vec = get(this->m_point_pmap,
                       *(this->m_first + indices[0])) - m_point_on_axis;
      FT v = vec * a;
      vec = vec - ((vec * a) * a);
      vec = vec * (1.0 / sqrt(vec.squared_length()));

      FT a1 = acos(vec * d1);
      FT a2 = acos(vec * d2);

      FT u = ((a2 < M_PI_2) ? 2 * M_PI - a1 : a1) * c;

      parameterSpace[0] = std::pair<FT, FT>(u, v);

      min[0] = max[0] = u;
      min[1] = max[1] = v;

      for (size_t i = 0;i<indices.size();i++) {
        Vector vec = get(this->m_point_pmap,
                         *(this->m_first + indices[i])) - m_point_on_axis;
        FT v = vec * a;
        vec = vec - ((vec * a) * a);
        vec = vec * (1.0 / sqrt(vec.squared_length()));

        FT a1 = acos(vec * d1);
        FT a2 = acos(vec * d2);

        FT u = ((a2 < M_PI_2) ? 2 * M_PI - a1 : a1) * c;

        min[0] = (std::min<FT>)(min[0], u);
        max[0] = (std::max<FT>)(max[0], u);
        min[1] = (std::min<FT>)(min[1], v);
        max[1] = (std::max<FT>)(max[1], v);

        parameterSpace[i] = std::pair<FT, FT>(u, v);
      }
    }


    void squared_distance(std::vector<FT> &dists,
                          const std::vector<int> &shapeIndex,
                          const std::vector<size_t> &indices) {
      Vector a = m_axis.to_vector();
      a = a * (1.0 / sqrt(a.squared_length()));
      for (size_t i = 0;i<indices.size();i++) {

        if (shapeIndex[indices[i]] == -1) {
          Vector v = get(this->m_point_pmap,
                         *(this->m_first + indices[i]))
                         - m_point_on_axis;

          v = v - ((v * a) * a);
          dists[i] = sqrt(v.squared_length()) - m_radius;
          dists[i] = dists[i] * dists[i];
        }
      }
    }

    void cos_to_normal(std::vector<FT> &angles,
                       const std::vector<int> &shapeIndex, 
                       const std::vector<size_t> &indices) const {
      Vector a = m_axis.to_vector();
      a = a * (1.0 / sqrt(a.squared_length()));
      for (size_t i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {

          Vector v = get(this->m_point_pmap,
                         *(this->m_first + indices[i])) 
                         - m_point_on_axis;

          v = v - ((v * a) * a);
          v = v * (1.0 / sqrt(v.squared_length()));
          angles[i] = abs(v * get(this->m_normal_pmap,
                                  *(this->m_first + indices[i])));
        }
      }
    }

    FT cos_to_normal(const Point &_p, const Vector &_n) const {
      Vector a = m_axis.to_vector();
      a = a * (1.0 / sqrt(a.squared_length()));
      Vector v = _p - m_point_on_axis;
      v = v - ((v * a) * a);
      v = v * (1.0 / sqrt(v.squared_length()));
      return abs(v * _n);
    }

    virtual size_t required_samples() const {
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
    Line m_axis;
    Point m_point_on_axis;
      
      /// \endcond
  };
}
#endif
