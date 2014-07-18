#ifndef CGAL_SHAPE_DETECTION_3_TORUS_SHAPE_H
#define CGAL_SHAPE_DETECTION_3_TORUS_SHAPE_H

#include "Shape_base.h"
#include <set>
#include <math.h>
#include <cmath>

/*!
 \file Torus_shape.h
 */

namespace CGAL {
    /*!
     \brief Torus_shape implements Shape_base. The torus is parameterized by the symmetry axis, its center on the axis and the major and minor radii.     */
  template <class Sd_traits>
  class Torus_shape : public Shape_base<Sd_traits> {
  public:
    typedef typename Sd_traits::Input_iterator Input_iterator; ///< random access iterator for input data.
    typedef typename Sd_traits::Geom_traits::FT FT; ///< number type.
    typedef typename Sd_traits::Geom_traits::Point_3 Point;///< point type.
    typedef typename Sd_traits::Geom_traits::Vector_3 Vector;///< vector type.
      typedef typename Sd_traits::Point_pmap Point_pmap;  ///< property map to access the location of an input point.
      typedef typename Sd_traits::Normal_pmap Normal_pmap;  ///< property map to access the unoriented normal of an input point.

    Torus_shape() : Shape_base<Sd_traits>() {}
      
      /*!
       Direction of the symmetry axis.
       */
    Vector axis() const {
      return m_axis;
    }
      /*!
       Center point on the symmetry axis.
       */
    Point center() const {
      return m_center;
    }
      /*!
       Helper function to write center point, symmetry axis and the two radii into a string.
       */
    std::string info() const {
      std::stringstream sstr;
      sstr << "Type: torus center(" << m_center.x() << ", " << m_center.y() << ", " << m_center.z() << ") axis(" << m_axis.x() << ", " << m_axis.y() << ", " << m_axis.z() << ") major radius = " << m_majorRad << " minor radius = " << m_minorRad << " #Pts: " << this->m_indices.size();

      return sstr.str();
    }
      
      /*!
       Major radius of the torus.
       */

    FT major_radius() const {
      return m_majorRad;
    }
      
      /*!
            Minor radius of the torus.
            */
    FT minor_radius() const {
      return m_minorRad;
    }

  protected:
      /// \cond SKIP_IN_MANUAL
      void create_shape(const std::vector<int> &indices) {
      Point p1 = get(this->m_pointPMap, *(this->m_first + indices[0]));
      Point p2 = get(this->m_pointPMap, *(this->m_first + indices[1]));
      Point p3 = get(this->m_pointPMap, *(this->m_first + indices[2]));
      Point p4 = get(this->m_pointPMap, *(this->m_first + indices[3]));

      Vector n1 = get(this->m_normalPMap, *(this->m_first + indices[0]));
      Vector n2 = get(this->m_normalPMap, *(this->m_first + indices[1]));
      Vector n3 = get(this->m_normalPMap, *(this->m_first + indices[2]));
      Vector n4 = get(this->m_normalPMap, *(this->m_first + indices[3]));

      // Implemented method from 'Geometric least-squares fitting of spheres, cylinders, cones and tori' by G. Lukacs,A.D. Marshall, R. R. Martin
      double a01 = CGAL::cross_product(n1, n2) * n3;
      double b01 = CGAL::cross_product(n1, n2) * n4;
      double a0 = CGAL::cross_product(p3 - p2, n1) * n3;
      double b0 = CGAL::cross_product(p4 - p2, n1) * n4;
      double a1 = CGAL::cross_product(p1 - p3, n2) * n3;
      double b1 = CGAL::cross_product(p1 - p4, n2) * n4;
      double a = CGAL::cross_product(p1 - p3, p2 - p1) * n3;
      double b = CGAL::cross_product(p1 - p4, p2 - p1) * n4;

      double div = 1.0 / (b1 * a01 - b01 * a1);
      double p = ((a01 * b + b1 * a0 - b0 * a1 - b01 * a)) * div * 0.5;
      double q = (b * a0 - b0 * a) * div;

      FT root = p * p - q;
      if (p * p - q < 0)
        root = 0;

      double y1 = -p - sqrt(root);
      double y2 = -p + sqrt(root);
      double x1 = -(a1 * y1 + a) / (a01 * y1 + a0);
      double x2 = -(a1 * y2 + a) / (a01 * y2 + a0);

      // 1. center + axis
      FT majorRad1, minorRad1, dist1;
      Point c1;
      Vector axis1;
      if (is_finite(x1) && is_finite(y1)) {
        c1 = p1 + n1 * x1;
        axis1 = c1 - (p2 + n2 * y1);
        axis1 = axis1 / sqrt(axis1.squared_length());
        //dist1 = getCircle(first, c1, axis1, majorRad1, minorRad1);
      }
      else dist1 = FLT_MAX;

      // 2. center + axis
      FT majorRad2 = 0, minorRad2 = 0, dist2;
      Point c2;
      Vector axis2;
      if (is_finite(x2) && is_finite(y2)) {
        c2 = p1 + n1 * x2;
        axis2 = c2 - (p2 + n2 * y2);
        axis2 = axis2 / sqrt(axis2.squared_length());
        //dist2 = getCircle(first, c2, axis2, majorRad2, minorRad2);
      }
      else dist2 = FLT_MAX;


      if (dist1 < dist2) {
        if (dist1 > this->m_epsilon)
          return;
        m_center = c1;
        m_axis = axis1;
        m_majorRad = majorRad1;
        m_minorRad = sqrt(minorRad1);
      }
      else {
        if (dist2 > this->m_epsilon)
          return;
        m_center = c2;
        m_axis = axis2;
        m_majorRad = majorRad2;
        m_minorRad = sqrt(minorRad2);
      }

      this->m_isValid = true;

      //create primitive
      //validate points and normals
    }

    void parameters(std::vector<std::pair<FT, FT> > &parameterSpace, const std::vector<int> &indices, FT min[2], FT max[2]) const {
      return;
    }

    void parameter_extend(const Point &center, FT width, FT min[2], FT max[2]) const {
      return;
    }

    FT squared_distance(const Point &_p) const {
      Vector d = _p - m_center;
      // height over symmetry plane
      FT p = d * m_axis;
      // distance from axis in plane
      FT l = sqrt(d * d - p * p);

      /*
      Vector inPlane = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, d));
      if (inPlane * d < 0)
      inPlane = -inPlane;*/

      // inPlane distance from circle
      FT l2 = m_majorRad - l;

      // distance from torus
      l = sqrt(p * p + l2 * l2) - m_minorRad;

      return l * l;
    }

    void squared_distance(std::vector<FT> &dists, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) {
      for (unsigned int i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {
          Point po = get(this->m_pointPMap, *(this->m_first + i));
          Vector d = po - m_center;
          // height over symmetry plane
          FT p = d * m_axis;
          // distance from axis in plane
          FT l = sqrt(d * d - p * p);

          // inPlane distance from circle
          FT l2 = m_majorRad - l;

          // distance from torus
          l = sqrt(p * p + l2 * l2) - m_minorRad;
          dists[i] = l * l;
        }
      }
    }

    void cos_to_normal(std::vector<FT> &angles, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) const {
      for (unsigned int i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {
          Vector d = get(this->m_pointPMap, *(this->m_first + i)) - m_center;
          // height over symmetry plane
          //FT p = d * m_axis;
          // distance from axis in plane
          //FT l = sqrt(d * d - p * p);

          Vector inPlane = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, d));
          if (inPlane * d < 0)
            inPlane = -inPlane;
          inPlane = inPlane / sqrt(inPlane.squared_length());

          d = get(this->m_pointPMap, *(this->m_first + i)) - (m_center + inPlane * m_majorRad);
          d = d / sqrt(d.squared_length());
          angles[i] = abs(d * get(this->m_normalPMap, *(this->m_first + i)));
        }
      }
    }

    FT cos_to_normal(const Point &_p, const Vector &_n) const {
      Vector d = _p - m_center;
      // height over symmetry plane
      //FT p = d * m_axis;
      // distance from axis in plane
      //FT l = sqrt(d * d - p * p);

      Vector inPlane = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, d));
      if (inPlane * d < 0)
        inPlane = -inPlane;
      inPlane = inPlane / sqrt(inPlane.squared_length());

      d = _p - (m_center + inPlane * m_majorRad);
      d = d / sqrt(d.squared_length());
      return abs(d * _n);
    }
      
      virtual int required_samples() const {
          return 4;
      }

    virtual bool supports_connected_component() const {
      return false;
    }

    virtual bool wraps_u() const {
      return false;
    }

    virtual bool wraps_v() const {
      return false;
    }

  private:
    /*
    FT getCircle(Point &center, const Vector &axis, FT &majorRad, FT &minorRad) const {
    // create spin image
    Geom_traits::Point_2 pts[4];
    for (unsigned int i = 0;i<4;i++) {
    Vector d = first[i].first - center;
    FT p = d * axis;
    pts[i] = Geom_traits::Point_2(p, sqrt(d * d - p * p));
    }
    Geom_traits::Circle_2 c(pts[0], pts[1], pts[2]);
    minorRad = c.squared_radius();
    majorRad = c.center().y();
    center = center + c.center().x() * axis;

    return abs((pts[3] - c.center()).squared_length() - c.squared_radius());
    }*/

    Point m_center;
    Vector m_axis;
    FT m_majorRad;
    FT m_minorRad;
      /// \endcond
  };
}
#endif