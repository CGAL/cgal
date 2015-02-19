#ifndef CGAL_SHAPE_DETECTION_3_CONE_SHAPE_H
#define CGAL_SHAPE_DETECTION_3_CONE_SHAPE_H

#include "Shape_base.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/*!
 \file Cone_shape.h
 */

// CODE REVIEW
// fix naming
// check degenerate cases, eg every time a division occurs
// is m_pointOnPrimitive used?
// 


namespace CGAL {
    /*!
     \brief Cone_shape implements Shape_base.
      The cone is parameterized by its apex, the axis and the opening angle.
      This representation models an infinite single-cone
       and does not consider a base plane.
     \ingroup PkgPointSetShapeDetection3
     */

  template <class Sd_traits>
  class Cone_shape : public Shape_base<Sd_traits> {
  public:
    /// \cond SKIP_IN_MANUAL
    typedef typename Sd_traits::Input_iterator Input_iterator; ///< random access iterator for input data.
    typedef typename Sd_traits::Point_pmap Point_pmap;  ///< property map to access the location of an input point.
    typedef typename Sd_traits::Normal_pmap Normal_pmap;  ///< property map to access the unoriented normal of an input point.
    /// \endcond

    typedef typename Sd_traits::Geom_traits::FT FT; ///< number type.
    typedef typename Sd_traits::Geom_traits::Point_3 Point;///< point type.
    typedef typename Sd_traits::Geom_traits::Vector_3 Vector;///< vector type.

      Cone_shape() : Shape_base<Sd_traits>() {}

	
      
      /*!
       Opening angle between the axis and the surface of the cone.
       */
      FT angle() const {
          return m_angle;
      }
      
      /*!
         The apex of the cone.
              */
      Point apex() const {
          return m_apex;
      }
      
      /*!
       The axis points from the apex into the cone.
       */
      Vector axis() const {
          return m_axis;
      }
      
      /*!
       Helper function to write apex, axis and angle of the cone and
       number of assigned points into a string.
       */
      std::string info() const {
          std::stringstream sstr;
          
          sstr << "Type: cone apex: (" << m_apex.x() << ", " << m_apex.y();
          sstr << ", " << m_apex.z() << ") axis: (" << m_axis.x() << ", ";
          sstr << m_axis.y() << ", " << m_axis.z() << ") angle:" << m_angle;
          sstr << " #Pts: " << this->m_indices.size();
          
          return sstr.str();
      }

      /*!
      Computes squared Euclidean distance from query point to the shape.
      */ 
      FT squared_distance(const Point &_p) const {
        Vector toApex = _p - m_apex;
        FT a = toApex.squared_length();

        // projection on axis
        FT b = toApex * m_axis;

        // distance to axis
        FT l = sqrt(a - b * b);
        FT c = m_cosAng * l;
        FT d = m_nSinAng * b;

        // far on other side?
        return (b < 0 && c - d < 0) ? a : abs(c + d) * abs(c + d);
      }

  protected:
      /// \cond SKIP_IN_MANUAL
      virtual void create_shape(const std::vector<size_t> &indices) {
      Point p1 = get(this->m_point_pmap, *(this->m_first + indices[0]));
      Point p2 = get(this->m_point_pmap, *(this->m_first + indices[1]));
      Point p3 = get(this->m_point_pmap, *(this->m_first + indices[2]));

      Vector n1 = get(this->m_normal_pmap, *(this->m_first + indices[0]));
      Vector n2 = get(this->m_normal_pmap, *(this->m_first + indices[1]));
      Vector n3 = get(this->m_normal_pmap, *(this->m_first + indices[2]));

      // first calculate intersection of three planes -> apex

      Vector lineDir = CGAL::cross_product(n1, n2);
      lineDir = lineDir * 1.0 / (sqrt(lineDir.squared_length()));

      // lineDir not normalized direction of intersection lines
      //  of two planes (p1, n1) and (p2, n2)
      // get point on line by moving point p1 onto line
      Vector orthLineInPlane = CGAL::cross_product(n1, lineDir);
      orthLineInPlane = orthLineInPlane * 
                        1.0 / (sqrt(orthLineInPlane.squared_length()));

      // distance of p1 to (p2, n2)
      FT d = (p1 - CGAL::ORIGIN) * n2 - (p2 - CGAL::ORIGIN) * n2;
      // projection of orthLineInPlane onto p2
      FT l = orthLineInPlane * n2;
      Point pointOnLine = p1 - (d/l) * orthLineInPlane;


      // checking
      d = (pointOnLine - p1) * n1; // should be 0
      d = (pointOnLine - p2) * n2; // should be 0


      // distance of pLineDir to (p3, n3)
      d = (pointOnLine - CGAL::ORIGIN) * n3 - (p3 - CGAL::ORIGIN) * n3;
      l = lineDir * n3;
      m_apex = pointOnLine - (d/l) * lineDir;


      // checking
      d = (m_apex - p1) * n1; // should be 0
      d = (m_apex - p2) * n2; // should be 0
      d = (m_apex - p3) * n3; // should be 0

      // 2. find axis
      Vector v1 = p1 - m_apex;
      v1 = v1 * 1.0 / (sqrt(v1.squared_length()));
      Point c1 = m_apex + v1;

      Vector v2 = p2 - m_apex;
      v2 = v2 * 1.0 / (sqrt(v2.squared_length()));
      Point c2 = m_apex + v2;

      Vector v3 = p3 - m_apex;
      v3 = v3 * 1.0 / (sqrt(v3.squared_length()));
      Point c3 = m_apex + v3;

      m_axis = CGAL::cross_product(c1 - c2, c1 - c3);
      m_axis = (orthLineInPlane * m_axis < 0) ? -m_axis : m_axis;
      m_axis = m_axis * 1.0 / sqrt(m_axis.squared_length());

      m_angle = acos(v1 * m_axis) + acos(v2 * m_axis) + acos(v3 * m_axis);
      m_angle /= 3;
      if (m_angle < 0 || m_angle > M_PI / 2.12)
        return;

      m_nSinAng = -sin(m_angle);
      m_cosAng = cos(m_angle);

      this->m_isValid = true;
    }

    void parameters(std::vector<std::pair<FT, FT> > &parameterSpace,
                    const std::vector<size_t> &indices,
                    FT min[2],
                    FT max[2]) const {
        Vector a = CGAL::cross_product(Vector(0, 1, 0), m_axis);
        if (a.squared_length() < 0.01)
          a = CGAL::cross_product(Vector(0, 0, 1), m_axis);
        a = a * 1.0 / sqrt(a.squared_length());
        Vector b = CGAL::cross_product(m_axis, a);
        b = b * 1.0 / sqrt(b.squared_length());
    }

    void squared_distance(std::vector<FT> &dists,
      const std::vector<int> &shapeIndex,
      const std::vector<size_t> &indices) {

      for (size_t i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {
          Vector toApex = get(this->m_point_pmap,
                              *(this->m_first + indices[i])) - m_apex;

          FT a = toApex.squared_length();

          // projection on axis
          FT b = toApex * m_axis;

          // distance to axis
          FT l = sqrt(a - b * b);
          FT c = m_cosAng * l;
          FT d = m_nSinAng * b;

          // far on other side?
          dists[i] = (b < 0 && c - d < 0) ? a : abs(c + d) * abs(c + d);
          if (dists[i] > 0.01) {
            int asd;
            asd = 2;
          }
        }
      }
    }

    void cos_to_normal(std::vector<FT> &angles,
                       const std::vector<int> &shapeIndex, 
                       const std::vector<size_t> &indices) const {
      for (size_t i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {
          // construct vector orthogonal to axis in direction of the point
          Vector a = get(this->m_point_pmap, 
                         *(this->m_first + indices[i])) - m_apex;

          Vector b = CGAL::cross_product(m_axis, 
                                         CGAL::cross_product(m_axis, a));
          b = (a * b < 0) ? -b : b;
          b = b * 1.0 / sqrt(b.squared_length());
          b = m_cosAng * b + m_nSinAng * m_axis;

          angles[i] = abs(get(this->m_normal_pmap, 
                              *(this->m_first + indices[i])) * b);
        }
      }
    }

    FT cos_to_normal(const Point &_p, const Vector &_n) const {
      // construct vector orthogonal to axis in direction of the point
      Vector a = _p - m_apex;
      Vector b = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, a));
      b = (a * b < 0) ? -b : b;
      b = b * 1.0 / sqrt(b.squared_length());
      b = m_cosAng * b + m_nSinAng * m_axis;

      return abs(_n * b);
    }
      
      virtual size_t required_samples() const {
          return 3;
      }

    virtual bool supports_connected_component() const {
      return false;
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
    FT m_angle;
    Point m_apex;
    Vector m_axis;
    Point m_pointOnPrimitive;
    FT m_nSinAng, m_cosAng;
      /// \endcond
  };
}
#endif
