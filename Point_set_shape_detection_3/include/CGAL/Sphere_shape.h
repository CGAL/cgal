#ifndef CGAL_SHAPE_DETECTION_3_SPHERE_SHAPE_H
#define CGAL_SHAPE_DETECTION_3_SPHERE_SHAPE_H

#include "Shape_base.h"
#include <set>

/*!
 \file Sphere_shape.h
 */

namespace CGAL {
    /*!
     \brief Sphere_shape implements Shape_base. The sphere is parameterized by its center and the radius.
     */
  template <class Sd_traits>
  class Sphere_shape : public Shape_base<Sd_traits> {
  public:
    typedef typename Sd_traits::Input_iterator Input_iterator;///< random access iterator for input data.
    typedef typename Sd_traits::Geom_traits::FT FT; ///< number type.
    //typedef typename Sd_traits::Geom_traits::Line_3 Line;
    typedef typename Sd_traits::Geom_traits::Point_3 Point;///< point type.
    typedef typename Sd_traits::Geom_traits::Vector_3 Vector;///< vector type.
    //typedef typename Sd_traits::Geom_traits::Plane_3 Plane;
    typedef typename Sd_traits::Geom_traits::Sphere_3 Sphere;///< sphere type.
    typedef typename Sd_traits::Point_pmap Point_pmap;   ///< property map to access the location of an input point.
    typedef typename Sd_traits::Normal_pmap Normal_pmap; ///< property map to access the unoriented normal of an input point.

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
      
      /*!
       Access to the radius of the sphere.
       */
    FT radius() const {
      return m_radius;
    }

  protected:
      /// \cond SKIP_IN_MANUAL
      void create_shape(const std::vector<size_t> &indices) {
      Point p1 = get(this->m_point_pmap, *(this->m_first + indices[0]));
      Point p2 = get(this->m_point_pmap, *(this->m_first + indices[1]));
      Point p3 = get(this->m_point_pmap, *(this->m_first + indices[2]));

      Vector n1 = get(this->m_normal_pmap, *(this->m_first + indices[0]));
      Vector n2 = get(this->m_normal_pmap, *(this->m_first + indices[1]));
      Vector n3 = get(this->m_normal_pmap, *(this->m_first + indices[2]));


      // Determine center: select midpoint of shortest line segment
      //  between p1 and p2
      // implemented from "3D game engine design" by Eberly 2001

      Vector diff = p1 - p2;
      FT a = n1 * n1;
      FT b = -(n1 * n2);
      FT c = n2 * n2;
      FT d = n1 * diff;

      FT det = abs(a * c - b * b);

      // parallel?
      if (det < 0.00001) {
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

      if (abs(d1-d2) > 2 * this->m_epsilon) {
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

      m_radius = (d1 + d2) * 0.5;

      if (abs(d3 - m_radius) > this->m_epsilon ||
          n3 * v3 < this->m_normal_threshold) {
        this->m_isValid = false;
        return;
      }

      m_sphere = Sphere(center, m_radius * m_radius);
    }

    void parameters(std::vector<std::pair<FT, FT> > &parameterSpace,
                    const std::vector<size_t> &indices, FT min[2], 
                    FT max[2]) const {
    }

    void parameter_extend(const Point &center, 
                          FT width, 
                          FT min[2],
                          FT max[2]) const {
    }

    FT squared_distance(const Point &_p) const {
      FT d = sqrt((m_sphere.center() - _p).squared_length()) - m_radius;
      return d*d;
    }

    void squared_distance(std::vector<FT> &dists,
      const std::vector<int> &shapeIndex,
      const std::vector<size_t> &indices) {

      for (size_t i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {
          dists[i] = sqrt((m_sphere.center()
            - get(this->m_point_pmap,
                  *(this->m_first + indices[i]))).squared_length())
            - m_radius;

          dists[i] = dists[i] * dists[i];
        }
      }
    }

    void cos_to_normal(std::vector<FT> &angles,
                       const std::vector<int> &shapeIndex, 
                       const std::vector<size_t> &indices) const {
      for (size_t i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {
          Vector n = m_sphere.center()
            - get(this->m_point_pmap, *(this->m_first + indices[i]));

          n = n * (1.0 / (sqrt(n.squared_length())));
          angles[i] = abs(get(this->m_normal_pmap,
                              *(this->m_first + indices[i])) * n);
        }
      }
    }

    FT cos_to_normal(const Point &_p, const Vector &_n) const {
      Vector n = m_sphere.center() - _p;
      n = n * (1.0 / (sqrt(n.squared_length())));
      return abs(_n * n);
    }
      
      virtual size_t required_samples() const {
          return 3;
      }

    // U is longitude
    virtual bool supports_connected_component() const {
      return false;
    }

    virtual bool wraps_u() const {
      return true;
    }

    virtual bool wraps_v() const {
      return false;
    }

  private:
    Sphere m_sphere;
    FT m_radius;
/// \endcond
  };
}
#endif
