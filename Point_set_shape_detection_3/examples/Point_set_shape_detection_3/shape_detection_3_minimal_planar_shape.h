#ifndef MY_PLANE_SHAPE_H
#define MY_PLANE_SHAPE_H

#include "Shape_base.h"
#include <set>

/*!
 \file Plane_shape.h
 */
namespace CGAL {
    /*!
     \brief Plane_shape implements Shape_base. The plane is parameterized by the normal vector and the distance to the origin.
     */
  template <class Sd_traits>
  class Plane_shape : public Shape_base<Sd_traits> {
  public:
    typedef typename Sd_traits::Geom_traits::FT FT;///< number type.
    typedef typename Sd_traits::Geom_traits::Point_3 Point;///< point type.
    typedef typename Sd_traits::Geom_traits::Plane_3 Plane_3;///< plane type for conversion operator.

  public:
    Plane_shape() : Shape_base<Sd_traits>() {}
            
    /*!
      Provides the squared Euclidean distance of the point to the shape.
      */
    FT squared_distance(const Point &_p) const {
      FT d = (_p - m_point_on_primitive) * m_normal;
      return d * d;
    }

  protected:

    /*!
      Constructs the shape based on a minimal set of samples from the input data.
     */          
    virtual void create_shape(const std::vector<size_t> &indices) {
      Point p1 = get(this->m_point_pmap, *(this->m_first + indices[0]));
      Point p2 = get(this->m_point_pmap, *(this->m_first + indices[1]));
      Point p3 = get(this->m_point_pmap, *(this->m_first + indices[2]));

      m_normal = CGAL::cross_product(p1 - p2, p1 - p3);

      m_normal = m_normal * (1.0 / sqrt(m_normal.squared_length()));
      m_d = -(p1[0] * m_normal[0] + p1[1] * m_normal[1] + p1[2] * m_normal[2]);
	  
	  m_isValid = true;
    }

    /*!
      Provides the squared Euclidean distance of a set of points.
     */    
    void squared_distance(std::vector<FT> &dists, const std::vector<int> &shapeIndex, const std::vector<size_t> &indices) {
      for (size_t i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1) {
          FT d = (get(this->m_point_pmap, *(this->m_first + indices[i])) - m_point_on_primitive) * m_normal;
          dists[i] = d * d;
        }
      }
    }

    void cos_to_normal(std::vector<FT> &angles, const std::vector<int> &shapeIndex, const std::vector<size_t> &indices) const {
      for (size_t i = 0;i<indices.size();i++) {
        if (shapeIndex[indices[i]] == -1)
          angles[i] = abs(get(this->m_normal_pmap, *(this->m_first + indices[i])) * m_normal);
      }
    }
    
    virtual size_t required_samples() const {
      return 3;
    } 
    
    std::string info() const {
      std::stringstream sstr;
      sstr << "Type: plane (" << m_normal.x() << ", " << m_normal.y() << ", " << m_normal.z() << ")x - " << m_d << " = 0" << " #Pts: " << this->m_indices.size();

      return sstr.str();
    }

  private:
    Point m_point_on_primitive;
    Vector m_normal;
    FT m_d;
  };
}
#endif
