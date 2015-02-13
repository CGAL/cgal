#ifndef MY_PLANE_SHAPE_H
#define MY_PLANE_SHAPE_H

#include "Shape_base.h"

/*!
 \file My_Plane.h
 */
namespace CGAL {
    /*!
     \brief My_Plane derives from Shape_base. The plane is parameterized by the normal vector and the distance to the origin.
     */
  template <class Sd_traits>
  class My_Plane : public Shape_base<Sd_traits> {
  public:
    typedef typename Sd_traits::Geom_traits::FT FT;///< number type.
    typedef typename Sd_traits::Geom_traits::Point_3 Point;///< point type.

  public:
    My_Plane() : Shape_base<Sd_traits>() {}
            
    /*!
      Computes the squared Euclidean distance from the query point to the shape.
      */
    FT squared_distance(const Point &p) const {
      const FT sd = (p - m_point_on_primitive) * m_normal;
      return sd * sd;
    }

  protected:

    /*!
      Constructs the shape based on a minimal set of samples from the input data.
     */          
    virtual void create_shape(const std::vector<size_t> &indices) {
      const Point p1 = this->get_point(indices[0]);
      const Point p2 = this->get_point(indices[1]);
      const Point p3 = this->get_point(indices[2]);

      m_normal = CGAL::cross_product(p1 - p2, p1 - p3);

      m_normal = m_normal * (1.0 / sqrt(m_normal.squared_length()));
      m_d = -(p1[0] * m_normal[0] + p1[1] * m_normal[1] + p1[2] * m_normal[2]);
	  
	  m_isValid = true;
    }

    /*!
      Computes squared Euclidean distance from a set of points.
     */    
    void squared_distance(std::vector<FT> &dists, const std::vector<size_t> &indices) {
      for (size_t i = 0; i < indices.size(); i++) {
        const FT sd = (this->get_point(indices[i]) - m_point_on_primitive) * m_normal;
        dists[i] = sd * sd;
      }
    }

    void cos_to_normal(std::vector<FT> &angles, const std::vector<size_t> &indices) const {
      for (size_t i = 0; i < indices.size(); i++)
        angles[i] = abs(this->get_normal(indices[i]) * m_normal);
    }

    virtual size_t required_samples() const {
      return 3;
    } 
    
    std::string info() const {
      std::stringstream sstr;
      sstr << "Type: plane (" << m_normal.x() << ", " 
		              << m_normal.y() << ", " 
			      << m_normal.z() << ")x - " << m_d << " = 0" << " #Pts: " << this->m_indices.size();

      return sstr.str();
    }

  private:
    Point m_point_on_primitive;
    Vector m_normal;
    FT m_d;
  };
}
#endif
