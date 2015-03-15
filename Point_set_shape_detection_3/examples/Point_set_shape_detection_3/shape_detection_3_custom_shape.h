#ifndef MY_PLANE_SHAPE_H
#define MY_PLANE_SHAPE_H

#include "Shape_base.h"


namespace CGAL {
    /*
     My_Plane derives from Shape_base. The plane is defined by
	 its normal vector and distance to the origin.
     */
  template <class ERTraits>
  class My_Plane : public Shape_base<ERTraits> {
  public:
    typedef typename ERTraits::Geom_traits::FT FT;///< number type.
    typedef typename ERTraits::Geom_traits::Point_3 Point;///< point type.

  public:
    My_Plane() : Shape_base<ERTraits>() {}

    //  Computes squared Euclidean distance from query point to the shape.
    virtual FT squared_distance(const Point &p) const {
      const FT sd = (p - m_point_on_primitive) * m_normal;
      return sd * sd;
    }

  protected:
    // Constructs shape based on minimal set of samples from the input data.    
    virtual void create_shape(const std::vector<size_t> &indices) {
      const Point p1 = this->get_point(indices[0]);
      const Point p2 = this->get_point(indices[1]);
      const Point p3 = this->get_point(indices[2]);

      m_normal = CGAL::cross_product(p1 - p2, p1 - p3);

      m_normal = m_normal * (1.0 / sqrt(m_normal.squared_length()));
      m_d = -(p1[0] * m_normal[0] + p1[1] * m_normal[1] + p1[2] * m_normal[2]);
	  
	  m_isValid = true;
    }

    // Computes squared Euclidean distance from a set of points.
	virtual void squared_distance(std::vector<FT> &dists,
                          const std::vector<size_t> &indices) {
      for (size_t i = 0; i < indices.size(); i++) {
        const FT sd = (this->get_point(indices[i])
                       - m_point_on_primitive) * m_normal;
        dists[i] = sd * sd;
      }
    }
    
    /*
      Computes the normal deviation between shape and
	  a set of points with normals.
     */
    virtual void cos_to_normal(std::vector<FT> &angles,
                       const std::vector<size_t> &indices) const {
      for (size_t i = 0; i < indices.size(); i++)
        angles[i] = abs(this->get_normal(indices[i]) * m_normal);
    }
    
    // Returns the number of required samples for construction.
    virtual size_t required_samples() const {
      return 3;
    }

    // Returns a string with shape parameters.
    virtual std::string info() const {
      std::stringstream sstr;
      sstr << "Type: plane (" << m_normal.x() << ", " 
        << m_normal.y() << ", " << m_normal.z() << ")x - " <<
        m_d << " = 0" << " #Pts: " << this->m_indices.size();

      return sstr.str();
    }

  private:
    Point m_point_on_primitive;
    Vector m_normal;
    FT m_d;
  };
}
#endif
