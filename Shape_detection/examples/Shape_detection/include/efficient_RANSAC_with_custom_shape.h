#ifndef MY_PLANE_SHAPE_H
#define MY_PLANE_SHAPE_H

#include <CGAL/number_utils.h>
#include <CGAL/Shape_detection/Efficient_RANSAC.h>

// My_Plane is derived from Shape_base. The plane is represented by
// its normal vector and distance to the origin.
template <class Traits>
class My_Plane : public CGAL::Shape_detection::Shape_base<Traits> {

public:

  typedef typename Traits::FT FT;           // number type
  typedef typename Traits::Point_3 Point;   // point type
  typedef typename Traits::Vector_3 Vector; // vector type

  My_Plane() :
  CGAL::Shape_detection::Shape_base<Traits>()
  { }

  // Compute squared Euclidean distance from query point to the shape.
  virtual FT squared_distance(const Point& p) const {

    const FT sd = (this->constr_vec(m_point_on_primitive, p)) * m_normal;
    return sd * sd;
  }

  Vector plane_normal() const {
    return m_normal;
  }

  FT d() const {
    return m_d;
  }

  // Return a string with shape parameters.
  virtual std::string info() const {

    std::stringstream sstr;

    sstr << "Type: plane (" << this->get_x(m_normal) << ", "
    << this->get_y(m_normal) << ", " << this->get_z(m_normal) << ")x - " <<

    m_d << " = 0" << " #Pts: " << this->m_indices.size();

    return sstr.str();
  }

protected:

  // Construct shape base on a minimal set of samples from the input data.
  virtual void create_shape(const std::vector<std::size_t>& indices) {

    const Point p1 = this->point(indices[0]);
    const Point p2 = this->point(indices[1]);
    const Point p3 = this->point(indices[2]);

    m_normal = this->cross_pdct(p1 - p2, p1 - p3);

    m_normal = m_normal * (1.0 / sqrt(this->sqlen(m_normal)));
    m_d = -(p1[0] * m_normal[0] + p1[1] * m_normal[1] + p1[2] * m_normal[2]);

    m_point_on_primitive = p1;
    this->m_is_valid = true;
  }

  // Compute squared Euclidean distance from a set of points.
  virtual void squared_distance(
    const std::vector<std::size_t>& indices,
    std::vector<FT>& dists) const {

    for (std::size_t i = 0; i < indices.size(); ++i) {

      const FT sd = (this->point(indices[i]) - m_point_on_primitive) * m_normal;
      dists[i] = sd * sd;
    }
  }

  // Compute the normal deviation between a shape and
  // a set of points with normals.
  virtual void cos_to_normal(
    const std::vector<std::size_t>& indices,
    std::vector<FT>& angles) const {

    for (std::size_t i = 0; i < indices.size(); ++i)
      angles[i] = CGAL::abs(this->normal(indices[i]) * m_normal);
  }

  // Return the number of required samples for construction.
  virtual std::size_t minimum_sample_size() const {
    return 3;
  }

private:

  Point  m_point_on_primitive;
  Vector m_normal;
  FT     m_d;

};

#endif // MY_PLANE_SHAPE_H
