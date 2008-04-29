#ifndef POINT_SET_3_H
#define POINT_SET_3_H

#include <CGAL/Point_with_normal_3.h>
#include <Gyroviz_point_3.h>
#include <CGAL/surface_reconstruction_assertions.h>

#include <CGAL/boost/graph/properties.h>

#include <algorithm>
#include <vector>
#include <GL/gl.h>


/// The Point_set_3 class is array of points + normals of type PointWithNormal_3
/// (in fact Gyroviz_point_3 to support algorithms specific to Gyroviz). 
/// It provides:
/// - accessors: points and normals iterators, property maps
/// - OpenGL rendering
/// - bounding box
/// 
/// CAUTION: invalidate_bounding_box() must be called 
/// after modifying the points.
/// 
/// @heading Parameters:
/// @param Gt       Geometric traits class.

template <class Gt>
class Point_set_3 : public std::vector<Gyroviz_point_3<Gt> >
{
// Private types
private:

  // Base class 
  typedef std::vector<Gyroviz_point_3<Gt> > Base;

  // Auxiliary class to build a normals iterator
  template <class Node>
  struct Project_normal {
    typedef Node                  argument_type;
    typedef typename Node::Normal Normal;
    typedef Normal                result_type;
    typedef CGAL::Arity_tag<1>    Arity;
    Normal&       operator()(Node& x)       const { return x.normal(); }
    const Normal& operator()(const Node& x) const { return x.normal(); }
  };

// Public types
public:

  // Repeat base class' types
  /// @cond SKIP_IN_MANUAL
  using Base::iterator;
  using Base::const_iterator;
  /// @endcond

  // Geometric types
  typedef Gt  Geom_traits; ///<Geometric traits class.
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Vector_3 Vector;
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
  typedef typename Geom_traits::Sphere_3 Sphere;
  
  typedef Gyroviz_point_3<Gt> Point_with_normal; ///<Model of PointWithNormal_3
  typedef Gyroviz_point_3<Gt> Gyroviz_point;     ///<Model of PointWithNormal_3 + cameras
  typedef typename Point_with_normal::Normal Normal; ///<Model of OrientedNormal_3 concept.

  // Iterator over points
  typedef std::vector<Point_with_normal>::iterator 
                                    Point_iterator;      
  typedef std::vector<Point_with_normal>::const_iterator 
                                    Point_const_iterator;      

  // Iterator over normals
  typedef CGAL::Iterator_project<iterator, 
                                  Project_normal<Point_with_normal> >  
                                    Normal_iterator;      
  typedef CGAL::Iterator_project<const_iterator, 
                                  Project_normal<Point_with_normal> >  
                                    Normal_const_iterator;      

// Data members
private:

    // Indicate if m_barycenter, m_bounding_box and m_diameter_standard_deviation below are valid
    mutable bool m_bounding_box_is_valid;
    mutable Iso_cuboid_3 m_bounding_box; // m_points's bounding box
    mutable Point m_barycenter; // m_points's barycenter
    mutable FT m_diameter_standard_deviation; // m_points's standard deviation

// Public methods
public:

  /// Default constructor.
  Point_set_3()
  {
    m_bounding_box_is_valid = false;
  }

  // Default copy constructor and operator =() are fine.

  // Repeat base class' public methods used below
  /// @cond SKIP_IN_MANUAL
  Base::begin;
  Base::end;
  Base::size;
  /// @endcond

  // Get first/last iterators over normals.
  Normal_iterator normals_begin()             { return Normal_iterator(begin()); }
  Normal_const_iterator normals_begin() const { return Normal_iterator(begin()); }
  Normal_iterator normals_end()               { return Normal_iterator(end()); }
  Normal_const_iterator normals_end() const   { return Normal_iterator(end()); }

  /// Get the bounding box.
  Iso_cuboid_3 bounding_box() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_bounding_box;
  }

  /// Get bounding sphere.
  Sphere bounding_sphere() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    // Center point
    FT mx = 0.5 * (m_bounding_box.xmax() + m_bounding_box.xmin());
    FT my = 0.5 * (m_bounding_box.ymax() + m_bounding_box.ymin());
    FT mz = 0.5 * (m_bounding_box.zmax() + m_bounding_box.zmin());
    Point center(mx,my,mz);

    // Squared radius
    FT dx = m_bounding_box.xmax() - m_bounding_box.xmin();
    FT dy = m_bounding_box.ymax() - m_bounding_box.ymin();
    FT dz = m_bounding_box.zmax() - m_bounding_box.zmin();
    FT squared_radius = dx*dx + dy*dy + dz*dz;

    return Sphere(center, squared_radius);
  }

  /// Get points barycenter.
  Point barycenter() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_barycenter;
  }

  /// Get the standard deviation of the distance to barycenter.
  FT diameter_standard_deviation() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_diameter_standard_deviation;
  }

  // Get the region of interest, ignoring the outliers.
  // This method is used to define the OpenGL arcball sphere.
  Sphere region_of_interest() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    // A good candidate is a sphere containing the dense region of the point cloud:
    // - center point is barycenter
    // - Radius is 2 * standard deviation
    float radius = 2.f * (float)m_diameter_standard_deviation;
    return Sphere(m_barycenter, radius*radius);
  }

  /// Update barycenter, bounding box, bounding sphere and standard deviation.
  /// Owner is responsible to call this function after modifying the triangulation.
  void invalidate_bounding_box()
  {
    m_bounding_box_is_valid = false;
  }

  // Draw m_points[] points using OpenGL calls.
  void gl_draw_vertices(unsigned char r, unsigned char g, unsigned char b,
                        float size) const
  {
    ::glPointSize(size);
    ::glColor3ub(r,g,b);
    ::glBegin(GL_POINTS);
    for (Point_const_iterator it = begin(); it != end(); it++)
    {
      const Point& p = *it;
      ::glVertex3d(p.x(),p.y(),p.z());
    }
    ::glEnd();
  }

  // Draw m_points[] normals using OpenGL calls.
  void gl_draw_normals(unsigned char r, unsigned char g, unsigned char b,
                       FT scale = 1.0) const
  {
    // Draw *oriented* normals
    ::glColor3ub(r,g,b);
    ::glBegin(GL_LINES);
    for (Point_const_iterator it = begin(); it != end(); it++)
    {
      const Point& p = *it;
      Normal n = it->normal();
      if ( n.is_normal_oriented() && n.get_vector() != CGAL::NULL_VECTOR )
      {
        Point q = p + scale * n.get_vector();
        glVertex3d(p.x(),p.y(),p.z());
        glVertex3d(q.x(),q.y(),q.z());
      }
    }
    ::glEnd();

    // Draw *non-oriented* normals
    ::glColor3ub(255,0,0);
    ::glBegin(GL_LINES);
    for (Point_const_iterator it = begin(); it != end(); it++)
    {
      const Point& p = *it;
      Normal n = it->normal();
      if ( !n.is_normal_oriented() && n.get_vector() != CGAL::NULL_VECTOR )
      {
        Point q = p + scale * n.get_vector();
        glVertex3d(p.x(),p.y(),p.z());
        glVertex3d(q.x(),q.y(),q.z());
      }
    }
    ::glEnd();
  }

// Private methods:
private:

  // Recompute barycenter, bounding box, bounding sphere and standard deviation.
  void update_bounding_box() const
  {
    if (begin() == end())
      return;

    // Update bounding box and barycenter.
    // TODO: we should use the functions in PCA component instead.
    FT xmin,xmax,ymin,ymax,zmin,zmax;
    xmin = ymin = zmin =  1e38;
    xmax = ymax = zmax = -1e38;
    Vector v = CGAL::NULL_VECTOR;
    FT norm = 0;
    for (Point_const_iterator it = begin(); it != end(); it++)
    {
      const Point& p = *it;

      // update bbox
      xmin = (std::min)(p.x(),xmin);
      ymin = (std::min)(p.y(),ymin);
      zmin = (std::min)(p.z(),zmin);
      xmax = (std::max)(p.x(),xmax);
      ymax = (std::max)(p.y(),ymax);
      zmax = (std::max)(p.z(),zmax);

      // update barycenter
      v = v + (p - CGAL::ORIGIN);
      norm += 1;
    }
    //
    Point p(xmin,ymin,zmin);
    Point q(xmax,ymax,zmax);
    m_bounding_box = Iso_cuboid_3(p,q);
    //
    m_barycenter = CGAL::ORIGIN + v / norm;

    // Compute standard deviation of the distance to barycenter
    typename Geom_traits::Compute_squared_distance_3 sqd;
    FT sq_radius = 0;
    for (Point_const_iterator it = begin(); it != end(); it++)
    {
      sq_radius += sqd(*it, m_barycenter);
    }
    sq_radius /= size();
    m_diameter_standard_deviation = CGAL::sqrt(sq_radius);

    m_bounding_box_is_valid = true;
  }

}; // end of class Point_set_3


/// Helper class: type of the "vertex_point" property map
/// of an Point_set_3 object.
template <class Gt>
class Point_set_vertex_point_const_map 
{
public:
    typedef Point_set_3<Gt> Point_set;
    typedef typename Gt::Point_3 Point_3;  

    // Property maps required types
    typedef boost::readable_property_map_tag          category;
    typedef Point_3                                   value_type;
    typedef value_type                                reference;
    typedef typename Point_set::Point_const_iterator  key_type;

    Point_set_vertex_point_const_map(const Point_set&) {}

    /// Free function to access the map elements.
    friend inline 
    reference get(const Point_set_vertex_point_const_map&, key_type p)
    {
      return *p;
    }
};

/// Free function to get the "vertex_point" property map
/// of an Point_set_3 object.
template <class Gt>
inline
Point_set_vertex_point_const_map<Gt> 
get(CGAL::vertex_point_t, const Point_set_3<Gt>& points) 
{
  Point_set_vertex_point_const_map<Gt> aMap(points);
  return aMap;
}


/// Helper class: type of the "vertex_normal" property map
/// of an Point_set_3 object.
template <class Gt>
class Point_set_vertex_normal_map 
  : public boost::put_get_helper<typename Point_set_3<Gt>::Point_with_normal::Normal&, 
                                  Point_set_vertex_normal_map<Gt> >
{
public:
    typedef Point_set_3<Gt> Point_set;
    typedef typename Point_set::Point_with_normal::Normal Normal;  

    // Property maps required types
    typedef boost::lvalue_property_map_tag            category;
    typedef Normal                                    value_type;
    typedef Normal&                                   reference;
    typedef typename Point_set::Point_iterator        key_type;

    Point_set_vertex_normal_map(const Point_set&) {}

    /// Access the map elements.
    reference operator[](key_type p) const { return p->normal(); }
};

/// Free function to get the "vertex_normal" property map
/// of an Point_set_3 object.
template <class Gt>
inline
Point_set_vertex_normal_map<Gt> 
get(boost::vertex_normal_t, const Point_set_3<Gt>& points) 
{
  Point_set_vertex_normal_map<Gt> aMap(points);
  return aMap;
}


/// Helper class: type of the "vertex_index" property map
/// of an Point_set_3 object.
template <class Gt>
class Point_set_vertex_index_const_map 
{
public:
    typedef Point_set_3<Gt> Point_set;

    // Property maps required types
    typedef boost::readable_property_map_tag          category;
    typedef unsigned int                              value_type;
    typedef value_type                                reference;
    typedef typename Point_set::Point_const_iterator  key_type;

    Point_set_vertex_index_const_map(const Point_set& points) 
    : m_points(points) 
    {}

    /// Free function to access the map elements.
    friend inline
    reference get(const Point_set_vertex_index_const_map& map, key_type p)
    {
      return std::distance(map.m_points.begin(), p);
    }
    
private:
  const Point_set& m_points;
};

/// Free function to get the "vertex_index" property map
/// of an Point_set_3 object.
template <class Gt>
inline
Point_set_vertex_index_const_map<Gt> 
get(boost::vertex_index_t, const Point_set_3<Gt>& points) 
{
  Point_set_vertex_index_const_map<Gt> aMap(points);
  return aMap;
}


/// Helper class: type of the "vertex_cameras" property map
/// of a Point_set_3 object.
template <class Gt>
class Point_set_vertex_cameras_const_map 
{
public:
    typedef Point_set_3<Gt> Point_set;
    typedef typename Point_set::Gyroviz_point::Camera_const_iterator 
                                                      Camera_const_iterator;  

    // Property maps required types
    typedef boost::readable_property_map_tag          category;
    typedef std::pair<Camera_const_iterator,Camera_const_iterator>  
                                                      value_type;
    typedef value_type                                reference;
    typedef typename Point_set::Point_const_iterator  key_type;

    Point_set_vertex_cameras_const_map(const Point_set&) {}

    /// Free function to access the map elements.
    friend inline
    reference get(const Point_set_vertex_cameras_const_map&, key_type p)
    {
      return std::make_pair(p->cameras_begin(), p->cameras_end());
    }
};

/// Free function to get the "vertex_cameras" property map
/// of a Point_set_3 object.
template <class Gt>
inline
Point_set_vertex_cameras_const_map<Gt> 
get(boost::vertex_cameras_t, const Point_set_3<Gt>& points) 
{
    Point_set_vertex_cameras_const_map<Gt> aMap(points);
    return aMap;
}


#endif // POINT_SET_3_H
