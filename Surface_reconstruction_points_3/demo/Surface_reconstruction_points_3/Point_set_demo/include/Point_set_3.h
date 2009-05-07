// Author: Laurent Saboret

#ifndef POINT_SET_3_H
#define POINT_SET_3_H

#include <CGAL/Point_with_normal_3.h>
#include <CGAL/point_set_property_map.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>

#include <UI_point_3.h>

#include <algorithm>
#include <deque>

#ifdef CGAL_GLEW_ENABLED
# include <GL/glew.h>
#else
# include <CGAL/gl.h>
#endif


/// The Point_set_3 class is array of points + normals of type
/// Point_with_normal_3<Gt> (in fact
/// UI_point_3 to support a selection flag and an optional radius).
/// It provides:
/// - accessors: points and normals iterators, property maps
/// - OpenGL rendering
/// - bounding box
///
/// CAUTION:
/// - User is responsible to call invalidate_bounds() after adding, moving or removing points.
/// - User is responsible to partition the point set as oriented/unoriented normals 
///   and to set unoriented_points_begin() appropriately.
///
/// @heading Parameters:
/// @param Gt       Geometric traits class.

template <class Gt>
class Point_set_3 : public std::deque<UI_point_3<Gt> >
{
// Private types
private:

  // Base class
  typedef std::deque<UI_point_3<Gt> > Base;

  // Auxiliary class to build a normals iterator
  template <class Node>
  struct Project_normal {
    typedef Node                  argument_type;
    typedef typename Node::Vector Vector;
    typedef Vector                result_type;
    Vector&       operator()(Node& x)       const { return x.normal(); }
    const Vector& operator()(const Node& x) const { return x.normal(); }
  };

// Public types
public:

  // Repeat base class' types
  /// @cond SKIP_IN_MANUAL
  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator const_iterator;
  /// @endcond

  // Classic CGAL geometric types
  typedef Gt  Geom_traits; ///<Geometric traits class.
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;  ///< == Point_3<Geom_traits>
  typedef typename Geom_traits::Vector_3 Vector; ///< == Vector_3<Geom_traits>
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;
  typedef typename Geom_traits::Sphere_3 Sphere;

  /// Type of points in Point_set_3
  typedef UI_point_3<Gt> UI_point; ///< Position + normal + selection flag
  // Its superclass:
  typedef typename UI_point::Point_with_normal Point_with_normal; ///< Position + normal

  // Iterator over Point_3 points
  typedef typename std::deque<UI_point>::iterator        Point_iterator;
  typedef typename std::deque<UI_point>::const_iterator  Point_const_iterator;

  // Iterator over normals
  typedef CGAL::Iterator_project<iterator,
                                 Project_normal<UI_point> >
                                                      Normal_iterator;
  typedef CGAL::Iterator_project<const_iterator,
                                 Project_normal<UI_point> >
                                                      Normal_const_iterator;

// Data members
private:

  // Indicate if m_barycenter, m_bounding_box, m_bounding_sphere and
  // m_diameter_standard_deviation below are valid.
  mutable bool m_bounding_box_is_valid;

  mutable Iso_cuboid m_bounding_box; // point set's bounding box
  mutable Sphere m_bounding_sphere; // point set's bounding sphere
  mutable Point m_barycenter; // point set's barycenter
  mutable FT m_diameter_standard_deviation; // point set's standard deviation

  unsigned int m_nb_selected_points; // number of selected points

  // Iterator over the first point with an unoriented normal.
  Point_iterator m_unoriented_points_begin;

// Public methods
public:

  /// Default constructor.
  Point_set_3()
  {
    m_nb_selected_points = 0;
    m_bounding_box_is_valid = false;
    m_unoriented_points_begin = end();
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
  Normal_const_iterator normals_begin() const { return Normal_const_iterator(begin()); }
  Normal_iterator normals_end()               { return Normal_iterator(end()); }
  Normal_const_iterator normals_end() const   { return Normal_const_iterator(end()); }

  // Get/set the iterator over the first point with an unoriented normal.
  /// User is responsible to partition the point set as oriented/unoriented normals 
  /// and to set unoriented_points_begin() appropriately.
  Point_iterator&      unoriented_points_begin()       { return m_unoriented_points_begin; }
  Point_const_iterator unoriented_points_begin() const { return m_unoriented_points_begin; }

  /// Get the number of selected points.
  unsigned int nb_selected_points() const { return m_nb_selected_points; }

  /// Mark a point as selected/not selected.
  void select(UI_point* point, bool is_selected = true)
  {
    if (point->is_selected() != is_selected)
    {
      point->select(is_selected);
      m_nb_selected_points += (is_selected ? 1 : -1);
    }
  }

  /// Mark a range of points as selected/not selected.
  ///
  /// @param first Iterator over first point to select/unselect.
  /// @param beyond Past-the-end iterator.
  void select(iterator first, iterator beyond,
              bool is_selected = true)
  {
    for (iterator it = first; it != beyond; it++)
      it->select(is_selected);

    m_nb_selected_points = std::count_if(begin(), end(),
                                         std::mem_fun_ref(&UI_point::is_selected));
  }

  /// Delete selected points.
  // Note: this call resets unoriented_points_begin().
  void delete_selection()
  {
    // erase-remove idiom
    erase(std::remove_if(begin(), end(), std::mem_fun_ref(&UI_point::is_selected)),
          end());

    // TODO: Update m_unoriented_points_begin. 
    //       This is tricky as erase() invalidates iterators.
    m_unoriented_points_begin = end();

    m_nb_selected_points = 0;
    invalidate_bounds();
  }

  /// Get the bounding box.
  Iso_cuboid bounding_box() const
  {
    if (!m_bounding_box_is_valid)
      update_bounds();

    return m_bounding_box;
  }

  /// Get bounding sphere.
  Sphere bounding_sphere() const
  {
    if (!m_bounding_box_is_valid)
      update_bounds();

    return m_bounding_sphere;
  }

  /// Get points barycenter.
  Point barycenter() const
  {
    if (!m_bounding_box_is_valid)
      update_bounds();

    return m_barycenter;
  }

  /// Get the standard deviation of the distance to barycenter.
  FT diameter_standard_deviation() const
  {
    if (!m_bounding_box_is_valid)
      update_bounds();

    return m_diameter_standard_deviation;
  }

  // Get the region of interest, ignoring the outliers.
  // This method is used to define the OpenGL arcball sphere.
  Sphere region_of_interest() const
  {
    if (!m_bounding_box_is_valid)
      update_bounds();

    // A good candidate is a sphere containing the dense region of the point cloud:
    // - center point is barycenter
    // - Radius is 2 * standard deviation
    float radius = 2.f * (float)m_diameter_standard_deviation;
    return Sphere(m_barycenter, radius*radius);
  }

  /// Update barycenter, bounding box, bounding sphere and standard deviation.
  /// User is responsible to call invalidate_bounds() after adding, moving or removing points.
  void invalidate_bounds()
  {
    m_bounding_box_is_valid = false;
  }

  // Draw points using OpenGL calls.
  // Preconditions: OpenGL point size and color must be set.
  void gl_draw_vertices() const
  {
    // Draw *non-selected* points
    if (m_nb_selected_points < size())
    {
      ::glBegin(GL_POINTS);
      for (const_iterator it = begin(); it != end(); it++)
      {
        const UI_point& p = *it;
        if ( ! p.is_selected() )
          ::glVertex3dv(&p.x());
      }
      ::glEnd();
    }

    // Draw *selected* points
    if (m_nb_selected_points > 0)
    {
      ::glPointSize(4.f);    // selected => bigger
      ::glColor3ub(255,0,0); // selected => red
      ::glBegin(GL_POINTS);
      for (const_iterator it = begin(); it != end(); it++)
      {
        const UI_point& p = *it;
        if (p.is_selected())
          ::glVertex3dv(&p.x());
      }
      ::glEnd();
    }
  }

  // Draw normals using OpenGL calls.
  // Preconditions: OpenGL line width and color must be set.
  void gl_draw_normals(float scale = 1.0) const // scale applied to normal length
  {
    // Draw normals of *non-selected* points
    if (m_nb_selected_points < size())
    {
      // Draw *oriented* normals
      ::glBegin(GL_LINES);
      for (const_iterator it = begin(); it != unoriented_points_begin(); it++)
      {
        const UI_point& p = *it;
        const Vector& n = p.normal();
        if (!p.is_selected())
        {
          Point q = p + scale * n;
          ::glVertex3d(p.x(),p.y(),p.z());
          ::glVertex3d(q.x(),q.y(),q.z());
        }
      }
      ::glEnd();

      // Draw *non-oriented* normals
      ::glColor3ub(245,184,0);       // non oriented => orange
      ::glBegin(GL_LINES);
      for (const_iterator it = unoriented_points_begin(); it != end(); it++)
      {
        const UI_point& p = *it;
        const Vector& n = p.normal();
        if (!p.is_selected())
        {
          Point q = p + scale * n;
          ::glVertex3d(p.x(),p.y(),p.z());
          ::glVertex3d(q.x(),q.y(),q.z());
        }
      }
      ::glEnd();
    }

    // Draw normals of *selected* points
    if (m_nb_selected_points > 0)
    {
      ::glColor3ub(255,0,0); // selected => red
      ::glBegin(GL_LINES);
      for (const_iterator it = begin(); it != end(); it++)
      {
        const UI_point& p = *it;
        const Vector& n = p.normal();
        if (p.is_selected())
        {
          Point q = p + scale * n;
          ::glVertex3d(p.x(),p.y(),p.z());
          ::glVertex3d(q.x(),q.y(),q.z());
        }
      }
      ::glEnd();
    }
  }

  // Draw oriented points with radius using OpenGL calls.
  // Preconditions: must be used inbetween calls to GlSplat library
  void gl_draw_splats() const
  {
    // TODO add support for selection
    ::glBegin(GL_POINTS);
    for (const_iterator it = begin(); it != end(); it++)
    {
      const UI_point& p = *it;
      //::glColor4f(0.5,0.6,0.7,1.);
      ::glNormal3dv(&p.normal().x());
#ifdef CGAL_GLEW_ENABLED
      ::glMultiTexCoord1d(GL_TEXTURE2, p.radius());
#endif
      ::glVertex3dv(&p.x());
    }
    ::glEnd();
  }

// Private methods:
private:

  /// Recompute barycenter, bounding box, bounding sphere and standard deviation.
  void update_bounds() const
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
    m_bounding_box = Iso_cuboid(p,q);
    //
    m_barycenter = CGAL::ORIGIN + v / norm;

    // bounding sphere
    CGAL::Min_sphere_d< CGAL::Optimisation_d_traits_3<Gt> > ms3(begin(), end());
    m_bounding_sphere = Sphere(ms3.center(), ms3.squared_radius());

    // Compute standard deviation of the distance to barycenter
    typename Geom_traits::Compute_squared_distance_3 sqd;
    FT sq_radius = 0;
    for (Point_const_iterator it = begin(); it != end(); it++)
        sq_radius += sqd(*it, m_barycenter);
    sq_radius /= size();
    m_diameter_standard_deviation = CGAL::sqrt(sq_radius);

    m_bounding_box_is_valid = true;
  }

}; // end of class Point_set_3


///// Helper class: type of the "vertex_point" property map
///// of an Point_set_3 object.
//template <class Gt>
//class Point_set_vertex_point_const_map
//{
//public:
//  typedef Point_set_3<Gt> Point_set;
//  typedef typename Gt::Point_3 Point_3;
//
//  // Property maps required types
//  typedef boost::readable_property_map_tag          category;
//  typedef Point_3                                   value_type;
//  typedef value_type                                reference;
//  typedef typename Point_set::Point_const_iterator  key_type;
//
//  Point_set_vertex_point_const_map(const Point_set&) {}
//
//  /// Free function to access the map elements.
//  friend inline
//  reference get(const Point_set_vertex_point_const_map&, key_type p)
//  {
//    return *p;
//  }
//};
//
///// Free function to get the "vertex_point" property map
///// of an Point_set_3 object.
//template <class Gt>
//inline
//Point_set_vertex_point_const_map<Gt>
//get(CGAL::vertex_point_t, const Point_set_3<Gt>& points)
//{
//  Point_set_vertex_point_const_map<Gt> aMap(points);
//  return aMap;
//}
//
//
//namespace boost {
//
///// Helper type and constant to get a "vertex_normal" property map.
//enum vertex_normal_t { vertex_normal } ;
//BOOST_INSTALL_PROPERTY(vertex, normal);
//
//} // namespace boost
//
//
///// Helper class: type of the "vertex_normal" property map
///// of an Point_set_3 object.
//template <class Gt>
//class Point_set_vertex_normal_map
//  : public boost::put_get_helper<typename Gt::Vector_3&,
//                                 Point_set_vertex_normal_map<Gt> >
//{
//public:
//  typedef Point_set_3<Gt> Point_set;
//  typedef typename Gt::Vector_3 Vector_3;
//
//  // Property maps required types
//  typedef boost::lvalue_property_map_tag            category;
//  typedef Vector_3                                  value_type;
//  typedef Vector_3&                                 reference;
//  typedef typename Point_set::iterator              key_type;
//
//  Point_set_vertex_normal_map(const Point_set&) {}
//
//  /// Access the map elements.
//  reference operator[](key_type p) const { return p->normal(); }
//};
//
///// Free function to get the "vertex_normal" property map
///// of an Point_set_3 object.
//template <class Gt>
//inline
//Point_set_vertex_normal_map<Gt>
//get(boost::vertex_normal_t, const Point_set_3<Gt>& points)
//{
//  Point_set_vertex_normal_map<Gt> aMap(points);
//  return aMap;
//}


#endif // POINT_SET_3_H
