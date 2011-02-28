// Author: Laurent Saboret, Nader Salman, Gael Guennebaud

#ifndef POINT_SET_3_H
#define POINT_SET_3_H

#include <CGAL/property_map.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>

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

// Public types
public:

  // Repeat base class' types
  /// @cond SKIP_IN_MANUAL
  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator const_iterator;

  using Base::erase;

  /// @endcond

  // Classic CGAL geometric types
  typedef Gt  Geom_traits; ///< Geometric traits class.
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;  ///< typedef to Geom_traits::Point_3
  typedef typename Geom_traits::Vector_3 Vector; ///< typedef to Geom_traits::Vector_3
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;
  typedef typename Geom_traits::Sphere_3 Sphere;

  /// Type of points in Point_set_3
  typedef UI_point_3<Gt> UI_point; ///< Position + normal + selection flag
  // Its superclass:
  typedef typename UI_point::Point_with_normal Point_with_normal; ///< Position + normal

  // Iterator over Point_3 points
  typedef typename std::deque<UI_point>::iterator        Point_iterator;
  typedef typename std::deque<UI_point>::const_iterator  Point_const_iterator;

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

  bool m_radii_are_uptodate;

// Public methods
public:

  /// Default constructor.
  Point_set_3()
  {
    m_nb_selected_points = 0;
    m_bounding_box_is_valid = false;
    m_radii_are_uptodate = false;
  }

  // Default copy constructor and operator =() are fine.

  // Repeat base class' public methods used below
  /// @cond SKIP_IN_MANUAL
  Base::begin;
  Base::end;
  Base::size;
  /// @endcond

  /// Gets the number of selected points.
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

  /// Deletes selected points.
  void delete_selection()
  {
    // Deletes selected points using erase-remove idiom
    erase(std::remove_if(begin(), end(), std::mem_fun_ref(&UI_point::is_selected)),
          end());

    // after erase(), use Scott Meyer's "swap trick" to trim excess capacity
    Point_set_3(*this).swap(*this);

    m_nb_selected_points = 0;
    invalidate_bounds();
  }

  /// Gets the bounding box.
  Iso_cuboid bounding_box() const
  {
    if (!m_bounding_box_is_valid)
      update_bounds();

    return m_bounding_box;
  }

  /// Gets bounding sphere.
  Sphere bounding_sphere() const
  {
    if (!m_bounding_box_is_valid)
      update_bounds();

    return m_bounding_sphere;
  }

  /// Gets points barycenter.
  Point barycenter() const
  {
    if (!m_bounding_box_is_valid)
      update_bounds();

    return m_barycenter;
  }

  /// Gets the standard deviation of the distance to barycenter.
  FT diameter_standard_deviation() const
  {
    if (!m_bounding_box_is_valid)
      update_bounds();

    return m_diameter_standard_deviation;
  }

  // Gets the region of interest, ignoring the outliers.
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
      // Draw normals
      ::glBegin(GL_LINES);
      for (const_iterator it = begin(); it != end(); it++)
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
      ::glNormal3dv(&p.normal().x());
#ifdef CGAL_GLEW_ENABLED
      ::glMultiTexCoord1d(GL_TEXTURE2, p.radius());
#endif
      ::glVertex3dv(&p.x());
    }
    ::glEnd();
  }


  bool are_radii_uptodate() const { return m_radii_are_uptodate; }
  void set_radii_uptodate(bool /*on*/) { m_radii_are_uptodate = false; }

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

    // Computes bounding sphere
    typedef CGAL::Min_sphere_of_spheres_d_traits_3<Gt,FT> Traits;
    typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
    typedef typename Traits::Sphere Traits_sphere;
    //
    // Represents points by a set of spheres with 0 radius
    std::vector<Traits_sphere> spheres;
    for (Point_const_iterator it = begin(); it != end(); it++)
      spheres.push_back(Traits_sphere(*it,0));
    //
    // Computes min sphere
    Min_sphere ms(spheres.begin(),spheres.end());
    typename Min_sphere::Cartesian_const_iterator coord = ms.center_cartesian_begin();
    FT cx = *coord++;
    FT cy = *coord++;
    FT cz = *coord++;
    m_bounding_sphere = Sphere(Point(cx,cy,cz), ms.radius()*ms.radius());

    // Computes standard deviation of the distance to barycenter
    typename Geom_traits::Compute_squared_distance_3 sqd;
    FT sq_radius = 0;
    for (Point_const_iterator it = begin(); it != end(); it++)
        sq_radius += sqd(*it, m_barycenter);
    sq_radius /= size();
    m_diameter_standard_deviation = CGAL::sqrt(sq_radius);

    m_bounding_box_is_valid = true;
  }

}; // end of class Point_set_3


#endif // POINT_SET_3_H
