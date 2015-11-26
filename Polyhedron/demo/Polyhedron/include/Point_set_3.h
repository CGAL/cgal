// Author: Laurent Saboret, Nader Salman, Gael Guennebaud

#ifndef POINT_SET_3_H
#define POINT_SET_3_H

#include <CGAL/property_map.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <UI_point_3.h>
#include <algorithm>
#include <vector>
# include <CGAL/gl.h>

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
/// - Selecting points changes the order of the points in the
///   container. If selection is *not* empty, it becomes invalid after
///   adding, moving or removing points, user is reponsible to call
///   unselect_all() in those cases.
///
/// @heading Parameters:
/// @param Gt       Geometric traits class.

template <class Gt>
class Point_set_3 : public std::vector<UI_point_3<Gt> >
{
// Private types
private:

  // Base class
  typedef std::vector<UI_point_3<Gt> > Base;

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

// Data members
private:

  // Indicate if m_barycenter, m_bounding_box, m_bounding_sphere and
  // m_diameter_standard_deviation below are valid.
  mutable bool m_bounding_box_is_valid;

  mutable Iso_cuboid m_bounding_box; // point set's bounding box
  mutable Sphere m_bounding_sphere; // point set's bounding sphere
  mutable Point m_barycenter; // point set's barycenter
  mutable FT m_diameter_standard_deviation; // point set's standard deviation

  std::size_t m_nb_selected; // handle selection

  bool m_radii_are_uptodate;

  // Assignment operator not implemented and declared private to make
  // sure nobody uses the default one without knowing it
  Point_set_3& operator= (const Point_set_3&)
  {
    return *this;
  }

  
// Public methods
public:

  /// Default constructor.
  Point_set_3()
  {
    m_nb_selected = 0;
    m_bounding_box_is_valid = false;
    m_radii_are_uptodate = false;
  }

  // copy constructor 
  Point_set_3 (const Point_set_3& p) : Base (p)
  {
    m_bounding_box_is_valid = p.m_bounding_box_is_valid;
    m_bounding_box = p.m_bounding_box;
    m_barycenter = p.m_barycenter;
    m_diameter_standard_deviation = p.m_diameter_standard_deviation;

    m_nb_selected = p.nb_selected_points ();
    
    m_radii_are_uptodate = p.m_radii_are_uptodate;
  }

  // Repeat base class' public methods used below
  /// @cond SKIP_IN_MANUAL
  using Base::begin;
  using Base::end;
  using Base::size;
  /// @endcond

  iterator first_selected() { return end() - m_nb_selected; }
  const_iterator first_selected() const { return end () - m_nb_selected; }
  void set_first_selected(iterator it)
  {
    m_nb_selected = static_cast<std::size_t>(std::distance (it, end()));
  }

  // Test if point is selected
  bool is_selected(const_iterator it) const
  {
    return static_cast<std::size_t>(std::distance (it, end())) <= m_nb_selected;
  }

  /// Gets the number of selected points.
  std::size_t nb_selected_points() const
  {
    return m_nb_selected;
  }

  /// Mark a point as selected/not selected.
  void select(iterator it, bool selected = true)
  {
    bool currently = is_selected (it);
    iterator first = first_selected();
    if (currently && !selected)
      {
        std::swap (*it, *first);
        -- m_nb_selected;
      }
    else if (!currently && selected)
      {
        std::swap (*it, *first);
        ++ m_nb_selected;
      }
  }

  void select_all()
  {
    m_nb_selected = size ();
  }
  void unselect_all()
  {
    m_nb_selected = 0;
  }


  // Invert selection
  void invert_selection()
  {
    iterator sel = end() - 1;
    iterator unsel = begin();

    iterator first = first_selected();

    while (sel != first - 1 && unsel != first)
      std::swap (*(sel --), *(unsel ++));
    
    m_nb_selected = size() - m_nb_selected;
  }

  /// Deletes selected points.
  void delete_selection()
  {
    // Deletes selected points using erase-remove idiom
    erase (first_selected(), end ());

    // after erase(), use Scott Meyer's "swap trick" to trim excess capacity
    Point_set_3(*this).swap(*this);
    m_nb_selected = 0;
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
    for (const_iterator it = begin(); it != end(); it++)
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
    typedef CGAL::Min_sphere_of_points_d_traits_3<Gt,FT> Traits;
    typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;

    Min_sphere ms(begin(),end());

    typename Min_sphere::Cartesian_const_iterator coord = ms.center_cartesian_begin();
    FT cx = *coord++;
    FT cy = *coord++;
    FT cz = *coord++;
    m_bounding_sphere = Sphere(Point(cx,cy,cz), ms.radius()*ms.radius());

    // Computes standard deviation of the distance to barycenter
    typename Geom_traits::Compute_squared_distance_3 sqd;
    FT sq_radius = 0;
    for (const_iterator it = begin(); it != end(); it++)
        sq_radius += sqd(*it, m_barycenter);
    sq_radius /= FT(size());
    m_diameter_standard_deviation = CGAL::sqrt(sq_radius);

    m_bounding_box_is_valid = true;
  }

}; // end of class Point_set_3


#endif // POINT_SET_3_H
