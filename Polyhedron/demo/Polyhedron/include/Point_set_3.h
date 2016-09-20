// Copyright (c) 2007-2016  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Nader Salman, Gael Guennebaud

#ifndef POINT_SET_3_H
#define POINT_SET_3_H

#include <CGAL/property_map.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/Point_set_3.h>

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
class Point_set_3 : public CGAL::Point_set_3<Gt>
{
private:

  // Base class
  typedef CGAL::Point_set_3<Gt> Base;
  
public:

  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator const_iterator;
  typedef typename Base::Item Item;
  
  // Classic CGAL geometric types
  typedef Gt  Geom_traits; ///< Geometric traits class.
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;  ///< typedef to Geom_traits::Point_3
  typedef typename Geom_traits::Vector_3 Vector; ///< typedef to Geom_traits::Vector_3
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;
  typedef typename Geom_traits::Sphere_3 Sphere;

  typedef typename Base::template Property_map<double>::type Double_prop;
  typedef typename Base::template Property_map<unsigned char>::type Byte_prop;
  
private:
  
  // Indicate if m_barycenter, m_bounding_box, m_bounding_sphere and
  // m_diameter_standard_deviation below are valid.
  mutable bool m_bounding_box_is_valid;
  mutable Iso_cuboid m_bounding_box; // point set's bounding box
  mutable Sphere m_bounding_sphere; // point set's bounding sphere
  mutable Point m_barycenter; // point set's barycenter
  mutable FT m_diameter_standard_deviation; // point set's standard deviation

  bool m_radii_are_uptodate;

  Double_prop m_radius;
  Byte_prop m_red;
  Byte_prop m_green;
  Byte_prop m_blue;

  // Assignment operator not implemented and declared private to make
  // sure nobody uses the default one without knowing it
  Point_set_3& operator= (const Point_set_3&)
  {
    return *this;
  }

  
public:
  Point_set_3 ()
  {
    m_bounding_box_is_valid = false;
    m_radii_are_uptodate = false;
  }

  // copy constructor 
  Point_set_3 (const Point_set_3& p) : Base ()
  {
    m_bounding_box_is_valid = p.m_bounding_box_is_valid;
    m_bounding_box = p.m_bounding_box;
    m_barycenter = p.m_barycenter;
    m_diameter_standard_deviation = p.m_diameter_standard_deviation;
    m_radii_are_uptodate = p.m_radii_are_uptodate;
  }

  iterator begin() { return Base::begin(); }
  iterator end() { return Base::removed_end(); }
  const_iterator begin() const { return Base::begin(); }
  const_iterator end() const { return Base::removed_end(); }
  std::size_t size() const { return this->m_base.size(); }

  bool add_radius()
  {
    bool out = false;
    boost::tie (m_radius, out) = this->template add_property<double> ("radius");
    return out;
  }
  double& radius (Item index) { return m_radius[this->m_indices[index]]; }
  const double& radius (Item index) const { return m_radius[this->m_indices[index]]; }
  double& radius (iterator it) { return m_radius[*it]; }
  const double& radius (const_iterator it) const { return m_radius[*it]; }

  void test () { }
  
  bool check_colors()
  {
    bool found = false;
    
    boost::tie (m_red, found) = this->template property<unsigned char>("red");
    if (!found)
      {
        boost::tie (m_red, found) = this->template property<unsigned char>("r");
        if (!found)
          return false;
      }

    boost::tie (m_green, found) = this->template property<unsigned char>("green");
    if (!found)
      {
        boost::tie (m_green, found) = this->template property<unsigned char>("g");
        if (!found)
          return false;
      }

    boost::tie (m_blue, found) = this->template property<unsigned char>("blue");
    if (!found)
      {
        boost::tie (m_blue, found) = this->template property<unsigned char>("b");
        if (!found)
          return false;
      }

    return true;
  }

  bool has_colors() const
  {
    return m_blue != Byte_prop();
  }
    
  const unsigned char& red (const_iterator it) { return m_red[*it]; }
  const unsigned char& green (const_iterator it) { return m_green[*it]; }
  const unsigned char& blue (const_iterator it) { return m_blue[*it]; }
  
  iterator first_selected() { return this->removed_begin(); }
  const_iterator first_selected() const { return this->removed_begin(); }
  void set_first_selected(iterator it)
  {
    this->remove_from (it);
  }

  const_iterator begin_or_selection_begin() const
  {
    return (this->m_nb_removed == 0 ? begin() : first_selected());
  }
  iterator begin_or_selection_begin()
  {
    return (this->m_nb_removed == 0 ? begin() : first_selected());
  }


  // Test if point is selected
  bool is_selected(const_iterator it) const
  {
    return static_cast<std::size_t>(std::distance (it, end())) <= this->m_nb_removed;
  }

  /// Gets the number of selected points.
  std::size_t nb_selected_points() const
  {
    return this->m_nb_removed;
  }

  /// Mark a point as selected/not selected.
  void select(iterator it, bool selected = true)
  {
    bool currently = is_selected (it);
    iterator first = this->removed_begin();
    if (currently && !selected)
      {
        std::swap (*it, *first);
        -- this->m_nb_removed;
      }
    else if (!currently && selected)
      {
        std::swap (*it, *first);
        ++ this->m_nb_removed;
      }
  }

  void select_all()
  {
    this->m_nb_removed = size ();
  }
  void unselect_all()
  {
    this->m_nb_removed = 0;
  }

  // Invert selection
  void invert_selection()
  {
    iterator sel = end() - 1;
    iterator unsel = begin();

    iterator first = this->removed_begin();

    while (sel != first - 1 && unsel != first)
      std::swap (*(sel --), *(unsel ++));
    
    this->m_nb_removed = size() - this->m_nb_removed;
  }

  /// Deletes selected points.
  void delete_selection()
  {
    this->collect_garbage();
    invalidate_bounds();
  }

  void merge_with (const Point_set_3& other)
  {
    if (!(other.has_normals()) && other.has_normals())
      this->add_normal_property();

    this->m_base.transfer (other.m_base);

    unselect_all();

    // Reset indices
    for (std::size_t i = 0; i < this->m_base.size(); ++ i)
      this->m_indices[i] = i;
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
      const Point& p = this->point(*it);
      
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

    Min_sphere ms(this->m_points.begin(), this->m_points.end());

    typename Min_sphere::Cartesian_const_iterator coord = ms.center_cartesian_begin();
    FT cx = *coord++;
    FT cy = *coord++;
    FT cz = *coord++;
    m_bounding_sphere = Sphere(Point(cx,cy,cz), ms.radius()*ms.radius());

    // Computes standard deviation of the distance to barycenter
    typename Gt::Compute_squared_distance_3 sqd;
    FT sq_radius = 0;
    for (const_iterator it = begin(); it != end(); it++)
      sq_radius += sqd(this->point(*it), m_barycenter);
    sq_radius /= FT(size());
    m_diameter_standard_deviation = CGAL::sqrt(sq_radius);

    m_bounding_box_is_valid = true;
  }
  
}; // end of class Point_set_3


#endif // POINT_SET_3_H
