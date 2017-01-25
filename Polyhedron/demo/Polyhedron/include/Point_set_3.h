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
// Author(s)     : Laurent Saboret, Nader Salman, Gael Guennebaud, Simon Giraudot

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
class Point_set_3 : public CGAL::Point_set_3<typename Gt::Point_3,
                                             typename Gt::Vector_3>
{
  
public:

  typedef Gt  Geom_traits; ///< Geometric traits class.
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;  ///< typedef to Geom_traits::Point_3
  typedef typename Geom_traits::Vector_3 Vector; ///< typedef to Geom_traits::Vector_3
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;
  typedef typename Geom_traits::Sphere_3 Sphere;

  // Base class
  typedef CGAL::Point_set_3<Point, Vector> Base;

  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator const_iterator;
  typedef typename Base::Index Index;
  
  // Classic CGAL geometric types

  typedef typename Base::template Property_map<double> Double_map;
  typedef typename Base::template Property_map<unsigned char> Byte_map;
  
private:
  
  // Indicate if m_barycenter, m_bounding_box, m_bounding_sphere and
  // m_diameter_standard_deviation below are valid.
  mutable bool m_bounding_box_is_valid;
  mutable Iso_cuboid m_bounding_box; // point set's bounding box
  mutable Sphere m_bounding_sphere; // point set's bounding sphere
  mutable Point m_barycenter; // point set's barycenter
  mutable FT m_diameter_standard_deviation; // point set's standard deviation

  bool m_radii_are_uptodate;

  Double_map m_radius;
  Byte_map m_red;
  Byte_map m_green;
  Byte_map m_blue;
  Double_map m_fred;
  Double_map m_fgreen;
  Double_map m_fblue;

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
  Point_set_3 (const Point_set_3& p) : Base (p)
  {
    check_colors();
    m_bounding_box_is_valid = p.m_bounding_box_is_valid;
    m_bounding_box = p.m_bounding_box;
    m_barycenter = p.m_barycenter;
    m_diameter_standard_deviation = p.m_diameter_standard_deviation;
    m_radii_are_uptodate = p.m_radii_are_uptodate;
  }

  iterator begin() { return this->m_indices.begin(); }
  iterator end() { return this->m_indices.end(); }
  const_iterator begin() const { return this->m_indices.begin(); }
  const_iterator end() const { return this->m_indices.end(); }
  std::size_t size() const { return this->m_base.size(); }

  bool add_radius()
  {
    bool out = false;
    boost::tie (m_radius, out) = this->template add_property_map<double> ("radius", 0.);
    return out;
  }
  double& radius (const Index& index) { return m_radius[index]; }
  const double& radius (const Index& index) const { return m_radius[index]; }

  bool check_colors()
  {
    bool found = false;
    
    boost::tie (m_red, found) = this->template property_map<unsigned char>("red");
    if (!found)
      {
        boost::tie (m_red, found) = this->template property_map<unsigned char>("r");
        if (!found)
          return get_float_colors();
      }

    boost::tie (m_green, found) = this->template property_map<unsigned char>("green");
    if (!found)
      {
        boost::tie (m_green, found) = this->template property_map<unsigned char>("g");
        if (!found)
          return false;
      }

    boost::tie (m_blue, found) = this->template property_map<unsigned char>("blue");
    if (!found)
      {
        boost::tie (m_blue, found) = this->template property_map<unsigned char>("b");
        if (!found)
          return false;
      }

    return true;
  }

  bool get_float_colors()
  {
    bool found = false;

    boost::tie (m_fred, found) = this->template property_map<double>("red");
    if (!found)
      {
        boost::tie (m_fred, found) = this->template property_map<double>("r");
        if (!found)
          return get_las_colors();
      }

    boost::tie (m_fgreen, found) = this->template property_map<double>("green");
    if (!found)
      {
        boost::tie (m_fgreen, found) = this->template property_map<double>("g");
        if (!found)
          return false;
      }

    boost::tie (m_fblue, found) = this->template property_map<double>("blue");
    if (!found)
      {
        boost::tie (m_fblue, found) = this->template property_map<double>("b");
        if (!found)
          return false;
      }
    return true;
  }

  bool get_las_colors()
  {
    bool found = false;
    
    typedef typename Base::template Property_map<unsigned short> Ushort_map;
    Ushort_map red, green, blue;
    
    boost::tie (red, found) = this->template property_map<unsigned short>("R");
    if (!found)
      return false;

    boost::tie (green, found) = this->template property_map<unsigned short>("G");
    if (!found)
      return false;

    boost::tie (blue, found) = this->template property_map<unsigned short>("B");
    if (!found)
      return false;

    unsigned int bit_short_to_char = 0;
    for (iterator it = begin(); it != end(); ++ it)
      if (get(red, *it) > 255
          || get(green, *it) > 255
          || get(blue, *it) > 255)
        {
          bit_short_to_char = 8;
          break;
        }

    m_red = this->template add_property_map<unsigned char>("r").first;
    m_green = this->template add_property_map<unsigned char>("g").first;
    m_blue = this->template add_property_map<unsigned char>("b").first;
    for (iterator it = begin(); it != end(); ++ it)
      {
        put (m_red, *it, (unsigned char)((get(red, *it) >> bit_short_to_char)));
        put (m_green, *it, (unsigned char)((get(green, *it) >> bit_short_to_char)));
        put (m_blue, *it, (unsigned char)((get(blue, *it) >> bit_short_to_char)));
      }
    this->remove_property_map(red);
    this->remove_property_map(green);
    this->remove_property_map(blue);
    
    return true;
  }

  bool has_colors() const
  {
    return (m_blue != Byte_map() || m_fblue != Double_map());
  }

  bool has_byte_colors() const
  {
    return (m_blue != Byte_map());
  }
    
  double red (const Index& index) const
  { return (m_red == Byte_map()) ? m_fred[index]  : double(m_red[index]) / 255.; }
  double green (const Index& index) const
  { return (m_green == Byte_map()) ? m_fgreen[index]  : double(m_green[index]) / 255.; }
  double blue (const Index& index) const
  { return (m_blue == Byte_map()) ? m_fblue[index]  : double(m_blue[index]) / 255.; }

  
  iterator first_selected() { return this->m_indices.end() - this->m_nb_removed; }
  const_iterator first_selected() const { return this->m_indices.end() - this->m_nb_removed; }
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
    return this->is_removed (it);
  }

  /// Gets the number of selected points.
  std::size_t nb_selected_points() const
  {
    return this->number_of_removed_points();
  }

  /// Mark a point as selected/not selected.
  void select(iterator it, bool selected = true)
  {
    bool currently = is_selected (it);
    iterator first = this->first_selected();
    --first;
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

    iterator first = this->first_selected();

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

  void merge_with (Point_set_3& other)
  {
    if (!(this->has_normal_map()) && other.has_normal_map())
      this->add_normal_map();
    if (!(this->has_colors()) && other.has_colors())
      {
        if (other.template has_property_map<unsigned char>("red"))
          {
            boost::tie (m_red, boost::tuples::ignore)
              = this->template add_property_map<unsigned char>("red", 0);
            boost::tie (m_green, boost::tuples::ignore)
              = this->template add_property_map<unsigned char>("green", 0);
            boost::tie (m_blue, boost::tuples::ignore)
              = this->template add_property_map<unsigned char>("blue", 0);
          }
        else
          {
            boost::tie (m_red, boost::tuples::ignore)
              = this->template add_property_map<unsigned char>("r", 0);
            boost::tie (m_green, boost::tuples::ignore)
              = this->template add_property_map<unsigned char>("g", 0);
            boost::tie (m_blue, boost::tuples::ignore)
              = this->template add_property_map<unsigned char>("b", 0);
          }
      }

    this->join (other);
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
