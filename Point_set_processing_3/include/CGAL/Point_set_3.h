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

#ifndef CGAL_POINT_SET_3_H
#define CGAL_POINT_SET_3_H

#include <CGAL/property_map.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/Surface_mesh/Properties.h>

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

namespace CGAL {

template <class Gt>
class Point_set_3
{
public:

  typedef Point_set_3<Gt> Point_set;
  typedef typename Gt::FT FT;
  typedef typename Gt::Point_3 Point;
  typedef typename Gt::Vector_3 Vector;
  typedef typename Gt::Iso_cuboid_3 Iso_cuboid;
  typedef typename Gt::Sphere_3 Sphere;
  

  typedef typename std::size_t Item;

  typedef typename Properties::Property_container<Item> Base;
  typedef typename Properties::Property_map<Item, std::size_t> Index_pmap;
  typedef typename Properties::Property_map<Item, Point> Point_pmap;
  typedef typename Properties::Property_map<Item, Vector> Vector_pmap;

  typedef typename Index_pmap::Array::vector_type::iterator iterator;
  typedef typename Index_pmap::Array::vector_type::const_iterator const_iterator;

protected:

  
  struct Index_back_inserter {

    typedef std::output_iterator_tag iterator_category;
    typedef std::size_t              value_type;
    typedef std::ptrdiff_t           difference_type;
    typedef void                     pointer;
    typedef void                     reference;

  private:

    Point_set& ps;
    std::size_t ind;
  
  public:
  
    Index_back_inserter(Point_set& ps, std::size_t ind=0) : ps(ps), ind(ind) {}
    Index_back_inserter& operator++() { return *this; }
    Index_back_inserter& operator++(int) { return *this; }
    Index_back_inserter& operator*() { return *this; }
    Index_back_inserter& operator= (std::size_t& ind)
    {
      if(ps.size() <= (typename Point_set::Item(ind)))
        ps.add_item();
      put(ps.indices(), Point_set::Item(ind),ind);
      ++ ind;
      return *this;
    }
                                  
  };

  template <typename Property>
  struct Property_push_pmap
  {
    typedef typename Property::value_type value_type;

    Point_set& ps;
    Property& prop;
    std::size_t ind;

    Property_push_pmap(Point_set& ps, Property& prop, std::size_t ind=0) : ps(ps), prop(prop), ind(ind) {}
    inline friend void put(Property_push_pmap& pm, std::size_t& i, typename Property::value_type& t)
    {
      if(pm.ps.size() <= (pm.ind))
        pm.ps.add_item();
      put(pm.prop, Point_set::Item(pm.ind), t);
      i = pm.ind;
      ++pm.ind;
    }
  };
  
  typedef Property_push_pmap<Point_pmap> Point_push_pmap;
  typedef Property_push_pmap<Vector_pmap> Normal_push_pmap;

  Base m_base;
  Point_pmap m_points;
  Index_pmap m_indices;
  Vector_pmap m_normals;

  // Assignment operator not implemented and declared private to make
  // sure nobody uses the default one without knowing it
  Point_set_3& operator= (const Point_set_3&)
  {
    return *this;
  }

  // Copy constructor not implemented and declared private to make
  // sure nobody uses the default one without knowing it
  Point_set_3 (const Point_set_3& p)
  {
  }

  
public:

  Point_set_3 () : m_base()
  {
    m_indices = m_base.template add<std::size_t> ("index").first;
    m_points = m_base.template add<Point> ("point").first;
  }

  void push_back (const Point& p)
  {
    add_item();
    m_indices[size()-1] = size()-1;
    m_points[size()-1] = p;
  }
  void push_back (const Point& p, const Vector& n)
  {
    push_back (p);
    assert (has_normals());
    m_normals[size()-1] = n;
  }

  Index_pmap& indices()
  {
    return m_indices;
  }

  Point_pmap& points()
  {
    return m_points;
  }

  Vector_pmap& normals()
  {
    return m_normals;
  }

  iterator begin() { return m_indices.array().begin(); }
  iterator end() { return m_indices.array().end(); }
  const_iterator begin() const { return m_indices.array().begin(); }
  const_iterator end() const { return m_indices.array().end(); }
  bool empty() const { return (m_base.size() == 0); }
  std::size_t size () const { return m_base.size(); }
  void clear()
  {
    m_base.clear();
    m_indices = m_base.template add<std::size_t> ("index").first;
    m_points = m_base.template add<Point> ("point").first;
  }
  Point& operator[] (Item index) { return m_points[m_indices[index]]; }
  const Point& operator[] (Item index) const { return m_points[m_indices[index]]; }
  Point& point (Item index) { return m_points[m_indices[index]]; }
  const Point& point (Item index) const { return m_points[m_indices[index]]; }
  std::size_t& item (Item index) { return m_indices[index]; }
  const std::size_t& item (Item index) const { return m_indices[index]; }

  void erase (iterator first, iterator beyond)
  {
    if (beyond != end())
      {
        // TODO
      }
    
    if (are_indices_up_to_date())
      {
        m_base.erase (*first, *beyond);
      }
    else
      {
        std::size_t size = (beyond - first);
        apply_indices_change();
        m_base.resize (m_base.size() - size);
      }
  }

  void add_item ()
  {
    m_base.push_back();
  }
  
  void apply_indices_change()
  {
    for (std::size_t i = 0; i < size(); ++ i)
      if (i != m_indices[i])
        {
          m_base.swap (i, m_indices[i]);
          -- i;
        }
  }
  
  void reset_indices()
  {
    std::size_t i = 0;
    for (iterator it = begin(); it != end(); ++ it, ++ i)
      *it = i;
  }

  bool are_indices_up_to_date()
  {
    std::size_t i = 0;
    for (iterator it = begin(); it != end(); ++ it, ++ i)
      if (*it != i)
        return false;
    return true;
  }
  
  Index_back_inserter index_back_inserter ()
  {
    return Index_back_inserter (*this, size());
  }
  Point_push_pmap point_push_pmap ()
  {
    return Property_push_pmap<Point_pmap> (*this, m_points, size());
  }
  Normal_push_pmap normal_push_pmap ()
  {
    return Property_push_pmap<Vector_pmap> (*this, m_normals, size());
  }

    
  bool has_normals() const
  {
    std::pair<Vector_pmap, bool> pm = m_base.template get<Vector> ("normal");
    return pm.second;
  }
  bool add_normal_property()
  {
    bool out = false;
    boost::tie (m_normals, out) = m_base.template add<Vector> ("normal");
    return out;
  }
  void remove_normal_property()
  {
    m_base.remove (m_normals);
  }
  Vector& normal (Item index) { return m_normals[m_indices[index]]; }
  const Vector& normal (Item index) const { return m_normals[m_indices[index]]; }

  template <typename T>
  bool has_property (const std::string& name) const
  {
    std::pair<typename Properties::template Property_map<Item, T>, bool>
      pm = m_base.template get<T> (name);
    return pm.second;
  }
  template <typename T>
  bool add_property (const std::string& name)
  {
    std::pair<typename Properties::template Property_map<Item, T>, bool>
      pm = m_base.template add<T> (name);
    return pm.second;
  }
  
  template <typename PMap>
  void remove_property (PMap& prop)
  {
    m_base.remove (prop);
  }
  
  template <typename T>
  bool remove_property (const std::string& name)
  {
    std::pair<typename Properties::template Property_map<Item, T>, bool>
      pm = m_base.template get<T> (name);
    if (!(pm.second))
      return false;
    remove_property (pm.first);
    return true;
  }

  template <typename T>
  T& property (typename Properties::template Property_map<Item, T>& pmap, std::size_t index)
  {
    return pmap[index];
  }
  template <typename T>
  const T& property (typename Properties::template Property_map<Item, T>& pmap, std::size_t index) const
  {
    return property (pmap, index);
  }
  
  template <typename T>
  T& property (const std::string& name, std::size_t index)
  {
    std::pair<typename Properties::template Property_map<Item, T>, bool>
      pm = m_base.template get<T> (name);
    return property (pm.first, index);
  }
  template <typename T>
  const T& property (const std::string& name, std::size_t index) const
  {
    return property (name, index);
  }
  
}; // end of class Point_set_3

} // namespace CGAL

#endif // CGAL_POINT_SET_3_H
