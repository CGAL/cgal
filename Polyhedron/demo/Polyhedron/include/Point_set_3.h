// Author: Laurent Saboret, Nader Salman, Gael Guennebaud

#ifndef POINT_SET_3_H
#define POINT_SET_3_H

#include <CGAL/property_map.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/Surface_mesh.h>

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
class Point_set_3
{
public:
  
  typedef typename Gt::FT FT;
  typedef typename Gt::Point_3 Point;
  typedef typename Gt::Vector_3 Vector;
  typedef CGAL::cpp11::array<FT, 3> Color;
  
  typedef CGAL::Surface_mesh<Point> Base;
  typedef typename Base::Vertex_index Index;
  typedef typename Base::template Property_map<Index, Vector> Vector_pmap;
  typedef typename Base::template Property_map<Index, Color> Color_pmap;
  
private:
  Base m_base;
  Vector_pmap m_normals;
  Color_pmap m_colors;
  
public:

  Point_set_3 () : m_base() { }

  void push_back (const Point& p)
  {
    m_base.add_vertex (p);
  }

  bool empty() const { return m_base.is_empty(); }
  std::size_t size () const { return m_base.number_of_vertices(); }
  Point& operator[] (std::size_t index) { return m_base.point (Index (index)); }
  const Point& operator[] (std::size_t index) const { return (*this)[index]; }

  
  bool has_normals() const
  {
    std::pair<Vector_pmap, bool> pm = m_base.template property_map<Index, Vector> ("normal");
    return pm.second;
  }
  bool add_normal_property()
  {
    bool out = false;
    boost::tie (m_normals, out) = m_base.template add_property_map<Index, Vector> ("normal");
    return out;
  }
  void remove_normal_property()
  {
    m_base.remove_property_map (m_normals);
  }
  Vector& normal (std::size_t index) { return m_normals[Index (index)]; }
  const Vector& normal (std::size_t index) const { return this->normal(index); }

  bool has_colors() const
  {
    std::pair<Color_pmap, bool> pm = m_base.template property_map<Index, Color> ("color");
    return pm.second;
  }
  bool add_color_property()
  {
    bool out = false;
    boost::tie (m_colors, out) = m_base.template add_property_map<Index, Color> ("color");
    return out;
  }
  void remove_color_property()
  {
    m_base.remove_property_map (m_colors);
  }
  Color& color (std::size_t index) { return m_colors[Index (index)]; }
  const Color& color (std::size_t index) const { return this->color(index); }

  

private:


  
}; // end of class Point_set_3


#endif // POINT_SET_3_H
