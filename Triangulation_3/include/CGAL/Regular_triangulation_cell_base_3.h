// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Christophe Delage <christophe.delage@sophia.inria.fr>

// cell of a triangulation of any dimension <=3
// with hidden points (for the regular triangulation)

#ifndef CGAL_REGULAR_TRIANGULATION_CELL_BASE_3_H
#define CGAL_REGULAR_TRIANGULATION_CELL_BASE_3_H

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/assertions.h>
#include <CGAL/Hidden_point_memory_policy.h>
#include <CGAL/Triangulation_cell_base_3.h>

#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <list>

namespace CGAL {

template < typename GT,
           typename Cb = Triangulation_cell_base_3<GT >,
           typename Memory_policy = Keep_hidden_points,
           typename C = std::list<typename GT::Weighted_point_3> >
class Regular_triangulation_cell_base_3
  : public Cb
{
public:
  typedef typename Cb::Vertex_handle                   Vertex_handle;
  typedef typename Cb::Cell_handle                     Cell_handle;

  typedef GT                                           Geom_traits;
  typedef typename Geom_traits::Point_3                Point_3;
  typedef typename Geom_traits::Weighted_point_3       Point;

  typedef C                                            Point_container;
  typedef typename Point_container::iterator           Point_iterator;
  typedef typename Point_container::const_iterator     Point_const_iterator;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other                 Cb2;
    typedef Regular_triangulation_cell_base_3<GT, Cb2, Memory_policy, C>  Other;
  };

  Regular_triangulation_cell_base_3()
    : Cb() {}

  Regular_triangulation_cell_base_3(Vertex_handle v0,
                                    Vertex_handle v1,
				    Vertex_handle v2,
                                    Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {}

  Regular_triangulation_cell_base_3(Vertex_handle v0,
                                    Vertex_handle v1,
				    Vertex_handle v2,
                                    Vertex_handle v3,
				    Cell_handle   n0,
                                    Cell_handle   n1,
				    Cell_handle   n2,
                                    Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {}

  // because we can't use default templates in ansi...
  Point_iterator hidden_points_begin()
  { return hidden_points_begin_internal<Memory_policy>(); }
  Point_iterator hidden_points_end()
  { return hidden_points_end_internal<Memory_policy>(); }

  // const versions
  Point_const_iterator hidden_points_begin() const
  { return hidden_points_begin_internal<Memory_policy>(); }
  Point_const_iterator hidden_points_end() const
  { return hidden_points_end_internal<Memory_policy>(); }

  void hide_point(const Point& p)
  { hide_point_internal<Memory_policy>(p); }
  void unhide_point(const Point_iterator pit)
  { unhide_point_internal<Memory_policy>(pit); }

  // Memory_policy is Tag_true -------------------------------------------------
  template<typename Tag>
  Point_iterator hidden_points_begin_internal(typename boost::enable_if_c<Tag::value>::type* = NULL)
  {  return _hidden.begin(); }
  template<typename Tag>
  Point_iterator hidden_points_end_internal(typename boost::enable_if_c<Tag::value>::type* = NULL)
  { return _hidden.end(); }

  template<typename Tag>
  Point_const_iterator hidden_points_begin_internal(typename boost::enable_if_c<Tag::value>::type* = NULL) const
  { return _hidden.begin(); }
  template<typename Tag>
  Point_const_iterator hidden_points_end_internal(typename boost::enable_if_c<Tag::value>::type* = NULL) const
  { return _hidden.end(); }

  template<typename Tag>
  void hide_point_internal(const Point& p, typename boost::enable_if_c<Tag::value>::type* = NULL)
  { _hidden.push_back(p); }
  template<typename Tag>
  void unhide_point_internal(const Point_iterator pit, typename boost::enable_if_c<Tag::value>::type* = NULL)
  { _hidden.erase(pit); }

  // Memory_policy is Tag_false ------------------------------------------------
  template<typename Tag>
  Point_iterator hidden_points_begin_internal(typename boost::disable_if_c<Tag::value>::type* = NULL)
  { return hidden_points_end(); }
  template<typename Tag>
  Point_iterator hidden_points_end_internal(typename boost::disable_if_c<Tag::value>::type* = NULL)
  { return _hidden.end(); }

    // const versions
  template<typename Tag>
  Point_const_iterator hidden_points_begin_internal(typename boost::disable_if_c<Tag::value>::type* = NULL) const
  { return hidden_points_end(); }
  template<typename Tag>
  Point_const_iterator hidden_points_end_internal(typename boost::disable_if_c<Tag::value>::type* = NULL) const
  { return _hidden.end(); }

  template<typename Tag>
  void hide_point_internal(const Point&, typename boost::disable_if_c<Tag::value>::type* = NULL)
  { }
  template<typename Tag>
  void unhide_point_internal(const Point_iterator, typename boost::disable_if_c<Tag::value>::type* = NULL)
  { }

  //note this function is not requested by the RegularTriangulationCellBase_3
  //it should be replaced everywhere by weighted_circumcenter()
  // but remains here for backward compatibility
  template<typename GT_>
  Point_3 circumcenter(const GT_& gt) const
  {
    CGAL_static_assertion((boost::is_same<Point_3,
      typename GT_::Construct_weighted_circumcenter_3::result_type>::value));
      return gt.construct_weighted_circumcenter_3_object()
        (this->vertex(0)->point(),
         this->vertex(1)->point(),
         this->vertex(2)->point(),
         this->vertex(3)->point());
  }

  Point_3 circumcenter() const
  {
    return circumcenter(Geom_traits());
  }

  template<typename GT_>
  Point_3 weighted_circumcenter(const GT_& gt) const
  {
    CGAL_static_assertion((boost::is_same<Point_3,
      typename GT_::Construct_weighted_circumcenter_3::result_type>::value));
      return gt.construct_weighted_circumcenter_3_object()
        (this->vertex(0)->point(),
         this->vertex(1)->point(),
         this->vertex(2)->point(),
         this->vertex(3)->point());
  }

  Point_3 weighted_circumcenter() const
  {
    return weighted_circumcenter(Geom_traits());
  }

private:
  Point_container _hidden;
};

} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_CELL_BASE_3_H
