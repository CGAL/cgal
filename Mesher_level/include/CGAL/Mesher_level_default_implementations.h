// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESHER_LEVEL_DEFAULT_IMPLEMENTATIONS_H
#define CGAL_MESHER_LEVEL_DEFAULT_IMPLEMENTATIONS_H

#include <CGAL/Mesher_level.h>

namespace CGAL {

/** This class implements the two get_triangulation_ref() functions.
    \param Tr The triangulation type */
template <class Tr>
class Triangulation_ref_impl
{
  Tr& tr;
public:
  Triangulation_ref_impl(Tr& t) : tr(t) 
  {
  }
  
  Tr& triangulation_ref_impl()
  {
    return tr;
  }
  const Tr& triangulation_ref_impl() const
  {
    return tr;
  }

}; // end class Triangulation_ref_impl<Tr>

/** This struct implements an empty private_test_point_conflict_impl()
    function. */
struct No_private_test_point_conflict 
{
  template <typename Point, typename Zone>
  Mesher_level_conflict_status
  private_test_point_conflict_impl(const Point&, const Zone&) const 
  {
    return NO_CONFLICT;
  }
}; // end No_private_test_point_conflict

/** This struct implements an empty test_point_conflict_from_superior_impl()
    function. */
struct No_test_point_conflict_from_superior
{
  template <typename Point, typename Zone>
  Mesher_level_conflict_status
  test_point_conflict_from_superior_impl(const Point&, const Zone&) const 
  {
    return NO_CONFLICT;
  }
}; // end No_test_point_conflict_from_superior

/** This struct implements empty functions:
      - private_test_point_conflict_impl() and
      - test_point_conflict_from_superior_impl().
*/
struct No_test_point_conflict : 
  public No_private_test_point_conflict,
  public No_test_point_conflict_from_superior
{
};

/** This struct implements an empty before_insertion_impl()
    function. */
struct No_before_insertion
{
  template <typename Cell_handle, typename Point, typename Zone>
  void before_insertion_impl(const Cell_handle&, const Point&,
			     Zone& )
  {
  }
}; // end No_before_insertion

/** This struct implements an empty after_insertion_impl()
    function. */
struct No_after_insertion
{
  template <typename Vertex_handle>
  void after_insertion_impl(const Vertex_handle&)
  {
  }
}; // end No_after_insertion

/** This struct implements an empty after_insertion_impl()
    function. */
struct No_after_no_insertion
{
  template <typename Cell_handle, typename Point, typename Zone>
  void after_no_insertion_impl(const Cell_handle&, const Point&,
			       const Zone& )
  {
  }
}; // end No_after_no_insertion

/** This struct implements empty functions:
      - before_insertion_impl(),
      - after_insertion_impl(),
      - after_no_insertion_impl()
*/
struct No_before_after_insertion :
  public No_after_insertion,
  public No_before_insertion,
  public No_after_no_insertion
{
};

/** This struct implements an empty before_conflicts_impl() function. */
struct No_before_conflicts {
  template <typename Face_handle, typename Point>
  void before_conflicts_impl(const Face_handle&, const Point&)
  {
  }
};

}  // end namespace CGAL

#endif // CGAL_MESHER_LEVEL_DEFAULT_IMPLEMENTATIONS_H
