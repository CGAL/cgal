// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_ITEMS_2_H
#define CGAL_STRAIGHT_SKELETON_ITEMS_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/Straight_skeleton_vertex_base_2.h>
#include <CGAL/Straight_skeleton_halfedge_base_2.h>
#include <CGAL/Straight_skeleton_face_base_2.h>

namespace CGAL {

class Straight_skeleton_items_2
{
public:

  template<class Refs, class Traits>
  struct Vertex_wrapper
  {
    typedef typename Traits::RT RT ;
    typedef typename Traits::Point_2 Point ;
    typedef Straight_skeleton_vertex_base_2 < Refs, Point, RT > Vertex; 
  };
  
  template<class Refs, class Traits> 
  struct Halfedge_wrapper
  {
    typedef typename Traits::Segment_2 Segment ;
    typedef Straight_skeleton_halfedge_base_2 < Refs, Segment > Halfedge; 
  };
  
  template<class Refs, class Traits> 
  struct Face_wrapper 
  {
    typedef Straight_skeleton_face_base_2 < Refs > Face; 
  } ;
};

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_ITEMS_2_H //
// EOF //
 
