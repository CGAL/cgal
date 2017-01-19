// Copyright (c) 2015  Università della Svizzera italiana.
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
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_HIERARCHY_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_HIERARCHY_2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>

namespace CGAL {

template < class Gt,
	   class ST = Segment_Delaunay_graph_storage_traits_2<Gt>,
	   class STag = Tag_false,
	   class D_S = Triangulation_data_structure_2<
              Segment_Delaunay_graph_hierarchy_vertex_base_2<
		Segment_Delaunay_graph_vertex_base_2<ST> >,
              Segment_Delaunay_graph_face_base_2<Gt> >,
	   class LTag = Tag_false>
class Segment_Delaunay_graph_Linf_hierarchy_2
  : public Segment_Delaunay_graph_hierarchy_2<Gt,ST,STag,D_S,LTag,
             Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag> >
{
  protected:
  typedef Segment_Delaunay_graph_hierarchy_2<Gt,ST,STag,D_S,LTag,
             Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag> >
          Base;
  typedef Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>
          Self;

  public:
  // CONSTRUCTORS
  //-------------
  Segment_Delaunay_graph_Linf_hierarchy_2(const Gt& gt = Gt())
    : Base(gt) {}

  template<class Input_iterator>
  Segment_Delaunay_graph_Linf_hierarchy_2(Input_iterator first,
				     Input_iterator beyond,
				     const Gt& gt=Gt())
    : Base(first, beyond, gt)
  {
  }

  Segment_Delaunay_graph_Linf_hierarchy_2(const Self& other)
    : Base(other) {}

  Self& operator=(const Self& other)
  {
    if ( this != &other ) {
      Base::operator=(other);
    }
    return *this;
  }

};


} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_HIERARCHY_2_H
