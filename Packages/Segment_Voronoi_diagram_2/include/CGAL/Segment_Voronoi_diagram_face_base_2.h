// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>




#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H

#include <CGAL/Triangulation_ds_face_base_2.h>

#ifdef USE_INPLACE_EDGE_LIST
#  include <CGAL/Apollonius_graph_face_base_2.h>
#endif

CGAL_BEGIN_NAMESPACE


#ifdef USE_INPLACE_EDGE_LIST
template <class Gt,
	  class Fb = Triangulation_ds_face_base_2<> >
class Segment_Voronoi_diagram_face_base_2
  : public Apollonius_graph_face_base_2<Gt,Fb>
{};

#else

template <class Gt,
	  class Fb = Triangulation_ds_face_base_2<> >
class Segment_Voronoi_diagram_face_base_2
  : public Fb
{
protected:
  // local types
  typedef typename Fb::Triangulation_data_structure    SVDDS;

public:
  // TYPES
  //------
  typedef Gt                            Geom_traits;
  typedef Fb                            Base;

  typedef SVDDS         Segment_Voronoi_diagram_data_structure_2;

  typedef typename SVDDS::Vertex_handle  Vertex_handle;
  typedef typename SVDDS::Face_handle    Face_handle;
  typedef typename SVDDS::Edge           Edge;


  template <typename SVDDS2>
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<SVDDS2>::Other  Vb2;
    typedef Segment_Voronoi_diagram_face_base_2<Gt,Vb2>      Other;
  }; 


public:
  // CREATION
  //---------
  Segment_Voronoi_diagram_face_base_2() : Base() {}

  Segment_Voronoi_diagram_face_base_2(Vertex_handle v0,
				      Vertex_handle v1,
				      Vertex_handle v2)
    : Base(v0,v1,v2) {}

  Segment_Voronoi_diagram_face_base_2(Vertex_handle v0,
				      Vertex_handle v1,
				      Vertex_handle v2,
				      Face_handle n0,
				      Face_handle n1,
				      Face_handle n2)
    : Base(v0,v1,v2,n0,n1,n2) {}
};

#endif

CGAL_END_NAMESPACE 

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H
