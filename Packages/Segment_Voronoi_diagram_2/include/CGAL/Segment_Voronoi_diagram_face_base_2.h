// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Segment_Voronoi_diagram_face_base_2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H

#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Apollonius_graph_face_base_2.h>

CGAL_BEGIN_NAMESPACE


template <class Gt,
	  class Fb = Triangulation_ds_face_base_2<> >
class Segment_Voronoi_diagram_face_base_2
  : public Apollonius_graph_face_base_2<Gt,Fb>
{};


CGAL_END_NAMESPACE 

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H
