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
#include <CGAL/Apollonius_graph_face_base_2.h>

CGAL_BEGIN_NAMESPACE


template <class Gt,
	  class Fb = Triangulation_ds_face_base_2<> >
class Segment_Voronoi_diagram_face_base_2
  : public Apollonius_graph_face_base_2<Gt,Fb>
{};


CGAL_END_NAMESPACE 

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_FACE_BASE_2_H
