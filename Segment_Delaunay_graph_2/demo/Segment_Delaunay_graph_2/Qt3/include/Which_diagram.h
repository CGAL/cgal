// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL$
// $Id$
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef WHICH_DIAGRAM_H
#define WHICH_DIAGRAM_H

#include <CGAL/basic.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>

namespace CGAL {

template<class Matching_class> struct Which_diagram;

template<class Gt, class DS_, class LTag>
struct Which_diagram< Segment_Delaunay_graph_2<Gt,DS_,LTag> >
{
  typedef Tag_false Is_hierarchy;
};

template<class Gt, class STag, class DS_, class LTag>
struct Which_diagram< Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS_,LTag> >
{
  typedef Tag_true  Is_hierarchy;
};

} //namespace CGAL


#endif // WHICH_DIAGRAM_H
