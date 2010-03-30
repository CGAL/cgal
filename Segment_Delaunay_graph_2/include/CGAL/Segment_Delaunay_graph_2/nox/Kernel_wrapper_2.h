// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France) and
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
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_KERNEL_WRAPPER_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_KERNEL_WRAPPER_2_H

#include <CGAL/Segment_Delaunay_graph_2/basic.h>

#include <CGAL/Segment_Delaunay_graph_2/nox/Constructions_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/nox/Segment_Delaunay_graph_site_2.h>


CGAL_BEGIN_NAMESPACE

CGAL_SEGMENT_DELAUNAY_GRAPH_2_BEGIN_NAMESPACE


template<class Kernel_base_2, class ITag>
class Kernel_wrapper_2 : public Kernel_base_2
{
public:
  typedef Kernel_base_2    Kernel_base;
  //  typedef ITag             Intersections_tag;
  typedef CGAL::Tag_false  Intersections_tag;

  typedef CGAL::Segment_Delaunay_graph_site_2<Kernel_base> Site_2;

  typedef Construct_sdg_site_2<Site_2,ITag>  Construct_site_2;
};


CGAL_SEGMENT_DELAUNAY_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_KERNEL_WRAPPER_2_H
