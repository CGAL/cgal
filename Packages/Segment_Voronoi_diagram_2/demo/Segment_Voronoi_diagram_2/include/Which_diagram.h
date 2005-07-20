// Copyright (c) 2005  INRIA Sophia-Antipolis (France) and
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

#ifndef WHICH_DIAGRAM_H
#define WHICH_DIAGRAM_H

#include <CGAL/basic.h>
#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>

CGAL_BEGIN_NAMESPACE

template<class Matching_class> struct Which_diagram;

template<class Gt, class DS, class LTag>
struct Which_diagram< Segment_Voronoi_diagram_2<Gt,DS,LTag> >
{
  typedef Tag_false Is_hierarchy;
};

template<class Gt, class STag, class DS, class LTag>
struct Which_diagram< Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> >
{
  typedef Tag_true  Is_hierarchy;
};

CGAL_END_NAMESPACE


#endif // WHICH_DIAGRAM_H
