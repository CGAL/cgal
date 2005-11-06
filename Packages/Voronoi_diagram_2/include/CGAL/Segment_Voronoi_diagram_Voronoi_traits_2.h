// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Segment_Voronoi_diagram_nearest_site_2.h>
#include <CGAL/Voronoi_diagram_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>
#include <CGAL/Voronoi_diagram_2/Construct_dual_points.h>
#include <CGAL/Voronoi_diagram_2/Voronoi_traits_functors.h>


CGAL_BEGIN_NAMESPACE

template<class SVD2>
struct Segment_Voronoi_diagram_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Voronoi_traits_base_2
  <SVD2,
   CGAL_VORONOI_DIAGRAM_2_INS::Site_accessor<typename SVD2::Site_2,
					     SVD2,Tag_false>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_Voronoi_point_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_nearest_site_2<SVD2> >
{
  typedef typename SVD2::Point_2                   Point_2;
  typedef typename SVD2::Site_2                    Site_2;

  typedef Tag_false                                Has_remove;
};


CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H
