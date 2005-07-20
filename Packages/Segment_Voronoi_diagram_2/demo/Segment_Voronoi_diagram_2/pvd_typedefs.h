// Copyright (c) 2004,2005  INRIA Sophia-Antipolis (France) and
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

#ifndef PVD_TYPEDEFS_H
#define PVD_TYPEDEFS_H

#include <CGAL/basic.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>
#include <CGAL/Segment_Voronoi_diagram_filtered_traits_2.h>

struct Rep : public CGAL::Simple_cartesian<double> {};
struct ERep : public CGAL::Simple_cartesian<CGAL::Gmpq> {};

#if 0
namespace CGAL {

  CGAL::Gmpq sqrt(const CGAL::Gmpq& x)
  {
    return CGAL::Gmpq(  sqrt( CGAL::to_double(x) )  );
  }

}
#endif

typedef CGAL::Sqrt_field_tag  MTag;
typedef CGAL::Ring_tag        EMTag;

typedef CGAL::Tag_false      ITag;
typedef CGAL::Tag_true       STag;


typedef
CGAL::Segment_Voronoi_diagram_filtered_traits_without_intersections_2<Rep,
								      MTag,
								      ERep,
								      EMTag>
Gt;

#include <CGAL/Segment_Voronoi_diagram_vertex_base_with_info_2.h>

typedef Gt::Point_2            Point_2;
typedef Gt::Segment_2          Segment;
typedef CGAL::Polygon_2<Rep>   Polygon_2;
typedef Gt::Site_2             Site;

typedef CGAL::Segment_Voronoi_diagram_vertex_base_2<Gt,ITag>          Vb;
typedef CGAL::Segment_Voronoi_diagram_vertex_base_with_info_2<Vb,int> Vbi;
typedef CGAL::Segment_Voronoi_diagram_hierarchy_vertex_base_2<Vbi>    Vbh;
typedef CGAL::Triangulation_face_base_2<Gt>                           Fb;
typedef CGAL::Triangulation_data_structure_2<Vbh,Fb>                  DS;




typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS>   SVD_2;
//typedef CGAL::Segment_Voronoi_diagram_2<Gt,DS>          SVD_2;

#endif  // PVD_TYPEDEFS_H
