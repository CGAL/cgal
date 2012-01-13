// Copyright (c) 2004,2005,2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef PDG_TYPEDEFS_H
#define PDG_TYPEDEFS_H

#include <CGAL/basic.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

struct Rep : public CGAL::Simple_cartesian<double> {};
#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
struct ERep : public CGAL::Simple_cartesian<CORE::Expr> {};
#else
struct ERep : public CGAL::Simple_cartesian<CGAL::Gmpq> {};
#endif

#if 0
namespace CGAL {

  CGAL::Gmpq sqrt(const CGAL::Gmpq& x)
  {
    return CGAL::Gmpq(  sqrt( CGAL::to_double(x) )  );
  }

}
#endif

typedef CGAL::Field_with_sqrt_tag  MTag;
#ifdef CGAL_USE_CORE
typedef CGAL::Field_with_sqrt_tag  EMTag;
#else
typedef CGAL::Integral_domain_without_division_tag        EMTag;
#endif

typedef CGAL::Tag_false      ITag;
typedef CGAL::Tag_true       STag;


typedef
CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<Rep,
								     MTag,
								     ERep,
								     EMTag>
Gt;

#include <CGAL/Segment_Delaunay_graph_vertex_base_with_info_2.h>

typedef Gt::Point_2            Point_2;
typedef Gt::Segment_2          Segment;
typedef CGAL::Polygon_2<Rep>   Polygon_2;
typedef Gt::Site_2             Site;

typedef CGAL::Segment_Delaunay_graph_storage_traits_2<Gt>             ST;
typedef CGAL::Segment_Delaunay_graph_vertex_base_2<ST>                Vb;
typedef CGAL::Segment_Delaunay_graph_vertex_base_with_info_2<Vb,int>  Vbi;
typedef CGAL::Segment_Delaunay_graph_hierarchy_vertex_base_2<Vbi>     Vbh;
typedef CGAL::Triangulation_face_base_2<Gt>                           Fb;
typedef CGAL::Triangulation_data_structure_2<Vbh,Fb>                  DS_;



typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,ST,STag,DS_>   SDG_2;
//typedef CGAL::Segment_Delaunay_graph_2<Gt,ST,DS_>          SDG_2;

#endif  // PDG_TYPEDEFS_H
