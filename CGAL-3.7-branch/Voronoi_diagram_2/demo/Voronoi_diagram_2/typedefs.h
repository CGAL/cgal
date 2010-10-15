// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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

#ifndef VD_TYPEDEFS_H
#define VD_TYPEDEFS_H

#include <CGAL/basic.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/IO/Qt_widget_Apollonius_site_2.h>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>

#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Regular_triangulation_adaptation_traits_2.h>
#include <CGAL/Apollonius_graph_adaptation_traits_2.h>

struct CK : public CGAL::Simple_cartesian<double> {};
typedef CK Rep;

typedef CGAL::Filtered_kernel<CK>                              DT_GT;
struct RT_GT : public CGAL::Regular_triangulation_filtered_traits_2<DT_GT> {} ;
struct AG_GT : public CGAL::Apollonius_graph_filtered_traits_2<CK>     {};

typedef DT_GT::Point_2            Point_2;
typedef RT_GT::Weighted_point_2   Weighted_point_2;
typedef AG_GT::Site_2             Site_2;
typedef Rep::Circle_2             Circle_2;

typedef CGAL::Delaunay_triangulation_2<DT_GT>        DT2;
typedef CGAL::Regular_triangulation_2<RT_GT>         RT2;
typedef CGAL::Apollonius_graph_2<AG_GT>              AG2;

typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT2>  DT_AT2;
typedef CGAL::Regular_triangulation_adaptation_traits_2<RT2>   RT_AT2;
typedef CGAL::Apollonius_graph_adaptation_traits_2<AG2>        AG_AT2;

typedef CGAL::Voronoi_diagram_2<DT2,DT_AT2>          VD2;
typedef CGAL::Voronoi_diagram_2<RT2,RT_AT2>          PD2;
typedef CGAL::Voronoi_diagram_2<AG2,AG_AT2>          AD2;

#include "Virtual_Voronoi_diagram_2.h"

typedef CGAL::Virtual_Voronoi_diagram_2   VVD2;

#endif  // VD_TYPEDEFS_H
