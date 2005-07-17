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
#include <CGAL/Delaunay_triangulation_Voronoi_traits_2.h>
#include <CGAL/Regular_triangulation_Voronoi_traits_2.h>
#include <CGAL/Apollonius_graph_Voronoi_traits_2.h>

struct CK : public CGAL::Simple_cartesian<double> {};
typedef CK Rep;

typedef CGAL::Filtered_kernel<CK>                              DT_GT;
typedef CGAL::Regular_triangulation_filtered_traits_2<DT_GT>   RT_GT;
typedef CGAL::Apollonius_graph_filtered_traits_2<CK>           AG_GT;

typedef DT_GT::Point_2            Point_2;
typedef RT_GT::Weighted_point_2   Weighted_point_2;
typedef AG_GT::Site_2             Site_2;
typedef Rep::Circle_2             Circle_2;

typedef CGAL::Delaunay_triangulation_2<DT_GT>        DT2;
typedef CGAL::Regular_triangulation_2<RT_GT>         RT2;
typedef CGAL::Apollonius_graph_2<AG_GT>              AG2;

typedef CGAL::Delaunay_triangulation_Voronoi_traits_2<DT2>  DT_VT2;
typedef CGAL::Regular_triangulation_Voronoi_traits_2<RT2>   RT_VT2;
typedef CGAL::Apollonius_graph_Voronoi_traits_2<AG2>        AG_VT2;

typedef CGAL::Voronoi_diagram_2<DT2,DT_VT2>          VD2;
typedef CGAL::Voronoi_diagram_2<RT2,RT_VT2>          PD2;
typedef CGAL::Voronoi_diagram_2<AG2,AG_VT2>          AD2;

#include "Virtual_Voronoi_diagram_2.h"

typedef CGAL::Virtual_Voronoi_diagram_2   VVD2;

#endif  // VD_TYPEDEFS_H
