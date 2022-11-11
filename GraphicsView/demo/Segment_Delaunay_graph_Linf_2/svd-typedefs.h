#ifndef SVD_TYPEDEFS_H
#define SVD_TYPEDEFS_H

#include <CGAL/basic.h>

#define USE_FILTERED_TRAITS

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Segment_Delaunay_graph_Linf_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>

#ifdef CGAL_USE_CORE
#  include <CGAL/CORE_Expr.h>
#endif

#if defined(USE_FILTERED_TRAITS) || !defined(CGAL_USE_CORE)
typedef  CGAL::Simple_cartesian<double> Rep;
#else
typedef CGAL::Simple_cartesian<CORE::Expr> Rep;
#endif

#ifdef USE_FILTERED_TRAITS
#ifdef CGAL_USE_CORE
typedef CGAL::Field_with_sqrt_tag MTag;
typedef CGAL::Field_with_sqrt_tag EMTag;
typedef CGAL::Simple_cartesian<CORE::Expr> ERep;
typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<Rep,
                                                          MTag,
                                                          ERep,
                                                          EMTag>
Gt;
#else
typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<Rep> Gt;
#endif
#else
typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Rep,CGAL::Field_tag> Gt;
#endif

typedef Gt::Point_2            Point_2;
typedef Gt::Segment_2          Segment_2;
//typedef Gt::Direction_2        Direction_2;

typedef CGAL::Polygon_2<Rep>   Polygon_2;
typedef Gt::Site_2             Site;

typedef CGAL::Tag_true         STag;

typedef CGAL::Segment_Delaunay_graph_storage_traits_2<Gt>       ST;

typedef CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag>    SDG_2;
//typedef CGAL::Segment_Delaunay_graph_2<Gt>          SDG_2;

#endif  // SVD_TYPEDEFS_H
