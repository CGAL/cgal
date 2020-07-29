// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>

#ifndef SVD_TYPEDEFS_H
#define SVD_TYPEDEFS_H

#include <CGAL/basic.h>

#define USE_FILTERED_TRAITS

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

#ifdef CGAL_USE_CORE
#  include <CGAL/CORE_Expr.h>
#endif

#if defined(USE_FILTERED_TRAITS) || !defined(CGAL_USE_CORE)
typedef  CGAL::Simple_cartesian<double> Rep;
#else
typedef CGAL::Simple_cartesian<CORE::Expr> Ref;
#endif

#ifdef USE_FILTERED_TRAITS
#ifdef CGAL_USE_CORE
typedef CGAL::Field_with_sqrt_tag MTag;
typedef CGAL::Field_with_sqrt_tag EMTag;
typedef CGAL::Simple_cartesian<CORE::Expr> ERep;
struct Gt
  : public CGAL::Segment_Delaunay_graph_filtered_traits_2<Rep,
                                                          MTag,
                                                          ERep,
                                                          EMTag>
{};
#else
struct Gt
  : public CGAL::Segment_Delaunay_graph_filtered_traits_2<Rep> {};
#endif
#else
struct Gt
  : public CGAL::Segment_Delaunay_graph_traits_2<Rep,CGAL::Field_tag> {};
#endif

typedef Gt::Point_2            Point_2;
typedef Gt::Segment_2          Segment_2;

typedef CGAL::Polygon_2<Rep>   Polygon_2;
typedef Gt::Site_2             Site;

typedef CGAL::Tag_true         STag;

typedef CGAL::Segment_Delaunay_graph_storage_traits_2<Gt>       ST;

typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,ST,STag>    SDG_2;
//typedef CGAL::Segment_Delaunay_graph_2<Gt>          SDG_2;

#endif  // SVD_TYPEDEFS_H
