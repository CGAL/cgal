// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-11 $
// release_date  : $CGAL_Date: 2002/08/04 $
//
// file          : include/CGAL/Polygons_bops_2.h
// package       : Map_overlay (1.12)
// maintainer    : Efi Fogel <efif@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra          <estere@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_POLYGON_BOPS_2_H
#define CGAL_POLYGON_BOPS_2_H

#ifndef CGAL_PM_WALK_ALONG_LINE_POINT_LOCATION_H
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif

#ifndef PLANAR_MAP_2
#include <CGAL/Planar_map_2.h>
#endif

#ifndef CGAL_BOP_DEFAULT_DCEL_H
#include <CGAL/Bop_default_dcel.h>
#endif

#ifndef CGAL_MAP_OVERLAY_DEFAULT_NOTIFIER_H
#include <CGAL/Map_overlay_default_notifier.h>
#endif

#ifndef CGAL_MAP_OVERLAY_H
#include <CGAL/Map_overlay.h>
#endif

#ifndef CGAL_BOOLEAN_OPERATIONS_2_H
#include <CGAL/Boolean_operations_2.h>
#endif

#ifndef CGAL_SWEEP_TO_CONSTRUCT_PLANAR_MAP_2_H
#include <CGAL/sweep_to_construct_planar_map_2.h>
#endif

#ifndef CGAL_POLYGON_2_H
#include <CGAL/Polygon_2.h>
#endif

#ifndef CGAL_HOLES_SPLIT_DCEL_H
#include <CGAL/Bops/Holes_split_dcel.h>
#endif

#ifndef CGAL_HOLES_SPLIT_NOTIFIER_H
#include <CGAL/Bops/Holes_split_notifier.h>
#endif

#ifndef CGAL_HOLES_SPLIT_H
#include <CGAL/Bops/Holes_split.h>
#endif

#ifndef CGAL_POLYGONS_FROM_FACES_H
#include <CGAL/Bops/Polygons_from_faces.h>
#endif

#ifndef CGAL_POLYGONS_DO_INTERSECT_2_H
#include <CGAL/Bops/Polygons_do_intersect_2.h>
#endif
#ifndef CGAL_POLYGONS_INTERSECTION_2_H
#include <CGAL/Bops/Polygons_intersection_2.h>
#endif
#ifndef CGAL_POLYGONS_UNION_2_H
#include <CGAL/Bops/Polygons_union_2.h>
#endif
#ifndef CGAL_POLYGONS_DIFFERENCE_2_H
#include <CGAL/Bops/Polygons_difference_2.h>
#endif
#ifndef CGAL_POLYGONS_SYMMETRIC_DIFFERENCE_2_H
#include <CGAL/Bops/Polygons_symmetric_difference_2.h>
#endif

#endif // CGAL_POLYGONS_BOPS2_H










