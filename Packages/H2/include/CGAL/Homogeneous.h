// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : Homogeneous.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_HOMOGENEOUS_H
#define CGAL_HOMOGENEOUS_H

#define CGAL_REP_CLASS_DEFINED

#include <CGAL/basic.h>
#include <CGAL/Quotient.h>
#include <CGAL/user_classes.h>
#include <CGAL/basic_classes.h>

#include <CGAL/Homogeneous/Aff_transformationH2.h>
#include <CGAL/Homogeneous/CircleH2.h>
#include <CGAL/Homogeneous/DirectionH2.h>
#include <CGAL/Homogeneous/Iso_rectangleH2.h>
#include <CGAL/Homogeneous/LineH2.h>
#include <CGAL/Homogeneous/PointH2.h>
#include <CGAL/Homogeneous/RayH2.h>
#include <CGAL/Homogeneous/SegmentH2.h>
#include <CGAL/Homogeneous/TriangleH2.h>
#include <CGAL/Homogeneous/VectorH2.h>
#include <CGAL/Homogeneous/Data_accessorH2.h>

#include <CGAL/Homogeneous/Aff_transformationH3.h>
#include <CGAL/Homogeneous/DirectionH3.h>
#include <CGAL/Homogeneous/Iso_cuboidH3.h>
#include <CGAL/Homogeneous/LineH3.h>
#include <CGAL/Homogeneous/PlaneH3.h>
#include <CGAL/Homogeneous/PointH3.h>
#include <CGAL/Homogeneous/RayH3.h>
#include <CGAL/Homogeneous/SegmentH3.h>
#include <CGAL/Homogeneous/SphereH3.h>
#include <CGAL/Homogeneous/TetrahedronH3.h>
#include <CGAL/Homogeneous/TriangleH3.h>
#include <CGAL/Homogeneous/VectorH3.h>

#include <CGAL/Homogeneous/basic_constructionsH2.h>
#include <CGAL/Homogeneous/distance_predicatesH2.h>
#include <CGAL/Homogeneous/predicates_on_directionsH2.h>
#include <CGAL/Homogeneous/predicates_on_linesH2.h>
#include <CGAL/Homogeneous/predicates_on_pointsH2.h>
#include <CGAL/Homogeneous/predicates_on_segmentsH2.h>
#include <CGAL/Homogeneous/predicates_on_rtH2.h>

#include <CGAL/Homogeneous/basic_constructionsH3.h>
#include <CGAL/Homogeneous/distance_predicatesH3.h>
#include <CGAL/Homogeneous/orientation_predicatesH3.h>
#include <CGAL/Homogeneous/predicates_on_pointsH3.h>
#include <CGAL/Homogeneous/predicates_on_pointsH2.h>

#include <CGAL/Homogeneous/homogeneous_rep.h>

#include <CGAL/iterator_traits_pointer_specs_for_homogeneous_kernel.h>

#endif // CGAL_HOMOGENEOUS_H
