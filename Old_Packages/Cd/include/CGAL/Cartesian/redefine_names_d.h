// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/redefine_names_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#define CGAL_CARTESIAN_REDEFINE_NAMES_D_H

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL

// The following defines allow to keep a uniform presentation in
// Point_d.h for both advanced and non-advanced kernels. This is a
// temporary situation. When the non-advanced kernel disappears, all
// classes ...Cd should be renamed ..._d and this will be it for this
// file. Meanwhile, this hack allows to have both versions at once.

// In Lutz' and Michael's design, PointCd<R> becomes
// Point_d<R,Cartesian_tag>, hence the following defines.

#define PointCd Point_d
#define VectorCd Vector_d
#define DirectionCd Direction_d
#define PlaneCd Plane_d
#define LineCd Line_d
#define RayCd Ray_d
#define SegmentCd Segment_d
#define TriangleCd Triangle_d
#define TetrahedronCd Tetrahedron_d
#define SimplexCd Simplex_d
#define Aff_transformationCd Aff_transformation_d
#define Data_accessorCd Data_accessor_d

// There is one more problem in unifying the two designs
// We must also take care of the declarations in .C files
// Specifically, in Lutz' and Michael's design, functions are defined as
// PointCd<R,Cartesian_tag>::function() {}
// and in Stefan's design, their definition looks like
// PointCd<R>::function() {}
// We unify this with PointCd<R CGAL_CTAG >::function
#define CGAL_CTAG   , Cartesian_tag

#define CGAL_ADVANCED_KERNEL_PARTIAL_SPEC <R_,Cartesian_tag>

#else

#define CGAL_CTAG

#define CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
 
// Note that it was not possible to keep the opposite (changing the
// Point_d into PointCd for Stefan's design) because it would also
// change "typename R::Point_d" into "typename R::PointCd" inside the
// classes definition. A successful hack would then have to analyze
// which Point_d is intended, hence too complicated.

#endif

#endif // CGAL_CARTESIAN_REDEFINE_NAMES_D_H
