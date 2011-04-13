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
// file          : include/CGAL/Cartesian/redefine_names_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CARTESIAN_REDEFINE_NAMES_3_H

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL

// The following defines allow to keep a uniform presentation in
// Point_3.h for both advanced and non-advanced kernels. This is a
// temporary situation. When the non-advanced kernel disappears, all
// classes ...C3 should be renamed ..._3 and this will be it for this
// file. Meanwhile, this hack allows to have both versions at once.

// In Lutz' and Michael's design, PointC3<R> becomes
// Point_3<R,Cartesian_tag>, hence the following defines.

#define PointC3 Point_3
#define VectorC3 Vector_3
#define DirectionC3 Direction_3
#define LineC3 Line_3
#define PlaneC3 Plane_3
#define RayC3 Ray_3
#define SegmentC3 Segment_3
#define SphereC3 Sphere_3
#define TriangleC3 Triangle_3
#define TetrahedronC3 Tetrahedron_3
#define Iso_cuboidC3 Iso_cuboid_3
#define Aff_transformationC3 Aff_transformation_3
#define Data_accessorC3 Data_accessor_3

// There is one more problem in unifying the two designs
// We must also take care of the declarations in .C files
// Specifically, in Lutz' and Michael's design, functions are defined as
// PointC3<R,Cartesian_tag>::function() {}
// and in Stefan's design, their definition looks like
// PointC3<R>::function() {}
// We unify this with PointC3<R CGAL_CTAG >::function
#define CGAL_CTAG   , Cartesian_tag

// This is the mark of a partial specialization.  Used by all kernel classes.
#define CGAL_ADVANCED_KERNEL_PARTIAL_SPEC <R_,Cartesian_tag>

#else

#define CGAL_CTAG
#define CGAL_ADVANCED_KERNEL_PARTIAL_SPEC

// Note that it was not possible to keep the opposite (changing the
// Point_3 into PointC3 for Stefan's design) because it would also
// change "typename R::Point_3" into "typename R::PointC3" inside the
// classes definition. A successful hack would then have to analyze
// which Point_3 is intended, hence too complicated.

#endif

#endif // CGAL_CARTESIAN_REDEFINE_NAMES_3_H
