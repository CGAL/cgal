// ======================================================================
//
// Copyright (c) 1999,2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : Kernel_archetype.h
// package       : Kernel_23
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Matthias Baesken
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_KERNEL_ARCHETYPE
#define CGAL_KERNEL_ARCHETYPE

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/functional_base.h>
#include <CGAL/Quotient.h>
#include <CGAL/concept_archetype_2.h>
#include <CGAL/concept_archetype_3.h>

#include <CGAL/Kernel/concept_archetype_functors.h>


CGAL_BEGIN_NAMESPACE

class Kernel_archetype {
public:

  // 2d

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_FT)
  typedef double                     FT;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_RT)
  typedef double                     RT;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_POINT_2)
  typedef Test_point_2               Point_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_VECTOR_2)
  typedef Test_vector_2              Vector_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_DIRECTION_2)
  typedef Test_direction_2           Direction_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LINE_2)
  typedef Test_line_2                Line_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_RAY_2)
  typedef Test_ray_2                 Ray_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_SEGMENT_2)
  typedef Test_segment_2             Segment_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_TRIANGLE_2)
  typedef Test_triangle_2            Triangle_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ISO_RECTANGLE_2)
  typedef Test_iso_rectangle_2       Iso_rectangle_2;
#endif
   
#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CIRCLE_2)
  typedef Test_circle_2              Circle_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_OBJECT_2)
  typedef CGAL::Object               Object_2;
#endif  

  // 3d
#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_POINT_3)
  typedef Test_point_3               Point_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_VECTOR_3)
  typedef Test_vector_3              Vector_3;  
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_DIRECTION_3)
  typedef Test_direction_3           Direction_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ISO_CUBOID_3)
  typedef Test_iso_cuboid_3          Iso_cuboid_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LINE_3)
  typedef Test_line_3                Line_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_RAY_3)
  typedef Test_ray_3                 Ray_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_SEGMENT_3)
  typedef Test_segment_3             Segment_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_SPHERE_3)
  typedef Test_sphere_3              Sphere_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_PLANE_3)
  typedef Test_plane_3               Plane_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_TRIANGLE_3)
  typedef Test_triangle_3            Triangle_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_TETRAHEDRON_3)
  typedef Test_tetrahedron_3         Tetrahedron_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_OBJECT_3)
  typedef CGAL::Object               Object_3;  
#endif
  
  // functors and access functions ...
  
// predicate ...  
#define CGAL_Kernel_pred(Y,Z) typedef CGALca::Y<Kernel_archetype> Y; \
Y Z() const {return Y();}
// accessor function ...
#define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

#include <CGAL/Kernel/concept_archetype_interface_macros.h>  
  
};

CGAL_END_NAMESPACE

#endif
