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
  typedef Point_2_archetype               Point_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_VECTOR_2)
  typedef Vector_2_archetype              Vector_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_DIRECTION_2)
  typedef Direction_2_archetype           Direction_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LINE_2)
  typedef Line_2_archetype                Line_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_RAY_2)
  typedef Ray_2_archetype                 Ray_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_SEGMENT_2)
  typedef Segment_2_archetype             Segment_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_TRIANGLE_2)
  typedef Triangle_2_archetype            Triangle_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ISO_RECTANGLE_2)
  typedef Iso_rectangle_2_archetype       Iso_rectangle_2;
#endif
   
#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CIRCLE_2)
  typedef Circle_2_archetype              Circle_2;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_OBJECT_2)
  typedef CGAL::Object                    Object_2;
#endif  

  // 3d
#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_POINT_3)
  typedef Point_3_archetype               Point_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_VECTOR_3)
  typedef Vector_3_archetype              Vector_3;  
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_DIRECTION_3)
  typedef Direction_3_archetype           Direction_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ISO_CUBOID_3)
  typedef Iso_cuboid_3_archetype          Iso_cuboid_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LINE_3)
  typedef Line_3_archetype                Line_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_RAY_3)
  typedef Ray_3_archetype                 Ray_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_SEGMENT_3)
  typedef Segment_3_archetype             Segment_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_SPHERE_3)
  typedef Sphere_3_archetype              Sphere_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_PLANE_3)
  typedef Plane_3_archetype               Plane_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_TRIANGLE_3)
  typedef Triangle_3_archetype            Triangle_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_TETRAHEDRON_3)
  typedef Tetrahedron_3_archetype         Tetrahedron_3;
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_OBJECT_3)
  typedef CGAL::Object                    Object_3;  
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
