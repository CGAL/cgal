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
// file          : concept_archetype_3.h
// package       : Kernel_23
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Matthias Baesken
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_CONCEPT_ARCHETYPE_3
#define CGAL_CONCEPT_ARCHETYPE_3

// 3d CA types ...

#define CGAL_concept_archetype_constructors(T) \
template<class T1> T(const T1&) { } \
template<class T1,class T2> T(const T1&,const T2&) { } \
template<class T1,class T2,class T3> T(const T1&,const T2&,const T3&) { } \
template<class T1,class T2,class T3,class T4> \
T(const T1&,const T2&,const T3&,const T4&) { } \
template<class T1,class T2,class T3,class T4,class T5> \
T(const T1&,const T2&,const T3&,const T4&,const T5&) { }


CGAL_BEGIN_NAMESPACE

struct Test_point_3 {
  Test_point_3() {  }
  Test_point_3(const Test_point_3& t) { }
  
  Test_point_3& operator=(const Test_point_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_point_3)
#endif  
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_point_3& obj1, const Test_point_3& obj2)
{ return true; }

inline bool operator!=(const Test_point_3& obj1, const Test_point_3& obj2)
{ return true; }
#endif

struct Test_segment_3 {
  Test_segment_3() {  }
  Test_segment_3(const Test_segment_3& t) { }
  
  Test_segment_3& operator=(const Test_segment_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_segment_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_segment_3& obj1, const Test_segment_3& obj2)
{ return true; }

inline bool operator!=(const Test_segment_3& obj1, const Test_segment_3& obj2)
{ return true; }
#endif

struct Test_line_3 {
  Test_line_3() {  }
  Test_line_3(const Test_line_3& t) { }
  
  Test_line_3& operator=(const Test_line_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_line_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_line_3& obj1, const Test_line_3& obj2)
{ return true; }

inline bool operator!=(const Test_line_3& obj1, const Test_line_3& obj2)
{ return true; }
#endif

struct Test_ray_3 {
  Test_ray_3() {  }
  Test_ray_3(const Test_ray_3& t) { }
  
  Test_ray_3& operator=(const Test_ray_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_ray_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_ray_3& obj1, const Test_ray_3& obj2)
{ return true; }

inline bool operator!=(const Test_ray_3& obj1, const Test_ray_3& obj2)
{ return true; }
#endif

struct Test_vector_3 {
  Test_vector_3() {  }
  Test_vector_3(const Test_vector_3& t) { }
  
  Test_vector_3& operator=(const Test_vector_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_vector_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_vector_3& obj1, const Test_vector_3& obj2)
{ return true; }

inline bool operator!=(const Test_vector_3& obj1, const Test_vector_3& obj2)
{ return true; }
#endif

struct Test_direction_3 {
  Test_direction_3() {  }
  Test_direction_3(const Test_direction_3& t) { }
  
  Test_direction_3& operator=(const Test_direction_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_direction_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_direction_3& obj1, 
                       const Test_direction_3& obj2)
{ return true; }

inline bool operator!=(const Test_direction_3& obj1, 
                       const Test_direction_3& obj2)
{ return true; }
#endif

struct Test_plane_3 {
  Test_plane_3() {  }
  Test_plane_3(const Test_plane_3& t) { }
  
  Test_plane_3& operator=(const Test_plane_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_plane_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_plane_3& obj1, const Test_plane_3& obj2)
{ return true; }

inline bool operator!=(const Test_plane_3& obj1, const Test_plane_3& obj2)
{ return true; }
#endif

struct Test_iso_cuboid_3 {
  Test_iso_cuboid_3() {  }
  Test_iso_cuboid_3(const Test_iso_cuboid_3& t) { }
  
  Test_iso_cuboid_3& operator=(const Test_iso_cuboid_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_iso_cuboid_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_iso_cuboid_3& obj1, 
                       const Test_iso_cuboid_3& obj2)
{ return true; }

inline bool operator!=(const Test_iso_cuboid_3& obj1, 
                       const Test_iso_cuboid_3& obj2)
{ return true; }
#endif

struct Test_sphere_3 {
  Test_sphere_3() {  }
  Test_sphere_3(const Test_sphere_3& t) { }
  
  Test_sphere_3& operator=(const Test_sphere_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_sphere_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_sphere_3& obj1, const Test_sphere_3& obj2)
{ return true; }

inline bool operator!=(const Test_sphere_3& obj1, const Test_sphere_3& obj2)
{ return true; }
#endif

struct Test_triangle_3 {
  Test_triangle_3() {  }
  Test_triangle_3(const Test_triangle_3& t) { }
  
  Test_triangle_3& operator=(const Test_triangle_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_triangle_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_triangle_3& obj1, const Test_triangle_3& obj2)
{ return true; }

inline bool operator!=(const Test_triangle_3& obj1, const Test_triangle_3& obj2)
{ return true; }
#endif

struct Test_tetrahedron_3 {
  Test_tetrahedron_3() {  }
  Test_tetrahedron_3(const Test_tetrahedron_3& t) { }
  
  Test_tetrahedron_3& operator=(const Test_tetrahedron_3& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_tetrahedron_3)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_tetrahedron_3& obj1, 
                       const Test_tetrahedron_3& obj2)
{ return true; }

inline bool operator!=(const Test_tetrahedron_3& obj1, 
                       const Test_tetrahedron_3& obj2)
{ return true; }
#endif

CGAL_END_NAMESPACE

#undef CGAL_concept_archetype_constructors

#endif
