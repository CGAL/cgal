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
// file          : concept_archetype_2.h
// package       : Kernel_23
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Matthias Baesken
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_CONCEPT_ARCHETYPE_2
#define CGAL_CONCEPT_ARCHETYPE_2

// 2d CA types (and number types)...

#define CGAL_concept_archetype_constructors(T) \
template<class T1> T(const T1&) { } \
template<class T1,class T2> T(const T1&,const T2&) { } \
template<class T1,class T2,class T3> T(const T1&,const T2&,const T3&) { } \
template<class T1,class T2,class T3,class T4> \
T(const T1&,const T2&,const T3&,const T4&) { } \
template<class T1,class T2,class T3,class T4,class T5> \
T(const T1&,const T2&,const T3&,const T4&,const T5&) { }


CGAL_BEGIN_NAMESPACE

/* this was replaced by a C++ built-in NT ...
struct Test_ft {
  Test_ft() {  }
  Test_ft(const Test_ft& t) { }
  
  Test_ft& operator=(const Test_ft& t) { return *this; }
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_ft& obj1, const Test_ft& obj2)
{ return true; }

inline bool operator!=(const Test_ft& obj1, const Test_ft& obj2)
{ return true; }
#endif

struct Test_rt {
  Test_rt() {  }
  Test_rt(const Test_rt& t) { }
  
  Test_rt& operator=(const Test_rt& t) { return *this; }
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_rt& obj1, const Test_rt& obj2)
{ return true; }

inline bool operator!=(const Test_rt& obj1, const Test_rt& obj2)
{ return true; }
#endif
*/


struct Test_point_2 {
  Test_point_2() {  }
  Test_point_2(const Test_point_2& t) { }
  
  Test_point_2& operator=(const Test_point_2& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_point_2)
#endif  
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_point_2& obj1, 
                       const Test_point_2& obj2)
{ return true; }

inline bool operator!=(const Test_point_2& obj1, 
                       const Test_point_2& obj2)
{ return true; }
#endif

struct Test_segment_2 {
  Test_segment_2() {  }
  Test_segment_2(const Test_segment_2& t) { }
  
  Test_segment_2& operator=(const Test_segment_2& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_segment_2)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_segment_2& obj1, const Test_segment_2& obj2)
{ return true; }

inline bool operator!=(const Test_segment_2& obj1, const Test_segment_2& obj2)
{ return true; }
#endif

struct Test_line_2 {
  Test_line_2() {  }
  Test_line_2(const Test_line_2& t) { }
  
  Test_line_2& operator=(const Test_line_2& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_line_2)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_line_2& obj1, const Test_line_2& obj2)
{ return true; }

inline bool operator!=(const Test_line_2& obj1, const Test_line_2& obj2)
{ return true; }
#endif

struct Test_ray_2 {
  Test_ray_2() {  }
  Test_ray_2(const Test_ray_2& t) { }
  
  Test_ray_2& operator=(const Test_ray_2& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_ray_2)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_ray_2& obj1, const Test_ray_2& obj2)
{ return true; }

inline bool operator!=(const Test_ray_2& obj1, const Test_ray_2& obj2)
{ return true; }
#endif

struct Test_vector_2 {
  Test_vector_2() {  }
  Test_vector_2(const Test_vector_2& t) { }
  
  Test_vector_2& operator=(const Test_vector_2& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_vector_2)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_vector_2& obj1, const Test_vector_2& obj2)
{ return true; }

inline bool operator!=(const Test_vector_2& obj1, const Test_vector_2& obj2)
{ return true; }
#endif

struct Test_direction_2 {
  Test_direction_2() {  }
  Test_direction_2(const Test_direction_2& t) { }
  
  Test_direction_2& operator=(const Test_direction_2& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_direction_2)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_direction_2& obj1, 
                       const Test_direction_2& obj2)
{ return true; }

inline bool operator!=(const Test_direction_2& obj1, 
                       const Test_direction_2& obj2)
{ return true; }
#endif

struct Test_triangle_2 {
  Test_triangle_2() {  }
  Test_triangle_2(const Test_triangle_2& t) { }
  
  Test_triangle_2& operator=(const Test_triangle_2& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_triangle_2)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_triangle_2& obj1, 
                       const Test_triangle_2& obj2)
{ return true; }

inline bool operator!=(const Test_triangle_2& obj1, 
                       const Test_triangle_2& obj2)
{ return true; }
#endif

struct Test_circle_2 {
  Test_circle_2() {  }
  Test_circle_2(const Test_circle_2& t) { }
  
  Test_circle_2& operator=(const Test_circle_2& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_circle_2)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_circle_2& obj1, const Test_circle_2& obj2)
{ return true; }

inline bool operator!=(const Test_circle_2& obj1, const Test_circle_2& obj2)
{ return true; }
#endif

struct Test_iso_rectangle_2 {
  Test_iso_rectangle_2() {  }
  Test_iso_rectangle_2(const Test_iso_rectangle_2& t) { }
  
  Test_iso_rectangle_2& operator=(const Test_iso_rectangle_2& t) 
  { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Test_iso_rectangle_2)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Test_iso_rectangle_2& obj1, 
                       const Test_iso_rectangle_2& obj2)
{ return true; }

inline bool operator!=(const Test_iso_rectangle_2& obj1, 
                       const Test_iso_rectangle_2& obj2)
{ return true; }
#endif

CGAL_END_NAMESPACE

#undef CGAL_concept_archetype_constructors

#endif
