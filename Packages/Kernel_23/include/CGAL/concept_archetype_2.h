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



struct Point_2_archetype {
  Point_2_archetype() {  }
  Point_2_archetype(const Point_2_archetype& t) { }
  
  Point_2_archetype& operator=(const Point_2_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Point_2_archetype)
#endif  
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Point_2_archetype& obj1, 
                       const Point_2_archetype& obj2)
{ return true; }

inline bool operator!=(const Point_2_archetype& obj1, 
                       const Point_2_archetype& obj2)
{ return true; }
#endif

struct Segment_2_archetype {
  Segment_2_archetype() {  }
  Segment_2_archetype(const Segment_2_archetype& t) { }
  
  Segment_2_archetype& operator=(const Segment_2_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Segment_2_archetype)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Segment_2_archetype& obj1, const Segment_2_archetype& obj2)
{ return true; }

inline bool operator!=(const Segment_2_archetype& obj1, const Segment_2_archetype& obj2)
{ return true; }
#endif

struct Line_2_archetype {
  Line_2_archetype() {  }
  Line_2_archetype(const Line_2_archetype& t) { }
  
  Line_2_archetype& operator=(const Line_2_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Line_2_archetype)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Line_2_archetype& obj1, const Line_2_archetype& obj2)
{ return true; }

inline bool operator!=(const Line_2_archetype& obj1, const Line_2_archetype& obj2)
{ return true; }
#endif

struct Ray_2_archetype {
  Ray_2_archetype() {  }
  Ray_2_archetype(const Ray_2_archetype& t) { }
  
  Ray_2_archetype& operator=(const Ray_2_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Ray_2_archetype)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Ray_2_archetype& obj1, const Ray_2_archetype& obj2)
{ return true; }

inline bool operator!=(const Ray_2_archetype& obj1, const Ray_2_archetype& obj2)
{ return true; }
#endif

struct Vector_2_archetype {
  Vector_2_archetype() {  }
  Vector_2_archetype(const Vector_2_archetype& t) { }
  
  Vector_2_archetype& operator=(const Vector_2_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Vector_2_archetype)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Vector_2_archetype& obj1, const Vector_2_archetype& obj2)
{ return true; }

inline bool operator!=(const Vector_2_archetype& obj1, const Vector_2_archetype& obj2)
{ return true; }
#endif

struct Direction_2_archetype {
  Direction_2_archetype() {  }
  Direction_2_archetype(const Direction_2_archetype& t) { }
  
  Direction_2_archetype& operator=(const Direction_2_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Direction_2_archetype)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Direction_2_archetype& obj1, 
                       const Direction_2_archetype& obj2)
{ return true; }

inline bool operator!=(const Direction_2_archetype& obj1, 
                       const Direction_2_archetype& obj2)
{ return true; }
#endif

struct Triangle_2_archetype {
  Triangle_2_archetype() {  }
  Triangle_2_archetype(const Triangle_2_archetype& t) { }
  
  Triangle_2_archetype& operator=(const Triangle_2_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Triangle_2_archetype)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Triangle_2_archetype& obj1, 
                       const Triangle_2_archetype& obj2)
{ return true; }

inline bool operator!=(const Triangle_2_archetype& obj1, 
                       const Triangle_2_archetype& obj2)
{ return true; }
#endif

struct Circle_2_archetype {
  Circle_2_archetype() {  }
  Circle_2_archetype(const Circle_2_archetype& t) { }
  
  Circle_2_archetype& operator=(const Circle_2_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Circle_2_archetype)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Circle_2_archetype& obj1, const Circle_2_archetype& obj2)
{ return true; }

inline bool operator!=(const Circle_2_archetype& obj1, const Circle_2_archetype& obj2)
{ return true; }
#endif

struct Iso_rectangle_2_archetype {
  Iso_rectangle_2_archetype() {  }
  Iso_rectangle_2_archetype(const Iso_rectangle_2_archetype& t) { }
  
  Iso_rectangle_2_archetype& operator=(const Iso_rectangle_2_archetype& t) 
  { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Iso_rectangle_2_archetype)
#endif   
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Iso_rectangle_2_archetype& obj1, 
                       const Iso_rectangle_2_archetype& obj2)
{ return true; }

inline bool operator!=(const Iso_rectangle_2_archetype& obj1, 
                       const Iso_rectangle_2_archetype& obj2)
{ return true; }
#endif

CGAL_END_NAMESPACE

#undef CGAL_concept_archetype_constructors

#endif
