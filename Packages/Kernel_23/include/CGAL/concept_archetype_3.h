// Copyright (c) 1999,2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Matthias Baesken

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

struct Point_3_archetype {
  Point_3_archetype() {  }
  Point_3_archetype(const Point_3_archetype& t) { }
  
  Point_3_archetype& operator=(const Point_3_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Point_3_archetype)
#endif  
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Point_3_archetype& obj1, 
                       const Point_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Point_3_archetype& obj1, 
                       const Point_3_archetype& obj2)
{ return true; }
#endif

struct Segment_3_archetype {
  Segment_3_archetype() {  }
  Segment_3_archetype(const Segment_3_archetype& t) { }
  
  Segment_3_archetype& operator=(const Segment_3_archetype& t){return *this;}
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Segment_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Segment_3_archetype& obj1, 
                       const Segment_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Segment_3_archetype& obj1, 
                       const Segment_3_archetype& obj2)
{ return true; }
#endif

struct Line_3_archetype {
  Line_3_archetype() {  }
  Line_3_archetype(const Line_3_archetype& t) { }
  
  Line_3_archetype& operator=(const Line_3_archetype& t){return *this;}
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Line_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Line_3_archetype& obj1, 
                       const Line_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Line_3_archetype& obj1, 
                       const Line_3_archetype& obj2)
{ return true; }
#endif

struct Ray_3_archetype {
  Ray_3_archetype() {  }
  Ray_3_archetype(const Ray_3_archetype& t) { }
  
  Ray_3_archetype& operator=(const Ray_3_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Ray_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Ray_3_archetype& obj1, 
                       const Ray_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Ray_3_archetype& obj1, 
                       const Ray_3_archetype& obj2)
{ return true; }
#endif

struct Vector_3_archetype {
  Vector_3_archetype() {  }
  Vector_3_archetype(const Vector_3_archetype& t) { }
  
  Vector_3_archetype& operator=(const Vector_3_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Vector_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Vector_3_archetype& obj1, 
                       const Vector_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Vector_3_archetype& obj1, 
                       const Vector_3_archetype& obj2)
{ return true; }
#endif

struct Direction_3_archetype {
  Direction_3_archetype() {  }
  Direction_3_archetype(const Direction_3_archetype& t) { }
  
  Direction_3_archetype& operator=(const Direction_3_archetype& t)
  { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Direction_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Direction_3_archetype& obj1, 
                       const Direction_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Direction_3_archetype& obj1, 
                       const Direction_3_archetype& obj2)
{ return true; }
#endif

struct Plane_3_archetype {
  Plane_3_archetype() {  }
  Plane_3_archetype(const Plane_3_archetype& t) { }
  
  Plane_3_archetype& operator=(const Plane_3_archetype& t) { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Plane_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Plane_3_archetype& obj1, 
                       const Plane_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Plane_3_archetype& obj1, 
                       const Plane_3_archetype& obj2)
{ return true; }
#endif

struct Iso_cuboid_3_archetype {
  Iso_cuboid_3_archetype() {  }
  Iso_cuboid_3_archetype(const Iso_cuboid_3_archetype& t) { }
  
  Iso_cuboid_3_archetype& operator=(const Iso_cuboid_3_archetype& t) 
  { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Iso_cuboid_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Iso_cuboid_3_archetype& obj1, 
                       const Iso_cuboid_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Iso_cuboid_3_archetype& obj1, 
                       const Iso_cuboid_3_archetype& obj2)
{ return true; }
#endif

struct Sphere_3_archetype {
  Sphere_3_archetype() {  }
  Sphere_3_archetype(const Sphere_3_archetype& t) { }
  
  Sphere_3_archetype& operator=(const Sphere_3_archetype& t) 
  { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Sphere_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Sphere_3_archetype& obj1, 
                       const Sphere_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Sphere_3_archetype& obj1, 
                       const Sphere_3_archetype& obj2)
{ return true; }
#endif

struct Triangle_3_archetype {
  Triangle_3_archetype() {  }
  Triangle_3_archetype(const Triangle_3_archetype& t) { }
  
  Triangle_3_archetype& operator=(const Triangle_3_archetype& t) 
  { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Triangle_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Triangle_3_archetype& obj1, 
                       const Triangle_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Triangle_3_archetype& obj1, 
                       const Triangle_3_archetype& obj2)
{ return true; }
#endif

struct Tetrahedron_3_archetype {
  Tetrahedron_3_archetype() {  }
  Tetrahedron_3_archetype(const Tetrahedron_3_archetype& t) { }
  
  Tetrahedron_3_archetype& operator=(const Tetrahedron_3_archetype& t) 
  { return *this; }
  
#if defined(CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS)
  CGAL_concept_archetype_constructors(Tetrahedron_3_archetype)
#endif    
};

#if defined(CGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS)
inline bool operator==(const Tetrahedron_3_archetype& obj1, 
                       const Tetrahedron_3_archetype& obj2)
{ return true; }

inline bool operator!=(const Tetrahedron_3_archetype& obj1, 
                       const Tetrahedron_3_archetype& obj2)
{ return true; }
#endif

struct Cartesian_coordinate_const_iterator_3_archetype {
  Cartesian_coordinate_const_iterator_3_archetype() { }
  Cartesian_coordinate_const_iterator_3_archetype& 
  operator=(const Cartesian_coordinate_const_iterator_3_archetype&) 
  { return *this; }
};

CGAL_END_NAMESPACE

#undef CGAL_concept_archetype_constructors

#endif
