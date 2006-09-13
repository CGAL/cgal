#ifndef COORDINATE_SPHERE_3_H
#define COORDINATE_SPHERE_3_H

#include <CGAL/Tools/Coordinate_index.h>
/* Coordiante */

#if 0
inline Coordinate_index sweep_coordinate() {
  return Coordinate_index::Z();
}
inline Coordinate_index plane_coordinate(int i) {
  if (i ==0) return Coordinate_index::X();
  else return Coordinate_index::Y();
}
inline unsigned int project(Coordinate_index i) {
  if (i == Coordinate_index::X()) return 0;
  else return 1;
}

inline Coordinate_index other_plane_coordinate(Coordinate_index i) {
  if (i == Coordinate_index::Y()) return Coordinate_index::X();
  else {
    CGAL_assertion(i== Coordinate_index::X());
    return Coordinate_index::Y();
    
  }
}
#else 
inline Coordinate_index sweep_coordinate() {
  return Coordinate_index::X();
}
inline Coordinate_index plane_coordinate(int i) {
  if (i ==0) return Coordinate_index::Y();
  else return Coordinate_index::Z();
}
inline unsigned int project(Coordinate_index i) {
  if (i == Coordinate_index::Y()) return 0;
  else {
    CGAL_assertion(i== Coordinate_index::Z());
    return 1;
  }
}

inline Coordinate_index other_plane_coordinate(Coordinate_index i) {
  if (i == Coordinate_index::Y()) return Coordinate_index::Z();
  else {
    CGAL_assertion(i== Coordinate_index::Z());
    return Coordinate_index::Y();
  }
}

template <class Vector>
inline Vector sweep_vector() {
  int v[3]={0,0,0};
  v[sweep_coordinate().index()]=1;
  return Vector(v[0], v[1], v[2]);
}

template <class Point, class NT>
inline Point sweep_point(NT n) {
  NT v[3]={0,0,0};
  v[sweep_coordinate().index()]=n;
  return Point(v[0], v[1], v[2]);
}

#endif
#endif
