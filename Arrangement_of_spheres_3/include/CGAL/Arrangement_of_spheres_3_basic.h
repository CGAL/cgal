#ifndef CGAL_ARRANGEMENT_OF_SPHERES_3_BASIC_H
#define CGAL_ARRANGEMENT_OF_SPHERES_3_BASIC_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Simple_cartesian.h>

//#define CGAL_ARRANGEMENT_OF_SPHERES_3_USE_TEMPLATES





#define CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE CGAL_BEGIN_NAMESPACE namespace CGAL_Arrangement_of_spheres_3_internal {
#define CGAL_AOS3_END_INTERNAL_NAMESPACE CGAL_END_NAMESPACE }

#define CGAL_AOS3_INTERNAL_NS CGAL::CGAL_Arrangement_of_spheres_3_internal

#include <CGAL/Arrangement_of_spheres_3/Coordinate_index.h>
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE
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
CGAL_AOS3_END_INTERNAL_NAMESPACE




#ifndef CGAL_AOS3_USE_TEMPLATES
#define CGAL_AOS3_TYPENAME

#define CGAL_AOS3_TEMPLATE
#define CGAL_AOS3_TARG
#define CGAL_AOS3_TRAITS typedef Arrangement_of_spheres_traits_3 Traits;
namespace CGAL{
typedef struct Geom_traits: Simple_cartesian<double>{} Arrangement_of_spheres_3_geom_traits;
}

#include <CGAL/Arrangement_of_spheres_traits_3.h>





#else
#define CGAL_AOS3_TYPENAME typename
#define CGAL_AOS3_TEMPLATE template <class Traits_t>
#define CGAL_AOS3_TARG <Traits_t>
#define CGAL_AOS3_TRAITS typedef Traits_t Traits;
#define CGAL_AOS3_BEGIN_IMPL
#define CGAL_AOS3_END_IMPL


#endif



#endif
