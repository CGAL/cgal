#ifndef CGAL_NEF3_DELEGATED_RAY_SHOOTER_H
#define CGAL_NEF3_DELEGATED_RAY_SHOOTER_H

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Nef_3/SNC_point_locator.h>

namespace CGAL {

template<typename SNC_>
class Ray_shooter : 
public Modifier_base<CGAL::SNC_point_locator<CGAL::SNC_decorator<SNC_> > {

  typedef SNC_                                   SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>     SNC_decorator;
  typedef CGAL::SNC_point_locator<SNC_decorator> SNC_point_locator;
  
  Object_handle& ores;
  Vector_3 dir;

  Ray_shooter(const Vector_3 vec, Object_handle& o) : ores(o), dir(vec);

  void operator()(SNC_point_locator& pl) {
    o_res = pl.shoot(dir);
  }

};

} //namespace CGAL
#endif // CGAL_NEF3_DELEGATED_RAY_SHOOTER_H
