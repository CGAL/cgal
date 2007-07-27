#ifndef CGAL_ARRANGEMENT_OF_SPHERES_3_H
#define CGAL_ARRANGEMENT_OF_SPHERES_3_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>

#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>

CGAL_BEGIN_NAMESPACE
CGAL_AOS3_TEMPLATE
class Arrangement_of_spheres_3 {
  CGAL_AOS3_TRAITS;

    Traits tr_;
public:
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Combinatorial_cross_section CGAL_AOS3_TARG Cross_section;


  Arrangement_of_spheres_3(const Traits& tr): tr_(tr){
  }
  
  
  
  void sweep_to(CGAL_AOS3_TYPENAME Traits::FT max, Cross_section &cs);
  void initialize_at(CGAL_AOS3_TYPENAME Traits::FT z, Cross_section &cs);
  
};


CGAL_END_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Arrangement_of_spheres_3_impl.h>
#endif

#endif
