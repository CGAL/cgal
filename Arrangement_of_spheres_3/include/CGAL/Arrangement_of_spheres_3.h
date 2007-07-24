#ifndef CGAL_ARRANGEMENT_OF_SPHERES_3_H
#define CGAL_ARRANGEMENT_OF_SPHERES_3_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>

#include <CGAL/Arrangement_of_spheres_3/Cross_section.h>

CGAL_BEGIN_NAMESPACE
CGAL_AOS3_TEMPLATE
class Arrangement_of_spheres_3 {
  CGAL_AOS3_TRAITS;

  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Cross_section CGAL_AOS3_TARG Cross_section;
  Traits tr_;
  Cross_section cs_;
public:
  Arrangement_of_spheres_3(const Traits& tr): tr_(tr), cs_(tr_.number_of_sphere_3s()){
  }
  
  
  
  void sweep_to(CGAL_AOS3_TYPENAME Traits::FT max);
  void initialize_at(CGAL_AOS3_TYPENAME Traits::FT z);
  
  CGAL_GET(Cross_section, cross_section, return cs_);
};


CGAL_END_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Arrangement_of_spheres_3_impl.h>
#endif

#endif
