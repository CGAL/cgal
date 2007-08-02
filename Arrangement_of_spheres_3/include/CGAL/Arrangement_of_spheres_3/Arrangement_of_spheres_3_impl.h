#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_initializer.h>


CGAL_BEGIN_NAMESPACE

CGAL_AOS3_TEMPLATE
void
Arrangement_of_spheres_3 CGAL_AOS3_TARG::initialize_at(CGAL_AOS3_TYPENAME Traits::FT end,
						       Cross_section &cs) {
  cs.set_number_of_spheres(tr_.number_of_sphere_3s());
  std::cout << "Initialize at " << end << std::endl;
  CGAL_AOS3_INTERNAL_NS::Cross_section_initializer CGAL_AOS3_TARG csi(cs, tr_);
  csi(end);
  
}


CGAL_AOS3_TEMPLATE
void
Arrangement_of_spheres_3 CGAL_AOS3_TARG::sweep_to(CGAL_AOS3_TYPENAME Traits::FT end,
						  Cross_section &cs) {
  std::cout << "Sweeping to " << end << std::endl;
  
}


CGAL_END_NAMESPACE
