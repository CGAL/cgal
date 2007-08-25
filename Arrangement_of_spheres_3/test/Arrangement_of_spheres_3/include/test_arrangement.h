#ifndef CGAL_TEST_ARRANGEMENT_H
#define CGAL_TEST_ARRANGEMENT_H
#include <CGAL/IO/qt_debug_examiner_viewer_2.h>
#include <vector>


template <class A> 
void test_arrangement() {
  typedef typename A::Traits::Geom_traits K;
  typedef typename A::Traits Traits_t;
  typedef typename K::Point_3 P;
  typedef typename K::Sphere_3 S;
  std::vector<S> spheres;
  spheres.push_back(S(P(10,10,10),1));
  typename A::Traits tr(spheres.begin(), spheres.end());
  A arr(tr);
  typename A::Cross_section cs;
  //cs.set_number_of_spheres(tr.number_of_sphere_3s());
  arr.initialize_at(10,cs);

  //typedef typename CGAL_AOS3_INTERNAL_NS::Cross_section_qt_viewer CGAL_AOS3_TARG CSV;
  //CSV csv(tr, cs, CGAL::qt_debug_examiner_viewer_2__);
  //csv(10);
}
 
#endif
