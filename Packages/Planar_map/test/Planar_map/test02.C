#include <CGAL/basic.h>

#ifndef VOID_TEST
#include "numrep1.h"
#include <cassert>
#endif //VOID_TEST

#include <iostream>

#ifndef VOID_TEST
#include <CGAL/Pm_default_dcel.h>

#include <CGAL/Pm_segment_traits_2.h>

#if STRATEGY == 2
#include <CGAL/Pm_naive_point_location.h>
#else
#if STRATEGY == 3
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif
#endif

#include <CGAL/Planar_map_2.h>

#include "numrep2.h"


typedef CGAL::Pm_default_dcel<Traits>   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits> Planar_map;
typedef Traits::Point_2                 Point_2;
typedef Traits::X_curve                 X_curve_2;

typedef Planar_map::Halfedge_handle Halfedge_handle;
typedef Planar_map::Face_iterator Face_iterator;
#endif //VOID_TEST

int main()
{
#ifndef VOID_TEST
#if STRATEGY == 3
  CGAL::Pm_walk_along_line_point_location<Planar_map> pl;
  Planar_map Pm(&pl);
#else
#if STRATEGY == 2
  CGAL::Pm_naive_point_location<Planar_map> pl;
  Planar_map Pm(&pl);
#else
  Planar_map Pm;
#endif
#endif
 
  int n; std::cin >> n;
  while (n--) {
    inputt x1, y1, x2, y2;
    std::cin >> x1 >> y1 >> x2 >> y2;
    Halfedge_handle h = Pm.insert(X_curve_2(Point_2(x1,y1),Point_2(x2,y2)));
    std::cout << "Inserted (" << h->curve() << ")" << std::endl;
  }

  assert(Pm.is_valid());

  X_curve_2 seg = (Pm.halfedges_begin())->curve();

  std::cout << "\nseg = " << seg << std::endl;

  Point_2 s = seg.source();
  Point_2 t = seg.target();
  Point_2 mid_point = s + (t-s)/TestR::RT(2);
  
  std::cout << "\nmid_point = " << mid_point << std::endl;

  std::cout << "\nsplitting edge " << (Pm.halfedges_begin())->curve();
  std::cout << " ... " ; 

  Halfedge_handle h = Pm.split_edge(Pm.halfedges_begin(),
                                    X_curve_2(s, mid_point),
                                    X_curve_2(mid_point, t));
  assert(Pm.is_valid());
  std::cout << "map valid" << std::endl;


  std::cout << "degree of h->target is " << h->target()->degree() << std::endl;
  std::cout << "h->next is " << h->next()->curve() << std::endl  ;

  std::cout << "merging... " ; 
  Pm.merge_edge(h,h->next_halfedge(),seg);
  assert(Pm.is_valid());
  std::cout << "map valid" << std::endl;  

  Face_iterator fit;
  Planar_map::Size i=0;
  for (fit=Pm.faces_begin(); fit != Pm.faces_end(); ++fit)
    ++i;
  assert(Pm.number_of_faces() == i);
#else //VOID_TEST
  usage();
#endif //VOID_TEST

  return 0;
}
