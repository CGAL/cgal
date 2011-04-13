#ifndef VOID_TEST
#include "numrep1.h"
#include <cassert>
#endif //VOID_TEST

#include <iostream>

#ifndef VOID_TEST
#include <CGAL/Pm_default_dcel.h>

#if TESTR == 1
#include <CGAL/Pm_segment_epsilon_traits.h>
#else
#include <CGAL/Pm_segment_exact_traits.h>
#endif

#if STRATEGY == 2
#include <CGAL/Pm_naive_point_location.h>
#else
#if STRATEGY == 3
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif
#endif

#include <CGAL/Planar_map_2.h>

#include "numrep2.h"


typedef CGAL::Pm_default_dcel< Traits > Dcel;
typedef CGAL::Planar_map_2< Dcel, Traits >  Planar_map;
typedef Traits::Point Point;
typedef Traits::X_curve Curve;

typedef Planar_map::Halfedge_handle Halfedge_handle;
typedef Planar_map::Face_iterator Face_iterator;
#endif //VOID_TEST

int main(int argc, char *argv[])
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
    std::cout << "Inserted (";
    std::cout << (Pm.insert(Curve(Point(x1,y1),Point(x2,y2))))->curve();
    std::cout <<")" << std::endl;
  }

  assert(Pm.is_valid());

  Curve seg = (Pm.halfedges_begin())->curve();

  std::cout << "\nseg = " << seg << std::endl;

  Point s=seg.source();
  Point t=seg.target();
  Point mid_point=s+(t-s)/TestR::RT(2);
  
  std::cout << "\nmid_point = " << mid_point << std::endl;

  std::cout << "\nsplitting edge " << (Pm.halfedges_begin())->curve();
  std::cout << " ... " ; 

  Halfedge_handle h = Pm.split_edge(Pm.halfedges_begin(),
                                    Curve(s,mid_point),Curve(mid_point,t));
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



