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

typedef CGAL::Pm_default_dcel< Traits >     Dcel;
typedef CGAL::Planar_map_2< Dcel, Traits >  Planar_map;
typedef Traits::Point                       Point;
typedef Traits::X_curve                     Curve;

typedef Planar_map::Vertex_handle           Vertex_handle;
typedef Planar_map::Halfedge_handle         Halfedge_handle;
typedef Planar_map::Face_handle             Face_handle;

typedef Planar_map::Face_iterator           Face_iterator;
typedef Planar_map::Halfedge_iterator       Halfedge_iterator;
typedef Planar_map::Vertex_iterator         Vertex_iterator;
typedef Planar_map::Holes_iterator          Holes_iterator;
typedef Planar_map::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;
typedef Planar_map::Ccb_halfedge_circulator Ccb_halfedge_circulator;
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

  Vertex_handle   v;
  Halfedge_handle e;
  
  int n; std::cin >> n;
  while (n--) {
    inputt x1, y1, x2, y2;
    std::cin >> x1 >> y1 >> x2 >> y2;
    e = Pm.insert(Curve(Point(x1,y1),Point(x2,y2)));
    v = e->source();
    std::cout << "Inserted "<< e->curve() << std::endl;
  }

  assert(Pm.is_valid());

  Halfedge_iterator eit;
  for (eit=Pm.halfedges_begin(); eit != Pm.halfedges_end(); eit++)
    std::cout << (*eit).curve() << std::endl;
  
  Halfedge_around_vertex_circulator eavcirc=v->incident_halfedges(), 
                                    eavcirc_done(eavcirc);

  std::cout << "Halfedges incident to " << v->point() << std::endl;
  do {
    std::cout << (*eavcirc).curve() << std::endl;
  } while (eavcirc != eavcirc_done);

  // Vertices test
  // This used to be tst26, (Shai, Aug. 06, 2000)
  Vertex_iterator vit;
  std::cout << "All vertices" << std::endl;
  for (vit=Pm.vertices_begin(); vit != Pm.vertices_end(); vit++)
    std::cout << (*vit).point() << std::endl;

  // Edges on ccb test
  // This used to be tst27, (Shai, Aug. 06, 2000)
  Ccb_halfedge_circulator hcirc = e->ccb(), hcirc_done(hcirc);
  std::cout << "All edges on ccb of edge " << e->curve() << std::endl;
  do {
    std::cout << (*hcirc).curve() << std::endl;
    ++hcirc;
  } while (hcirc != hcirc_done);

  // Holes of unbounded face test
  // This used to be tst28, (Shai, Aug. 06, 2000)
  Face_handle f = Pm.unbounded_face();
  Holes_iterator hit;
  std::cout << "All holes of unbounded face" << std::endl;
  for (hit=f->holes_begin(); hit != f->holes_end(); hit++)
    std::cout << (*hit)->curve() << std::endl;
#else //VOID_TEST
  usage();
#endif //VOID_TEST

  return 0;
}








