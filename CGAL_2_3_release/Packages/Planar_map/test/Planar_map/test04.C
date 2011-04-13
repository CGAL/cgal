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

typedef Planar_map::Halfedge_iterator       Halfedge_iterator;
typedef Planar_map::Locate_type             Locate_type;
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

  int n, count_n; 
  std::cin >> n;
  count_n = n;
  while (n--) {
    inputt x1, y1, x2, y2;
    std::cin >> x1 >> y1 >> x2 >> y2;
    e = Pm.insert(Curve(Point(x1,y1),Point(x2,y2)));
    v = e->source();
    std::cout << "Inserted "<< e->curve() << std::endl;
  }

  assert(Pm.is_valid());

  // Check remove function of the Planar Map
  Halfedge_iterator eit, end_eit;
  Point             mid_point, s, t;
  Locate_type       lt;
  Curve             seg;
  
  eit     = Pm.halfedges_begin();
  end_eit = Pm.halfedges_end();
  while ( eit != end_eit && 
          count_n > 0 ) // to avoid infinite loops in case of error
    {
      seg = eit->curve();
      std::cout << "Removing "<< seg << std::endl;

      s         = seg.source();
      t         = seg.target();
      mid_point = s+(t-s)/TestR::RT(2);

      Pm.remove_edge(eit);
      assert(Pm.is_valid());

      //verify that edge is not in planar map

      e = Pm.locate(mid_point,lt);
      CGAL_assertion(lt==Planar_map::FACE ||
                     lt==Planar_map::UNBOUNDED_FACE);
      std::cout << "middle point of removed edge, " << mid_point << " - ";
      if (lt==Planar_map::UNBOUNDED_FACE)
	std::cout << "UNBOUNDED_FACE" << std::endl;
      else if (lt==Planar_map::FACE)
	std::cout << "FACE" << std::endl;

      eit     = Pm.halfedges_begin();
      end_eit = Pm.halfedges_end();
      count_n--;
    }

  // verify that exactly n edges were removed.
  CGAL_assertion(count_n == 0);

  // verify that Planar Map is empty
  std::cout << "Number of faces   : " << Pm.number_of_faces() << std::endl;
  std::cout << "Number of halfedges " << Pm.number_of_halfedges() << std::endl;
  std::cout << "Number of vertices: " << Pm.number_of_vertices() << std::endl;
  CGAL_assertion(Pm.number_of_faces() == 1);
  CGAL_assertion(Pm.number_of_halfedges() == 0);
  CGAL_assertion(Pm.number_of_vertices() == 0);
#else //VOID_TEST
  usage();
#endif //VOID_TEST

  return 0;
}








