#ifndef TRIANGULATE_POLYGON_WITH_HOLES_2_H
#define TRIANGULATE_POLYGON_WITH_HOLES_2_H

#include <CGAL/pm_from_polygon_with_holes_2.h>
#include <CGAL/partition_y_monotone_2.h>
#include <CGAL/triangulate_monotone_polygon_2.h>
#include <CGAL/Unique_hash_map.h>

#ifdef _DEBUG_WINDOW
#include <CGAL/IO/Pm_Window_stream.h>
extern CGAL::Window_stream W;
#endif

namespace CGAL {

template <class Iterator, class Pm>
inline void divide_pm_by_diagonals( Iterator begin, Iterator beyond, 
				    Pm& pm) {
  typedef typename Pm::X_curve X_curve;
  typedef typename Pm::Vertex_handle Vertex_handle;

  for( Iterator i = begin; i != beyond; ++i) {
    Vertex_handle v1(i->first.current_circulator()->source());
    Vertex_handle v2(i->second.current_circulator()->source());
    X_curve segment( v1->point(), v2->point());
    pm.insert_at_vertices( segment, v1, v2);
  }
}

template <class BoundaryInputIterator, class TriangleOutputIterator, 
	  class Traits>
void triangulate_polygon_with_holes_2( BoundaryInputIterator begin,
				    BoundaryInputIterator end,
				    TriangleOutputIterator triangles,
				    const Traits& traits) {

  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Triangle_2 Triangle_2;

  typedef typename Traits::Planar_map_2 Planar_map_2;
  typedef typename Traits::Partition_traits_2 Partition_traits_2;
  typedef typename Traits::Monotone_polygon_triangulation_traits_2
    Monotone_polygon_triangulation_traits_2;

  typedef typename Traits::Circulator_project Circulator_project;
  typedef std::list<Circulator_project> Cycles_list;
  typedef typename Cycles_list::const_iterator Cycle_iterator;

  typedef typename Traits::Diagonal Diagonal;
  typedef std::list<Diagonal> Diagonal_list;
  typedef typename Diagonal_list::const_iterator Diagonal_iterator;

#define USING(x) typedef typename Planar_map_2::x x;
  USING(Vertex_handle);
  USING(Face_handle);
  USING(Face_iterator);
  USING(Holes_iterator);
  USING(Ccb_halfedge_circulator);
#undef USING

  // create a planar map to support the face creation 

  Planar_map_2 pm;
  Face_handle polygon = pm_from_polygon_with_holes_2( begin, end, pm);
  //Face_handle polygon = pm_from_nef_facet( begin, end, pm);

#ifdef _DEBUG_WINDOW
  W.clear();
  W << BLUE << pm;
  Point_2 pause; W >> pause;
#endif
  
  // mark the faces that correspond to holes

  Unique_hash_map< Face_handle, bool> hole_mark(false);
  for( Holes_iterator hi = polygon->holes_begin(); 
       hi != polygon->holes_end();
       ++hi) {
    Face_handle hole_face = (*hi)->twin()->face();
    hole_mark[hole_face] = true;
  }

  // get a list of the cycles bounding the input polygon on the planar map

  Cycles_list cycles;  
  cycles.push_back(Circulator_project(polygon->outer_ccb()));
  for( Holes_iterator hi = polygon->holes_begin(); 
       hi != polygon->holes_end(); ++hi)
    cycles.push_back(Circulator_project((*hi)->ccb()));

  // divide the polygon into y-monotone pieces

  Diagonal_list diagonals;
  partition_y_monotone_2( cycles.begin(), cycles.end(), 
			  std::back_inserter(diagonals),
			  Partition_traits_2());
  divide_pm_by_diagonals( diagonals.begin(), diagonals.end(), pm);

  // mark the faces corresponding to y-monoton pieces
  
  Unique_hash_map< Face_handle, bool> piece_mark(false);
  for( Face_iterator fi = pm.faces_begin(); fi != pm.faces_end(); ++fi) {
    if( fi->is_unbounded()) 
      continue;
    if( hole_mark[fi])
      continue;
    piece_mark[fi] = true;
  }

#ifdef _DEBUG_WINDOW
  W.clear();
  W << GREEN << pm;
  W >> pause;
#endif

  // triangulate each monotone part ...

  for( Face_iterator fi = pm.faces_begin(); fi != pm.faces_end(); ++fi) {
    if( !piece_mark[fi])
      continue;
    Face_handle piece = fi;
    CGAL_assertion( std::distance( piece->holes_begin(),
				   piece->holes_end()) == 0); // no holes
    Diagonal_list diagonals;
    Circulator_project circulator(piece->outer_ccb());    
    triangulate_monotone_polygon_2( circulator, std::back_inserter(diagonals), 
				    Monotone_polygon_triangulation_traits_2());
    divide_pm_by_diagonals( diagonals.begin(), diagonals.end(), pm);
  }

#ifdef _DEBUG_WINDOW
  W.clear();
  W << YELLOW << pm;
#endif

  // and extract triangles from the topological map 
  
  for( Face_iterator fi = pm.faces_begin(); fi != pm.faces_end(); ++fi) {
    if( fi->is_unbounded()) 
      continue;
    if( hole_mark[fi])
      continue;
    Face_handle triangle = fi;
    CGAL_assertion( std::distance( triangle->holes_begin(),
				   triangle->holes_end()) == 0); // no holes
    Ccb_halfedge_circulator c(triangle->outer_ccb());
    CGAL_assertion_code( Ccb_halfedge_circulator cend(c));
    Point_2 p, q, r;
    p = c->source()->point();
    c++;
    q = c->source()->point();
    c++;
    r = c->source()->point();
    CGAL_assertion( ++c == cend);
    *triangles++ = Triangle_2( p, q, r);
  }

  return;
}

}

#endif // TRIANGULATE_POLYGON_WITH_HOLES_2_H

