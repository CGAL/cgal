#ifndef PM_FROM_POLYGON_WITH_HOLES_2_H
#define PM_FROM_POLYGON_WITH_HOLES_2_H

#undef _DEBUG 
#define _DEBUG 5
#include <CGAL/Nef_2/debug.h>

// Construction of a planar subdivision from a polygon with holes
//
// Requirements.
//
// Input: A set of circulators over the boundaries of a simple polygon
//        where the first circulator is the outer boundary and the rest
//        are inner boundaries.
//        Outer boundary is ccw oriented, inner boundaries are cw oriented.
// Ouput: A planar map divided by the input polygon

namespace CGAL {

template <class VertexCirculator, class PM>
typename PM::Face_handle
insert_cycle_in_face_interior( VertexCirculator vertices, 
			       typename PM::Face_handle bounding_face,
			       PM& pm) 
{
  typedef typename PM::X_curve X_curve;
  typedef typename PM::Halfedge_handle Halfedge_handle;
  typedef typename PM::Face_handle Face_handle;

  // edge 1
  VertexCirculator curr(vertices), prev(curr);
  curr++;
  CGAL_assertion( curr != prev);
  X_curve first_segment( *prev, *curr);
  Halfedge_handle first_edge = pm.insert_in_face_interior( first_segment, 
							   bounding_face);

  // edge 2 ... n-1
  VertexCirculator done(vertices);
  done--;
  CGAL_assertion( curr != done);
  Halfedge_handle edge = first_edge;
  do {
    prev++;
    curr++;
    X_curve segment( *prev, *curr);
    edge = pm.insert_from_vertex( segment, edge);
  }
  while( curr != done);

  // edge n
  X_curve last_segment( *curr, first_edge->source()->point());
  edge = pm.insert_at_vertices( last_segment, 
				edge->target(), first_edge->source());
  return edge->face();
}

template <class BoundaryIterator, class PM>
typename PM::Face_handle
pm_from_polygon_with_holes_2( BoundaryIterator begin, BoundaryIterator end, 
			      PM& pm) 
{
  typedef typename BoundaryIterator::value_type Circulator;
  typedef typename PM::Face_handle Face_handle;
 
  BoundaryIterator bi = begin;
  Circulator outer = *bi;
  // create outer cycle
  Face_handle polygon;
  polygon = insert_cycle_in_face_interior( outer, pm.unbounded_face(), pm);
  // create inner cycle(s)
  for( ++bi; bi != end; ++bi) {
    Circulator inner = *bi;
    insert_cycle_in_face_interior( inner, polygon, pm);
  }
  return polygon;
}

}

#endif // PM_FROM_POLYGON_WITH_HOLES_2_H
