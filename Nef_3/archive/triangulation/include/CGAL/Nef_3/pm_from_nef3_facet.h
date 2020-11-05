#ifndef PM_FROM_NEF3_FACET_2_H
#define PM_FROM_NEF3_FACET_2_H

#undef _DEBUG
#define _DEBUG 5
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <typename Cycle_circulator, typename Planar_map>
typename Planar_map::Face_handle
add_cycle( Cycle_circulator c, typename Planar_map::Face_handle outer_cycle, Planar_map& pm) {

  typedef typename Planar_map::X_curve X_curve;
  typedef typename Planar_map::Halfedge_handle Halfedge_handle;
  typedef typename Planar_map::Face_handle Face_handle;
  typedef typename Planar_map::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;

  // edge 1
  Cycle_circulator curr(c), prev(curr);
  curr++;
  CGAL_assertion( curr != prev);
  X_curve first_segment( *prev, *curr);
  CGAL_NEF_TRACEN("inserting first segment "<<first_segment);
  Halfedge_handle first_edge =
    pm.insert_in_face_interior( first_segment, outer_cycle);

  // edge 2 ... n-1
  Cycle_circulator done(c);
  done--;
  Halfedge_handle edge = first_edge;
  while( curr != done) {
    prev++;
    curr++;
    X_curve segment( *prev, *curr);
    Halfedge_around_vertex_circulator hav(edge->target()->incident_halfedges()), hend(hav);
    CGAL_For_all( hav, hend)
      if( hav->source()->point() == segment.target())
        break;
    if( hav != hend) { // segment already inserted
      CGAL_NEF_TRACEN("segment already inserted (bubble): "<<segment);
      edge = hav->twin();
    }
    else if( edge->source()->point() == segment.target()) {
      CGAL_NEF_TRACEN("segment already inserted (dagger): "<<segment);
      edge = edge->twin();
    }
    else {
      CGAL_NEF_TRACEN("interseting segment "<<segment);
      edge = pm.insert_from_vertex( segment, edge);
    }
  }

  // edge n
  CGAL_assertion( *c == first_edge->source()->point());
  X_curve last_segment( *curr, *c);
  Halfedge_around_vertex_circulator hav(edge->target()->incident_halfedges()), hend(hav);
  CGAL_For_all( hav, hend)
    if( hav->source()->point() == last_segment.target())
      break;
  Halfedge_handle last_edge;
  if( hav != hend) { // segment already inserted
    CGAL_NEF_TRACEN("last segment already inserted (bubble): "<<last_segment);
    last_edge = hav->twin();
  }
    else if( edge->source()->point() == last_segment.target()) {
      CGAL_NEF_TRACEN("last segment already inserted (dagger): "<<last_segment);
      last_edge = edge->twin();
    }
  else {
    CGAL_NEF_TRACEN("inserting last segment "<<last_segment);
    last_edge = pm.insert_at_vertices( last_segment, edge->target(), first_edge->source());
  }

  Face_handle face = last_edge->face();
  if( face->is_unbounded())
    face = last_edge->twin()->face();

  return face;
}


template <typename SNC_structure, typename Planar_map, typename Projector>
typename Planar_map::Face_handle
add_hole( typename SNC_structure::Halffacet_cycle_const_iterator ci,
          typename Planar_map::Face_handle face, Planar_map& pm,
          Projector projector) {

  typedef typename Planar_map::Face_handle Face_handle;
  typedef typename SNC_structure::SNC_decorator SNC_decorator;
  typedef typename SNC_structure::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator
    SHalfedge_around_facet_const_circulator;

  typedef typename SNC_structure::Kernel Kernel;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Point_3 Point_3;

  typedef std::vector<Point_2> Point_2_vector;
  typedef typename Point_2_vector::const_iterator Point_2_vector_iterator;
  typedef Circulator_from_iterator<Point_2_vector_iterator> Circulator;

  SNC_decorator D;
  Face_handle hole;
  if( ci.is_shalfedge()) {
    SHalfedge_const_handle se(ci);
    SHalfedge_around_facet_const_circulator sc(se), send(sc);
    Point_2_vector cycle;
    cycle.reserve(circulator_distance(sc, send));
    CGAL_For_all( sc, send)
      cycle.push_back(projector(sc->source()->source()->point()));
    CGAL_NEF_TRACEN("number of vertices on cycle: "<<circulator_distance( sc, send));
#ifdef _DEBUG
    //CGAL_NEF_TRACEN("vertices on cycle:");
    //std::copy( cycle.begin(), cycle.end(), std::ostream_iterator<Point_2>( std::cerr, "\n"));
#endif
    hole = add_cycle( Circulator(cycle.begin(), cycle.end()), face, pm);
  }
  else if(ci.is_shalfloop()) {
    CGAL_warning_msg( 0, "isolated points are not supported on planar maps");
  }
  else
    CGAL_error_msg( "wrong handle");
  return hole;
}

template <typename SNC_structure, typename Planar_map, typename Projector>
typename Planar_map::Face_handle
add_face( typename SNC_structure::Halffacet_cycle_const_iterator ci,
           typename Planar_map::Face_handle face, Planar_map& pm,
           Projector projector) {
  typedef typename SNC_structure::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator
    SHalfedge_around_facet_const_circulator;
  CGAL_precondition( ci.is_shalfedge());
  CGAL_assertion_code(SHalfedge_const_handle se(ci));
  CGAL_assertion_code(SHalfedge_around_facet_const_circulator sc(se));
  CGAL_assertion_code(SHalfedge_around_facet_const_circulator send(sc));
  CGAL_precondition( circulator_distance(sc, send) > 2);
  return add_hole<SNC_structure>( ci, face, pm, projector);
}

template <typename SNC_structure, typename Planar_map, typename Projector>
typename Planar_map::Face_handle
pm_from_nef3_facet( typename SNC_structure::Halffacet_const_handle f,
                    Planar_map& pm,
                    Projector projector) {
  typedef typename SNC_structure::Halffacet_cycle_const_iterator
    Halffacet_cycle_const_iterator;
  typedef typename Planar_map::Face_handle Face_handle;

  Halffacet_cycle_const_iterator hci = f->facet_cycles_begin();
  CGAL_NEF_TRACEN("adding outer boundary...");
  Face_handle pm_face = add_face<SNC_structure>
    ( hci, pm.unbounded_face(), pm, projector);
  ++hci;
  CGAL_assertion( !pm_face->is_unbounded());
  for( ; hci != f->facet_cycles_end(); ++hci) {
    CGAL_NEF_TRACEN("adding hole...");
    add_hole<SNC_structure>( hci, pm_face, pm, projector);
  }
  return pm_face;
}

} //namespace CGAL

#endif // PM_FROM_NEF3_FACET_2_H
