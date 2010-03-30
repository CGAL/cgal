// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>


// class implementation continued
//=================================

CGAL_BEGIN_NAMESPACE

//====================================================================
//====================================================================
//                   CONSTRUCTORS
//====================================================================
//====================================================================

// copy constructor
template<class Gt, class D_S, class LTag>
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
Segment_Delaunay_graph_nox_2(const Segment_Delaunay_graph_nox_2& other)
  : DG(other.geom_traits())
{
  Segment_Delaunay_graph_nox_2&
    non_const_other = const_cast<Segment_Delaunay_graph_nox_2&>(other);
  copy(non_const_other);
  CGAL_postcondition( is_valid() );
}

// assignment operator
template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Self&
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
operator=(const Self& other)
{
  if ( this != &other ) {
    Segment_Delaunay_graph_nox_2&
      non_const_other = const_cast<Segment_Delaunay_graph_nox_2&>(other);
    copy(non_const_other);
  }
  return (*this);
}

//====================================================================
//====================================================================
//                   METHODS FOR INSERTION
//====================================================================
//====================================================================

//--------------------------------------------------------------------
// insertion of first three sites
//--------------------------------------------------------------------

template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
insert_first(const Site_2& s)
{
  CGAL_precondition( number_of_vertices() == 0 );

  Vertex_handle v = this->_tds.insert_second();
  v->set_site(s);
  return v;
}

template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
insert_second(const Site_2& s)
{
  CGAL_precondition( number_of_vertices() == 1 );
  CGAL_assertion( s.is_point() );

  Vertex_handle v0(finite_vertices_begin());

  // MK: change the equality test between points by the functor in
  // geometric traits

  // v0->site() is actually a point
  if ( same_points(s, v0->site()) ) {
#if 0
    // merge info of identical sites
    merge_info(v0, s);
#endif
    return v0;
  }

  return create_vertex_dim_up(s);
}


template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
insert_third(const Site_2& t)
{
  CGAL_precondition( number_of_vertices() == 2 );
  CGAL_precondition( t.is_point() );

  if ( only_points ) {
    return DG::insert(t.point());
  }

  // p0 and p1 are actually points
  Vertex_handle v0 = finite_vertices_begin();
  Vertex_handle v1 = ++finite_vertices_begin();
  const Site_2& t0 = v0->site();
  const Site_2& t1 = v1->site();

  if ( same_points(t, t0) ) {
#if 0
    // merge info of identical sites
    merge_info(v0, ss);
#endif
    return v0;
  }
  if ( same_points(t, t1) ) {
#if 0
    // merge info of identical sites
    merge_info(v1, ss);
#endif
    return v1;
  }

  Vertex_handle v = create_vertex_dim_up(t);

  Face_handle f(finite_faces_begin());

  const Site_2& s1 = f->vertex(0)->site();
  const Site_2& s2 = f->vertex(1)->site();
  const Site_2& s3 = f->vertex(2)->site();

  Orientation o =
    geom_traits().orientation_2_object()(s1, s2, s3);

  if ( o != COLLINEAR ) {
    if ( o == RIGHT_TURN ) {
      f->reorient();
      for (int i = 0; i < 3; i++) {
	f->neighbor(i)->reorient();
      }
    }
  } else {
    typename Geom_traits::Compare_x_2 compare_x =
      geom_traits().compare_x_2_object();

    Comparison_result xcmp12 = compare_x(s1, s2);
    if ( xcmp12 == SMALLER ) {        // x1 < x2
      Comparison_result xcmp23 = compare_x(s2, s3);
      if ( xcmp23 == SMALLER ) {            // x2 < x3
	flip(f, f->index(v1));
      } else {
	Comparison_result xcmp31 = compare_x(s3, s1);
	if ( xcmp31 == SMALLER ) {          // x3 < x1
	  flip(f, f->index(v0));
	} else {                            // x1 < x3 < x2
	  flip(f, f->index(v)); 
	}
      }
    } else if ( xcmp12 == LARGER ) {  // x1 > x2
      Comparison_result xcmp32 = compare_x(s3, s2);
      if ( xcmp32 == SMALLER ) {            // x3 < x2
	flip(f, f->index(v1));
      } else {
	Comparison_result xcmp13 = compare_x(s1, s3);
	if ( xcmp13 == SMALLER ) {          // x1 < x3
	  flip(f, f->index(v0));
	} else {                            // x2 < x3 < x1
	  flip(f, f->index(v));
	}
      }
    } else {                          // x1 == x2
      typename Geom_traits::Compare_y_2 compare_y =
	geom_traits().compare_y_2_object();

      Comparison_result ycmp12 = compare_y(s1, s2);
      if ( ycmp12 == SMALLER ) {      // y1 < y2
	Comparison_result ycmp23 = compare_y(s2, s3);
	if ( ycmp23 == SMALLER ) {          // y2 < y3
	  flip(f, f->index(v1));
	} else {
	  Comparison_result ycmp31 = compare_y(s3, s1);
	  if ( ycmp31 == SMALLER ) {        // y3 < y1
	    flip(f, f->index(v0));
	  } else {                          // y1 < y3 < y2
	    flip(f, f->index(v));
	  }
	}
      } else if ( ycmp12 == LARGER ) { // y1 > y2
	Comparison_result ycmp32 = compare_y(s3, s2);
	if ( ycmp32 == SMALLER ) {           // y3 < y2
	  flip(f, f->index(v1));
	} else {
	  Comparison_result ycmp13 = compare_y(s1, s3);
	  if ( ycmp13 == SMALLER ) {         // y1 < y3
	    flip(f, f->index(v0));
	  } else {                           // y2 < y3 < y1
	    flip(f, f->index(v));
	  }
	}
      } else {
	// this line should never have been reached
	CGAL_error();
      }
    }
  }

  return v;
}


template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
insert_third(const Site_2& s, Vertex_handle , Vertex_handle )
{
  CGAL_precondition( number_of_vertices() == 2 );

  //  this can only be the case if the first site is a segment
  CGAL_precondition( dimension() == 1 );

  only_points = false;

  Vertex_handle v = create_vertex_dim_up(s);

  Face_circulator fc = incident_faces(v);

  while ( true ) {
    Face_handle f(fc);
    if ( !is_infinite(f) ) {
      flip(f, f->index(v));
      break;
    }
    ++fc;
  }
  
  return v;
}

//--------------------------------------------------------------------
// insertion of a point
//--------------------------------------------------------------------



template<class Gt, class D_S, class LTag>
bool
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
is_flippable(Face_handle f, int i)
{
  // an edge is considered flippable if both its endpoints have degree
  // at least 3
  Vertex_handle v1 = f->vertex(ccw(i));
  {
    int count = 0;
    Vertex_circulator vc = incident_vertices(v1), done(vc);
    do {
      ++count;
    } while ( count < 3 || ++vc != done );
    if ( count < 3 ) return false;
  }

  Vertex_handle v2 = f->vertex(cw(i));
  int count = 0;
  Vertex_circulator vc = incident_vertices(v2), done(vc);
  do {
    ++count;
  } while ( count < 3 || ++vc != done );
  return count == 3;
}


template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
propagating_flip(Face_handle f, int i, const Site_2& t)
{
  if ( !is_flippable(f, i) ) return;

  Face_handle n = f->neighbor(i);

  //  Edge esym = sym_edge(f, i);
  //  if (  is_infinite( esym.first->vertex(esym.second) )   ) { return; }

  Sign s = incircle(n, t);

  if ( s != NEGATIVE ) { return; }

  bool interior_in_conflict = edge_interior(f, i, t, s);
  
  if ( interior_in_conflict ) {
    flip(f, i);
    return;
  }

  flip(f, i);
  propagating_flip(f, i, t);
  i = n->index( f->vertex(i) );
  propagating_flip(n, i, t);
}


template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
insert_point(const Site_2& t, Vertex_handle vnear)
{
  CGAL_precondition( t.is_point() );

  if ( only_points ) {
    if ( vnear == Vertex_handle() ) {
      return DG::insert(t.point());
    } else {
      return DG::insert(t.point(), vnear->face());
    }
  }

  int n = number_of_vertices();
  if ( n == 0 ) {
    return insert_first(t);
  } else if ( n == 1 ) {
    return insert_second(t);
  } else if ( n == 2 ) {
    return insert_third(t);
  }

  CGAL_assertion( number_of_vertices() > 2 );

  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  // MK::ERROR: I need to write a insert_point_no_search method that
  // does not search for the nearest neighbor; this should be used by
  // insert_point. Below the first version of the code is correct. The
  // second is what the insert_point method should do before calling
  // insert_point_no_search.

  // first find the nearest neighbor
#if 1
  Vertex_handle  vnearest = nearest_neighbor( t, vnear );
#else
  Vertex_handle vnearest;
  if ( vnear == Vertex_handle() ) {
    vnearest = nearest_neighbor( t, vnear );
  } else {
    vnearest = vnear;
  }
#endif

  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  //*********************************************************************


  // check is the point has already been inserted or lies on
  // a segment already inserted

  //  Arrangement_type at_res = arrangement_type(t, vnearest);
  if ( vnearest->site().is_point() ) {
    //    if ( at_res == AT2::IDENTICAL ) {
    if ( same_points(vnearest->site(), t) ) {
#if 0
      // merge info of identical sites
      merge_info(vnearest, ss);
#endif
      return vnearest;
    }
  } else {
#if !defined(CGAL_NO_ASSERTIONS) && !defined(NDEBUG)
    Arrangement_type at_res = arrangement_type(t, vnearest);
    CGAL_assertion( vnearest->site().is_segment() );
    CGAL_assertion( at_res != AT2::TOUCH_1 );
    CGAL_assertion( at_res != AT2::TOUCH_2 );
    CGAL_assertion( at_res == AT2::DISJOINT );
#endif
  }

  return insert_point2(t, vnearest);
}


template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
insert_point2(const Site_2& t, Vertex_handle vnearest)
{
  CGAL_precondition( t.is_point() );
  CGAL_assertion( number_of_vertices() > 2 );

  CGAL_expensive_precondition
    ( nearest_neighbor(t, vnearest) == vnearest );

  // find the first conflict

#if !defined(CGAL_NO_ASSERTIONS) && !defined(NDEBUG)
  // verify that there are no intersections...
  Vertex_circulator vc = incident_vertices(vnearest);
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    Arrangement_type at_res = arrangement_type(t, vv);

    CGAL_assertion( at_res == AT2::DISJOINT );
    ++vc;
  } while ( vc != vc_start );
#endif

  // first look for conflict with vertex
  Face_circulator fc_start = incident_faces(vnearest);
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;

  do {
    Face_handle f(fc);

    s = incircle(f, t);

    f->tds_data().set_incircle_sign(s);

    if ( s == NEGATIVE ) {
      start_f = f;
      break;
    }
    ++fc;
  } while ( fc != fc_start );

  // we are not in conflict with a Voronoi vertex, so we have to
  // be in conflict with the interior of a Voronoi edge
  if ( s != NEGATIVE ) {
    Edge_circulator ec_start = incident_edges(vnearest);
    Edge_circulator ec = ec_start;

    bool interior_in_conflict(false);
    Edge e;
    do {
      e = *ec;

      Sign s1 = e.first->tds_data().incircle_sign();
      Sign s2 = e.first->neighbor(e.second)->tds_data().incircle_sign();

      if ( s1 == s2 ) {
	interior_in_conflict = edge_interior(e, t, s1);
      } else {
	// It seems that there was a problem here when one of the
	// signs was positive and the other zero. In this case we
	// still check pretending that both signs where positive
	interior_in_conflict = edge_interior(e, t, POSITIVE);
      }

      if ( interior_in_conflict ) { break; }
      ++ec;
    } while ( ec != ec_start );

#if 0
    Face_circulator fc_start = incident_faces(vnearest);
    Face_circulator fc = fc_start;

    do {
      fc->tds_data().clear();
      ++fc;
    } while ( fc != fc_start );
#endif

    CGAL_assertion( interior_in_conflict );

    return insert_degree_2(e, t);
  }

#ifdef CGAL_SDG_INSERT_WITH_FLIPS
  // this is the new vertex
  Vertex_handle v = this->_tds.insert_in_face(start_f);
  v->set_site(t);
  Face_handle f = v->face();
  Face_handle next;
  int i;
  Face_handle start(f);
  do {
    i = f->index(v);
    next = f->neighbor(ccw(i));
    propagating_flip(f, i, t);
    f = next;
  } while ( next != start );
  return v;
#else
  // we are in conflict with a Voronoi vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;

  Triple<bool, Vertex_handle, Arrangement_type>
    vcross(false, Vertex_handle(), AT2::DISJOINT);

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.
  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, t, l, vcross);

  CGAL_assertion( !vcross.first );

  Vertex_handle v = create_vertex(t);

  retriangulate_conflict_region(v, l);

  return v;
#endif
}


//--------------------------------------------------------------------
// insertion of a segment
//--------------------------------------------------------------------
template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
insert_segment(const Site_2& t, Vertex_handle vnear)
{
  CGAL_precondition( t.is_segment() );

  if ( is_degenerate_segment(t) ) {
#if 0
    Storage_site_2 ss_src = ss.source_site();
    convert_info(ss_src, ss, true);
#endif
    return insert_point(t.source(), vnear);
  }

#if 0
  Storage_site_2 ss_src = ss.source_site();
  convert_info(ss_src, ss, true);
  Storage_site_2 ss_trg = ss.target_site();
  convert_info(ss_trg, ss, false);
#endif

  Vertex_handle v0 = insert_point( t.source(), vnear );
  CGAL_assertion( is_valid() );
  Vertex_handle v1 = insert_point( t.target(), v0 );
  CGAL_assertion( is_valid() );

  only_points = false;

  if ( number_of_vertices() == 2 ) {
    return insert_third(t, v0, v1);
  }

  return insert_segment_interior(t, v0);
}


template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
insert_segment_interior(const Site_2& t, Vertex_handle vnearest)
{
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( number_of_vertices() > 2 );

  CGAL_assertion( vnearest != Vertex_handle() );

  only_points = false;

  // find the first conflict

  // first check if the segment has already been inserted...
  //***************************************************************
  // MK:: What happens when I ask for the nearest neighbor of a
  // segment and the nearest neighbor is the segment itself??????
  //***************************************************************
  Vertex_circulator vc = incident_vertices(vnearest);
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    if ( is_infinite(vv) ) {
      vc++;
      continue;
    }

    if ( vv->site().is_segment() && same_segments(vv->site(), t) ) {
#if 0
      // merge info of identical items
      merge_info(vv, ss);
#endif
      return vv;
    }
    ++vc;
  } while ( vc != vc_start );

  // first look for conflict with vertex
  Face_circulator fc_start = incident_faces(vnearest);
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;

  do {
    Face_handle f(fc);

    s = incircle(f, t);

    f->tds_data().set_incircle_sign(s);

    if ( s == NEGATIVE ) {
      start_f = f;
      break;
    }
    ++fc;
  } while ( fc != fc_start );

  // segments must have a conflict with at least one vertex
  CGAL_assertion( s == NEGATIVE );

#ifdef CGAL_SDG_INSERT_WITH_FLIPS
  // this is the new vertex
  Vertex_handle v = this->_tds.insert_in_face(start_f);
  v->set_site(t);
  Face_handle f = v->face();
  Face_handle next;
  int i;
  Face_handle start(f);
  do {
    i = f->index(v);
    next = f->neighbor(ccw(i));
    propagating_flip(f, i, t);
    f = next;
  } while ( next != start );
  return v;
#else
  // we are in conflict with a Voronoi vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;

  Triple<bool, Vertex_handle, Arrangement_type>
    vcross(false, Vertex_handle(), AT2::DISJOINT);

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.
  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, t, l, vcross);

  CGAL_assertion( vcross.third == AT2::DISJOINT ||
		  vcross.third == AT2::CROSSING ||
		  vcross.third == AT2::INTERIOR );

  // no intersecting segment has been found; we insert the segment as
  // usual...
  Vertex_handle v = create_vertex(t);

  retriangulate_conflict_region(v, l);

  return v;
#endif
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------
// helper methods for insertion (find conflict region)
//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
initialize_conflict_region(const Face_handle& f, List& l)
{


  l.clear();
  for (int i = 0; i < 3; i++) {
    l.push_back(sym_edge(f, i));
  }
}


template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
expand_conflict_region(const Face_handle& f, const Site_2& t,
		       List& l,
		       Triple<bool,Vertex_handle,Arrangement_type>& vcross)
{
  if ( f->tds_data().is_in_conflict() ) { return; }

  // this is done to stop the recursion when intersecting segments
  // are found
  if ( vcross.first ) { return; }

  // setting fm[f] to true means that the face has been reached and
  // that the face is available for recycling. If we do not want the
  // face to be available for recycling we must set this flag to
  // false.
  f->tds_data().mark_in_conflict();
  fhc_.push_back(f);

  //  CGAL_assertion( fm.find(f) != fm.end() );

  for (int i = 0; i < 3; i++) {
    Face_handle n = f->neighbor(i);

    bool face_registered = n->tds_data().is_in_conflict();

    if ( !face_registered ) {
      for (int j = 0; j < 3; j++) {
	Vertex_handle vf = n->vertex(j);

	if ( is_infinite(vf) ) { continue; }

#if 0
	Arrangement_type at_res = arrangement_type(t, vf);

	CGAL_assertion( vcross.third == AT2::DISJOINT ||
			vcross.third == AT2::CROSSING ||
			vcross.third == AT2::INTERIOR );

	if ( vf->is_segment() ) {
	  CGAL_assertion( at_res != AT2::IDENTICAL );
	  CGAL_assertion( at_res != AT2::TOUCH_11_INTERIOR_1 );
	  CGAL_assertion( at_res != AT2::TOUCH_12_INTERIOR_1 );

	  if ( at_res == AT2::CROSSING ) {
	    vcross.first = true;
	    vcross.second = vf;
	    vcross.third = AT2::CROSSING;
	    l.clear();
	    fhc_.clear();
	    return;
	  } else {
	    CGAL_assertion ( at_res == AT2::DISJOINT ||
			     at_res == AT2::TOUCH_1 ||
			     at_res == AT2::TOUCH_2 ||
			     at_res == AT2::TOUCH_11 ||
			     at_res == AT2::TOUCH_12 ||
			     at_res == AT2::TOUCH_21 ||
			     at_res == AT2::TOUCH_22 );
	    // we do nothing in these cases
	  }
	} else {
	  CGAL_assertion( vf->is_point() );
	  if ( at_res == AT2::INTERIOR ) {
	    vcross.first = true;
	    vcross.second = vf;
	    vcross.third = AT2::INTERIOR;
	    l.clear();
	    fhc_.clear();
	    return;
	  }
	}
#endif
      }
    }

    Sign s = incircle(n, t);

    n->tds_data().set_incircle_sign(s);

    Sign s_f = f->tds_data().incircle_sign();

    if ( s == POSITIVE ) { continue; }
    if ( s != s_f ) { continue; }

    bool interior_in_conflict = edge_interior(f, i, t, s);

    if ( !interior_in_conflict ) { continue; }

    if ( face_registered ) { continue; }

    Edge e = sym_edge(f, i);

    CGAL_assertion( l.is_in_list(e) );
    int j = this->_tds.mirror_index(f, i);
    Edge e_before = sym_edge(n, ccw(j));
    Edge e_after = sym_edge(n, cw(j));
    if ( !l.is_in_list(e_before) ) {
      l.insert_before(e, e_before);
    }
    if ( !l.is_in_list(e_after) ) {
      l.insert_after(e, e_after);
    }
    l.remove(e);

    expand_conflict_region(n, t, l, vcross);

    // this is done to stop the recursion when intersecting segments
    // are found
    //    if ( fm.size() == 0 && l.size() == 0 ) { return; }
    if ( vcross.first ) { return; }
  } // for-loop
}


//--------------------------------------------------------------------
// retriangulate conflict region
//--------------------------------------------------------------------

template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
add_bogus_vertex(Edge e, List& l)
{
  Edge esym = sym_edge(e);
  Face_handle g1 = e.first;
  Face_handle g2 = esym.first;

  Vertex_handle v = insert_degree_2(e);

  Face_circulator fc(v);
  Face_handle f1(fc);
  Face_handle f2(++fc);
  int i1 = f1->index(v);
  int i2 = f2->index(v);

  CGAL_assertion( ((f1->neighbor(i1) == g1) && (f2->neighbor(i2) == g2)) ||
		  ((f1->neighbor(i1) == g2) && (f2->neighbor(i2) == g1)) );

  Edge ee, eesym;
  if ( f1->neighbor(i1) == g1 ) {
    ee = Edge(f2, i2);
    eesym = Edge(f1, i1);
  } else {
    ee = Edge(f1, i1);
    eesym = Edge(f2, i2);
  }

  l.replace(e, ee);
  l.replace(esym, eesym);

  return v;
}


template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_list
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
add_bogus_vertices(List& l)
{
#if defined(USE_INPLACE_LIST)
  Vertex_list vertex_list;

  Edge e_start = l.front();
  Edge e = e_start;

  std::list<Edge> edge_list;

  do {
    Edge esym = sym_edge(e);
    if ( l.is_in_list(esym) ) {
      if ( !esym.first->tds_data().is_selected(esym.second) ) {
	e.first->tds_data().mark_selected(e.second);
	edge_list.push_back(e);
      }
    }
    e = l.next(e);
  } while ( e != e_start );

  e_start = l.front();
  e = e_start;
  do {
    if ( e.first->tds_data().is_selected(e.second) ) {
      e.first->tds_data().mark_unselected(e.second);
    }
  } while ( e != e_start );

  typename std::list<Edge>::iterator it;

  for (it = edge_list.begin();  it != edge_list.end(); ++it) {
    Vertex_handle v = add_bogus_vertex(*it, l);
    vertex_list.push_back(v);
  }

  return vertex_list;
#else
  Vertex_list vertex_list;

  std::set<Edge> edge_list;

  edge_list.clear();

  Edge e_start = l.front();
  Edge e = e_start;

  do {
    Edge esym = sym_edge(e);
    if ( l.is_in_list(esym) &&
	 edge_list.find(esym) == edge_list.end() ) {
      edge_list.insert(e);
    }
    e = l.next(e);
  } while ( e != e_start );

  typename std::set<Edge>::iterator it;

  for (it = edge_list.begin();  it != edge_list.end(); ++it) {
    Vertex_handle v = add_bogus_vertex(*it, l);
    vertex_list.push_back(v);
  }

  return vertex_list;
#endif
}

template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
remove_bogus_vertices(Vertex_list& vl)
{
  while ( vl.size() > 0 ) {
    Vertex_handle v = vl.front();
    vl.pop_front();
    remove_degree_2(v);
  }
}


template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
retriangulate_conflict_region(Vertex_handle v, List& l)
{
  // 1. add the bogus vetrices
  Vertex_list dummy_vertices = add_bogus_vertices(l);

  // 2. repair the face pointers...
  Edge e_start = l.front();
  Edge eit = e_start;
  do {
    Edge esym = sym_edge(eit);
    Face_handle f = eit.first;
    int k = eit.second;
    CGAL_assertion( !l.is_in_list(esym) );
    CGAL_assertion( !f->tds_data().is_in_conflict() );

    f->vertex(ccw(k))->set_face(f);
    f->vertex( cw(k))->set_face(f);
    eit = l.next(eit);
  } while ( eit != e_start );

  // 3. copy the edge list to a vector of edges and clear the edge list
  // MK:: here I actually need to copy the edges to an std::list<Edge>, or
  // even better add iterators to the list of type List
#if 0
  std::vector<Edge> ve(l.size());

  Edge efront = l.front();
  Edge e = efront;
  unsigned int k = 0;
  do {
    ve[k] = e;
    ++k;
    e = l.next(e);
  } while ( e != efront );

  l.clear();

  // 4. retriangulate the hole
  this->_tds.star_hole(v, ve.begin(), ve.end());
#else
  this->_tds.star_hole(v, l.begin(), l.end());
  l.clear();
#endif

  // 5. remove the bogus vertices
  remove_bogus_vertices(dummy_vertices);

  // 6. remove the unused faces
  typename std::vector<Face_handle>::iterator it;
  for (it = fhc_.begin(); it != fhc_.end(); ++it) {
    (*it)->tds_data().clear();
    this->_tds.delete_face( *it );
  }

  fhc_.clear();

  // 7. DONE!!!!
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------
// combinatorial operations
//--------------------------------------------------------------------
//--------------------------------------------------------------------

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// point location
//--------------------------------------------------------------------
//--------------------------------------------------------------------
template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
nearest_neighbor(const Site_2& p,
		 Vertex_handle start_vertex) const
{
  CGAL_precondition( p.is_point() );

  if ( number_of_vertices() == 0 ) {
    return Vertex_handle();
  }

  if ( start_vertex == Vertex_handle() ) {
    start_vertex = finite_vertex();
  }

  //  if ( start_vertex == NULL ) { return start_vertex; }

  Vertex_handle vclosest;
  Vertex_handle v = start_vertex;

  if ( number_of_vertices() < 3 ) {
    vclosest = v;
    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      Vertex_handle v1(vit);
      if ( v1 != vclosest /*&& !is_infinite(v1)*/ ) {
	Site_2 t0 = vclosest->site();
	Site_2 t1 = v1->site();
	if ( side_of_bisector(t0, t1, p) == ON_NEGATIVE_SIDE ) {
	  vclosest = v1;
	}
      }
    }
    return vclosest;
  }

  do {
    vclosest = v;
    Site_2 t0 = v->site();
    //    if ( t0.is_point() && same_points(p, t0) ) {
    //      return vclosest;
    //    }
    Vertex_circulator vc_start = incident_vertices(v);
    Vertex_circulator vc = vc_start;
    do {
      if ( !is_infinite(vc) ) {
	Vertex_handle v1(vc);
	Site_2 t1 = v1->site();
	Oriented_side os = side_of_bisector(t0, t1, p);

	if ( os == ON_NEGATIVE_SIDE ) {
	  v = v1;
	  break;
	}
      }
      ++vc;
    } while ( vc != vc_start );
  } while ( vclosest != v );

  return vclosest;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// methods for the predicates
//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<class Gt, class D_S, class LTag>
Sign
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
incircle(const Face_handle& f, const Site_2& q) const
{
  if ( !is_infinite(f) ) {
    return incircle(f->vertex(0)->site(),
		    f->vertex(1)->site(),
		    f->vertex(2)->site(), q);
  }

  int inf_i(-1); // to avoid compiler warning
  for (int i = 0; i < 3; i++) {
    if ( is_infinite(f->vertex(i)) ) {
      inf_i = i;
      break;
    }
  }
  return incircle( f->vertex( ccw(inf_i) )->site(),
		   f->vertex(  cw(inf_i) )->site(), q );
}


template<class Gt, class D_S, class LTag>
Sign
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
incircle(const Vertex_handle& v0, const Vertex_handle& v1,
	      const Vertex_handle& v2, const Vertex_handle& v) const
{
  CGAL_precondition( !is_infinite(v) );

  if ( !is_infinite(v0) && !is_infinite(v1) &&
       !is_infinite(v2) ) {
    return incircle(v0->site(), v1->site(),
		    v2->site(), v->site());
  }

  if ( is_infinite(v0) ) {
    CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) );
    return incircle( v1->site(), v2->site(), v->site());
  }
  if ( is_infinite(v1) ) {
    CGAL_precondition( !is_infinite(v0) && !is_infinite(v2) );
    return incircle( v2->site(), v0->site(), v->site());
  }

  CGAL_assertion( is_infinite(v2) );
  CGAL_precondition( !is_infinite(v0) && !is_infinite(v1) );
  return incircle( v0->site(), v1->site(), v->site());
}


// this the finite edge interior predicate for a degenerate edge
template<class Gt, class D_S, class LTag>
bool
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
finite_edge_interior(const Face_handle& f, int i, const Site_2& q,
		     Sign sgn, int) const
{
  if ( !is_infinite( this->_tds.mirror_vertex(f, i) ) ) {
    CGAL_precondition( is_infinite(f->vertex(i)) );

    Face_handle g = f->neighbor(i);
    int j = this->_tds.mirror_index(f, i);

    return finite_edge_interior(g, j, q, sgn, 0 /* degenerate */);
  }

  CGAL_precondition( is_infinite( this->_tds.mirror_vertex(f, i) ) );

  const Site_2& t1 = f->vertex( ccw(i) )->site();
  const Site_2& t2 = f->vertex(  cw(i) )->site();

  if ( is_infinite(f->vertex(i)) ) {
    return finite_edge_interior(t1, t2, q, sgn);
  }

  const Site_2& t3 = f->vertex(i)->site();
  return finite_edge_interior(t1, t2, t3, q, sgn);
}

template<class Gt, class D_S, class LTag>
bool
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
finite_edge_interior(const Vertex_handle& v1, const Vertex_handle& v2,
		     const Vertex_handle& v3, const Vertex_handle& v4,
		     const Vertex_handle& v, Sign sgn, int) const
{
  CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) && 
		     !is_infinite(v) );
  if ( !is_infinite( v4 ) ) {
    CGAL_precondition( is_infinite(v3) );

    return
      finite_edge_interior(v2, v1, v4, v3, v, sgn, 0 /* degenerate */);
  }

  CGAL_precondition( is_infinite( v4 ) );

  const Site_2& t1 = v1->site();
  const Site_2& t2 = v2->site();
  const Site_2& q = v->site();

  if ( is_infinite(v3) ) {
    return finite_edge_interior(t1, t2, q, sgn);
  }

  const Site_2& t3 = v3->site();
  return finite_edge_interior(t1, t2, t3, q, sgn);
}

template<class Gt, class D_S, class LTag>
bool
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
infinite_edge_interior(const Face_handle& f, int i,
		       const Site_2& q, Sign sgn) const
{
  if ( !is_infinite( f->vertex(ccw(i)) ) ) {
    CGAL_precondition( is_infinite( f->vertex(cw(i)) ) );
    Face_handle g = f->neighbor(i);
    int j = this->_tds.mirror_index(f, i);

    return infinite_edge_interior(g, j, q, sgn);
  }

  CGAL_precondition( is_infinite( f->vertex(ccw(i)) ) );

  const Site_2& t2 = f->vertex(  cw(i) )->site();
  const Site_2& t3 = f->vertex(     i  )->site();
  const Site_2& t4 = this->_tds.mirror_vertex(f, i)->site();

  return infinite_edge_interior(t2, t3, t4, q, sgn);
}


template<class Gt, class D_S, class LTag>
bool
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
infinite_edge_interior(const Vertex_handle& v1,
		       const Vertex_handle& v2,
		       const Vertex_handle& v3,
		       const Vertex_handle& v4,
		       const Vertex_handle& v, Sign sgn) const
{
  CGAL_precondition( !is_infinite(v3) && !is_infinite(v4) && 
		     !is_infinite(v) );

  if ( !is_infinite( v1 ) ) {
    CGAL_precondition( is_infinite( v2 ) );

    return infinite_edge_interior(v2, v1, v4, v3, v, sgn);
  }

  CGAL_precondition( is_infinite( v1 ) );

  const Site_2& t2 = v2->site();
  const Site_2& t3 = v3->site();
  const Site_2& t4 = v4->site();
  const Site_2& q = v->site();

  return infinite_edge_interior(t2, t3, t4, q, sgn);
}




template<class Gt, class D_S, class LTag>
bool
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
edge_interior(const Vertex_handle& v1,
	      const Vertex_handle& v2,
	      const Vertex_handle& v3,
	      const Vertex_handle& v4,
	      const Vertex_handle& v, Sign sgn) const
{
  CGAL_precondition( !is_infinite(v) );

  bool is_inf_v1 = is_infinite(v1);
  bool is_inf_v2 = is_infinite(v2);
  bool is_inf_v3 = is_infinite(v3);
  bool is_inf_v4 = is_infinite(v4);

  bool result;

  if ( !is_inf_v1 && !is_inf_v2 && !is_inf_v3 && !is_inf_v4 ) {
    result = finite_edge_interior(v1, v2, v3, v4, v, sgn);
  } else if ( is_inf_v3 || is_inf_v4 ) {
    result = finite_edge_interior(v1, v2, v3, v4, v, sgn, 0/* degenerate */);
  } else {
    result = infinite_edge_interior(v1, v2, v3, v4, v, sgn);
  }

  return result;
}


template<class Gt, class D_S, class LTag>
bool
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
edge_interior(const Face_handle& f, int i,
	      const Site_2& q, Sign sgn) const
{
  Face_handle g = f->neighbor(i);

  bool is_inf_f = is_infinite(f);
  bool is_inf_g = is_infinite(g);

  bool result;

  if ( !is_inf_f && !is_inf_g ) {
    result = finite_edge_interior(f, i, q, sgn);
  } else if ( !is_inf_f || !is_inf_g ) {
    result = finite_edge_interior(f, i, q, sgn, 0 /* denegerate */);
  } else {
    if ( !is_infinite(f, i) ) {
      result = finite_edge_interior(f, i, q, sgn, 0 /* degenerate */);
    } else {
      result = infinite_edge_interior(f, i, q, sgn);
    }
  }

  return result;
}


template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Arrangement_type
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
arrangement_type(const Site_2& p, const Site_2& q) const
{
  typedef typename Geom_traits::Arrangement_type_2  AT2;
  typedef typename AT2::result_type                 Arrangement_type;

  Arrangement_type res = geom_traits().arrangement_type_2_object()(p, q);

  // The valeus that have to be treated are the following:
  // DISJOINT, TOUCH_1, TOUCH_2, CROSSING, IDENTICAL, INTERIOR,
  // TOUCH_11_INTERIOR_1, TOUCH_12_INTERIOR_1, TOUCH_21_INTERIOR_1 and
  // TOUCH_22_INTERIOR_1.
  //
  // The remaining values will either never appear because of one of
  // the following reasons:
  // 1. we insert the endpoints of the segments first and then the
  //    interior (OVERLAPPING_*, INTERIOR_*, TOUCH_*_INTERIOR_2).
  // 2. the values have no meaning since we consider the segments to
  //    be open (TOUCH_INTERIOR_*). In this case, the conflict will
  //    appear when we test with the endpoint.
  // 3. a conflict will first happen with an endpoint before testing
  //    for the segment (TOUCH_2*_INTERIOR_1). In this case the
  //    segment to be inserted will first find an endpoint in its
  //    interior before actually finding that there is another segment
  //    it overlaps with.

  CGAL_assertion( res != AT2::INTERIOR_1 );
  CGAL_assertion( res != AT2::INTERIOR_2 );

  CGAL_assertion( res != AT2::OVERLAPPING_11 );
  CGAL_assertion( res != AT2::OVERLAPPING_12 );
  CGAL_assertion( res != AT2::OVERLAPPING_21 );
  CGAL_assertion( res != AT2::OVERLAPPING_22 );

  CGAL_assertion( res != AT2::TOUCH_11_INTERIOR_2 );
  CGAL_assertion( res != AT2::TOUCH_21_INTERIOR_2 );
  CGAL_assertion( res != AT2::TOUCH_12_INTERIOR_2 );
  CGAL_assertion( res != AT2::TOUCH_22_INTERIOR_2 );

  CGAL_assertion( res != AT2::TOUCH_21_INTERIOR_1 );
  CGAL_assertion( res != AT2::TOUCH_22_INTERIOR_1 );

  if ( res == AT2::TOUCH_INTERIOR_12 || res == AT2::TOUCH_INTERIOR_21 ||
       res == AT2::TOUCH_INTERIOR_11 || res == AT2::TOUCH_INTERIOR_22 ) {
    return AT2::DISJOINT;
  }
  if ( res == AT2::TOUCH_11 || res == AT2::TOUCH_12 ||
       res == AT2::TOUCH_21 || res == AT2::TOUCH_22 ) {
    return AT2::DISJOINT;
  }

  return res;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// embedding and visualization methods and constructions for primal
// and dual
//--------------------------------------------------------------------
//--------------------------------------------------------------------

// primal
template<class Gt, class D_S, class LTag>
Object
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
primal(const Edge e) const
{
  typedef typename Gt::Line_2   Line_2;
  typedef typename Gt::Ray_2    Ray_2;

  CGAL_precondition( !is_infinite(e) );

  if ( this->dimension() == 1 ) {
    const Site_2& p = (e.first)->vertex(cw(e.second))->site();
    const Site_2& q = (e.first)->vertex(ccw(e.second))->site();

    Line_2 l = construct_sdg_bisector_2_object()(p,q);
    return make_object(l);
  }

  // dimension == 2
  // none of the two adjacent faces is infinite
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    const Site_2& p = (e.first)->vertex( ccw(e.second) )->site();
    const Site_2& q = (e.first)->vertex(  cw(e.second) )->site();
    const Site_2& r = (e.first)->vertex(     e.second  )->site();
    const Site_2& s = this->_tds.mirror_vertex(e.first, e.second)->site();
    return construct_sdg_bisector_segment_2_object()(p,q,r,s);
  }

  // both of the adjacent faces are infinite
  if ( is_infinite(e.first) &&
       is_infinite(e.first->neighbor(e.second)) )  {
    const Site_2& p = (e.first)->vertex(cw(e.second))->site();
    const Site_2& q = (e.first)->vertex(ccw(e.second))->site();
    Line_2 l = construct_sdg_bisector_2_object()(p,q);
    return make_object(l);
  }

  // only one of the adjacent faces is infinite
  CGAL_assertion( is_infinite( e.first ) ||
		  is_infinite( e.first->neighbor(e.second) )
		  );

  CGAL_assertion( !(is_infinite( e.first ) &&
		    is_infinite( e.first->neighbor(e.second) )
		    )
		  );

  CGAL_assertion( is_infinite(e.first->vertex(e.second)) ||
		  is_infinite(this->_tds.mirror_vertex(e.first, e.second)) );

  Edge ee = e;
  if ( is_infinite( e.first->vertex(e.second) )  ) {
    ee = sym_edge(e);
  }
  const Site_2& p = ee.first->vertex( ccw(ee.second) )->site();
  const Site_2& q = ee.first->vertex(  cw(ee.second) )->site();
  const Site_2& r = ee.first->vertex(     ee.second  )->site();

  Ray_2 ray = construct_sdg_bisector_ray_2_object()(p,q,r);
  return make_object(ray);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// validity test method
//--------------------------------------------------------------------
//--------------------------------------------------------------------
template<class Gt, class D_S, class LTag>
bool Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
is_valid(bool verbose, int level) const
{
  if (level < 0) { return true; }

  if (number_of_vertices() <= 1) {
    if ( verbose && number_of_vertices() == 1 ) {
      std::cerr << "SDGDS is ok... " << std::flush;
    }
    return true;
  }

  // level 0 test: check the TDS
  bool result = data_structure().is_valid(verbose, level);

  if ( result && verbose ) {
    std::cerr << "SDGDS is ok... " << std::flush;
  }

  if (level == 0) { return result; }

  // level 1 test: do the incircle tests
  if (number_of_vertices() < 3)  { return true; }

  for (All_edges_iterator eit = all_edges_begin();
       eit != all_edges_end(); ++eit) {
    Edge e = *eit;
    Face_handle f = e.first;

    Vertex_handle v = this->_tds.mirror_vertex(f, e.second);

    if ( f->vertex(e.second) == v ) { continue; }
    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->site()) != NEGATIVE );
    }
    Edge sym_e = sym_edge(e);
    f = sym_e.first;
    v = this->_tds.mirror_vertex(f, sym_e.second);

    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->site()) != NEGATIVE );
    }
  }

  if ( result && verbose ) {
    std::cerr << "Segment Delaunay graph is ok..." << std::flush;
  }
  if ( !result && verbose ) {
    std::cerr << "Segment Delaunay graph is NOT valid..." << std::flush;
  }

  return result;
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------
// misc
//--------------------------------------------------------------------
//--------------------------------------------------------------------


template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
print_error_message() const
{
  std::cerr << std::endl;
  std::cerr << "WARNING:" << std::endl;
  std::cerr << "A segment-segment intersection was found."
	    << std::endl;
  std::cerr << "The Segment_Delaunay_graph_nox_2 class is not configured"
	    << " to handle this situation." << std::endl;
  std::cerr << "Please look at the documentation on how to handle"
	    << " this behavior." << std::endl;
  std::cerr << std::endl;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// the copy method
//--------------------------------------------------------------------
//--------------------------------------------------------------------


template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
copy(Segment_Delaunay_graph_nox_2& other)
{
#if 0
  // second, copy the site representation info for the input sites
  // using the correct handles (i.e., the handles from the new point
  // container
  isc_.clear();
  typename Input_sites_container::iterator iit_other = other.isc_.begin();
  for (; iit_other != other.isc_.end(); ++iit_other) {
    Site_rep_2 old_srep = *iit_other;
    Site_rep_2 new_srep( hm[boost::tuples::get<0>(old_srep)],
			 hm[boost::tuples::get<1>(old_srep)],
			 boost::tuples::get<2>(old_srep) );
    isc_.insert( new_srep );
  }
  
  CGAL_assertion( pc_.size() == other.pc_.size() );
  CGAL_assertion( isc_.size() == other.isc_.size() );

#ifndef CGAL_NO_ASSERTIONS
  {
    Point_handle it_other = other.pc_.begin();
    Point_handle it_this = pc_.begin();
    for (; it_other != other.pc_.end(); ++it_other, ++it_this) {
      CGAL_assertion( *it_other == *it_this );
    }
  }
#endif
#endif
  // copy the value of the boolean flag concerning whether the diagram
  // has points only.
  only_points = other.only_points;

  // then copy the diagram
  DG::operator=(other);
#if 0
  // now we have to update the sotrage sites in each vertex of the
  // diagram and also update the 

  // then update the storage sites for each vertex
  Intersections_tag itag;

  Finite_vertices_iterator vit_other = other.finite_vertices_begin();
  Finite_vertices_iterator vit_this = finite_vertices_begin();
  for (; vit_other != other.finite_vertices_end(); vit_other++,
	 vit_this++) {
    Storage_site_2 ss_other = vit_other->storage_site();

#ifndef CGAL_NO_ASSERTIONS
    Storage_site_2 ss_this = vit_this->storage_site();
    if ( ss_other.is_segment() ) {
      CGAL_assertion( ss_this.is_segment() );
      CGAL_assertion( same_segments(ss_this.site(), ss_other.site()) );
    } else {
      CGAL_assertion( ss_this.is_point() );
      CGAL_assertion( same_points(ss_this.site(), ss_other.site()) );
    }
#endif

    Storage_site_2 new_ss_this = copy_storage_site(ss_other, hm, itag);
    vit_this->set_site( new_ss_this );
  }
#endif
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// getting endpoints of segments
//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
first_endpoint_of_segment(const Vertex_handle& v) const
{
  CGAL_assertion( v->is_segment() );
  // MK:: here I should be able to find a way to get the point without
  // constructing the Site_2 fe
  Site_2 fe = v->site().source_site();
  Vertex_circulator vc_start = incident_vertices(v);
  Vertex_circulator vc = vc_start;
  do {
    // Vertex_handle vv(vc);
    if ( !is_infinite(vc) && vc->is_point() ) {
      if ( same_points(fe, vc->site()) ) {
	return Vertex_handle(vc);
      }
    }
    vc++;
  } while ( vc != vc_start );

  // we should never reach this point
  CGAL_error();
  return Vertex_handle();
}

template<class Gt, class D_S, class LTag>
typename Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
second_endpoint_of_segment(const Vertex_handle& v) const
{
  CGAL_assertion( v->is_segment() );
  Site_2 fe = v->site().target_site();
  Vertex_circulator vc_start = incident_vertices(v);
  Vertex_circulator vc = vc_start;
  do {
    //      Vertex_handle vv(vc);
    if ( !is_infinite(vc) && vc->is_point() ) {
      if ( same_points(fe, vc->site()) ) {
	return Vertex_handle(vc);
      }
    }
    vc++;
  } while ( vc != vc_start );

  // we should never reach this point
  CGAL_error();
  return Vertex_handle();
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// file I/O
//--------------------------------------------------------------------
//--------------------------------------------------------------------

#if 0
template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
file_output(std::ostream& os, const Storage_site_2& t,
	    Point_handle_mapper& P) const
{
  CGAL_precondition( t.is_defined() );

  if ( t.is_point() ) {
    // 0 for point
    os << 0;
    if ( is_ascii(os) ) { os << ' '; }
    // 0 for input
    os << 0;
    if ( is_ascii(os) ) { os << ' '; }
    os << P[t.point()];
  } else { // t is a segment
    // 1 for segment
    os << 1;
    if ( is_ascii(os) ) { os << ' '; }
    // 0 for input
    os << 0;
    if ( is_ascii(os) ) { os << ' '; }
    os << P[t.source_of_supporting_site()];
    if ( is_ascii(os) ) { os << ' '; }
    os << P[t.target_of_supporting_site()];
  }
}

template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
file_input(std::istream& is, Storage_site_2& t,
	   const Point_handle_vector& P, const Tag_false&) const
{
  int type, input;
  is >> type >> input;
  CGAL_assertion( type == 0 || type == 1 );
  CGAL_assertion( input == 0 );
  if ( type == 0 ) {
    // we have an input point
    size_type p;
    is >> p;
    t = st_.construct_storage_site_2_object()(P[p]);
  } else {
    CGAL_assertion( type == 1 );
    // we have an input segment
    size_type p1, p2;
    is >> p1 >> p2;
    t = st_.construct_storage_site_2_object()(P[p1], P[p2]);
  }
}

template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
file_input(std::istream& is, Storage_site_2& t,
	   const Point_handle_vector& P, const Tag_true&) const
{
  int type, input;
  is >> type >> input;
  CGAL_assertion( type == 0 || type == 1 );
  CGAL_assertion( input >= 0 && input <= 3 );
  if ( type == 0 ) {
    // we have a point
    if ( input == 0 ) {
      // we have an input point
      size_type p;
      is >> p;
      t = st_.construct_storage_site_2_object()(P[p]);
    } else {
      // we have a point that is the intersection of two segments
      CGAL_assertion( input == 1 );
      size_type p1, p2, q1, q2;
      is >> p1 >> p2 >> q1 >> q2;
      t = st_.construct_storage_site_2_object()(P[p1], P[p2], P[q1], P[q2]);
    }
  } else {
    // we have a segment
    CGAL_assertion( type == 1 );
    if ( input == 0 ) {
      // we have an input segment
      size_type p1, p2;
      is >> p1 >> p2;
      t = st_.construct_storage_site_2_object()(P[p1], P[p2]);
    } else if ( input < 3 ) {
      // we have a segment whose source or target is input but not both
      size_type p1, p2, q1, q2;
      is >> p1 >> p2 >> q1 >> q2;
      t = st_.construct_storage_site_2_object()(P[p1], P[p2],
						P[q1], P[q2], input == 1);
    } else {
      // we have a segment whose neither its source nor its target is input
      CGAL_assertion( input == 3 );
      size_type p1, p2, q1, q2, r1, r2;
      is >> p1 >> p2 >> q1 >> q2 >> r1 >> r2;
      t = st_.construct_storage_site_2_object()(P[p1], P[p2],
						P[q1], P[q2],
						P[r1], P[r2]);      
    }
  }
}

//--------------------------------------------------------------------

template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
file_output(std::ostream& os, Point_handle_mapper& P,
	    bool print_point_container) const
{
  // ouput to a file
  size_type n = this->_tds.number_of_vertices();
  size_type m = this->_tds.number_of_full_dim_faces();

  CGAL_assertion( n >= 1 );

  if( is_ascii(os) ) {
    os << n << ' ' << m << ' ' << dimension() << std::endl;
  } else {
    os << n << m << dimension();
  }

  // points in point container and input sites container
  if ( print_point_container ) {
    if ( is_ascii(os) ) { os << std::endl; }
    os << pc_.size();
    if ( is_ascii(os) ) { os << std::endl; }
    for (const_Point_handle ph = pc_.begin(); ph != pc_.end(); ++ph) {
      os << *ph;
      if ( is_ascii(os) ) { os << std::endl; }
    }

    // print the input sites container
    if ( is_ascii(os) ) { os << std::endl; }
    os << isc_.size();
    if ( is_ascii(os) ) { os << std::endl; }
    for (typename Input_sites_container::const_iterator it = isc_.begin();
	 it != isc_.end(); ++it) {
      os << P[boost::tuples::get<0>(*it)];
      if ( is_ascii(os) ) { os << " "; }
      os << P[boost::tuples::get<1>(*it)];
      if ( is_ascii(os) ) { os << std::endl; }
    }
  }

  std::map<Vertex_handle,int> V;
  std::map<Face_handle,int> F;

  // first vertex (infinite vertex) 
  size_type inum = 0;
  V[infinite_vertex()] = inum++;
  
  // finite vertices
  if (is_ascii(os)) os << std::endl;
  for (Finite_vertices_iterator vit = finite_vertices_begin();
       vit != finite_vertices_end(); ++vit) {
    V[vit] = inum++;
    //    os << vit->site();
    file_output(os, vit->storage_site(), P);
    // write non-combinatorial info of the vertex
    //    os << *vit ;
    if ( is_ascii(os) ) { os << std::endl; }
  }
  if ( is_ascii(os) ) { os << std::endl; }

  // vertices of the faces
  inum = 0;
  int dim = (dimension() == -1 ? 1 :  dimension() + 1);
  for(All_faces_iterator fit = all_faces_begin();
      fit != all_faces_end(); ++fit) {
    F[fit] = inum++;
    for(int j = 0; j < dim ; ++j) {
      os << V[ fit->vertex(j) ];
      if( is_ascii(os) ) { os << ' '; }
    }
    // write non-combinatorial info of the face
    //    os << *fit ;
    if( is_ascii(os) ) { os << std::endl; }
  }
  if( is_ascii(os) ) { os << std::endl; }
    
  // neighbor pointers of the  faces
  for( All_faces_iterator it = all_faces_begin();
       it != all_faces_end(); ++it) {
    for(int j = 0; j < dimension()+1; ++j){
      os << F[ it->neighbor(j) ];
      if( is_ascii(os) ) { os << ' '; }
    }
    if( is_ascii(os) ) { os << std::endl; }
  }

  if ( is_ascii(os) ) { os << std::endl; }
}



template<class Gt, class D_S, class LTag>
void
Segment_Delaunay_graph_nox_2<Gt,D_S,LTag>::
file_input(std::istream& is, bool read_handle_vector,
	   Point_handle_vector& P)
{
  //input from file
  size_type n, m;
  int d;
  is >> n >> m >> d;

  CGAL_assertion( n >= 1 );

  size_type i = 0;
  Storage_site_2 ss;

  if ( read_handle_vector ) {
    pc_.clear();
    size_type np;
    is >> np;
    for (; i < np; i++) {
      Point_2 p;
      is >> p;
      std::pair<Point_handle,bool> res = pc_.insert(p);
      P.push_back(res.first);
      CGAL_assertion( P[i] == res.first );
    }

    // now read the input sites container
    isc_.clear();
    size_type nisc;
    is >> nisc;
    int id1, id2;
    for (i = 0; i < nisc; i++) {
      is >> id1 >> id2;
      isc_.insert( Site_rep_2(P[id1], P[id2], id1 == id2) );
    }
  }

  if ( n == 1 ) {
    CGAL_assertion( d == -1 );
    if ( number_of_vertices() > 0 ) { clear(); }
    return;
  }
  if ( n == 2 ) {
    CGAL_assertion( d == 0 );
    if ( number_of_vertices() > 0 ) { clear(); }
    file_input(is, ss, P, Intersections_tag());
    insert_first(ss, *ss.point());
    return;
  }
  if ( n == 3 ) {
    CGAL_assertion( d == 1 );
    if ( number_of_vertices() > 0 ) { clear(); }
    file_input(is, ss, P, Intersections_tag());
    insert_first(ss, *ss.point());  
    file_input(is, ss, P, Intersections_tag());
    insert_second(ss, *ss.point());  
    return;
  }

  if (this->_tds.number_of_vertices() != 0) { this->_tds.clear(); }

  this->_tds.set_dimension(d);

  std::vector<Vertex_handle> V(n);
  std::vector<Face_handle> F(m);

  // first vertex (infinite vertex)
  V[0] = this->_tds.create_vertex();
  this->set_infinite_vertex(V[0]);
  i = 1;

  // read vertices
  for (; i < n; ++i) {
    V[i] = this->_tds.create_vertex();
    file_input(is, ss, P, Intersections_tag());
    V[i]->set_site(ss);
    // read non-combinatorial info of the vertex
    //    is >> *(V[i]);
  }
  
  // Creation of the faces
  int index;
  int dim = (dimension() == -1 ? 1 : dimension() + 1);

  for (i = 0; i < m; ++i) {
    F[i] = this->_tds.create_face();
    for (int j = 0; j < dim ; ++j){
      is >> index;
      F[i]->set_vertex(j, V[index]);
      // The face pointer of vertices is set too often,
      // but otherwise we had to use a further map
      V[index]->set_face(F[i]);
    }
    // read in non-combinatorial info of the face
    //      is >> *(F[i]) ;
  }

  // Setting the neighbor pointers 
  for (i = 0; i < m; ++i) {
    for (int j = 0; j < dimension()+1; ++j){
      is >> index;
      F[i]->set_neighbor(j, F[index]);
    }
  }
}
#endif

CGAL_END_NAMESPACE

// EOF
