// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
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
// $Source$
// $Revision$ $Date$
// $Name$
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
template<class Gt, class DS, class LTag>
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
Segment_Voronoi_diagram_2(const Segment_Voronoi_diagram_2& other)
  : DG(other.geom_traits())
{
  Segment_Voronoi_diagram_2&
    non_const_other = const_cast<Segment_Voronoi_diagram_2&>(other);
  copy(non_const_other);
  CGAL_postcondition( is_valid() );
}

// assignment operator
template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Self&
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
operator=(const Self& other)
{
  if ( this != &other ) {
    Segment_Voronoi_diagram_2&
      non_const_other = const_cast<Segment_Voronoi_diagram_2&>(other);
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

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_first(const Point_2& p)
{
  CGAL_precondition( number_of_vertices() == 0 );

  Storage_site_2 ss = create_storage_site(p);

  //  return create_vertex_dim_up(ss);
  Vertex_handle v = this->_tds.insert_second();
  v->set_site(ss);
  return v;
}

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_second(const Point_2& p)
{
  CGAL_precondition( number_of_vertices() == 1 );
  // p0 is actually a point
  Site_2 p0 = finite_vertices_begin()->site();
  // MK: change the equality test between points by the functor in
  // geometric traits
  Site_2 tp = Site_2::construct_site_2(p);
  if ( same_points(tp,p0) ) {
    return Vertex_handle(finite_vertices_begin());
  }

  Storage_site_2 ss = create_storage_site(p);
  return create_vertex_dim_up(ss);
}

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_third(const Point_2& p)
{
  Site_2 t = Site_2::construct_site_2(p);
  Storage_site_2 ss = create_storage_site(p);
  return insert_third(t, ss);
}

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_third(const Site_2& t, const Storage_site_2& ss)
{
  CGAL_precondition( number_of_vertices() == 2 );

  // p0 and p1 are actually points
  Vertex_handle v0 = finite_vertices_begin();
  Vertex_handle v1 = ++finite_vertices_begin();
  Site_2 t0 = v0->site();
  Site_2 t1 = v1->site();

  if ( same_points(t, t0) ) { return v0; }
  if ( same_points(t, t1) ) { return v1; }

  Vertex_handle v = create_vertex_dim_up(ss);

  Face_handle f(finite_faces_begin());

  Site_2 s1 = f->vertex(0)->site();
  Site_2 s2 = f->vertex(1)->site();
  Site_2 s3 = f->vertex(2)->site();

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
	CGAL_assertion( false );
      }
    }
  }

  return v;
}


template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_third(Vertex_handle v0, Vertex_handle v1)
{
  CGAL_precondition( number_of_vertices() == 2 );

  //  this can only be the case if the first site is a segment
  CGAL_precondition( dimension() == 1 );

  Storage_site_2 ss = create_storage_site(v0, v1);
  Vertex_handle v = create_vertex_dim_up(ss);

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

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_point(const Point_2& p, Vertex_handle vnear)
{
  int n = number_of_vertices();
  if ( n == 0 ) {
    return insert_first(p);
  } else if ( n == 1 ) {
    return insert_second(p);
  } else if ( n == 2 ) {
    return insert_third(p);
  }

  Site_2 t = Site_2::construct_site_2(p);
  Storage_site_2 ss = create_storage_site(p);
  return insert_point(ss, t, vnear);
}


template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_point(const Storage_site_2& ss, const Site_2& t,
	     Vertex_handle vnear)
{
  CGAL_precondition( t.is_point() );
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
  Arrangement_type at_res = arrangement_type(t, vnearest);
  if ( vnearest->is_point() ) {
    if ( at_res == AT2::IDENTICAL ) {
      return vnearest;
    }
  } else {
    CGAL_assertion( vnearest->is_segment() );
    CGAL_assertion( at_res != AT2::TOUCH_1 );
    CGAL_assertion( at_res != AT2::TOUCH_2 );
    CGAL_assertion( at_res == AT2::DISJOINT || at_res == AT2::INTERIOR );
    if ( at_res == AT2::INTERIOR ) {
      CGAL_assertion( t.is_exact() );

      Vertex_triple vt = insert_exact_point_on_segment(ss, t, vnearest);
      return vt.first;
    } else {
      // the point to be inserted does not belong to the interior of a
      // segment
      CGAL_assertion( at_res == AT2::DISJOINT );
    }
  }

  return insert_point2(ss, t, vnearest);
}

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_point2(const Storage_site_2& ss, const Site_2& t,
	      Vertex_handle vnearest)
{
  CGAL_precondition( t.is_point() );
  CGAL_assertion( number_of_vertices() > 2 );

  CGAL_expensive_precondition
    ( nearest_neighbor(t, vnearest) == vnearest );

  // find the first conflict

#ifndef NDEBUG
  // verify that there are no intersections...
  Vertex_circulator vc = vnearest->incident_vertices();
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    Arrangement_type at_res = arrangement_type(t, vv);

    CGAL_assertion( at_res == AT2::DISJOINT );
    ++vc;
  } while ( vc != vc_start );
#endif

  // first look for conflict with vertex
  Face_circulator fc_start = vnearest->incident_faces();
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;

  std::map<Face_handle,Sign> sign_map;

  do {
    Face_handle f(fc);

    s = incircle(f, t);

    sign_map[f] = s;

    if ( s == NEGATIVE ) {
      start_f = f;
      break;
    }
    ++fc;
  } while ( fc != fc_start );

  // we are not in conflict with a Voronoi vertex, so we have to
  // be in conflict with the interior of a Voronoi edge
  if ( s != NEGATIVE ) {
    Edge_circulator ec_start = vnearest->incident_edges();
    Edge_circulator ec = ec_start;

    bool interior_in_conflict(false);
    Edge e;
    do {
      e = *ec;

      Sign s1 = sign_map[e.first];
      Sign s2 = sign_map[e.first->neighbor(e.second)];

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

    sign_map.clear();

    CGAL_assertion( interior_in_conflict );

    return insert_degree_2(e, ss);
  }


  // we are in conflict with a Voronoi vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;
  Face_map fm;

  Triple<bool, Vertex_handle, Arrangement_type>
    vcross(false, Vertex_handle(), AT2::DISJOINT);

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.
  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, t, ss, l, fm, sign_map, vcross);

  CGAL_assertion( !vcross.first );

  Vertex_handle v = create_vertex(ss);

  retriangulate_conflict_region(v, l, fm);

  return v;
}

//--------------------------------------------------------------------
// insertion of a point that lies on a segment
//--------------------------------------------------------------------

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Face_pair
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
find_faces_to_split(const Vertex_handle& v, const Site_2& t) const
{
  CGAL_precondition( v->is_segment() );

#ifndef NDEBUG
  {
    // count number of adjacent infinite faces
    Face_circulator fc = incident_faces(v);
    Face_circulator fc_start = fc;
    int n_inf = 0;
    do {
      if ( is_infinite(fc) ) { n_inf++; }
      fc++;
    } while ( fc != fc_start );
    CGAL_assertion( n_inf == 0 || n_inf == 2 || n_inf == 4 );
  }
#endif

  Face_circulator fc1 = incident_faces(v);
  Face_circulator fc2 = fc1; ++fc2;
  Face_circulator fc_start = fc1;
  Face_handle f1, f2;
  bool found_f1 = false, found_f2 = false;
  Site_2 sitev_supp = v->site().supporting_site();
  do {
    Face_handle ff1(fc1), ff2(fc2);

    Oriented_side os1, os2;

    if ( is_infinite(ff1) ) {
      int id_v = ff1->index(v);
      int cw_v = this->cw( id_v );
      int ccw_v = this->ccw( id_v );

      Site_2 sv_ep;
      if ( is_infinite( ff1->vertex(cw_v) ) ) {
	CGAL_assertion(  !is_infinite( ff1->vertex(ccw_v) )  );
	CGAL_assertion( ff1->vertex(ccw_v)->site().is_point() );
	sv_ep = ff1->vertex(ccw_v)->site();
      } else {
	CGAL_assertion(  !is_infinite( ff1->vertex( cw_v) )  );
	CGAL_assertion( ff1->vertex( cw_v)->site().is_point() );
	sv_ep = ff1->vertex( cw_v)->site();
      }

      os1 = oriented_side(sv_ep, sitev_supp, t);
    } else {
      os1 = oriented_side(fc1->vertex(0)->site(),
			  fc1->vertex(1)->site(),
			  fc1->vertex(2)->site(),
			  sitev_supp, t);
    }

    if ( is_infinite(ff2) ) {
      int id_v = ff2->index(v);
      int cw_v = this->cw( id_v );
      int ccw_v = this->ccw( id_v );

      Site_2 sv_ep;
      if ( is_infinite( ff2->vertex(cw_v) ) ) {
	CGAL_assertion(  !is_infinite( ff2->vertex(ccw_v) )  );
	CGAL_assertion( ff2->vertex(ccw_v)->site().is_point() );
	sv_ep = ff2->vertex(ccw_v)->site();
      } else {
	CGAL_assertion(  !is_infinite( ff2->vertex( cw_v) )  );
	CGAL_assertion( ff2->vertex( cw_v)->site().is_point() );
	sv_ep = ff2->vertex( cw_v)->site();
      }

      os2 = oriented_side(sv_ep, sitev_supp, t);
    } else {
      os2 = oriented_side(fc2->vertex(0)->site(),
			  fc2->vertex(1)->site(),
			  fc2->vertex(2)->site(),
			  sitev_supp, t);
    }

    if ( !found_f1 &&
	 os1 != ON_POSITIVE_SIDE && os2 == ON_POSITIVE_SIDE ) {
      f1 = ff2;
      found_f1 = true;
    }

    if ( !found_f2 &&
	 os1 == ON_POSITIVE_SIDE && os2 != ON_POSITIVE_SIDE ) {
      f2 = ff2;
      found_f2 = true;
    }

    if ( found_f1 && found_f2 ) { break; }

    ++fc1, ++fc2;
  } while ( fc_start != fc1 ); 


  CGAL_assertion( found_f1 && found_f2 );
  CGAL_assertion( f1 != f2 );

  return Face_pair(f1, f2);
}

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_triple
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_exact_point_on_segment(const Storage_site_2& ss, const Site_2& t,
			      Vertex_handle v)
{
  // splits the segment site v->site() in two and inserts represented by t
  // on return the three vertices are, respectively, the vertex
  // corresponding to t and the two subsegments of v->site()

  CGAL_assertion( t.is_point() );
  CGAL_assertion( t.is_exact() );

  Storage_site_2 ssitev = v->storage_site();  

  CGAL_assertion( ssitev.is_segment() );

  Face_pair fpair = find_faces_to_split(v, t);

  Quadruple<Vertex_handle, Vertex_handle, Face_handle, Face_handle>
    qq = this->_tds.split_vertex(v, fpair.first, fpair.second);

  Intersections_tag itag;
  // now I need to update the sites for vertices v1 and v2
  Vertex_handle v1 = qq.first;
  Storage_site_2 ssv1 = split_storage_site(ssitev, ss, 0, itag);
  v1->set_site( ssv1 );

  Vertex_handle v2 = qq.second;
  Storage_site_2 ssv2 = split_storage_site(ssitev, ss, 1, itag);
  v2->set_site( ssv2 );

  Vertex_handle vsx =
    this->_tds.insert_in_edge(qq.third, cw(qq.third->index(v1)));

  vsx->set_site(ss);

  return Vertex_triple(vsx, v1, v2);
}

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_triple
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_point_on_segment(const Storage_site_2& ss, const Site_2& t,
			Vertex_handle v, const Tag_true&)
{
  // splits the segment site v->site() in two and inserts the point of
  // intersection of t and v->site()
  // on return the three vertices are, respectively, the point of
  // intersection and the two subsegments of v->site()

  Storage_site_2 ssitev = v->storage_site();
  Storage_site_2 ssx = create_storage_site(ss, ssitev);

  Site_2 sitev = ssitev.site();

  Face_pair fpair = find_faces_to_split(v, ssx.site());

  Quadruple<Vertex_handle, Vertex_handle, Face_handle, Face_handle>
    qq = this->_tds.split_vertex(v, fpair.first, fpair.second);

  // now I need to update the sites for vertices v1 and v2
  Vertex_handle v1 = qq.first;
  Storage_site_2 ssv1;
  Site_2 sv1;
  if ( sitev.is_exact(0) ) {
    ssv1 = create_storage_site(ssitev, ss, true);
  } else {
    ssv1 = create_storage_site_type1(ssitev, ssitev, ss);
  }
  sv1 = ssv1.site();
  v1->set_site( ssv1 );

  Vertex_handle v2 = qq.second;
  Storage_site_2 ssv2;
  Site_2 sv2;
  if ( sitev.is_exact(1) ) {
    ssv2 = create_storage_site(ssitev, ss, false);
  } else {
    ssv2 = create_storage_site_type2(ssitev, ss, ssitev);
  }
  sv2 = ssv2.site();
  v2->set_site( ssv2 );

  Vertex_handle vsx =
    this->_tds.insert_in_edge(qq.third, cw(qq.third->index(v1)));

  vsx->set_site(ssx);

  return Vertex_triple(vsx, v1, v2);
}

//--------------------------------------------------------------------
// insertion of a segment
//--------------------------------------------------------------------
template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_segment(const Site_2& t, Vertex_handle vnear)
{
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( t.is_exact() );

  if ( is_degenerate_segment(t) ) {
    return insert_point(t.source(), vnear);
  }

  Vertex_handle v0 = insert_point( t.source(), vnear );
  Vertex_handle v1 = insert_point( t.target(), v0 );

  if ( number_of_vertices() == 2 ) {
    return insert_third(v0, v1);
  }

  Storage_site_2 ss = create_storage_site(v0, v1);
  return insert_segment_interior(t, ss, v0);
}


template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_segment_interior(const Site_2& t, const Storage_site_2& ss,
			Vertex_handle vnearest)
{
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( number_of_vertices() > 2 );

  CGAL_assertion( vnearest != Vertex_handle() );

  // find the first conflict

  // first look if there are intersections...
  Vertex_circulator vc = vnearest->incident_vertices();
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    if ( is_infinite(vv) ) {
      vc++;
      continue;
    }

    Arrangement_type at_res = arrangement_type(t, vv);

    if ( vv->is_segment() ) {
      if ( at_res == AT2::DISJOINT || at_res == AT2::TOUCH_1 ||
	   at_res == AT2::TOUCH_2 || at_res == AT2::TOUCH_11 ||
	   at_res == AT2::TOUCH_12 || at_res == AT2::TOUCH_21 ||
	   at_res == AT2::TOUCH_22 ) {
	// do nothing
      } else if ( at_res == AT2::IDENTICAL ) {
	return vv;
      } else if ( at_res == AT2::CROSSING ) {
	Intersections_tag itag;
	return insert_intersecting_segment(ss, t, vv, itag);
      } else if ( at_res == AT2::TOUCH_11_INTERIOR_1 ) {
	Intersections_tag itag;

	Vertex_handle vp = second_endpoint_of_segment(vv);	
	Storage_site_2 ssvp = vp->storage_site();
	Storage_site_2 sss = split_storage_site(ss, ssvp, 1, itag);

	return insert_segment_interior(sss.site(), sss, vp);
      } else if ( at_res == AT2::TOUCH_12_INTERIOR_1 ) {
	Intersections_tag itag;

	Vertex_handle vp = first_endpoint_of_segment(vv);	
	Storage_site_2 ssvp = vp->storage_site();
	Storage_site_2 sss = split_storage_site(ss, ssvp, 0, itag);

	return insert_segment_interior(sss.site(), sss, vp);
      } else {
	// this should never be reached; the only possible values for
	// at_res are DISJOINT, CROSSING, TOUCH_11_INTERIOR_1
	// and TOUCH_12_INTERIOR_1
	CGAL_assertion( false );
      }
    } else {
      CGAL_assertion( vv->is_point() );
      if ( at_res == AT2::INTERIOR ) {
	Storage_site_2 ssvv = vv->storage_site();
	if ( ssvv.is_exact() ) {
	  Intersections_tag itag;
	  Storage_site_2 ss1 = split_storage_site(ss, ssvv, 0, itag);
	  Storage_site_2 ss2 = split_storage_site(ss, ssvv, 1, itag);
	  insert_segment_interior(ss1.site(), ss1, vv);
	  return insert_segment_interior(ss2.site(), ss2, vv);
	} else {
	  // this should never be reached; the only possible values for
	  // at_res are DISJOINT and INTERIOR
	  CGAL_assertion( false );
	}
      }
    }
    ++vc;
  } while ( vc != vc_start );

  // first look for conflict with vertex
  Face_circulator fc_start = vnearest->incident_faces();
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;

  std::map<Face_handle,Sign> sign_map;

  do {
    Face_handle f(fc);

    s = incircle(f, t);

    sign_map[f] = s;

    if ( s == NEGATIVE ) {
      start_f = f;
      break;
    }
    ++fc;
  } while ( fc != fc_start );

  // segments must have a conflict with at least one vertex
  CGAL_assertion( s == NEGATIVE );

  // we are in conflict with a Voronoi vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;
  Face_map fm;

  Triple<bool, Vertex_handle, Arrangement_type>
    vcross(false, Vertex_handle(), AT2::DISJOINT);

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.
  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, t, ss, l, fm, sign_map, vcross);

  CGAL_assertion( vcross.third == AT2::DISJOINT ||
		  vcross.third == AT2::CROSSING ||
		  vcross.third == AT2::INTERIOR );

  // the following condition becomes true only if intersecting
  // segments are found
  if ( vcross.first ) {
    if ( t.is_segment() ) {
      if ( vcross.third == AT2::CROSSING ) {
	Intersections_tag itag;
	return insert_intersecting_segment(ss, t, vcross.second, itag);
      } else if ( vcross.third == AT2::INTERIOR ) {
	Storage_site_2 ssvv = vcross.second->storage_site();
	Intersections_tag itag;
	Storage_site_2 ss1 = split_storage_site(ss, ssvv, 0, itag);
	Storage_site_2 ss2 = split_storage_site(ss, ssvv, 1, itag);
	insert_segment_interior(ss1.site(), ss1, vcross.second);
	return insert_segment_interior(ss2.site(), ss2, vcross.second);
      } else {
	// this should never be reached; the only possible values for
	// vcross.third are CROSSING, INTERIOR and DISJOINT
	CGAL_assertion( false );
      }
    }
  }

  // no intersecting segment has been found; we insert the segment as
  // usual...
  Vertex_handle v = create_vertex(ss);

  retriangulate_conflict_region(v, l, fm);

  return v;
}


//--------------------------------------------------------------------
// insertion of an intersecting segment
//--------------------------------------------------------------------
template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     Tag_false)
{
  static int i = 0;
  if ( i == 0 ) {
    i = 1;
    print_error_message();
  }
  return Vertex_handle();
}

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     Tag_true tag)
{
  CGAL_precondition( t.is_segment() && v->is_segment() );

  const Storage_site_2& ssitev = v->storage_site();
  Site_2 sitev = ssitev.site();

  if ( same_segments(t, sitev) ) {
    return v;
  }

  Vertex_triple vt = insert_point_on_segment(ss, t, v, tag);

  Vertex_handle vsx = vt.first;
  
  Storage_site_2 ss3, ss4;
  Site_2 s3, s4;
  if ( t.is_exact(0) ) {
    ss3 = create_storage_site(ss, ssitev, true);
  } else {
    ss3 = create_storage_site_type1(ss, ss, ssitev);
  }
  s3 = ss3.site();

  if ( t.is_exact(1) ) {
    ss4 = create_storage_site(ss, ssitev, false);
  } else {
    ss4 = create_storage_site_type2(ss, ssitev, ss);
  }
  s4 = ss4.site();

  insert_segment_interior(s3, ss3, vsx);
  insert_segment_interior(s4, ss4, vsx);
  return vsx;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// helper methods for insertion (find conflict region)
//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
initialize_conflict_region(const Face_handle& f, List& l)
{


  l.clear();
  for (int i = 0; i < 3; i++) {
    l.push_back(sym_edge(f, i));
  }
}


template<class Gt, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
expand_conflict_region(const Face_handle& f, const Site_2& t,
		       const Storage_site_2& ss,
		       List& l, Face_map& fm,
		       std::map<Face_handle,Sign>& sign_map,
		       Triple<bool,Vertex_handle,Arrangement_type>&
		       vcross)
{
  if ( fm.find(f) != fm.end() ) { return; }

  // this is done to stop the recursion when intersecting segments
  // are found
  if ( vcross.first ) { return; }

  // setting fm[f] to true means that the face has been reached and
  // that the face is available for recycling. If we do not want the
  // face to be available for recycling we must set this flag to
  // false.
  fm[f] = true;

  //  CGAL_assertion( fm.find(f) != fm.end() );

  for (int i = 0; i < 3; i++) {
    Face_handle n = f->neighbor(i);

    bool face_registered = (fm.find(n) != fm.end());

    if ( !face_registered ) {
      for (int j = 0; j < 3; j++) {
	Vertex_handle vf = n->vertex(j);

	if ( is_infinite(vf) ) { continue; }

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
	    fm.clear();
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
	    fm.clear();
	    return;
	  }
	}
      }
    }

    Sign s = incircle(n, t);

    sign_map[n] = s;

    Sign s_f = sign_map[f];

    if ( s == POSITIVE ) { continue; }
    if ( s != s_f ) { continue; }

    bool interior_in_conflict = edge_interior(f, i, t, s);

    if ( !interior_in_conflict ) { continue; }

    if ( face_registered ) { continue; }

    Edge e = sym_edge(f, i);

    CGAL_assertion( l.is_in_list(e) );
    int j = f->mirror_index(i);
    Edge e_before = sym_edge(n, ccw(j));
    Edge e_after = sym_edge(n, cw(j));
    if ( !l.is_in_list(e_before) ) {
      l.insert_before(e, e_before);
    }
    if ( !l.is_in_list(e_after) ) {
      l.insert_after(e, e_after);
    }
    l.remove(e);

    expand_conflict_region(n, t, ss, l, fm, sign_map, vcross);

    // this is done to stop the recursion when intersecting segments
    // are found
    //    if ( fm.size() == 0 && l.size() == 0 ) { return; }
    if ( vcross.first ) { return; }
  } // for-loop
}


//--------------------------------------------------------------------
// retriangulate conflict region
//--------------------------------------------------------------------

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
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


template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_list
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
add_bogus_vertices(List& l)
{
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
}

template<class Gt, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
remove_bogus_vertices(Vertex_list& vl)
{
  while ( vl.size() > 0 ) {
    Vertex_handle v = vl.front();
    vl.pop_front();
    remove_degree_2(v);
  }
}


template<class Gt, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
retriangulate_conflict_region(Vertex_handle v, List& l, 
			      Face_map& fm)
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
    CGAL_assertion( fm.find(f) == fm.end() );
    f->vertex(ccw(k))->set_face(f);
    f->vertex( cw(k))->set_face(f);
    eit = l.next(eit);
  } while ( eit != e_start );

  // 3. copy the edge list to a vector of edges and clear the edge list
  std::vector<Edge> ve;

  Edge efront = l.front();
  Edge e = efront;
  do {
    ve.push_back(e);
    e = l.next(e);
  } while ( e != efront );

  l.clear();

  // 4. retriangulate the hole
  this->_tds.star_hole(v, ve.begin(), ve.end());

  // 5. remove the bogus vertices
  remove_bogus_vertices(dummy_vertices);

  // 6. remove the unused faces
  typename Face_map::iterator it;
  for (it = fm.begin(); it != fm.end(); ++it) {
    Face_handle fh = (*it).first;
    this->_tds.delete_face(fh);
  }

  fm.clear();

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
template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
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

template<class Gt, class DS, class LTag>
Sign
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
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


template<class Gt, class DS, class LTag>
Sign
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
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
template<class Gt, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
finite_edge_interior(const Face_handle& f, int i, const Site_2& q,
		     Sign sgn, int) const
{
  if ( !is_infinite( f->mirror_vertex(i) ) ) {
    CGAL_precondition( is_infinite(f->vertex(i)) );

    Face_handle g = f->neighbor(i);
    int j = f->mirror_index(i);

    return finite_edge_interior(g, j, q, sgn, 0 /* degenerate */);
  }

  CGAL_precondition( is_infinite( f->mirror_vertex(i) ) );

  Site_2 t1 = f->vertex( ccw(i) )->site();
  Site_2 t2 = f->vertex(  cw(i) )->site();

  if ( is_infinite(f->vertex(i)) ) {
    return finite_edge_interior(t1, t2, q, sgn);
  }

  Site_2 t3 = f->vertex(i)->site();
  return finite_edge_interior(t1, t2, t3, q, sgn);
}

template<class Gt, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
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

  Site_2 t1 = v1->site();
  Site_2 t2 = v2->site();
  Site_2 q = v->site();

  if ( is_infinite(v3) ) {
    return finite_edge_interior(t1, t2, q, sgn);
  }

  Site_2 t3 = v3->site();
  return finite_edge_interior(t1, t2, t3, q, sgn);
}

template<class Gt, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
infinite_edge_interior(const Face_handle& f, int i,
		       const Site_2& q, Sign sgn) const
{
  if ( !is_infinite( f->vertex(ccw(i)) ) ) {
    CGAL_precondition( is_infinite( f->vertex(cw(i)) ) );
    Face_handle g = f->neighbor(i);
    int j = f->mirror_index(i);

    return infinite_edge_interior(g, j, q, sgn);
  }

  CGAL_precondition( is_infinite( f->vertex(ccw(i)) ) );

  Site_2 t2 = f->vertex(  cw(i) )->site();
  Site_2 t3 = f->vertex(     i  )->site();
  Site_2 t4 = f->mirror_vertex(i)->site();

  return infinite_edge_interior(t2, t3, t4, q, sgn);
}


template<class Gt, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
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

  Site_2 t2 = v2->site();
  Site_2 t3 = v3->site();
  Site_2 t4 = v4->site();
  Site_2 q = v->site();

  return infinite_edge_interior(t2, t3, t4, q, sgn);
}




template<class Gt, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
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


template<class Gt, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
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


template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Arrangement_type
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
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
  if ( res == res == AT2::TOUCH_11 || res == AT2::TOUCH_12 ||
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
template<class Gt, class DS, class LTag>
Object
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
primal(const Edge e) const
{
  typedef typename Gt::Line_2   Line_2;
  typedef typename Gt::Ray_2    Ray_2;

  CGAL_precondition( !is_infinite(e) );

  if ( this->dimension() == 1 ) {
    Site_2 p = (e.first)->vertex(cw(e.second))->site();
    Site_2 q = (e.first)->vertex(ccw(e.second))->site();

    Line_2 l = construct_svd_bisector_2_object()(p,q);
    return make_object(l);
  }

  // dimension == 2
  // none of the two adjacent faces is infinite
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    Site_2 p = (e.first)->vertex( ccw(e.second) )->site();
    Site_2 q = (e.first)->vertex(  cw(e.second) )->site();
    Site_2 r = (e.first)->vertex(     e.second  )->site();
    Site_2 s = (e.first)->mirror_vertex(e.second)->site();
    return construct_svd_bisector_segment_2_object()(p,q,r,s);
  }

  // both of the adjacent faces are infinite
  if ( is_infinite(e.first) &&
       is_infinite(e.first->neighbor(e.second)) )  {
    Site_2 p = (e.first)->vertex(cw(e.second))->site();
    Site_2 q = (e.first)->vertex(ccw(e.second))->site();
    Line_2 l = construct_svd_bisector_2_object()(p,q);
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

  CGAL_assertion(  is_infinite( e.first->vertex(e.second) ) ||
		   is_infinite( e.first->mirror_vertex(e.second) )  );

  Edge ee = e;
  if ( is_infinite( e.first->vertex(e.second) )  ) {
    ee = sym_edge(e);
  }
  Site_2 p = ee.first->vertex( ccw(ee.second) )->site();
  Site_2 q = ee.first->vertex(  cw(ee.second) )->site();
  Site_2 r = ee.first->vertex(     ee.second  )->site();

  Ray_2 ray = construct_svd_bisector_ray_2_object()(p,q,r);
  return make_object(ray);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// validity test method
//--------------------------------------------------------------------
//--------------------------------------------------------------------
template<class Gt, class DS, class LTag>
bool Segment_Voronoi_diagram_2<Gt,DS,LTag>::
is_valid(bool verbose, int level) const
{
  if (level < 0) { return true; }

  if (number_of_vertices() <= 1) { return true; }

  // level 0 test: check the TDS
  bool result = data_structure().is_valid(verbose, level);

  if ( result && verbose ) {
    std::cerr << "SVDDS is ok... " << std::flush;
  }

  if (level == 0) { return result; }

  // level 1 test: do the incircle tests
  if (number_of_vertices() < 3)  { return true; }

  for (All_edges_iterator eit = all_edges_begin();
       eit != all_edges_end(); ++eit) {
    Edge e = *eit;
    Face_handle f = e.first;

    Vertex_handle v = f->mirror_vertex(e.second);

    if ( f->vertex(e.second) == v ) { continue; }
    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->site()) != NEGATIVE );
    }
    Edge sym_e = sym_edge(e);
    f = sym_e.first;
    v = f->mirror_vertex(sym_e.second);

    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->site()) != NEGATIVE );
    }
  }

  if ( result && verbose ) {
    std::cerr << "Segment Voronoi diagram is ok..." << std::flush;
  }
  if ( !result && verbose ) {
    std::cerr << "Segment Voronoi diagram is NOT valid..." << std::flush;
  }

  return result;
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------
// misc
//--------------------------------------------------------------------
//--------------------------------------------------------------------


template<class Gt, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
print_error_message() const
{
  std::cerr << std::endl;
  std::cerr << "WARNING:" << std::endl;
  std::cerr << "A segment-segment intersection was found."
	    << std::endl;
  std::cerr << "The segment Voronoi diagram class is not configured"
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

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Storage_site_2
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
copy_storage_site(const Storage_site_2& ss_other, Handle_map& hm,
		  const Tag_false&)
{
  if ( ss_other.is_segment() ) {
    Point_handle p0 = hm[ ss_other.point_handle(0) ];
    Point_handle p1 = hm[ ss_other.point_handle(1) ];

    return Storage_site_2(p0, p1);
  } else {
    Point_handle p0 = hm[ ss_other.point_handle(0) ];

    return Storage_site_2(p0);
  }
}

template<class Gt, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,DS,LTag>::Storage_site_2
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
copy_storage_site(const Storage_site_2& ss_other, Handle_map& hm,
		  const Tag_true&)
{
  if ( ss_other.is_segment() ) {
    if ( ss_other.is_exact() ) {
      Point_handle p0 = hm[ ss_other.point_handle(0) ];
      Point_handle p1 = hm[ ss_other.point_handle(1) ];

      return Storage_site_2(p0, p1);
    } else if ( ss_other.is_exact(0) ) {
      Point_handle p0 = hm[ ss_other.point_handle(0) ];
      Point_handle p1 = hm[ ss_other.point_handle(1) ];
      Point_handle p4 = hm[ ss_other.point_handle(4) ];
      Point_handle p5 = hm[ ss_other.point_handle(5) ];

      return Storage_site_2(p0, p1, p4, p5, true);
    } else if ( ss_other.is_exact(1) ) {
      Point_handle p0 = hm[ ss_other.point_handle(0) ];
      Point_handle p1 = hm[ ss_other.point_handle(1) ];
      Point_handle p2 = hm[ ss_other.point_handle(2) ];
      Point_handle p3 = hm[ ss_other.point_handle(3) ];

      return Storage_site_2(p0, p1, p2, p3, false);
    } else {
      Point_handle p0 = hm[ ss_other.point_handle(0) ];
      Point_handle p1 = hm[ ss_other.point_handle(1) ];
      Point_handle p2 = hm[ ss_other.point_handle(2) ];
      Point_handle p3 = hm[ ss_other.point_handle(3) ];
      Point_handle p4 = hm[ ss_other.point_handle(4) ];
      Point_handle p5 = hm[ ss_other.point_handle(5) ];

      return Storage_site_2(p0, p1, p2, p3, p4, p5);
    }
  } else {
    if ( ss_other.is_exact() ) {
      Point_handle p0 = hm[ ss_other.point_handle(0) ];
      return Storage_site_2(p0);
    } else {
      Point_handle p2 = hm[ ss_other.point_handle(2) ];
      Point_handle p3 = hm[ ss_other.point_handle(3) ];
      Point_handle p4 = hm[ ss_other.point_handle(4) ];
      Point_handle p5 = hm[ ss_other.point_handle(5) ];
      return Storage_site_2(p2, p3, p4, p5);
    }
  }
}

template<class Gt, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,DS,LTag>::
copy(Segment_Voronoi_diagram_2& other)
{
  // first copy the point container and input point container
  pc_ = other.pc_;
  isc_ = other.isc_;

  // then copy the diagram
  DG::operator=(other);

  // now we have to update the sotrage sites in each vertex of the
  // diagram and also update the 

  // first create a map between the old point handles and the new ones
  Handle_map hm;

  Point_handle it_other = other.pc_.begin();
  Point_handle it_this = pc_.begin();
  for (; it_other != other.pc_.end(); ++it_other, ++it_this) {
    hm[it_other] = it_this;
  }

  // then update the storage sites for each vertex
  Intersections_tag itag;

  Finite_vertices_iterator vit_other = other.finite_vertices_begin();
  Finite_vertices_iterator vit_this = finite_vertices_begin();
  for (; vit_other != other.finite_vertices_end(); vit_other++,
	 vit_this++) {
    Storage_site_2 ss_other = vit_other->storage_site();

#ifndef NDEBUG
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

}

CGAL_END_NAMESPACE

// EOF
