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


//--------------------------------------------------------------------
// test method
//--------------------------------------------------------------------
template<class Gt, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
is_valid(bool verbose, int level) const
{
  if (level < 0) { return true; }

  if (number_of_vertices() <= 1) { return true; }

  // level 0 test: check the TDS
  bool result = ds().is_valid(verbose, level);
  //  bool result(true);

  //  CGAL_assertion( result );

  if ( result && verbose ) {
    std::cout << "SVDDS is ok... " << std::flush;
  }

  if (level == 0) { return result; }

  // level 1 test: do the incircle tests

  if (number_of_vertices() < 3)  return true;

  //  CGAL_assertion(result);

  for (All_edges_iterator eit = all_edges_begin();
       eit != all_edges_end(); ++eit) {
    Edge e = *eit;
    Face_handle f = e.first;

    Vertex_handle v = f->mirror_vertex(e.second);

    if ( f->vertex(e.second) == v ) { continue; }
    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->site()) != NEGATIVE );
      //    CGAL_assertion(result);
    }
    Edge sym_e = sym_edge(e);
    f = sym_e.first;
    v = f->mirror_vertex(sym_e.second);

    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->site()) != NEGATIVE );
      //    CGAL_assertion(result);
    }
  }

  if ( result && verbose ) {
    std::cerr << "Segment Voronoi diagram is ok..." << std::flush;
  }
  if ( !result && verbose ) {
    std::cerr << "Segment Voronoi diagram is NOT valid..." << std::flush;
  }

  //  CGAL_assertion(result);
  return result;
}



//--------------------------------------------------------------------
// embedding and visualization methods and constructions for primal
// and dual
//--------------------------------------------------------------------

// circumcenter
template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Point_2
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
circumcenter(const Face_handle& f) const
{
  CGAL_precondition( this->dimension()==2 || !is_infinite(f) );
  return circumcenter(f->vertex(0)->site(),
		      f->vertex(1)->site(),
		      f->vertex(2)->site());
}

template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Point_2
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
circumcenter(const Site_2& t0, const Site_2& t1, 
	     const Site_2& t2) const
{
  return
    geom_traits().construct_svd_vertex_2_object()(t0, t1, t2);
}

// circumcircle
template<class Gt, class PC, class DS, class LTag>
inline
typename Gt::Circle_2
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
circumcircle(const Face_handle& f) const
{
  CGAL_precondition( this->dimension()==2 || !is_infinite(f) );
  return circumcircle(f->vertex(0)->site(),
		      f->vertex(1)->site(),
		      f->vertex(2)->site());
}

template<class Gt, class PC, class DS, class LTag>
inline
typename Gt::Circle_2
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
circumcircle(const Site_2& t0, const Site_2& t1, const Site_2& t2) const
{
  return Construct_svd_circle_2()(t0, t1, t2);
}


template<class Gt, class PC, class DS, class LTag>
inline
typename Gt::Line_2
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
circumcircle(const Point_2& p0, const Point_2& p1) const
{
  return
    geom_traits().construct_line_2_object()(p0, p1);
}


// primal
template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Point_2
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
primal(const Face_handle& f) const
{
  return circumcenter(f);
}


template<class Gt, class PC, class DS, class LTag>
inline
Object
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
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
// combinatorial operations
//--------------------------------------------------------------------
template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Edge
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
flip(Face_handle& f, int i)
{
  CGAL_precondition ( f != Face_handle() );
  CGAL_precondition (i == 0 || i == 1 || i == 2);
  CGAL_precondition( this->dimension()==2 ); 

  CGAL_precondition( f->vertex(i) != f->mirror_vertex(i) );

  this->_tds.flip(f, i);

  return Edge(f, ccw(i));
}

template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Edge
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
flip(Edge e)
{
  return flip(e.first, e.second);
}

/*
template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_in_face(Face_handle& f, const Weighted_point& p)
{
  Vertex_handle v = ds().insert_in_face( f );

  v->set_point(p);
  return v;
}
*/

template<class Gt, class PC, class DS, class LTag>
inline
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
is_degree_2(const Vertex_handle& v) const
{
  Face_circulator fc = v->incident_faces();
  Face_circulator fc1 = fc;
  ++(++fc1);
  return ( fc == fc1 );
}

template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_degree_2(Edge e)
{
  return this->_tds.insert_degree_2(e.first,e.second);
}

template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_degree_2(Edge e,	const Storage_site_2& ss)
{
  Vertex_handle v = insert_degree_2(e);
  v->set_site(ss);
  return v;
}

template<class Gt, class PC, class DS, class LTag>
inline
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_degree_2(Vertex_handle v)
{
  CGAL_precondition( is_degree_2(v) );

  this->_tds.remove_degree_2(v);
}


#if 0
template<class Gt, class PC, class DS, class LTag>
inline
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_degree_3(Vertex_handle v)
{
  remove_degree_3(v, NULL);
}


template<class Gt, class PC, class DS, class LTag>
inline
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_degree_3(Vertex_handle v, Face* f)
{
  CGAL_precondition( v->degree() == 3 );
  this->_tds.remove_degree_3(v, f);
}
#endif

//--------------------------------------------------------------------
// insertion of site
//--------------------------------------------------------------------

template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_first(const Point_2& p)
{
  CGAL_precondition( number_of_vertices() == 0 );

  Storage_site_2 ss = create_storage_site(p);

  //  return create_vertex_dim_up(ss);
  Vertex_handle v = this->_tds.insert_second();
  v->set_site(ss);
  return v;
}

template<class Gt, class PC, class DS, class LTag>
inline
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_second(const Point_2& p)
{
  CGAL_precondition( number_of_vertices() == 1 );
  // p0 is actually a point
  Site_2 p0 = finite_vertices_begin()->site();
  // MK: change the equality test between points by the functor in
  // geometric traits
  if ( are_same_points(Site_2(p),p0) ) {
    return Vertex_handle(finite_vertices_begin());
  }

  Storage_site_2 ss = create_storage_site(p);
  return create_vertex_dim_up(ss);
}

template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_third(const Point_2& p)
{
  CGAL_precondition( number_of_vertices() == 2 );

  Site_2 t(p);

  // p0 and p1 are actually points
  Vertex_handle v0 = finite_vertices_begin();
  Vertex_handle v1 = ++finite_vertices_begin();
  Site_2 t0 = v0->site();
  Site_2 t1 = v1->site();

  // MK::ERROR: change the equality test between points by the functor
  // in geometric traits
  if ( are_same_points(t, t0) ) { return v0; }
  if ( are_same_points(t, t1) ) { return v1; }

  Storage_site_2 ss = create_storage_site(p);
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

template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
//insert_third(const Point_2& p0, const Point_2 & p1)
insert_third(Vertex_handle v0, Vertex_handle v1)
{
  CGAL_precondition( number_of_vertices() == 2 );

  //  this can only be the case if the first site is a segment
  CGAL_precondition( ds().dimension() == 1 );

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


template<class Gt, class PC, class DS, class LTag>
inline
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
do_intersect(const Site_2& t, Vertex_handle v) const
{
  if ( is_infinite(v) ) { return false; }
  // add here the cases where t is a segment and intersects a point
  // and t is a point and lies in a segment

  if ( !t.is_segment() || !v->is_segment() ) { return false; }

  if ( !intersection_flag ) {
    return same_segments(t, v->site());
  }

  if ( do_intersect(t, v->site()) ) {
    //	print_error_message();
    return true;
  }
  return false;
}

template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_point(const Point_2& p, Vertex_handle vnear)
{
  if ( number_of_vertices() == 0 ) {
    return insert_first(p);
  } else if ( number_of_vertices() == 1 ) {
    return insert_second(p);
  } else if ( number_of_vertices() == 2 ) {
    return insert_third(p);
  }

  // first find the nearest neighbor
  Site_2 t(p);

  Vertex_handle  vnearest = nearest_neighbor( p, vnear );

  CGAL_assertion( vnearest != Vertex_handle() );

  // find the first conflict
  // check if it is already inserted
  if ( vnearest->is_point() && are_same_points(t, vnearest->site()) ) {
    return vnearest;
  }

  // MK: add here code that checks if the inserted segment has already
  // been inserted

  // find the first conflict

  // first look if there are intersections...
  Vertex_circulator vc = vnearest->incident_vertices();
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    if ( do_intersect(t, vv) ) {
      // ADD HERE CODE THAT DOES THE APPROPRIATE INSERTION
      return Vertex_handle();
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
      }

      if ( interior_in_conflict ) { break; }
      ++ec;
    } while ( ec != ec_start );

    sign_map.clear();

    CGAL_assertion( interior_in_conflict );

    Storage_site_2 ss = create_storage_site(p);
    return insert_degree_2(e, ss);
  }


  // we are in conflict with a Voronoi vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;
  Face_map fm;

  std::pair<bool, Vertex_handle> vcross(false, Vertex_handle());

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.

  Storage_site_2 ss = create_storage_site(p);

  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, t, ss, l, fm, sign_map, vcross, NULL);

  CGAL_assertion( !vcross.first );

  Vertex_handle v = create_vertex(ss);

  retriangulate_conflict_region(v, l, fm);

  return v;
}


template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_point(const Storage_site_2& ss, const Site_2& t,
	     Vertex_handle vnear)
{
  CGAL_assertion( false );

  CGAL_precondition( t.is_point() );
  CGAL_precondition( !t.is_exact() );
  CGAL_assertion( number_of_vertices() > 2 );

#if 0
  if ( number_of_vertices() == 0 ) {
    return insert_first(t.point());
  } else if ( number_of_vertices() == 1 ) {
    return insert_second(t.point());
  } else if ( number_of_vertices() == 2 ) {
    return insert_third(t.point());
  }
#endif

  // first find the nearest neighbor
  Vertex_handle  vnearest = nearest_neighbor( t.point(), vnear );

  CGAL_assertion( vnearest != Vertex_handle() );

  // find the first conflict
  // check if it is already inserted
  if ( vnearest->is_point() && are_same_points(t, vnearest->site()) ) {
    return vnearest;
  }

  // MK: add here code that checks if the inserted segment has already
  // been inserted

  // find the first conflict

  // first look if there are intersections...
  Vertex_circulator vc = vnearest->incident_vertices();
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    if ( do_intersect(t, vv) ) {
      // ADD HERE CODE THAT DOES THE APPROPRIATE INSERTION
      return Vertex_handle();
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

  std::pair<bool, Vertex_handle> vcross(false, Vertex_handle());

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.
  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, t, l, fm, sign_map, vcross, NULL);

  CGAL_assertion( !vcross.first );

  Vertex_handle v = create_vertex(ss);

  retriangulate_conflict_region(v, l, fm);

  return v;
}


template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_segment(const Site_2& t, Vertex_handle vnear,
	       bool insert_endpoints)
{
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( t.is_exact() );
  CGAL_precondition( insert_endpoints );

  if ( is_degenerate_segment(t) ) {
    return insert_point(t.source(), vnear);
  }

  if ( number_of_vertices() == 0 ) {
    Vertex_handle v0 = insert_first( t.source() );
    Vertex_handle v1 = insert_second( t.target() );
    Segment_2 s = t.segment();
    return insert_third(v0, v1);
  } else {
    Vertex_handle v0, v1;
    if ( number_of_vertices() == 1 ) {
      v0 = insert_second( t.source() );
      v1 = insert_third( t.target() );
    } else if ( number_of_vertices() == 2 ) {
      v0 = insert_third( t.source() );
      v1 = insert_point( t.target(), v0 );
    } else {
      v0 = insert_point( t.source(), vnear );
      v1 = insert_point( t.target(), v0 );
    }

    Storage_site_2 ss = create_storage_site(v0, v1);
    return insert_segment2( t, ss, v0, false );
    //    insert_point( t.source_site(), Vertex_handle() );
    //    insert_point( t.target_site(), Vertex_handle() );
  }
}


template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_segment2(const Site_2& t, const Storage_site_2& ss,
		Vertex_handle vnearest, bool insert_endpoints)
{
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( !insert_endpoints );
  CGAL_precondition( number_of_vertices() >= 2 );

#if 0
  if ( is_degenerate_segment(t) ) {
    return insert_point(t.source_site(), vnear);
  }

  if ( number_of_vertices() <= 1 ) {
    //    insert_point( t.source_site(), Vertex_handle() );
    //    insert_point( t.target_site(), Vertex_handle() );
    Vertex_handle v0 = insert_point( t.source(), Vertex_handle() );
    insert_point( t.target(), Vertex_handle() );
    return insert_segment( t, v0, false );
  } else if ( number_of_vertices() == 2 ) {
    Segment_2 s = t.segment();
    return insert_third(s.source(), s.target());
  }

  // first find the nearest neighbor
  Vertex_handle vnearest;
  if ( insert_endpoints ) {
    // if we insert a segment, insert the endpoints first
    vnearest = insert_point( t.source_site(), vnear );
    insert_point( t.target_site(), vnearest );
  } else {
    vnearest = vnear;
  }

#endif
  CGAL_assertion( vnearest != Vertex_handle() );
  // MK: add here code that checks if the inserted segment has already
  // been inserted; MAYBE THIS IS NOT NEEDED; I ALREADY DO IT IN
  // do_intersect

  // find the first conflict

  // first look if there are intersections...
  Vertex_circulator vc = vnearest->incident_vertices();
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    if ( same_segments(t, vv) ) {
      return vv;
    }
    if ( do_intersect(t, vv) ) {
      if ( t.is_segment() ) {
	return insert_intersecting_segment(ss, t, vv,
					   Intersections_tag());
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

#if 0
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
      }

      if ( interior_in_conflict ) { break; }
      ++ec;
    } while ( ec != ec_start );

    sign_map.clear();

    CGAL_assertion( interior_in_conflict );

    return insert_degree_2(e, t);
  }
#endif

  // we are in conflict with a Voronoi vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;
  Face_map fm;

  std::pair<bool, Vertex_handle> vcross(false, Vertex_handle());

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.
  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, t, ss, l, fm, sign_map, vcross, NULL);

  // the following condition becomes true only if intersecting
  // segments are found
  if ( vcross.first ) {
    if ( t.is_segment() ) {
      return insert_intersecting_segment(ss, t, vcross.second,
					 Intersections_tag());
      //      return vcross.second;
    }
  }

  Vertex_handle v = create_vertex(ss);

  retriangulate_conflict_region(v, l, fm);

  return v;
}


//--------------------------------------------------------------------
// insertion of intersecting site
//--------------------------------------------------------------------
template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     Tag_false)
{
  static int i = 0;
  if ( i == 0 ) {
    i = 1;
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
  return Vertex_handle();
}

template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     Tag_true)
{
  CGAL_precondition( t.is_segment() && v->is_segment() );

  if ( same_segments(t, v->site()) ) {
    return v;
  }

  //  Storage_site_2 ssx(ss.supporting_segment_handle(),
  //		     v->storage_site().supporting_segment_handle());
  Storage_site_2 ssitev = v->storage_site();
  Storage_site_2 ssx( ss.point_handle(0), ss.point_handle(1),
		      ssitev.point_handle(0), ssitev.point_handle(1) );

  Site_2 sitev(v->site());
  Site_2 sx(t.point(0), t.point(1), sitev.point(0), sitev.point(1));

  Face_circulator fc1 = incident_faces(v);
  Face_circulator fc2 = fc1; ++fc2;
  Face_circulator fc_start = fc1;
  Face_handle f1, f2;
  bool found_f1 = false, found_f2 = false;
  Site_2 sitev_supp(sitev.point(0), sitev.point(1));
  do {
    Face_handle ff1(fc1), ff2(fc2);
    Oriented_side os1 = oriented_side(fc1->vertex(0)->site(),
				      fc1->vertex(1)->site(),
				      fc1->vertex(2)->site(),
				      sitev_supp, sx);
    Oriented_side os2 = oriented_side(fc2->vertex(0)->site(),
				      fc2->vertex(1)->site(),
				      fc2->vertex(2)->site(),
				      sitev_supp, sx);
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

  CGAL_assertion( f1 != f2 );

  Quadruple<Vertex_handle, Vertex_handle, Face_handle, Face_handle>
    qq = this->_tds.split_vertex(v, f1, f2);

  // now I need to update the sites for vertices v1 and v2
  Vertex_handle v1 = qq.first;
  Storage_site_2 ssv1;
  Site_2 sv1;
  if ( sitev.is_exact(0) ) {
    sv1.set_segment(sitev.point(0), sitev.point(1),
		    t.point(0), t.point(1), true);
    //    ssv1.set_segment(ssitev.supporting_segment_handle(),
    //		     ss.supporting_segment_handle(), true);
    ssv1.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ss.point_handle(0), ss.point_handle(1), true);
  } else {
    sv1.set_segment(sitev.point(0), sitev.point(1),
		    sitev.point(2), sitev.point(3),
		    t.point(0), t.point(1));
    //    ssv1.set_segment(ssitev.supporting_segment_handle(),
    //		     ssitev.crossing_segment_handle(0),
    //		     ss.supporting_segment_handle());
    ssv1.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ssitev.point_handle(2), ssitev.point_handle(3),
		     ss.point_handle(0), ss.point_handle(1));
  }
  v1->set_site( ssv1 );

  Vertex_handle v2 = qq.second;
  Storage_site_2 ssv2;
  Site_2 sv2;
  if ( sitev.is_exact(1) ) {
    sv2.set_segment(sitev.point(0), sitev.point(1),
		    t.point(0), t.point(1), false);
    //    ssv2.set_segment(ssitev.supporting_segment_handle(),
    //		     ss.supporting_segment_handle(), false);
    ssv2.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ss.point_handle(0), ss.point_handle(1), false);
  } else {
    sv2.set_segment(sitev.point(0), sitev.point(1),
		    t.point(0), t.point(1),
		    sitev.point(4), sitev.point(5));
    //    ssv2.set_segment(ssitev.supporting_segment_handle(),
    //		     ss.supporting_segment_handle(),
    //		     ssitev.crossing_segment_handle(1));
    ssv2.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ss.point_handle(0), ss.point_handle(1),
		     ssitev.point_handle(4), ssitev.point_handle(5));
  }
  v2->set_site( ssv2 );

  Vertex_handle vsx =
    this->_tds.insert_in_edge(qq.third, cw(qq.third->index(v1)));

  vsx->set_site(ssx);

  Storage_site_2 ss3, ss4;
  Site_2 s3, s4;
  if ( t.is_exact(0) ) {
    s3.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), true);
    //    ss3.set_segment(ss.supporting_segment_handle(),
    //		    ssitev.supporting_segment_handle(), true);
    ss3.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ssitev.point_handle(1), ssitev.point_handle(1), true);
  } else {
    s3.set_segment(t.point(0), t.point(1),
		   t.point(2), t.point(3),
		   sitev.point(0), sitev.point(1));
    //    ss3.set_segment(ss.supporting_segment_handle(),
    //		    ss.crossing_segment_handle(0),
    //		    ssitev.supporting_segment_handle());
    ss3.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ss.point_handle(2), ss.point_handle(3),
		    ssitev.point_handle(0), ssitev.point_handle(1));
  }

  if ( t.is_exact(1) ) {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), false);
    //    ss4.set_segment(ss.supporting_segment_handle(),
    //		    ssitev.supporting_segment_handle(), false);
    ss4.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ssitev.point_handle(0), ssitev.point_handle(1), false);
  } else {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1),
		   t.point(4), t.point(5));
    //    ss4.set_segment(ss.supporting_segment_handle(),
    //		    ssitev.supporting_segment_handle(),
    //		    ss.crossing_segment_handle(1));
    ss4.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ssitev.point_handle(0), ssitev.point_handle(1),
		    ss.point_handle(4), ss.point_handle(5));
  }

  insert_segment2(s3, ss3, vsx, false);
  insert_segment2(s4, ss4, vsx, false);
  return vsx;
}

//--------------------------------------------------------------------
// find conflict region
//--------------------------------------------------------------------

template<class Gt, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
initialize_conflict_region(const Face_handle& f, List& l)
{


  l.clear();
  for (int i = 0; i < 3; i++) {
    l.push_back(sym_edge(f, i));
  }
}


template<class Gt, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
expand_conflict_region(const Face_handle& f, const Site_2& t,
		       const Storage_site_2& ss,
		       List& l, Face_map& fm,
		       std::map<Face_handle,Sign>& sign_map,
		       std::pair<bool,Vertex_handle>& vcross,
		       std::vector<Vh_triple*>* fe)
{
  if ( fm.find(f) != fm.end() ) { return; }

  //  Site_2 t = ss.site();

  // this is done to stop the recursion when intersecting segments
  // are found
  if ( vcross.first ) { return; }

  for (int i = 0; i < 3; i++) {
    Vertex_handle vf = f->vertex(i);
    if ( do_intersect(t, vf) ) {
      vcross.first = true;
      vcross.second = vf;
      l.clear();
      fm.clear();
      return;
    }
  }

  // setting fm[f] to true means that the face has been reached and
  // that the face is available for recycling. If we do not want the
  // face to be available for recycling we must set this flag to
  // false.
  fm[f] = true;

  //  CGAL_assertion( fm.find(f) != fm.end() );

  for (int i = 0; i < 3; i++) {
    Face_handle n = f->neighbor(i);

    Sign s = incircle(n, t);

    sign_map[n] = s;

    Sign s_f = sign_map[f];

    if ( s == POSITIVE ) { continue; }
    if ( s != s_f ) { continue; }

    bool interior_in_conflict = edge_interior(f, i, t, s);

    if ( !interior_in_conflict ) { continue; }

    if ( fm.find(n) != fm.end() ) {
#if 0
      Edge e = sym_edge(f, i);
      if ( l.is_in_list(e) ||
	   l.is_in_list(sym_edge(e)) ) {
	l.remove(e);
	l.remove(sym_edge(e));

	// we should have never reached this point...
	// this check is done mainly for debugging; in the final
	// version these if-statements should be removed.
	bool loop_in_conflict_region_found(false);
	CGAL_assertion( loop_in_conflict_region_found );
      }
#endif
      continue;
    }

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

    if ( fe != NULL ) {
      Vh_triple* vhq = new Vh_triple[1];

      (*vhq)[0] = NULL;
      (*vhq)[1] = n->vertex(     j  );
      (*vhq)[2] = n->vertex( ccw(j) );

      fe->push_back(vhq);
    }

    expand_conflict_region(n, t, ss, l, fm, sign_map, vcross, fe);

    // this is done to stop the recursion when intersecting segments
    // are found
    //    if ( fm.size() == 0 && l.size() == 0 ) { return; }
    if ( vcross.first ) { return; }
  } // for-loop
}


//--------------------------------------------------------------------
// retriangulate conflict region
//--------------------------------------------------------------------

template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
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


template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_list
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
add_bogus_vertices(List& l)
{
  Vertex_list vertex_list;

  static std::set<Edge> edge_list;

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

template<class Gt, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_bogus_vertices(Vertex_list& vl)
{
  while ( vl.size() > 0 ) {
    Vertex_handle v = vl.front();
    vl.pop_front();
    remove_degree_2(v);
  }
}

#if 0
template<class Gt, class PC, class DS, class LTag>
std::vector<typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Face*>
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
get_faces_for_recycling(Face_map& fm, unsigned int n_wanted)
{
  std::vector<Face*> vf;

  typename Face_map::iterator fmit;
  for (fmit = fm.begin(); fmit != fm.end(); ++fmit) {
    Face_handle f = (*fmit).first;
    if ( fm[f] == true ) { vf.push_back(f); }
  }

  while ( vf.size() < n_wanted ) {
    Face* fp = static_cast<Face*>(_tds.create_face());
    vf.push_back(fp);
  }

  while ( vf.size() > n_wanted ) {
    Face* fp = vf.back();
    vf.pop_back();
    _tds.delete_face(fp);
  }
  
  return vf;
}
#endif

template<class Gt, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
retriangulate_conflict_region(Vertex_handle v, List& l, 
			      Face_map& fm)
{
  //  Vertex_handle v = create_vertex(ss);
  //  Vertex_handle v = _tds.create_vertex();
  //  v->set_site(t);

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

  //  std::vector<Face*> vf = get_faces_for_recycling(fm, l.size());
  //  std::list<Face*> vf;

  // 3. copy the edge list to a vector of edges and clear the edge
  //    list
  std::vector<Edge> ve;

  Edge efront = l.front();
  Edge e = efront;
  do {
    ve.push_back(e);
    e = l.next(e);
  } while ( e != efront );

  l.clear();

  // 4. retriangulate the hole
  //  _tds.star_hole(v, ve.begin(), ve.end(), vf.begin(), vf.end());
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
  //  return v;
}



//--------------------------------------------------------------------
// point location
//--------------------------------------------------------------------
template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
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
// methods for the predicates
//----------------------------------------------------------------------


template<class Gt, class PC, class DS, class LTag>
inline bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
are_same_points(const Site_2& p, const Site_2& q) const
{
  return geom_traits().are_same_points_2_object()(p, q);
}

template<class Gt, class PC, class DS, class LTag>
inline bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
same_segments(const Site_2& t, Vertex_handle v) const
{
  if ( is_infinite(v) ) { return false; }
  if ( t.is_point() || v->site().is_point() ) { return false; }
  return same_segments(t, v->site());
}

template<class Gt, class PC, class DS, class LTag>
inline bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
same_segments(const Site_2& p, const Site_2& q) const
{
  CGAL_precondition( p.is_segment() && q.is_segment() );

  return
    (are_same_points(p.source_site(), q.source_site()) &&
     are_same_points(p.target_site(), q.target_site())) ||
    (are_same_points(p.source_site(), q.target_site()) &&
     are_same_points(p.target_site(), q.source_site()));
}

template<class Gt, class PC, class DS, class LTag>
inline Oriented_side
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
side_of_bisector(const Site_2 &t1, const Site_2 &t2, const Site_2 &q) const
{
  return geom_traits().oriented_side_of_bisector_2_object()(t1, t2, q);
}


template<class Gt, class PC, class DS, class LTag>
inline Sign
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
incircle(const Site_2 &t1, const Site_2 &t2,
	 const Site_2 &t3, const Site_2 &q) const
{
  return geom_traits().vertex_conflict_2_object()(t1, t2, t3, q);
}

template<class Gt, class PC, class DS, class LTag>
inline Sign
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
incircle(const Site_2 &t1, const Site_2 &t2,
	 const Site_2 &q) const
{
  return
    geom_traits().vertex_conflict_2_object()(t1, t2, q);
}


template<class Gt, class PC, class DS, class LTag>
Sign
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
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


template<class Gt, class PC, class DS, class LTag>
inline Sign
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
incircle(const Vertex_handle& v0, const Vertex_handle& v1,
	 const Vertex_handle& v) const
{
  CGAL_precondition( !is_infinite(v0) && !is_infinite(v1)
		     && !is_infinite(v) );

  return incircle( v0->site(), v1->site(), v->site());
}

template<class Gt, class PC, class DS, class LTag>
Sign
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
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


template<class Gt, class PC, class DS, class LTag>
inline bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
finite_edge_interior(const Site_2& t1, const Site_2& t2,
		     const Site_2& t3, const Site_2& t4,
		     const Site_2& q, Sign sgn) const
{
  return
    geom_traits().finite_edge_interior_conflict_2_object()(t1,t2,t3,t4,q,sgn);
}

template<class Gt, class PC, class DS, class LTag>
inline bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
finite_edge_interior(const Face_handle& f, int i,
		     const Site_2& q, Sign sgn) const
{
  CGAL_precondition( !is_infinite(f) &&
		     !is_infinite(f->neighbor(i)) );
  return finite_edge_interior( f->vertex( ccw(i) )->site(),
			       f->vertex(  cw(i) )->site(),
			       f->vertex(     i  )->site(),
			       f->mirror_vertex(i)->site(), q, sgn);
}

template<class Gt, class PC, class DS, class LTag>
inline bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
finite_edge_interior(const Vertex_handle& v1,
		     const Vertex_handle& v2,
		     const Vertex_handle& v3,
		     const Vertex_handle& v4,
		     const Vertex_handle& v, Sign sgn) const
{
  CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) &&
		     !is_infinite(v3) && !is_infinite(v4) &&
		     !is_infinite(v) );
  return finite_edge_interior( v1->site(), v2->site(),
			       v3->site(), v4->site(),
			       v->site(), sgn);
}

template<class Gt, class PC, class DS, class LTag>
inline bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
finite_edge_interior_degenerated(const Site_2& t1, const Site_2& t2,
				 const Site_2& t3, const Site_2& q,
				 Sign sgn) const
{
  return
    geom_traits().finite_edge_interior_conflict_2_object()(t1,t2,t3,q,sgn);
}

template<class Gt, class PC, class DS, class LTag>
inline bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
finite_edge_interior_degenerated(const Site_2& t1, const Site_2& t2,
				 const Site_2& q, Sign sgn) const
{
  return
    geom_traits().finite_edge_interior_conflict_2_object()(t1, t2, q, sgn);
}


template<class Gt, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
finite_edge_interior_degenerated(const Face_handle& f, int i,
				 const Site_2& q, Sign sgn) const
{
  if ( !is_infinite( f->mirror_vertex(i) ) ) {
    CGAL_precondition( is_infinite(f->vertex(i)) );

    Face_handle g = f->neighbor(i);
    int j = f->mirror_index(i);

    return finite_edge_interior_degenerated(g, j, q, sgn);
  }

  CGAL_precondition( is_infinite( f->mirror_vertex(i) ) );

  Site_2 t1 = f->vertex( ccw(i) )->site();
  Site_2 t2 = f->vertex(  cw(i) )->site();

  if ( is_infinite(f->vertex(i)) ) {
    return finite_edge_interior_degenerated(t1, t2, q, sgn);
  }

  Site_2 t3 = f->vertex(i)->site();
  return finite_edge_interior_degenerated(t1, t2, t3, q, sgn);
}

template<class Gt, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
finite_edge_interior_degenerated(const Vertex_handle& v1,
				 const Vertex_handle& v2,
				 const Vertex_handle& v3,
				 const Vertex_handle& v4,
				 const Vertex_handle& v, Sign sgn) const
{
  CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) && 
		     !is_infinite(v) );
  if ( !is_infinite( v4 ) ) {
    CGAL_precondition( is_infinite(v3) );

    return
      finite_edge_interior_degenerated(v2, v1, v4, v3, v, sgn);
  }

  CGAL_precondition( is_infinite( v4 ) );

  Site_2 t1 = v1->site();
  Site_2 t2 = v2->site();
  Site_2 q = v->site();

  if ( is_infinite(v3) ) {
    return finite_edge_interior_degenerated(t1, t2, q, sgn);
  }

  Site_2 t3 = v3->site();
  return finite_edge_interior_degenerated(t1, t2, t3, q, sgn);
}

template<class Gt, class PC, class DS, class LTag>
inline bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
infinite_edge_interior(const Site_2& t2, const Site_2& t3,
		       const Site_2& t4, const Site_2& q, Sign sgn) const
{
  return
    geom_traits().infinite_edge_interior_conflict_2_object()(t2,t3,t4,q,sgn);
}

template<class Gt, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
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


template<class Gt, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
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




template<class Gt, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
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
    result = finite_edge_interior_degenerated(v1, v2, v3, v4, v, sgn);
  } else {
    result = infinite_edge_interior(v1, v2, v3, v4, v, sgn);
  }

  return result;
}


template<class Gt, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
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
    result = finite_edge_interior_degenerated(f, i, q, sgn);
  } else {
    //    Edge e(f, i);
    if ( !is_infinite(f, i) ) {
      result = finite_edge_interior_degenerated(f, i, q, sgn);
    } else {
      result = infinite_edge_interior(f, i, q, sgn);
    }
  }

  return result;
}

/*
template<class Gt, class PC, class DS, class LTag>
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Conflict_type
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
finite_edge_conflict_type_degenerated(const Weighted_point& p1,
				      const Weighted_point& p2,
				      const Weighted_point& q) const
{
  Sign i1 = incircle(p1, p2, q);
  Sign i2 = incircle(p2, p1, q);

  if ( i1 == NEGATIVE && i2 == POSITIVE ) {
    return LEFT_VERTEX;
  } else if ( i1 == POSITIVE && i2 == NEGATIVE ) {
    return RIGHT_VERTEX;
  } else if ( i1 == POSITIVE && i2 == POSITIVE ) {
    bool b = finite_edge_interior_degenerated(p1, p2, q, false);
    return (b ? INTERIOR : NO_CONFLICT);
  } else if ( i1 == NEGATIVE && i2 == NEGATIVE ) {
    bool b = finite_edge_interior_degenerated(p1, p2, q, true);
    return (b ? ENTIRE_EDGE : BOTH_VERTICES);
  } else {
    // this should never be reached; the degenerated incircle never
    // returns ZERO
    bool not_ready_yet(false);
    assert( not_ready_yet );
  }

  // to satisfy compiler
  return NO_CONFLICT;
}
*/



//----------------------------------------------------------------------
// methods for disk removal
//----------------------------------------------------------------------
#if 0
template<class Gt, class PC, class DS, class LTag>
std::pair<
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle,
typename Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::Vertex_handle >
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
endpoint_vertices(Vertex_handle v) const
{
  CGAL_precondition( v->is_segment() );

  Vertex_circulator vc_start = incident_vertices(v);
  Vertex_circulator vc = vc_start;

  std::pair<Vertex_handle, Vertex_handle> vertices;
  do {
    Vertex_handle u(vc);
    if ( u->is_point() &&
	 are_same_points(u->site(), v->source_site()) ) {
      vertices.first = u;
    }

    if ( u->is_point() &&
	 are_same_points(u->site(), v->target_site()) ) {
      vertices.second = u;
    }

    ++vc;
  } while ( vc != vc_start );

  return vertices;
}

template<class Gt, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
is_endpoint_of_segment(Vertex_handle v) const
{
  if ( v->is_segment() ) { return false; }

  bool is_endpoint(false);
  Vertex_circulator vc_start = incident_vertices(v);
  Vertex_circulator vc = vc_start;

  do {
    Vertex_handle u(vc);
    if ( u->is_segment() &&
	 is_endpoint_of_segment(v->site(), u->site()) ) {
      is_endpoint = true;
      break;
    }
    ++vc;
  } while ( vc != vc_start );

  return is_endpoint;
}


template<class Gt, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_first(Vertex_handle v)
{
  Delaunay_graph::remove_first(v);
}

template<class Gt, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_second(Vertex_handle v)
{
  Delaunay_graph::remove_second(v);
}

template<class Gt, class PC, class DS, class LTag>
unsigned int
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_third(Vertex_handle v, bool remove_endpoints)
{
  if ( is_endpoint_of_segment(v) )  { return 0; }

  if ( v->is_segment() && remove_endpoints ) {
    // since we have only three sites and one is a segment, the other
    // must be its endpoints.
    clear();
    return 3;
  }

  if ( is_degree_2(v) ) {
    Face_handle fh(v->incident_faces());
    int i = fh->index(v);
    flip(fh, i);
  } else if ( v->degree() == 4 ) {
    Edge_circulator ec = v->incident_edges();
    for (int i = 0; i < 4; i++) {
      Edge e = *ec;
      Edge sym = sym_edge(e);
      if ( e.first->vertex(e.second) !=	sym.first->vertex(sym.second) ) {
	flip(e);
	break;
      }
      ++ec;
    }
  }

  this->_tds.remove_dim_down(v);

  return 1;
}


template<class Gt, class PC, class DS, class LTag>
unsigned int
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove(Vertex_handle v, bool remove_endpoints)
{
  CGAL_precondition( v != Vertex_handle() );
  CGAL_precondition( !is_infinite(v) );

  int num_removed(0);
  int n = number_of_vertices();
  if ( n == 1 ) {
    remove_first(v);
    num_removed = 1;
  } else if ( n == 2 ) {
    remove_second(v);
    num_removed = 1;
  } else if ( n == 3 ) {
    num_removed = remove_third(v, remove_endpoints);
  } else {
    int degree = v->degree();
    if ( degree == 2 ) {
      num_removed = remove_degree_2(v, remove_endpoints);
    } else if ( degree == 3 ) {
      num_removed = remove_degree_3(v, remove_endpoints);
    } else {
      num_removed = remove_degree_d(v, remove_endpoints);
    }
  }

  //  CGAL_assertion( is_valid(false, 2) );

  return num_removed;
}

template<class Gt, class PC, class DS, class LTag>
unsigned int
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_degree_2(Vertex_handle v, bool remove_endpoints)
{
  // segments are at least of degree 4
  CGAL_assertion( !v->is_segment() );

  if ( is_endpoint_of_segment(v) )  { return 0; }

  remove_degree_2(v);
  return 1;
}


template<class Gt, class PC, class DS, class LTag>
unsigned int
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_degree_3(Vertex_handle v, bool remove_endpoints)
{
  // segments are at least of degree 4
  CGAL_assertion( !v->is_segment() );

  if ( is_endpoint_of_segment(v) )  { return 0; }

  remove_degree_3(v);
  return 1;
}

template<class Gt, class PC, class DS, class LTag>
unsigned int
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_degree_d(Vertex_handle v, bool remove_endpoints)
{
  if ( v->is_segment() && remove_endpoints ) {
    std::pair<Vertex_handle, Vertex_handle> vertices
      = endpoint_vertices(v);

    unsigned int n_seg = remove(v, false);

    CGAL_assertion( n_seg == 1 );

    unsigned int n_src = remove( vertices.first, false );
    unsigned int n_trg = remove( vertices.second, false );
    return (n_seg + n_src + n_trg);
  }

  if ( is_endpoint_of_segment(v) )  { return 0; }

  minimize_degree(v);
  int deg = v->degree();
  if ( deg == 3 ) {
    return remove_degree_3(v, remove_endpoints);
  }
  if ( deg == 2 ) {
    return remove_degree_2(v, remove_endpoints);
  }

  Segment_Voronoi_diagram_2<Gt,PC,DS,LTag> svd_small;

  std::map<Vertex_handle,Vertex_handle> vmap;

  Vertex_circulator vc_start = v->incident_vertices();
  Vertex_circulator vc = v->incident_vertices();
  Vertex_handle vh_large, vh_small;
  do {
    vh_large = Vertex_handle(vc);
    if ( is_infinite(vh_large) ) {
      vh_small = svd_small.infinite_vertex();
      vmap[vh_small] = vh_large;
    } else { 
      if ( vc->is_point() ) {
	vh_small = svd_small.insert( vc->point() );
      } else {
	std::pair<Vertex_handle, Vertex_handle> vertices
	  = endpoint_vertices(vh_large);

	Vertex_handle v_src = 
	  svd_small.insert( vertices.first->point() );
	vmap[v_src] = vertices.first;

	Vertex_handle v_trg = 
	  svd_small.insert( vertices.second->point() );
	vmap[v_trg] = vertices.second;

	vh_small = svd_small.insert( vc->segment(), v_src );
      }

      if ( vh_small != Vertex_handle() ) {
	vmap[vh_small] = vh_large;
      }
    }
    ++vc;
  } while ( vc != vc_start );


  /*
    // need to check the following portion
  if ( ag_small.number_of_vertices() == 2 ) {
    CGAL_assertion( deg == 4 );
    Edge_circulator ec = v->incident_edges();
    for (int i = 0; i < 4; i++) {
      Edge e = *ec;
      Edge sym = sym_edge(e);
      if ( e.first->vertex(e.second) !=	sym.first->vertex(sym.second) ) {
	flip(e);
	break;
      }
      ++ec;
    }
    remove_degree_3(v);
    return;
  }
  */

  Vertex_handle vn;
  if ( v->is_segment() ) {
    vn = svd_small.nearest_neighbor( v->segment().source() );
  } else {
    vn = svd_small.nearest_neighbor( v->point() );
  }

  List l;
  Face_map fm;
  std::vector<Vh_triple*> flipped_edges;  

  svd_small.find_conflict_region_remove(v, vn, l, fm,
					&flipped_edges);

  l.clear();
  fm.clear();

  Edge_circulator ec;

  unsigned int num_fe = flipped_edges.size();
  for (unsigned int i = 0; i < num_fe; i++) {
    Vh_triple *vhq = flipped_edges[num_fe - i - 1];

    bool found(false);
    ec = v->incident_edges();
    Edge_circulator ec_start = ec;
    do {
      Edge e = *ec;
      if ( (e.first->vertex(  cw(e.second) ) == vmap[(*vhq)[1]] &&
	    e.first->vertex(     e.second  ) == vmap[(*vhq)[2]]) ||
	   (e.first->vertex( ccw(e.second) ) == vmap[(*vhq)[1]] &&
	    e.first->mirror_vertex(e.second) == vmap[(*vhq)[2]]) ) {
	flip(e);
	found = true;
	break;
      }
      ++ec;
    } while ( ec != ec_start );

    CGAL_assertion( found );
  }

  CGAL_precondition( v->degree() == 3 );

  this->_tds.remove_degree_3(v, NULL);

  for (unsigned int i = 0; i < num_fe; i++) {
    delete flipped_edges[i];
  }

  return 1;
}

/*
template<class Gt, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
remove_degree_d_vertex(Vertex_handle v)
{
  minimize_degree(v);
  int deg = v->degree();
  if ( deg == 3 ) {
    remove_degree_3(v);
    return;
  }
  if ( deg == 2 ) {
    remove_degree_2(v);
    return;
  }
  
  Segment_Voronoi_diagram_2<Gt,PC,DS,LTag> ag_small;

  std::map<Vertex_handle,Vertex_handle> vmap;

  Vertex_circulator vc_start = v->incident_vertices();
  Vertex_circulator vc = v->incident_vertices();
  Vertex_handle vh_large, vh_small;
  do {
    vh_large = Vertex_handle(vc);
    if ( is_infinite(vh_large) ) {
      vh_small = ag_small.infinite_vertex();
      vmap[vh_small] = vh_large;
    } else { 
      vh_small = ag_small.insert(vc->point());
      if ( vh_small != Vertex_handle() ) {
	vmap[vh_small] = vh_large;
      }
    }
    ++vc;
  } while ( vc != vc_start );

  if ( ag_small.number_of_vertices() == 2 ) {
    CGAL_assertion( deg == 4 );
    Edge_circulator ec = v->incident_edges();
    for (int i = 0; i < 4; i++) {
      Edge e = *ec;
      Edge sym = sym_edge(e);
      if ( e.first->vertex(e.second) !=	sym.first->vertex(sym.second) ) {
	flip(e);
	break;
      }
      ++ec;
    }
    remove_degree_3(v);
    return;
  }


  Vertex_handle vn = ag_small.nearest_neighbor(v->point());

  List l;
  Face_map fm;
  Vertex_map vm;
  std::vector<Vh_triple*> flipped_edges;  

  ag_small.find_conflict_region_remove(v, vn, l, fm, vm,
				       &flipped_edges);

  l.clear();
  fm.clear();
  vm.clear();

  Edge_circulator ec;

  unsigned int num_fe = flipped_edges.size();
  for (unsigned int i = 0; i < num_fe; i++) {
    Vh_triple *vhq = flipped_edges[num_fe - i - 1];

    bool found(false);
    ec = v->incident_edges();
    Edge_circulator ec_start = ec;
    do {
      Edge e = *ec;
      if ( (e.first->vertex(  cw(e.second) ) == vmap[(*vhq)[1]] &&
	    e.first->vertex(     e.second  ) == vmap[(*vhq)[2]]) ||
	   (e.first->vertex( ccw(e.second) ) == vmap[(*vhq)[1]] &&
	    e.first->mirror_vertex(e.second) == vmap[(*vhq)[2]]) ) {
	flip(e);
	found = true;
	break;
      }
      ++ec;
    } while ( ec != ec_start );

    CGAL_assertion( found );
  }
  CGAL_precondition( v->degree() == 3 );

  this->_tds.remove_degree_3(v, NULL);

  for (unsigned int i = 0; i < num_fe; i++) {
    delete flipped_edges[i];
  }
}
*/

template<class Gt, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
minimize_degree(Vertex_handle v)
{
  CGAL_precondition ( v->degree() > 3 );

  Face_circulator fc_start = v->incident_faces();
  Face_circulator fc = v->incident_faces();
  bool found(false);
  do {
    Face_handle f = Face_handle(fc);
    int i = ccw( f->index(v) );

    CGAL_assertion( f->vertex( cw(i) ) == v );

    Vertex_handle v0 = f->vertex( i );
    Vertex_handle v1 = f->mirror_vertex( i );

    bool is_admissible = (v0 != v1) &&
      !is_infinite(f) && !is_infinite( f->neighbor(i) );

    if ( is_admissible && is_degenerate_edge(f, i) ) {
      Edge e = flip(f, i);
      f = e.first;

      if ( !f->has_vertex(v) ) {
	f = e.first->neighbor(e.second);
	CGAL_assertion( f->has_vertex(v) );
      }

      fc = --( v->incident_faces(f) );
      fc_start = fc;
      found = true;
    } else {
      ++fc;
      found = false;
    }
  } while ( found || fc != fc_start );
}

template<class Gt, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
find_conflict_region_remove(const Vertex_handle& v,
			    const Vertex_handle& vnearest,
			    List& l, Face_map& fm,
			    //			    std::map<Face_handle,Sign>& sign_map,
			    std::vector<Vh_triple*>* fe)
{
  Site_2 t = v->site();

  // check if it is already inserted
  if ( t.is_point() && vnearest->is_point() &&
       are_same_points(t, vnearest->site()) ) {
    return;
  }

  CGAL_precondition( vnearest != Vertex_handle() );


  // find the first conflict

  // first look for conflict with vertex
  Face_circulator fc_start = vnearest->incident_faces();
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;
  do {
    Face_handle f(fc);
    //    int id = f->mirror_indexf->index(vnearest)
    s = incircle(f, t);

    if ( s == NEGATIVE ) {
      start_f = f;
      break;
    }
    ++fc;
  } while ( fc != fc_start );

  CGAL_assertion( s == NEGATIVE );

  // we are not in conflict with an Apollonius vertex, so we have to
  // be in conflict with the interior of an Apollonius edge
  if ( s != NEGATIVE ) {
    Edge_circulator ec_start = vnearest->incident_edges();
    Edge_circulator ec = ec_start;

    bool interior_in_conflict(false);
    Edge e;
    do {
      e = *ec;
      interior_in_conflict = edge_interior(e, t, false);

      if ( interior_in_conflict ) { break; }
      ++ec;
    } while ( ec != ec_start );

    CGAL_assertion( interior_in_conflict );

    l.push_back(e);
    l.push_back(sym_edge(e));
    return;
  }

  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, v->site(), l, fm, fe);
}
#endif

//----------------------------------------------------------------------
// access methods
//----------------------------------------------------------------------


template<class Gt, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::size_type
Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>::
number_of_incident_segments(Vertex_handle v) const
{
  CGAL_precondition( v->is_point() );
  Vertex_circulator vc_start = incident_vertices(v);
  Vertex_circulator vc = vc_start;

  unsigned int counter(0);
  do {
    Vertex_handle vn(vc);
    if ( vn->is_segment() &&
	 is_endpoint_of_segment(v->site(), vn->site()) ) {
      counter++;
    }
    ++vc;
  } while ( vc != vc_start );

  return counter;
}

CGAL_END_NAMESPACE
