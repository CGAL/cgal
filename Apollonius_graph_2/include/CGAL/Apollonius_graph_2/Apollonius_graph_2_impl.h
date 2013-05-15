// Copyright (c) 2003,2004,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>


#ifndef CGAL_APOLLONIUS_GRAPH_2_IMPL_H
#define CGAL_APOLLONIUS_GRAPH_2_IMPL_H

#include <CGAL/use.h>

// class implementation continued
//=================================

namespace CGAL {


//--------------------------------------------------------------------
// test method
//--------------------------------------------------------------------
template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
is_valid(bool verbose, int level) const
{
  if (level < 0) { return true; }

  if (number_of_vertices() <= 1) { return true; }

  // level 0 test: check the TDS
  bool result = this->_tds.is_valid(verbose, level);
  //  bool result(true);

  //  CGAL_assertion( result );

  if ( result && verbose ) {
    std::cerr << "AGDS is ok... " << std::flush;
  }

  if (level == 0) { return result; }

  // level 1 test: do the incircle tests

  if (number_of_vertices() < 3)  return true;

  //  CGAL_triangulation_assertion(result);

  for (All_edges_iterator eit = all_edges_begin();
       eit != all_edges_end(); ++eit) {
    Edge e = *eit;
    Face_handle f = e.first;
    Vertex_handle v = tds().mirror_vertex(f, e.second);
    if ( f->vertex(e.second) == v ) { continue; }
    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->site()) != NEGATIVE );
      //    CGAL_triangulation_assertion(result);
    }
    Edge sym_e = sym_edge(e);
    f = sym_e.first;
    v = tds().mirror_vertex(f, sym_e.second);
    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->site()) != NEGATIVE );
      //    CGAL_triangulation_assertion(result);
    }
  }

  if ( result && verbose ) {
    std::cerr << "Apollonius diagram is ok..." << std::flush;
  }
  if ( !result && verbose ) {
    std::cerr << "Apollonius diagram is NOT valid..." << std::flush;
  }

  //  CGAL_triangulation_assertion(result);
  return result;
}



//--------------------------------------------------------------------
// embedding and visualization methods and constructions for primal
// and dual
//--------------------------------------------------------------------

// circumcenter
template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Point_2
Apollonius_graph_2<Gt,Agds,LTag>::
circumcenter(const Face_handle& f) const
{
  CGAL_triangulation_precondition (dimension()==2 || !is_infinite(f));
  return circumcenter(f->vertex(0)->site(),
		      f->vertex(1)->site(),
		      f->vertex(2)->site());
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Point_2
Apollonius_graph_2<Gt,Agds,LTag>::
circumcenter(const Site_2& p0, const Site_2& p1, 
	     const Site_2& p2) const
{
  return
    geom_traits().construct_Apollonius_vertex_2_object()(p0, p1, p2);
}

// circumcircle
template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Site_2
Apollonius_graph_2<Gt,Agds,LTag>::
circumcircle(const Face_handle& f) const
{
  CGAL_triangulation_precondition (dimension()==2 || !is_infinite(f));
  return circumcircle(f->vertex(0)->site(),
		      f->vertex(1)->site(),
		      f->vertex(2)->site());
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Site_2
Apollonius_graph_2<Gt,Agds,LTag>::
circumcircle(const Site_2& p0, const Site_2& p1, 
	     const Site_2& p2) const
{
  return
    geom_traits().construct_Apollonius_site_2_object()(p0, p1, p2);
}


template<class Gt, class Agds, class LTag>
typename Gt::Line_2
Apollonius_graph_2<Gt,Agds,LTag>::
circumcircle(const Site_2& p0, const Site_2& p1) const
{
  return
    geom_traits().construct_Apollonius_site_2_object()(p0, p1);
}


// dual
template<class Gt, class Agds, class LTag>
typename Gt::Object_2
Apollonius_graph_2<Gt,Agds,LTag>::
dual(const Face_handle& f) const
{
  if ( !is_infinite(f) ) {
    Site_2 cc = circumcircle(f);
    return geom_traits().construct_object_2_object()(cc);
  }

  int i_inf(-1);
  for (int i = 0; i < 3; i++) {
    if ( is_infinite( f->vertex(i) ) ) {
      i_inf = i;
      break;
    }
  }
  typename Gt::Line_2 ll = circumcircle(f->vertex((i_inf+1)%3)->site(),
					f->vertex((i_inf+2)%3)->site());
  return geom_traits().construct_object_2_object()(ll);
}


template<class Gt, class Agds, class LTag>
typename Gt::Object_2
Apollonius_graph_2<Gt,Agds,LTag>::
dual(const Edge e) const
{
  CGAL_triangulation_precondition( !is_infinite(e) );

  if ( dimension() == 1 ) {
    Site_2 p = (e.first)->vertex(cw(e.second))->site();
    Site_2 q = (e.first)->vertex(ccw(e.second))->site();

    return construct_Apollonius_bisector_2_object()(p,q);
  }

  // dimension == 2
  // none of the two adjacent faces is infinite
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    Site_2 p = (e.first)->vertex( ccw(e.second) )->site();
    Site_2 q = (e.first)->vertex(  cw(e.second) )->site();
    Site_2 r = (e.first)->vertex(     e.second  )->site();
    Site_2 s = tds().mirror_vertex(e.first, e.second)->site();
    return construct_Apollonius_bisector_segment_2_object()(p,q,r,s);
  }

  // both of the adjacent faces are infinite
  if ( is_infinite(e.first) &&
       is_infinite(e.first->neighbor(e.second)) )  {
    Site_2 p = (e.first)->vertex(cw(e.second))->site();
    Site_2 q = (e.first)->vertex(ccw(e.second))->site();
    return construct_Apollonius_bisector_2_object()(p,q);
  }

  // only one of the adjacent faces is infinite
  CGAL_triangulation_assertion( is_infinite( e.first ) ||
				is_infinite( e.first->neighbor(e.second) )
				);

  CGAL_triangulation_assertion( !(is_infinite( e.first ) &&
				  is_infinite( e.first->neighbor(e.second) )
				  )
				);

  CGAL_triangulation_assertion
    (  is_infinite( e.first->vertex(e.second) ) ||
       is_infinite( tds().mirror_vertex(e.first, e.second) )  );

  Edge ee = e;
  if ( is_infinite( e.first->vertex(e.second) )  ) {
    ee = sym_edge(e);
  }
  Site_2 p = ee.first->vertex( ccw(ee.second) )->site();
  Site_2 q = ee.first->vertex(  cw(ee.second) )->site();
  Site_2 r = ee.first->vertex(     ee.second  )->site();

  return construct_Apollonius_bisector_ray_2_object()(p,q,r);
}



// primal
template<class Gt, class Agds, class LTag>
typename Gt::Object_2
Apollonius_graph_2<Gt,Agds,LTag>::
primal(const Edge e) const
{
  typedef typename Geom_traits::Segment_2  Segment;
  typedef typename Geom_traits::Ray_2      Ray;
  // typedef CGAL::Hyperbola_segment_2<Gt>    Hyperbola_segment;
  typedef CGAL::Parabola_segment_2<Gt>     Parabola_segment;
  //  typedef typename Geom_traits::Hyperbola_segment_2  Hyperbola_segment;
  //  typedef typename Geom_traits::Parabola_segment_2   Parabola_segment;

  //  CGAL_triangulation_precondition( !is_infinite(e) );

  if ( number_of_vertices() != 2 ) {
    if ( is_infinite(e) ) {
      Ray ray;
      if (  is_infinite( e.first->vertex(cw(e.second)) )  ) {
	Site_2 p = e.first->vertex( ccw(e.second) )->site();
	Site_2 r = e.first->vertex( e.second )->site();
	Site_2 s = tds().mirror_vertex( e.first, e.second )->site();
	ray = construct_Apollonius_primal_ray_2_object()(p,r,s);
      } else {
	CGAL_triangulation_assertion
	  (   is_infinite( e.first->vertex(ccw(e.second)) )   );
	Site_2 q = e.first->vertex( cw(e.second) )->site();
	Site_2 r = e.first->vertex( e.second )->site();
	Site_2 s = tds().mirror_vertex( e.first, e.second )->site();
	ray = construct_Apollonius_primal_ray_2_object()(q,s,r);
      }
      return make_object(ray);
    }
  }

  if ( dimension() == 1 ) {
    Site_2 p = (e.first)->vertex(cw(e.second))->site();
    Site_2 q = (e.first)->vertex(ccw(e.second))->site();
    Segment seg = construct_Apollonius_primal_segment_2_object()(p, q);
    return make_object(seg);

  }

  // dimension == 2
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    Site_2 p = (e.first)->vertex( ccw(e.second) )->site();
    Site_2 q = (e.first)->vertex(  cw(e.second) )->site();
    Site_2 r = (e.first)->vertex(     e.second  )->site();
    Site_2 s = tds().mirror_vertex(e.first, e.second)->site();
    return construct_Apollonius_primal_segment_2_object()(p,q,r,s);
  }

  // both of the adjacent faces are infinite
  if ( is_infinite(e.first) &&
       is_infinite(e.first->neighbor(e.second)) )  {
    Site_2 p = (e.first)->vertex(cw(e.second))->site();
    Site_2 q = (e.first)->vertex(ccw(e.second))->site();
    Segment seg = construct_Apollonius_primal_segment_2_object()(p,q);
    return make_object(seg);
  }

  // only one of the adjacent faces is infinite
  Edge ee = e;
  if ( is_infinite(e.first) ) {
    ee = sym_edge(e);
  }
  Site_2 p = (ee.first)->vertex( ccw(ee.second) )->site();
  Site_2 q = (ee.first)->vertex(  cw(ee.second) )->site();
  Site_2 r = (ee.first)->vertex(     ee.second  )->site();
  Parabola_segment ps =
    construct_Apollonius_primal_segment_2_object()(p,q,r);
  return make_object(ps);
}


//--------------------------------------------------------------------
// combinatorial operations
//--------------------------------------------------------------------
template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Edge
Apollonius_graph_2<Gt,Agds,LTag>::
flip(Face_handle& f, int i)
{
  CGAL_triangulation_precondition ( f != Face_handle() );
  CGAL_triangulation_precondition (i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition( dimension()==2 ); 

  CGAL_triangulation_precondition( f->vertex(i) != tds().mirror_vertex(f,i) );

  this->_tds.flip(f, i);

  return Edge(f, ccw(i));
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Edge
Apollonius_graph_2<Gt,Agds,LTag>::
flip(Edge e)
{
  return flip(e.first, e.second);
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
insert_in_face(Face_handle& f, const Site_2& p)
{
  Vertex_handle v = this->_tds.insert_in_face( f );

  v->set_site(p);
  return v;
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
is_degree_2(const Vertex_handle& v) const
{
  Face_circulator fc = incident_faces(v);
  Face_circulator fc1 = fc;
  ++(++fc1);
  return ( fc == fc1 );
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
insert_degree_2(Edge e)
{
  return this->_tds.insert_degree_2(e.first, e.second);
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
insert_degree_2(Edge e, const Site_2& p)
{
  Vertex_handle v = insert_degree_2(e);

  v->set_site(p);
  return v;
}


template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove_degree_2(Vertex_handle v)
{
  CGAL_triangulation_precondition( is_degree_2(v) );

  this->_tds.remove_degree_2( v );
}


template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove_degree_3(Vertex_handle v)
{
  remove_degree_3(v, Face_handle());
}


template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove_degree_3(Vertex_handle v, Face_handle f)
{
  CGAL_triangulation_precondition( degree(v) == 3 );
  this->_tds.remove_degree_3(v, f);
}

//--------------------------------------------------------------------
// insertion of weighted point
//--------------------------------------------------------------------

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
insert_first(const Site_2& p)
{
  CGAL_triangulation_precondition(number_of_vertices() == 0);
  Vertex_handle v = this->_tds.insert_second();
  v->set_site(p);
  return v;
  //  return Delaunay_graph::insert_first(p);
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
insert_second(const Site_2& p)
{
  CGAL_triangulation_precondition( number_of_vertices() == 1 );
  Vertex_handle vnew;
  Vertex_handle v(finite_vertices_begin());
  if ( is_hidden(v->site(), p) ) {
    v->add_hidden_site(p);
    vnew = Vertex_handle();  
  } else if ( is_hidden(p, v->site()) ) {
    v->add_hidden_site(v->site());
    v->set_site(p);
    vnew = v;
  } else {
    CGAL_triangulation_precondition(number_of_vertices() == 1);
    vnew = this->_tds.insert_dim_up(infinite_vertex(), true);
    vnew->set_site(p);
    //    vnew = Delaunay_graph::insert_second(p);
  }

  return vnew;
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
insert_third(const Site_2& p)
{
  CGAL_triangulation_precondition( number_of_vertices() == 2 );

  Face_handle f((*finite_edges_begin()).first);
  Vertex_handle v1(f->vertex(0));
  Vertex_handle v2(f->vertex(1));

  if ( is_hidden(v1->site(), p) ) {
    v1->add_hidden_site(p);
    return Vertex_handle();
  }
  if ( is_hidden(v2->site(), p) ) {
    v2->add_hidden_site(p);
    return Vertex_handle();
  }

  bool t1 = is_hidden(p, v1->site());
  bool t2 = is_hidden(p, v2->site());

  if ( t1 && !t2 ) {
    v1->add_hidden_site(v1->site());
    v1->set_site(p);
    return v1;
  } else if ( !t1 && t2 ) {
    v2->add_hidden_site(v2->site());
    v2->set_site(p);
    return v2;
  } else if ( t1 && t2 ) {
    v1->add_hidden_site(v1->site());
    v1->add_hidden_site(v2->site());
    v1->set_site(p);
    remove_second(v2);
    return v1;
  }

  Conflict_type ct =
    finite_edge_conflict_type_degenerated(v1->site(), v2->site(), p);

  if ( ct == RIGHT_VERTEX || ct == LEFT_VERTEX ) {
    Vertex_handle v =
      this->_tds.insert_dim_up(infinite_vertex(), ct == LEFT_VERTEX);
    v->set_site(p);
    return v;
  }


  // SPECIAL CODE -- START
  // code that should work for pseudo-circles, but only consumes time
  // for circles
  if ( ct == NO_CONFLICT ) {
    Conflict_type ct1 =
      infinite_edge_conflict_type(v1->site(), v2->site(), v2->site(), p);

    Conflict_type ct2 =
      infinite_edge_conflict_type(v2->site(), v1->site(), v1->site(), p);

    if ( ct1 == NO_CONFLICT && ct2 == NO_CONFLICT ) {
      return Vertex_handle();
    }
  }

  if ( ct == ENTIRE_EDGE ) {
    Conflict_type ct1 =
      infinite_edge_conflict_type(v1->site(), v2->site(), v2->site(), p);

    if ( ct1 == ENTIRE_EDGE ) {
      v1->add_hidden_site(v1->site());
      v2->add_hidden_site(v1->site());
      v1->set_site(p);
      return v1;
    }

    Conflict_type ct2 =
      infinite_edge_conflict_type(v2->site(), v1->site(), v1->site(), p);

    if ( ct2 == ENTIRE_EDGE ) {
      v1->add_hidden_site(v2->site());
      v2->add_hidden_site(v2->site());
      v2->set_site(p);
      return v2;
    }
  }
  // SPECIAL CODE -- END


  Vertex_handle v = this->_tds.insert_dim_up(infinite_vertex());
  v->set_site(p);
  if ( ct == NO_CONFLICT ) {
    Conflict_type ct1 =
      finite_edge_conflict_type_degenerated(v1->site(), p, v2->site());

    CGAL_assertion( ct1 == NO_CONFLICT || ct1 == ENTIRE_EDGE );
    Vertex_handle vv = (ct1 == NO_CONFLICT) ? v2 : v1;

    Edge_circulator ec = incident_edges(v);
    while ( true ) {
      if ( ec->first->vertex( ccw(ec->second) ) == vv ) {
	flip(*ec);
	break;
      }
      ++ec;
    }
  } else if ( ct == INTERIOR ) {
    Edge_circulator ec = incident_edges(v);

    while ( true ) {
      if ( is_infinite(ec) ) {
	flip(*ec);
	break;
      }
      ec++;
    }
  } else if ( ct == ENTIRE_EDGE ) {
    Face_circulator fc = incident_faces(v);

    while ( true ) {
      Face_handle ff(fc);
      if ( !is_infinite(ff) ) {
	flip(ff, ff->index(v));
	break;
      }
      ++fc;
    }
  } else {
    CGAL_assertion( ct == BOTH_VERTICES );

    Conflict_type ct1 =
      finite_edge_conflict_type_degenerated(v1->site(), p, v2->site());

    Edge_circulator ec =
      ( ct1 == INTERIOR ) ? incident_edges(v2) : incident_edges(v1);
    while ( true ) {
      if ( is_infinite(ec) ) {
	flip(*ec);
	break;
      }
      ec++;
    }
  }

  //  CGAL_triangulation_assertion( is_valid() );

  return v;
}


template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
insert(const Site_2& p, Vertex_handle vnear)
{
  if ( number_of_vertices() == 0 ) {
    return insert_first(p);
  }
  if ( number_of_vertices() == 1 ) {
    return insert_second(p);
  }
  if ( number_of_vertices() == 2 ) {
    return insert_third(p);
  }

  // first find the nearest neighbor
  Vertex_handle vnearest = nearest_neighbor(p.point(), vnear);

  CGAL_assertion( vnearest != Vertex_handle() );


  // check if it is hidden
  Site_2 wp_nearest = vnearest->site();
  if ( is_hidden(wp_nearest, p) ) {
    vnearest->add_hidden_site(p);
    return Vertex_handle();
  }

  // find the first conflict

  // first look for conflict with vertex
  Face_circulator fc_start = incident_faces(vnearest);
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;
  do {
    Face_handle f(fc);
    s = incircle(f, p);

    if ( s == NEGATIVE ) {
      start_f = f;
      break;
    }
    ++fc;
  } while ( fc != fc_start );

  // we are not in conflict with an Apollonius vertex, so we have to
  // be in conflict with the interior of an Apollonius edge
  if ( s != NEGATIVE ) {
    Edge_circulator ec_start = incident_edges(vnearest);
    Edge_circulator ec = ec_start;

    bool interior_in_conflict(false);
    Edge e;
    do {
      e = *ec;

      interior_in_conflict = edge_interior(e, p, false);

      if ( interior_in_conflict ) { break; }
      ++ec;
    } while ( ec != ec_start );

    if ( !interior_in_conflict ) {
      return Vertex_handle();
    }

    return insert_degree_2(e, p);
  }


  // we are in conflict with an Apollonius vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;
  Face_map fm;
  Vertex_map vm;

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.
  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, p, l, fm, vm, NULL);

  //  retriangulate_conflict_region(v, l, fm, vm);
  Vertex_handle v = retriangulate_conflict_region(p, l, fm, vm);

  fm.clear();
  vm.clear();


  return v;
}


//--------------------------------------------------------------------
// find conflict region
//--------------------------------------------------------------------
template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
find_conflict_region_remove(const Vertex_handle& v,
			    const Vertex_handle& vnearest,
			    List& l, Face_map& fm,
			    Vertex_map& vm,
			    std::vector<Vh_triple*>* fe)
{
  Site_2 p = v->site();
  // check if it is hidden
  Site_2 wp_nearest = vnearest->site();
  if ( is_hidden(wp_nearest, p) ) {
    vnearest->add_hidden_site(p);
    return;
  }

  CGAL_precondition( vnearest != Vertex_handle() );


  // find the first conflict

  // first look for conflict with vertex
  Face_circulator fc_start = incident_faces(vnearest);
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;
  do {
    Face_handle f(fc);
    s = incircle(f, p);

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
    Edge_circulator ec_start = incident_edges(vnearest);
    Edge_circulator ec = ec_start;

    bool interior_in_conflict(false);
    Edge e;
    do {
      e = *ec;
      interior_in_conflict = edge_interior(e, p, false);

      if ( interior_in_conflict ) { break; }
      ++ec;
    } while ( ec != ec_start );

    CGAL_assertion( interior_in_conflict );

    l.push_back(e);
    l.push_back(sym_edge(e));
    return;
  }

  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, v->site(), l, fm, vm, fe);
}

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
initialize_conflict_region(const Face_handle& f, List& l) const
{
  l.clear();
  for (int i = 0; i < 3; i++) {
    l.push_back(sym_edge(f, i));
  }
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
check_edge_for_hidden_sites(const Face_handle& f, int i,
			    const Site_2& p, Vertex_map& vm) const
{
  bool found(false);

  Vertex_handle v1 = f->vertex(ccw(i));
  if ( vm.find(v1) == vm.end() ) {
    if ( !is_infinite(v1) && is_hidden(p, v1->site()) ) {
      vm[v1] = true;
      found = true;
    }
  } else {
    found = true;
  }

  Vertex_handle v2 = f->vertex(cw(i));
  if ( vm.find(v2) == vm.end() ) {
    if ( !is_infinite(v2) && is_hidden(p, v2->site()) ) {
      vm[v2] = true;
      found = true;
    }
  } else {
    found = true;
  }

  return found;
}

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
expand_conflict_region(const Face_handle& f, const Site_2& p,
		       List& l, Face_map& fm, Vertex_map& vm,
		       std::vector<Vh_triple*>* fe)
{
  // setting fm[f] to true means that the face has been reached and
  // that the face is available for recycling. If we do not want the
  // face to be available for recycling we must set this flag to
  // false.
  fm[f] = true;

  //  CGAL_assertion( fm.find(f) != fm.end() );
  for (int i = 0; i < 3; i++) {
    bool hidden_found =
      check_edge_for_hidden_sites(f, i, p, vm);

    Face_handle n = f->neighbor(i);

    if ( !hidden_found ) {
      Sign s = incircle(n, p);
      if ( s != NEGATIVE ) { continue; }

      bool interior_in_conflict = edge_interior(f, i, p, true);

      if ( !interior_in_conflict ) { continue; }
    }

    if ( fm.find(n) != fm.end() ) {
      Edge e = sym_edge(f, i);
      if ( l.is_in_list(e) ||
	   l.is_in_list(sym_edge(e)) ) {
	l.remove(e);
	l.remove(sym_edge(e));
      }
      continue;
    }

    Edge e = sym_edge(f, i);

    CGAL_assertion( l.is_in_list(e) );
    int j = tds().mirror_index(f, i);
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

      (*vhq)[0] = Vertex_handle();
      (*vhq)[1] = n->vertex(     j  );
      (*vhq)[2] = n->vertex( ccw(j) );

      fe->push_back(vhq);
    }


    expand_conflict_region(n, p, l, fm, vm, fe);
  } // for-loop
}


//--------------------------------------------------------------------
// retriangulate conflict region
//--------------------------------------------------------------------

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
add_bogus_vertex(Edge e, List& l)
{
  Edge esym = sym_edge(e);
  Face_handle g1 = e.first;
  CGAL_assertion_code(Face_handle g2 = esym.first);

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

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_list
Apollonius_graph_2<Gt,Agds,LTag>::
add_bogus_vertices(List& l)
{
  Vertex_list vertex_list;
  Edge_list edge_list;

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

  typename Edge_list::iterator it;

  for (it = edge_list.begin();  it != edge_list.end(); ++it) {
    e = *it;
    Vertex_handle v = add_bogus_vertex(e, l);
    vertex_list.push_back(v);
  }

  return vertex_list;
}

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove_bogus_vertices(Vertex_list& vl)
{
  while ( vl.size() > 0 ) {
    Vertex_handle v = vl.front();
    vl.pop_front();
    remove_degree_2(v);
  }
}


template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
move_hidden_sites(Vertex_handle& vold, Vertex_handle& vnew)
{
  typename Vertex::Hidden_sites_iterator wpit;

  for (wpit = vold->hidden_sites_begin();
       wpit != vold->hidden_sites_end(); ++wpit) {
    vnew->add_hidden_site(*wpit);
  }

  vold->clear_hidden_sites_container();
}




template<class Gt, class Agds, class LTag>
std::vector<typename Apollonius_graph_2<Gt,Agds,LTag>::Face*>
Apollonius_graph_2<Gt,Agds,LTag>::
get_faces_for_recycling(Face_map& fm, unsigned int n_wanted)
{
  std::vector<Face*> vf;

  typename Face_map::iterator fmit;
  for (fmit = fm.begin(); fmit != fm.end(); ++fmit) {
    Face_handle f = (*fmit).first;
    if ( fm[f] == true ) { vf.push_back( f ); }
  }

  while ( vf.size() < n_wanted ) {
    Face* fp = static_cast<Face*>(this->_tds.create_face());
    vf.push_back(fp);
  }

  while ( vf.size() > n_wanted ) {
    Face* fp = vf.back();
    vf.pop_back();
    this->_tds.delete_face(fp);
  }
  
  return vf;
}

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove_hidden_vertices(Vertex_map& vm)
{
  typename Vertex_map::iterator it;

  for (it = vm.begin(); it != vm.end(); ++it) {
    Vertex_handle vhidden = (*it).first;
    this->_tds.delete_vertex( vhidden );
  }
  vm.clear();
}


template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
retriangulate_conflict_region(const Site_2& p,	List& l, 
			      Face_map& fm, Vertex_map& vm)
{
  size_type vmsize = vm.size();
  size_type num_vert = number_of_vertices();

  if ( num_vert - vmsize == 0 ) {
    // 1. copy all hidden sites to a temporary list
    Site_list wp_list;
    typename Vertex_map::iterator vmit;
    for (vmit = vm.begin(); vmit != vm.end(); ++vmit) {
      Vertex_handle vhidden = (*vmit).first;

      wp_list.push_back(vhidden->site());
      typename Vertex::Hidden_sites_iterator it;
      for (it = vhidden->hidden_sites_begin();
	   it != vhidden->hidden_sites_end(); ++it) {
	wp_list.push_back(*it);
      }	
      vhidden->clear_hidden_sites_container();
    }

    // 2. clear the current Apollonius diagram
    clear();

    // 3. add a new vertex
    Vertex_handle v = insert_first(p);

    // 4. add all old sites to the hidden site list of the
    // new site
    Site_list_iterator wpit;
    for (wpit = wp_list.begin(); wpit != wp_list.end(); ++wpit) {
      v->add_hidden_site(*wpit);
    }

    return v;
  } else if ( num_vert - vmsize == 1 ) {
    // 1. copy all hidden sites to a temporary list
    Site_list wp_list;
    typename Vertex_map::iterator vmit;
    for (vmit = vm.begin(); vmit != vm.end(); ++vmit) {
      Vertex_handle vhidden = (*vmit).first;

      wp_list.push_back(vhidden->site());
      typename Vertex::Hidden_sites_iterator it;
      for (it = vhidden->hidden_sites_begin();
	   it != vhidden->hidden_sites_end(); ++it) {
	wp_list.push_back(*it);
      }	
      vhidden->clear_hidden_sites_container();
    }

    // 2. find which vertex remains non-hidden and copy its hidden
    // sites to a local container
    Vertex_handle non_hidden;
    Finite_vertices_iterator vit = finite_vertices_begin();
    do {
      non_hidden = Vertex_handle(vit);
      ++vit;
    } while ( vm.find(non_hidden) != vm.end() );

    Site_2 p1 = non_hidden->site();
    Site_list wp_list1;
    typename Vertex::Hidden_sites_iterator it;
    for (it = non_hidden->hidden_sites_begin();
	 it != non_hidden->hidden_sites_end(); ++it) {
      wp_list1.push_back(*it);
    }	
    non_hidden->clear_hidden_sites_container();

    // 3. clear the current Apollonius graph
    clear();

    // 4. insert the two non-hidden sites and copy the corresponding
    // hidden sites
    Vertex_handle v1 = insert_first(p1);
    for (Site_list_iterator it = wp_list1.begin();
	 it != wp_list1.end(); ++it) {
      v1->add_hidden_site(*it);
    }

    Vertex_handle v = insert_second(p);
    for (Site_list_iterator it = wp_list.begin();
	 it != wp_list.end(); ++it) {
      v->add_hidden_site(*it);
    }

    return v;
  }

  Vertex_handle v = this->_tds.create_vertex();
  v->set_site(p);

  // 1. move all the hidden sites to the new one
  typename Vertex_map::iterator vmit;
  for (vmit = vm.begin(); vmit != vm.end(); ++vmit) {
    Vertex_handle vhidden = (*vmit).first;
    move_hidden_sites(vhidden, v);
    v->add_hidden_site(vhidden->site());
  }

  CGAL_precondition( number_of_vertices() - vm.size() >= 2 );

  // 2. add the bogus vetrices
  Vertex_list dummy_vertices = add_bogus_vertices(l);

  // 3. repair the face pointers...
  Edge e_start = l.front();
  Edge eit = e_start;
  do {
    CGAL_assertion_code(Edge esym =)sym_edge(eit);
    Face_handle f = eit.first;
    int k = eit.second;
    CGAL_assertion( !l.is_in_list(esym) );
    CGAL_assertion( fm.find(f) == fm.end() );
    f->vertex(ccw(k))->set_face(f);
    f->vertex( cw(k))->set_face(f);
    eit = l.next(eit);
  } while ( eit != e_start );

  //  std::vector<Face*> vf = get_faces_for_recycling(fm, l.size());
  std::list<Face*> vf;

  // 4. copy the edge list to a vector of edges and clear the in place 
  //    list
  typedef typename Agds::Edge Agds_edge;
  std::vector<Agds_edge> ve;

  Edge efront = l.front();
  Edge e = efront;
  do {
    ve.push_back(Agds_edge(e.first, e.second));
    e = l.next(e);
  } while ( e != efront );

  l.clear();

  // 5. remove the hidden vertices
  remove_hidden_vertices(vm);

  // 6. retriangulate the hole
  //  _tds.star_hole( v, ve.begin(), ve.end(), vf.begin(), vf.end());
  this->_tds.star_hole(v, ve.begin(), ve.end());

  // 7. remove the bogus vertices
  remove_bogus_vertices(dummy_vertices);

  // 8. remove the unused faces
  typename Face_map::iterator it;
  for (it = fm.begin(); it != fm.end(); ++it) {
    Face_handle fh = (*it).first;
    this->_tds.delete_face( fh );
  }

  CGAL_assertion( number_of_vertices() + vmsize == num_vert + 1 );

  // 9. DONE!!!!
  return v;
}



//--------------------------------------------------------------------
// point location
//--------------------------------------------------------------------
template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
nearest_neighbor(const Point_2& p) const
{
  return nearest_neighbor(p, Vertex_handle());
}


template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Vertex_handle
Apollonius_graph_2<Gt,Agds,LTag>::
nearest_neighbor(const Point_2& p,
		 Vertex_handle start_vertex) const
{
  if ( number_of_vertices() == 0 ) {
    return Vertex_handle();
  }

  if ( start_vertex == Vertex_handle() ) {
    start_vertex = finite_vertex();
  }

  //  if ( start_vertex == Vertex_handle() ) { return start_vertex; }

  Vertex_handle vclosest;
  Vertex_handle v = start_vertex;

  if ( number_of_vertices() < 3 ) {
    vclosest = v;
    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      Vertex_handle v1(vit);
      if ( v1 != vclosest /*&& !is_infinite(v1)*/ ) {
	Site_2 p1 = vclosest->site();
	Site_2 p2 = v1->site();
	if ( side_of_bisector(p1, p2, p) == ON_NEGATIVE_SIDE ) {
	  vclosest = v1;
	}
      }
    }
    return vclosest;
  }

  do {
    vclosest = v;
    Site_2 p1 = v->site();
    Vertex_circulator vc_start = incident_vertices(v);
    Vertex_circulator vc = vc_start;
    do {
      if ( !is_infinite(vc) ) {
	Vertex_handle v1(vc);
	Site_2 p2 = v1->site();
	if ( side_of_bisector(p1, p2, p) == ON_NEGATIVE_SIDE ) {
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

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
is_hidden(const Site_2 &p, const Site_2 &q) const
{
  return geom_traits().is_hidden_2_object()(p, q);
}

template<class Gt, class Agds, class LTag>
Oriented_side
Apollonius_graph_2<Gt,Agds,LTag>::
side_of_bisector(const Site_2 &p1,
		 const Site_2 &p2,
		 const Point_2 &p) const
{
  return geom_traits().oriented_side_of_bisector_2_object()(p1, p2, p);
}


template<class Gt, class Agds, class LTag>
Sign
Apollonius_graph_2<Gt,Agds,LTag>::
incircle(const Site_2 &p1, const Site_2 &p2,
	 const Site_2 &p3, const Site_2 &q) const
{
  return geom_traits().vertex_conflict_2_object()(p1, p2, p3, q);
}

template<class Gt, class Agds, class LTag>
Sign
Apollonius_graph_2<Gt,Agds,LTag>::
incircle(const Site_2 &p1, const Site_2 &p2,
	 const Site_2 &q) const
{
  return
    geom_traits().vertex_conflict_2_object()(p1, p2, q);
}


template<class Gt, class Agds, class LTag>
Sign
Apollonius_graph_2<Gt,Agds,LTag>::
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


template<class Gt, class Agds, class LTag>
Sign
Apollonius_graph_2<Gt,Agds,LTag>::
incircle(const Vertex_handle& v0, const Vertex_handle& v1,
	 const Vertex_handle& v) const
{
  CGAL_precondition( !is_infinite(v0) && !is_infinite(v1)
		     && !is_infinite(v) );

  return incircle( v0->site(), v1->site(), v->site());
}

template<class Gt, class Agds, class LTag>
Sign
Apollonius_graph_2<Gt,Agds,LTag>::
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


template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
finite_edge_interior(const Site_2& p1,
		     const Site_2& p2,
		     const Site_2& p3,
		     const Site_2& p4,
		     const Site_2& q, bool b) const
{
  if ( is_hidden(q, p1) ) { return true; }
  if ( is_hidden(q, p2) ) { return true; }
  return
    geom_traits().finite_edge_interior_conflict_2_object()(p1,p2,p3,p4,q,b);
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
finite_edge_interior(const Face_handle& f, int i,
		     const Site_2& p, bool b) const
{
  CGAL_precondition( !is_infinite(f) &&
		     !is_infinite(f->neighbor(i)) );
  return finite_edge_interior( f->vertex( ccw(i) )->site(),
			       f->vertex(  cw(i) )->site(),
			       f->vertex(     i  )->site(),
			       tds().mirror_vertex(f, i)->site(), p, b);
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
finite_edge_interior(const Vertex_handle& v1,
		     const Vertex_handle& v2,
		     const Vertex_handle& v3,
		     const Vertex_handle& v4,
		     const Vertex_handle& v, bool b) const
{
  CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) &&
		     !is_infinite(v3) && !is_infinite(v4) &&
		     !is_infinite(v) );
  return finite_edge_interior( v1->site(), v2->site(),
			       v3->site(), v4->site(),
			       v->site(), b);
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
finite_edge_interior_degenerated(const Site_2& p1,
				 const Site_2& p2,
				 const Site_2& p3,
				 const Site_2& q,
				 bool b) const
{
  if ( is_hidden(q, p1) ) { return true; }
  if ( is_hidden(q, p2) ) { return true; }
  return
    geom_traits().finite_edge_interior_conflict_2_object()(p1,p2,p3,q,b);
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
finite_edge_interior_degenerated(const Site_2& p1,
				 const Site_2& p2,
				 const Site_2& q,
				 bool b) const
{
  if ( is_hidden(q, p1) ) { return true; }
  if ( is_hidden(q, p2) ) { return true; }
  return
    geom_traits().finite_edge_interior_conflict_2_object()(p1, p2, q, b);
}


template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
finite_edge_interior_degenerated(const Face_handle& f, int i,
				 const Site_2& p, bool b) const
{
  if ( !is_infinite( tds().mirror_vertex(f, i) ) ) {
    CGAL_precondition( is_infinite(f->vertex(i)) );

    Face_handle g = f->neighbor(i);
    int j = tds().mirror_index(f, i);

    return finite_edge_interior_degenerated(g, j, p, b);
  }

  CGAL_precondition( is_infinite( tds().mirror_vertex(f, i) ) );

  Site_2 p1 = f->vertex( ccw(i) )->site();
  Site_2 p2 = f->vertex(  cw(i) )->site();

  if ( is_infinite(f->vertex(i)) ) {
    return finite_edge_interior_degenerated(p1, p2, p, b);
  }

  Site_2 p3 = f->vertex(i)->site();
  return finite_edge_interior_degenerated(p1, p2, p3, p, b);
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
finite_edge_interior_degenerated(const Vertex_handle& v1,
				 const Vertex_handle& v2,
				 const Vertex_handle& v3,
				 const Vertex_handle& v4,
				 const Vertex_handle& v, bool b) const
{
  CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) && 
		     !is_infinite(v) );
  if ( !is_infinite( v4 ) ) {
    CGAL_precondition( is_infinite(v3) );

    return
      finite_edge_interior_degenerated(v2, v1, v4, v3, v, b);
  }

  CGAL_precondition( is_infinite( v4 ) );

  Site_2 p1 = v1->site();
  Site_2 p2 = v2->site();
  Site_2 p = v->site();

  if ( is_infinite(v3) ) {
    return finite_edge_interior_degenerated(p1, p2, p, b);
  }

  Site_2 p3 = v3->site();
  return finite_edge_interior_degenerated(p1, p2, p3, p, b);
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
infinite_edge_interior(const Site_2& p2,
		       const Site_2& p3,
		       const Site_2& p4,
		       const Site_2& q,
		       bool b) const
{
  if ( is_hidden(q, p2) ) { return true; }
  return
    geom_traits().infinite_edge_interior_conflict_2_object()(p2,p3,p4,q,b);
}

template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
infinite_edge_interior(const Face_handle& f, int i,
		       const Site_2& p, bool b) const
{
  if ( !is_infinite( f->vertex(ccw(i)) ) ) {
    CGAL_precondition( is_infinite( f->vertex(cw(i)) ) );
    Face_handle g = f->neighbor(i);
    int j = tds().mirror_index(f, i);

    return infinite_edge_interior(g, j, p, b);
  }

  CGAL_precondition( is_infinite( f->vertex(ccw(i)) ) );

  Site_2 p2 = f->vertex(  cw(i) )->site();
  Site_2 p3 = f->vertex(     i  )->site();
  Site_2 p4 = tds().mirror_vertex(f, i)->site();

  return infinite_edge_interior(p2, p3, p4, p, b);
}


template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
infinite_edge_interior(const Vertex_handle& v1,
		       const Vertex_handle& v2,
		       const Vertex_handle& v3,
		       const Vertex_handle& v4,
		       const Vertex_handle& v, bool b) const
{
  CGAL_precondition( !is_infinite(v3) && !is_infinite(v4) && 
		     !is_infinite(v) );

  if ( !is_infinite( v1 ) ) {
    CGAL_precondition( is_infinite( v2 ) );

    return infinite_edge_interior(v2, v1, v4, v3, v, b);
  }

  CGAL_precondition( is_infinite( v1 ) );

  Site_2 p2 = v2->site();
  Site_2 p3 = v3->site();
  Site_2 p4 = v4->site();
  Site_2 p = v->site();

  return infinite_edge_interior(p2, p3, p4, p, b);
}




template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
edge_interior(const Vertex_handle& v1,
	      const Vertex_handle& v2,
	      const Vertex_handle& v3,
	      const Vertex_handle& v4,
	      const Vertex_handle& v, bool b) const
{
  CGAL_precondition( !is_infinite(v) );

  bool is_inf_v1 = is_infinite(v1);
  bool is_inf_v2 = is_infinite(v2);
  bool is_inf_v3 = is_infinite(v3);
  bool is_inf_v4 = is_infinite(v4);

  bool result;

  if ( !is_inf_v1 && !is_inf_v2 && !is_inf_v3 && !is_inf_v4 ) {
    result = finite_edge_interior(v1, v2, v3, v4, v, b);
  } else if ( is_inf_v3 || is_inf_v4 ) {
    result = finite_edge_interior_degenerated(v1, v2, v3, v4, v, b);
  } else {
    result = infinite_edge_interior(v1, v2, v3, v4, v, b);
  }

  return result;
}


template<class Gt, class Agds, class LTag>
bool
Apollonius_graph_2<Gt,Agds,LTag>::
edge_interior(const Face_handle& f, int i,
	      const Site_2& p, bool b) const
{
  Face_handle g = f->neighbor(i);

  bool is_inf_f = is_infinite(f);
  bool is_inf_g = is_infinite(g);

  bool result;

  if ( !is_inf_f && !is_inf_g ) {
    result = finite_edge_interior(f, i, p, b);
  } else if ( !is_inf_f || !is_inf_g ) {
    result = finite_edge_interior_degenerated(f, i, p, b);
  } else {
    //    Edge e(f, i);
    if ( !is_infinite(f, i) ) {
      result = finite_edge_interior_degenerated(f, i, p, b);
    } else {
      result = infinite_edge_interior(f, i, p, b);
    }
  }

  return result;
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Conflict_type
Apollonius_graph_2<Gt,Agds,LTag>::
infinite_edge_conflict_type(const Site_2& p2,
			    const Site_2& p3,
			    const Site_2& p4,
			    const Site_2& q) const
{
  Sign i1 = incircle(p2, p3, q);
  Sign i2 = incircle(p4, p2, q);

  CGAL_assertion( i1 != ZERO && i2 != ZERO );

  if ( i1 == NEGATIVE && i2 == POSITIVE ) {
    return LEFT_VERTEX;
  } else if ( i1 == POSITIVE && i2 == NEGATIVE ) {
    return RIGHT_VERTEX;
  } else if ( i1 == POSITIVE && i2 == POSITIVE ) {
    bool b = infinite_edge_interior(p2, p3, p4, q, false);
    return (b ? INTERIOR : NO_CONFLICT);
  } else if ( i1 == NEGATIVE && i2 == NEGATIVE ) {
    bool b = infinite_edge_interior(p2, p3, p4, q, true);
    return (b ? ENTIRE_EDGE : BOTH_VERTICES);
  } else {
    // this should never be reached; the degenerated incircle never
    // returns ZERO
    CGAL_error();
  }

  // to satisfy compiler
  return NO_CONFLICT;
}

template<class Gt, class Agds, class LTag>
typename Apollonius_graph_2<Gt,Agds,LTag>::Conflict_type
Apollonius_graph_2<Gt,Agds,LTag>::
finite_edge_conflict_type_degenerated(const Site_2& p1,
				      const Site_2& p2,
				      const Site_2& q) const
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
    CGAL_error();
  }

  // to satisfy compiler
  return NO_CONFLICT;
}


//----------------------------------------------------------------------
// methods for disk removal
//----------------------------------------------------------------------


template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove_first(Vertex_handle v)
{
  Delaunay_graph::remove_first(v);
}

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove_second(Vertex_handle v)
{
  Delaunay_graph::remove_second(v);
}

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove_third(Vertex_handle v)
{
  if ( is_degree_2(v) ) {
    Face_handle fh( incident_faces(v) );
    int i = fh->index(v);
    flip(fh, i);
  } else if ( degree(v) == 4 ) {
    Edge_circulator ec = incident_edges(v);
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

  this->_tds.remove_dim_down( v );
}


template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle() );
  CGAL_triangulation_precondition( !is_infinite(v) );

  // find a neighbor of v to use for point location of hidden sites to
  // be re-inserted
  Vertex_handle vnear;
  if ( /*StoreHidden*/ true ) {
    if ( number_of_vertices() > 10 ) {
      Vertex_circulator vc_start = incident_vertices(v);
      Vertex_circulator vc = vc_start;
      do {
	if ( !is_infinite(vc) ) {
	  vnear = Vertex_handle(vc);
	  break;
	}
	++vc;
      } while ( vc != vc_start );
    }
  }

  Site_list wp_list;
  typename Vertex::Hidden_sites_iterator wpit;
  for (wpit = v->hidden_sites_begin();
       wpit != v->hidden_sites_end(); ++wpit) {
    wp_list.push_back(*wpit);
  }

  size_type n = number_of_vertices();
  if ( n == 1 ) {
    remove_first(v);
  } else if ( n == 2 ) {
    remove_second(v);
  } else if ( n == 3 ) {
    remove_third(v);
  } else {
    std::size_t deg = degree(v);
    if ( deg == 2 ) {
      remove_degree_2(v);
    } else if ( deg == 3 ) {
      remove_degree_3(v);
    } else {
      remove_degree_d_vertex(v);
    }
  }

  Site_less_than_comparator less_than(geom_traits());
  std::sort(wp_list.begin(), wp_list.end(), less_than);
  for (unsigned int i = 0; i < wp_list.size(); i++) {
    vnear = insert(wp_list[i], vnear);
  }
}


template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
remove_degree_d_vertex(Vertex_handle v)
{
  minimize_degree(v);
  std::size_t deg = degree(v);
  if ( deg == 3 ) {
    remove_degree_3(v);
    return;
  }
  if ( deg == 2 ) {
    remove_degree_2(v);
    return;
  }
  
  Apollonius_graph_2<Gt,Agds,LTag> ag_small;

  std::map<Vertex_handle,Vertex_handle> vmap;

  Vertex_circulator vc_start = incident_vertices(v);
  Vertex_circulator vc = incident_vertices(v);
  Vertex_handle vh_large, vh_small;
  do {
    vh_large = Vertex_handle(vc);
    if ( is_infinite(vh_large) ) {
      vh_small = ag_small.infinite_vertex();
      vmap[vh_small] = vh_large;
    } else { 
      vh_small = ag_small.insert(vc->site());
      if ( vh_small != Vertex_handle() ) {
	vmap[vh_small] = vh_large;
      }
    }
    ++vc;
  } while ( vc != vc_start );

  if ( ag_small.number_of_vertices() == 2 ) {
    CGAL_assertion( deg == 4 );
    Edge_circulator ec = incident_edges(v);
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


  Vertex_handle vn = ag_small.nearest_neighbor(v->site().point());

  CGAL_assertion( vn != Vertex_handle() );

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

  std::size_t num_fe = flipped_edges.size();
  for (std::size_t i = 0; i < num_fe; i++) {
    Vh_triple *vhq = flipped_edges[num_fe - i - 1];

    bool found(false);
    ec = incident_edges(v);
    Edge_circulator ec_start = ec;
    do {
      Edge e = *ec;
      if ( (e.first->vertex(  cw(e.second) ) == vmap[(*vhq)[1]] &&
	    e.first->vertex(     e.second  ) == vmap[(*vhq)[2]]) ||
	   (e.first->vertex( ccw(e.second) ) == vmap[(*vhq)[1]] &&
	    tds().mirror_vertex(e.first,e.second) == vmap[(*vhq)[2]]) ) {
	flip(e);
	found = true;
	break;
      }
      ++ec;
    } while ( ec != ec_start );

    CGAL_assertion( found );
    CGAL_USE(found);
  }
  CGAL_triangulation_precondition( degree(v) == 3 );

  this->_tds.remove_degree_3( v, Face_handle() );

  for (unsigned int i = 0; i < num_fe; i++) {
    delete[] flipped_edges[i];
  }
}


template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::
minimize_degree(Vertex_handle v)
{
  CGAL_precondition ( degree(v) > 3 );

  Face_circulator fc_start = incident_faces(v);
  Face_circulator fc = incident_faces(v);
  bool found(false);
  do {
    Face_handle f = Face_handle(fc);
    int i = ccw( f->index(v) );

    CGAL_assertion( f->vertex( cw(i) ) == v );

    Vertex_handle v0 = f->vertex( i );
    Vertex_handle v1 = tds().mirror_vertex( f, i );

    bool is_admissible = (v0 != v1) &&
      !is_infinite(f) && !is_infinite( f->neighbor(i) );

    if ( is_admissible && is_degenerate_edge(f, i) ) {
      Edge e = flip(f, i);
      f = e.first;

      if ( !f->has_vertex(v) ) {
	f = e.first->neighbor(e.second);
	CGAL_assertion( f->has_vertex(v) );
      }

      fc = --( incident_faces(v,f) );
      fc_start = fc;
      found = true;
    } else {
      ++fc;
      found = false;
    }
  } while ( found || fc != fc_start );
}


//----------------------------------------------------------------------
// methods for I/O
//----------------------------------------------------------------------

template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::file_output(std::ostream& os) const
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

  std::map<Vertex_handle,int> V;
  std::map<Face_handle,int> F;

  // first vertex (infinite vertex) 
  int inum = 0;
  V[infinite_vertex()] = inum++;
  
  // finite vertices
  if (is_ascii(os)) os << std::endl;
  for (Finite_vertices_iterator vit = finite_vertices_begin();
       vit != finite_vertices_end(); ++vit) {
    V[vit] = inum++;
    os << vit->site();
    if ( is_ascii(os) ) { os << ' '; }
    os << vit->number_of_hidden_sites();
    typename Vertex::Hidden_sites_iterator hit;
    for (hit = vit->hidden_sites_begin(); hit != vit->hidden_sites_end();
	 ++hit) {
      if ( is_ascii(os) ) { os << ' '; }
      os << *hit;
    }
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


template<class Gt, class Agds, class LTag>
void
Apollonius_graph_2<Gt,Agds,LTag>::file_input(std::istream& is)
{
  //input from file
  size_type n, m;
  int d;
  is >> n >> m >> d;

  CGAL_assertion( n >= 1 );

  if ( n == 1 ) {
    CGAL_assertion( d == -1 );
    if ( number_of_vertices() > 0 ) { clear(); }
    return;
  }
  if ( n == 2 ) {
    CGAL_assertion( d == 0 );
    if ( number_of_vertices() > 0 ) { clear(); }
    Site_2 s;
    is >> s;
    Vertex_handle v = insert(s);
    unsigned int n_hidden;
    is >> n_hidden;
    for (unsigned int i = 0; i < n_hidden; i++) {
      is >> s;
      v->add_hidden_site(s);
    }
    return;
  }
  if ( n == 3 ) {
    CGAL_assertion( d == 1 );
    if ( number_of_vertices() > 0 ) { clear(); }
    for (int j = 0; j < 2; j++) {
      Site_2 s;
      is >> s;
      Vertex_handle v = insert(s);
      unsigned int n_hidden;
      is >> n_hidden;
      for (unsigned int i = 0; i < n_hidden; i++) {
	is >> s;
	v->add_hidden_site(s);
      }
    }
    return;
  }

  if (this->_tds.number_of_vertices() != 0) { this->_tds.clear(); }

  this->_tds.set_dimension(d);

  std::vector<Vertex_handle> V(n);
  std::vector<Face_handle> F(m);

  size_type i = 0;

  // first vertex (infinite vertex)
  V[0] = create_vertex();
  this->set_infinite_vertex(V[0]);
  i++;


  // read vertices
  for (; i < n; ++i) {
    V[i] = create_vertex();
    Site_2 s;
    is >> s;
    V[i]->set_site(s);
    unsigned int n_hidden;
    is >> n_hidden;
    for (unsigned int j = 0; j < n_hidden; j++) {
      is >> s;
      V[i]->add_hidden_site(s);
    }
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


} //namespace CGAL


#endif // CGAL_APOLLONIUS_GRAPH_2_IMPL_H
