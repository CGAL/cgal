// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Apollonius_graph_2.C
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================




// class implementation continued
//=================================

CGAL_BEGIN_NAMESPACE


//--------------------------------------------------------------------
// test method
//--------------------------------------------------------------------
template< class Gt, bool StoreHidden, class Agds >
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
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
    Vertex_handle v = f->mirror_vertex(e.second);
    if ( f->vertex(e.second) == v ) { continue; }
    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->point()) != NEGATIVE );
      //    CGAL_triangulation_assertion(result);
    }
    Edge sym_e = sym_edge(e);
    f = sym_e.first;
    v = f->mirror_vertex(sym_e.second);
    if ( !is_infinite(v) ) {
      result = result &&
	( incircle(f, v->point()) != NEGATIVE );
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
template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Point_2
Apollonius_graph_2<Gt,StoreHidden,Agds>::
circumcenter(const Face_handle& f) const
{
  CGAL_triangulation_precondition (dimension()==2 || !is_infinite(f));
  return circumcenter(f->vertex(0)->point(),
		      f->vertex(1)->point(),
		      f->vertex(2)->point());
}

template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Point_2
Apollonius_graph_2<Gt,StoreHidden,Agds>::
circumcenter(const Weighted_point_2& p0, const Weighted_point_2& p1, 
	     const Weighted_point_2& p2) const
{
  return
    geom_traits().construct_Apollonius_vertex_2_object()(p0, p1, p2);
}

// circumcircle
template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Weighted_point_2
Apollonius_graph_2<Gt,StoreHidden,Agds>::
circumcircle(const Face_handle& f) const
{
  CGAL_triangulation_precondition (dimension()==2 || !is_infinite(f));
  return circumcircle(f->vertex(0)->point(),
		      f->vertex(1)->point(),
		      f->vertex(2)->point());
}

template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Weighted_point_2
Apollonius_graph_2<Gt,StoreHidden,Agds>::
circumcircle(const Weighted_point_2& p0, const Weighted_point_2& p1, 
	     const Weighted_point_2& p2) const
{
  return
    geom_traits().construct_Apollonius_weighted_point_2_object()(p0, p1, p2);
}


template< class Gt, bool StoreHidden, class Agds >
inline
typename Gt::Line_2
Apollonius_graph_2<Gt,StoreHidden,Agds>::
circumcircle(const Weighted_point_2& p0, const Weighted_point_2& p1) const
{
  return
    geom_traits().construct_Apollonius_weighted_point_2_object()(p0, p1);
}


// dual
template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Point_2
Apollonius_graph_2<Gt,StoreHidden,Agds>::
dual(const Face_handle& f) const
{
  return circumcenter(f);
}


template< class Gt, bool StoreHidden, class Agds >
inline
Object
Apollonius_graph_2<Gt,StoreHidden,Agds>::
dual(const Edge e) const
{
  CGAL_triangulation_precondition( !is_infinite(e) );

  if ( dimension() == 1 ) {
    Weighted_point_2 p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point_2 q = (e.first)->vertex(ccw(e.second))->point();

    return geom_traits().construct_Apollonius_bisector_2_object()(p,q);
  }

  // dimension == 2
  // none of the two adjacent faces is infinite
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    Weighted_point_2 p = (e.first)->vertex( ccw(e.second) )->point();
    Weighted_point_2 q = (e.first)->vertex(  cw(e.second) )->point();
    Weighted_point_2 r = (e.first)->vertex(     e.second  )->point();
    Weighted_point_2 s = (e.first)->mirror_vertex(e.second)->point();
    return
      geom_traits().construct_Apollonius_bisector_segment_2_object()(p,q,r,s);
  }

  // both of the adjacent faces are infinite
  if ( is_infinite(e.first) &&
       is_infinite(e.first->neighbor(e.second)) )  {
    Weighted_point_2 p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point_2 q = (e.first)->vertex(ccw(e.second))->point();
    return geom_traits().construct_Apollonius_bisector_2_object()(p,q);
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
       is_infinite( e.first->mirror_vertex(e.second) )  );

  Edge ee = e;
  if ( is_infinite( e.first->vertex(e.second) )  ) {
    ee = sym_edge(e);
  }
  Weighted_point_2 p = ee.first->vertex( ccw(ee.second) )->point();
  Weighted_point_2 q = ee.first->vertex(  cw(ee.second) )->point();
  Weighted_point_2 r = ee.first->vertex(     ee.second  )->point();

  return geom_traits().construct_Apollonius_bisector_ray_2_object()(p,q,r);
}



// primal
template< class Gt, bool StoreHidden, class Agds >
inline
Object
Apollonius_graph_2<Gt,StoreHidden,Agds>::
primal(const Edge e) const
{
  typedef typename Geom_traits::Segment_2  Segment;
  typedef typename Geom_traits::Ray_2      Ray;
  typedef typename Geom_traits::Hyperbola_segment_2  Hyperbola_segment;
  typedef typename Geom_traits::Parabola_segment_2   Parabola_segment;

  //  CGAL_triangulation_precondition( !is_infinite(e) );

  if ( number_of_vertices() != 2 ) {
    if ( is_infinite(e) ) {
      Ray ray;
      if (  is_infinite( e.first->vertex(cw(e.second)) )  ) {
	Weighted_point_2 p = e.first->vertex( ccw(e.second) )->point();
	Weighted_point_2 r = e.first->vertex( e.second )->point();
	Weighted_point_2 s = e.first->mirror_vertex( e.second )->point();
	ray = geom_traits().construct_Apollonius_primal_ray_2_object()(p,r,s);
      } else {
	CGAL_triangulation_assertion
	  (   is_infinite( e.first->vertex(ccw(e.second)) )   );
	Weighted_point_2 q = e.first->vertex( cw(e.second) )->point();
	Weighted_point_2 r = e.first->vertex( e.second )->point();
	Weighted_point_2 s = e.first->mirror_vertex( e.second )->point();
	ray = geom_traits().construct_Apollonius_primal_ray_2_object()(q,s,r);
      }
      return make_object(ray);
    }
  }

  if ( dimension() == 1 ) {
    Weighted_point_2 p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point_2 q = (e.first)->vertex(ccw(e.second))->point();
    Segment seg =
      geom_traits().construct_Apollonius_primal_segment_2_object()(p, q);
    return make_object(seg);

  }

  // dimension == 2
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    Weighted_point_2 p = (e.first)->vertex( ccw(e.second) )->point();
    Weighted_point_2 q = (e.first)->vertex(  cw(e.second) )->point();
    Weighted_point_2 r = (e.first)->vertex(     e.second  )->point();
    Weighted_point_2 s = (e.first)->mirror_vertex(e.second)->point();
    return
      geom_traits().construct_Apollonius_primal_segment_2_object()(p,q,r,s);
  }

  // both of the adjacent faces are infinite
  if ( is_infinite(e.first) &&
       is_infinite(e.first->neighbor(e.second)) )  {
    Weighted_point_2 p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point_2 q = (e.first)->vertex(ccw(e.second))->point();
    Segment seg =
      geom_traits().construct_Apollonius_primal_segment_2_object()(p,q);
    return make_object(seg);
  }

  // only one of the adjacent faces is infinite
  Edge ee = e;
  if ( is_infinite(e.first) ) {
    ee = sym_edge(e);
  }
  Weighted_point_2 p = (ee.first)->vertex( ccw(ee.second) )->point();
  Weighted_point_2 q = (ee.first)->vertex(  cw(ee.second) )->point();
  Weighted_point_2 r = (ee.first)->vertex(     ee.second  )->point();
  Parabola_segment ps =
    geom_traits().construct_Apollonius_primal_segment_2_object()(p,q,r);
  return make_object(ps);
}


//--------------------------------------------------------------------
// combinatorial operations
//--------------------------------------------------------------------
template< class Gt, bool StoreHidden, class Agds >
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Edge
Apollonius_graph_2<Gt,StoreHidden,Agds>::
flip(Face_handle& f, int i)
{
  CGAL_triangulation_precondition ( &(*f) != NULL );
  CGAL_triangulation_precondition (i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition( dimension()==2 ); 

  CGAL_triangulation_precondition( f->vertex(i) != f->mirror_vertex(i) );

  this->_tds.flip( &(*f), i);

  return Edge(f, ccw(i));
}

template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Edge
Apollonius_graph_2<Gt,StoreHidden,Agds>::
flip(Edge e)
{
  return flip(e.first, e.second);
}

template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
insert_in_face(Face_handle& f, const Weighted_point_2& p)
{
  Vertex_handle v = static_cast<Vertex*>(this->_tds.insert_in_face( &(*f) ));

  v->set_point(p);
  return v;
}

template< class Gt, bool StoreHidden, class Agds >
inline
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
is_degree_2(const Vertex_handle& v) const
{
  Face_circulator fc = v->incident_faces();
  Face_circulator fc1 = fc;
  ++(++fc1);
  return ( fc == fc1 );
}

template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
insert_degree_2(Edge e)
{
  return this->_tds.insert_degree_2(&(*e.first),e.second);
}

template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
insert_degree_2(Edge e, const Weighted_point_2& p)
{
  Vertex_handle v = insert_degree_2(e);

  v->set_point(p);
  return v;
}


template< class Gt, bool StoreHidden, class Agds >
inline
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
remove_degree_2(Vertex_handle v)
{
  CGAL_triangulation_precondition( is_degree_2(v) );

  this->_tds.remove_degree_2( &(*v) );
}


template< class Gt, bool StoreHidden, class Agds >
inline
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
remove_degree_3(Vertex_handle v)
{
  remove_degree_3(v, NULL);
}


template< class Gt, bool StoreHidden, class Agds >
inline
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
remove_degree_3(Vertex_handle v, Face* f)
{
  CGAL_triangulation_precondition( v->degree() == 3 );
  this->_tds.remove_degree_3( &(*v), f);
}

//--------------------------------------------------------------------
// insertion of weighted point
//--------------------------------------------------------------------

template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
insert_first(const Weighted_point_2& p)
{
  return Delaunay_graph::insert_first(p);
}

template< class Gt, bool StoreHidden, class Agds >
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
insert_second(const Weighted_point_2& p)
{
  CGAL_triangulation_precondition( number_of_vertices() == 1 );
  Vertex_handle vnew;
  Vertex_handle v(finite_vertices_begin());
  if ( is_hidden(v->point(), p) ) {
    v->add_hidden_weighted_point(p);
    vnew = Vertex_handle(NULL);  
  } else if ( is_hidden(p, v->point()) ) {
    v->add_hidden_weighted_point(v->point());
    v->set_point(p);
    vnew = v;
  } else {
    vnew = Delaunay_graph::insert_second(p);
  }

  return vnew;
}

template< class Gt, bool StoreHidden, class Agds >
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
insert_third(const Weighted_point_2& p)
{
  CGAL_triangulation_precondition( number_of_vertices() == 2 );

  Vertex_handle v1(vertices_begin());
  Vertex_handle v2(++vertices_begin());

  if ( is_hidden(v1->point(), p) ) {
    v1->add_hidden_weighted_point(p);
    return Vertex_handle(NULL);
  }
  if ( is_hidden(v2->point(), p) ) {
    v2->add_hidden_weighted_point(p);
    return Vertex_handle(NULL);
  }

  bool t1 = is_hidden(p, v1->point());
  bool t2 = is_hidden(p, v2->point());

  if ( t1 && !t2 ) {
    v1->add_hidden_weighted_point(v1->point());
    v1->set_point(p);
    return v1;
  } else if ( !t1 && t2 ) {
    v2->add_hidden_weighted_point(v2->point());
    v2->set_point(p);
    return v2;
  } else if ( t1 && t2 ) {
    v1->add_hidden_weighted_point(v1->point());
    v1->add_hidden_weighted_point(v2->point());
    v1->set_point(p);
    remove_second(v2);
    return v1;
  }

  Vertex_handle v = this->_tds.insert_dim_up(infinite_vertex());
  v->set_point(p);

  Face_handle f(finite_faces_begin());

  Point_2 p1 = f->vertex(0)->point().point();
  Point_2 p2 = f->vertex(1)->point().point();
  Point_2 p3 = f->vertex(2)->point().point();

  Orientation o =
    geom_traits().orientation_2_object()(p1, p2, p3);

  if ( o != LEFT_TURN ) {
    f->reorient();
    for (int i = 0; i < 3; i++) {
      f->neighbor(i)->reorient();
    }
  }

  Conflict_type ct =
    finite_edge_conflict_type_degenerated(v1->point(), v2->point(), p);

  if ( ct == NO_CONFLICT ) {
    Oriented_side os =
      side_of_bisector(v1->point(), v2->point(), p);

    CGAL_assertion( os != ON_ORIENTED_BOUNDARY );
    Vertex_handle vv = ( os == ON_NEGATIVE_SIDE ) ? v1 : v2;

    Face_circulator fc = incident_faces(v);
    while ( true ) {
      Face_handle f(fc);
      int k = f->index(v);
      Vertex_handle vh = f->vertex(ccw(k));
      if ( vh == vv ) {
	flip(f, cw(k));
	break;
      }
      ++fc;
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
      Face_handle f(fc);
      if ( !is_infinite(f) ) {
	flip(f, f->index(v));
	break;
      }
      ++fc;
    }
  } else if ( ct == BOTH_VERTICES ) {


    Conflict_type ct1 =
      finite_edge_conflict_type_degenerated(v1->point(), p, v2->point());

    Edge_circulator ec;
    ec = ( ct1 == INTERIOR ) ? incident_edges(v2) : incident_edges(v1);
    while ( true ) {
      if ( is_infinite(ec) ) {
	flip(*ec);
	break;
      }
      ec++;
    }
  } else {
    CGAL_assertion( ct == RIGHT_VERTEX || ct == LEFT_VERTEX );
    // do nothing here
  }

  //  CGAL_triangulation_assertion( is_valid() );

  return v;
}



template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
insert(const Weighted_point_2& p)
{
  return insert(p, NULL);
}


template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
insert(const Weighted_point_2& p, Vertex_handle vnear)
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

  CGAL_assertion( &(*vnearest) != NULL );


  // check if it is hidden
  Weighted_point_2 wp_nearest = vnearest->point();
  if ( is_hidden(wp_nearest, p) ) {
    vnearest->add_hidden_weighted_point(p);
    return Vertex_handle(NULL);
  }

  // find the first conflict

  // first look for conflict with vertex
  Face_circulator fc_start = vnearest->incident_faces();
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
    Edge_circulator ec_start = vnearest->incident_edges();
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

    return insert_degree_2(e, p);
  }


  // we are in conflict with an Apollonius vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;
  Face_map fm;
  Vertex_map vm;

  //  Vertex_handle v = _tds.create_vertex();
  //  v->set_point(p);

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
template< class Gt, bool StoreHidden, class Agds >
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
find_conflict_region_remove(const Vertex_handle& v,
			    const Vertex_handle& vnearest,
			    List& l, Face_map& fm,
			    Vertex_map& vm,
			    std::vector<Vh_triple*>* fe)
{
  Weighted_point_2 p = v->point();
  // check if it is hidden
  Weighted_point_2 wp_nearest = vnearest->point();
  if ( is_hidden(wp_nearest, p) ) {
    vnearest->add_hidden_weighted_point(p);
    return;
  }

  CGAL_precondition( &(*vnearest) != NULL );


  // find the first conflict

  // first look for conflict with vertex
  Face_circulator fc_start = vnearest->incident_faces();
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;
  do {
    Face_handle f(fc);
    //    int id = f->mirror_indexf->index(vnearest)
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
    Edge_circulator ec_start = vnearest->incident_edges();
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

    l.push(e);
    l.push(sym_edge(e));
    return;
  }

  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, v->point(), l, fm, vm, fe);
}

template< class Gt, bool StoreHidden, class Agds >
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
initialize_conflict_region(const Face_handle& f, List& l)
{
  l.clear();
  for (int i = 0; i < 3; i++) {
    l.push(sym_edge(f, i));
  }
}

template< class Gt, bool StoreHidden, class Agds >
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
check_edge_for_hidden_weighted_points(const Face_handle& f, int i,
				      const Weighted_point_2& p,
				      Vertex_map& vm)
{
  bool found(false);

  Vertex_handle v1 = f->vertex(ccw(i));
  if ( vm.find(v1) == vm.end() ) {
    if ( !is_infinite(v1) && is_hidden(p, v1->point()) ) {
      vm[v1] = true;
      found = true;
    }
  } else {
    found = true;
  }

  Vertex_handle v2 = f->vertex(cw(i));
  if ( vm.find(v2) == vm.end() ) {
    if ( !is_infinite(v2) && is_hidden(p, v2->point()) ) {
      vm[v2] = true;
      found = true;
    }
  } else {
    found = true;
  }

  return found;
}

template< class Gt, bool StoreHidden, class Agds >
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
expand_conflict_region(const Face_handle& f, const Weighted_point_2& p,
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
      check_edge_for_hidden_weighted_points(f, i, p, vm);

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


    expand_conflict_region(n, p, l, fm, vm, fe);
  } // for-loop
}


//--------------------------------------------------------------------
// retriangulate conflict region
//--------------------------------------------------------------------

template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
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

template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_list
Apollonius_graph_2<Gt,StoreHidden,Agds>::
add_bogus_vertices(List& l)
{
  Vertex_list vertex_list;
  Edge_list edge_list;

  Edge e_start = l.front();
  Edge e = e_start;
  do {
    Edge esym = sym_edge(e);
    if ( l.is_in_list(esym) ) {
      bool found(false);
      typename Edge_list::iterator it;
      for (it = edge_list.begin(); it != edge_list.end(); ++it) {
	if ( *it == esym ) {
	  found = true;
	  break;
	}
      }
      if ( !found ) { edge_list.push_back(e); }
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

template< class Gt, bool StoreHidden, class Agds >
inline
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
remove_bogus_vertices(Vertex_list& vl)
{
  while ( vl.size() > 0 ) {
    Vertex_handle v = vl.front();
    vl.pop_front();
    remove_degree_2(v);
  }
}


template< class Gt, bool StoreHidden, class Agds >
inline
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
move_hidden_weighted_points(Vertex_handle& vold, Vertex_handle& vnew)
{
  typename Vertex_base::Hidden_weighted_point_iterator wpit;

  for (wpit = vold->hidden_weighted_points_begin();
       wpit != vold->hidden_weighted_points_end(); ++wpit) {
    vnew->add_hidden_weighted_point(*wpit);
  }

  vold->clear_hidden_weighted_point_container();
}




template< class Gt, bool StoreHidden, class Agds >
std::vector<typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Face*>
Apollonius_graph_2<Gt,StoreHidden,Agds>::
get_faces_for_recycling(Face_map& fm, unsigned int n_wanted)
{
  std::vector<Face*> vf;

  typename Face_map::iterator fmit;
  for (fmit = fm.begin(); fmit != fm.end(); ++fmit) {
    Face_handle f = (*fmit).first;
    if ( fm[f] == true ) { vf.push_back( &(*f) ); }
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

template< class Gt, bool StoreHidden, class Agds >
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
remove_hidden_vertices(Vertex_handle&v, Vertex_map& vm, Face_map& fm)
{
  typename Vertex_map::iterator it;

  for (it = vm.begin(); it != vm.end(); ++it) {
    Vertex_handle vhidden = (*it).first;
    this->_tds.delete_vertex( &(*vhidden) );
  }
  vm.clear();
}


template< class Gt, bool StoreHidden, class Agds >
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
retriangulate_conflict_region(const Weighted_point_2& p,	List& l, 
			      Face_map& fm, Vertex_map& vm)
{
  int vmsize = vm.size();
  int num_vert = number_of_vertices();

  if ( num_vert - vmsize == 0 ) {
    // 1. copy all hidden sites to a temporary list
    Weighted_point_list wp_list;
    typename Vertex_map::iterator vmit;
    for (vmit = vm.begin(); vmit != vm.end(); ++vmit) {
      Vertex_handle vhidden = (*vmit).first;

      wp_list.push_back(vhidden->point());
      typename Vertex_base::Hidden_weighted_point_iterator it;
      for (it = vhidden->hidden_weighted_points_begin();
	   it != vhidden->hidden_weighted_points_end(); ++it) {
	wp_list.push_back(*it);
      }	
      vhidden->clear_hidden_weighted_point_container();
    }

    // 2. clear the current Apollonius diagram
    clear();

    // 3. add a new vertex
    Vertex_handle v = Delaunay_graph::insert_first(p);

    // 4. add all old sites to the hidden weighted point list of the
    // new site
    Weighted_point_list_iterator wpit;
    for (wpit = wp_list.begin(); wpit != wp_list.end(); ++wpit) {
      v->add_hidden_weighted_point(*wpit);
    }

    return v;
  } else if ( num_vert - vmsize == 1 ) {
    // 1. copy all hidden sites to a temporary list
    Weighted_point_list wp_list;
    typename Vertex_map::iterator vmit;
    for (vmit = vm.begin(); vmit != vm.end(); ++vmit) {
      Vertex_handle vhidden = (*vmit).first;

      wp_list.push_back(vhidden->point());
      typename Vertex_base::Hidden_weighted_point_iterator it;
      for (it = vhidden->hidden_weighted_points_begin();
	   it != vhidden->hidden_weighted_points_end(); ++it) {
	wp_list.push_back(*it);
      }	
      vhidden->clear_hidden_weighted_point_container();
    }

    // 2. find which vertex remains non-hidden and copy its hidden
    // sites to a local container
    Vertex_handle non_hidden;
    Finite_vertices_iterator vit = finite_vertices_begin();
    do {
      non_hidden = Vertex_handle(vit);
      ++vit;
    } while ( vm.find(non_hidden) != vm.end() );

    Weighted_point_2 p1 = non_hidden->point();
    Weighted_point_list wp_list1;
    typename Vertex_base::Hidden_weighted_point_iterator it;
    for (it = non_hidden->hidden_weighted_points_begin();
	 it != non_hidden->hidden_weighted_points_end(); ++it) {
      wp_list1.push_back(*it);
    }	
    non_hidden->clear_hidden_weighted_point_container();

    // 3. clear the current Apollonius graph
    clear();

    // 4. insert the two non-hidden sites and copy the corresponding
    // hidden sites
    Vertex_handle v1 = Delaunay_graph::insert_first(p1);
    for (Weighted_point_list_iterator it = wp_list1.begin();
	 it != wp_list1.end(); ++it) {
      v1->add_hidden_weighted_point(*it);
    }

    Vertex_handle v = Delaunay_graph::insert_second(p);
    for (Weighted_point_list_iterator it = wp_list.begin();
	 it != wp_list.end(); ++it) {
      v->add_hidden_weighted_point(*it);
    }

    return v;
  }

  Vertex_handle v = this->_tds.create_vertex();
  v->set_point(p);

  // 1. move all the hidden weighted points to the new one
  typename Vertex_map::iterator vmit;
  for (vmit = vm.begin(); vmit != vm.end(); ++vmit) {
    Vertex_handle vhidden = (*vmit).first;
    move_hidden_weighted_points(vhidden, v);
    v->add_hidden_weighted_point(vhidden->point());
  }

  CGAL_precondition( number_of_vertices() - vm.size() >= 2 );

  // 2. add the bogus vetrices
  Vertex_list dummy_vertices = add_bogus_vertices(l);

  // 3. repair the face pointers...
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
  std::list<Face*> vf;

  // 4. copy the edge list to a vector of edges and clear the in place 
  //    list
  typedef typename Agds::Edge Agds_edge;
  std::vector<Agds_edge> ve;

  Edge efront = l.front();
  Edge e = efront;
  do {
    ve.push_back(Agds_edge(&(*e.first), e.second));
    e = l.next(e);
  } while ( e != efront );

  l.clear();

  // 5. remove the hidden vertices
  remove_hidden_vertices(v, vm, fm);

  // 6. retriangulate the hole
  //  _tds.star_hole( &(*v), ve.begin(), ve.end(), vf.begin(), vf.end());
  this->_tds.star_hole(&(*v), ve.begin(), ve.end());

  // 7. remove the bogus vertices
  remove_bogus_vertices(dummy_vertices);

  // 8. remove the unused faces
  typename Face_map::iterator it;
  for (it = fm.begin(); it != fm.end(); ++it) {
    Face_handle fh = (*it).first;
    this->_tds.delete_face( &(*fh) );
  }

  CGAL_assertion( number_of_vertices() == num_vert - vmsize + 1 );

  // 9. DONE!!!!
  return v;
}



//--------------------------------------------------------------------
// point location
//--------------------------------------------------------------------
template< class Gt, bool StoreHidden, class Agds >
inline
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
nearest_neighbor(const Point_2& p) const
{
  return nearest_neighbor(p, NULL);
}


template< class Gt, bool StoreHidden, class Agds >
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Vertex_handle
Apollonius_graph_2<Gt,StoreHidden,Agds>::
nearest_neighbor(const Point_2& p,
		 Vertex_handle start_vertex) const
{
  if ( number_of_vertices() == 0 ) {
    return Vertex_handle(NULL);
  }

  if ( &(*start_vertex) == NULL ) {
    start_vertex = finite_vertex();
  }

  //  if ( &(*start_vertex) == NULL ) { return start_vertex; }

  Vertex_handle vclosest;
  Vertex_handle v = start_vertex;

  if ( number_of_vertices() < 3 ) {
    vclosest = v;
    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      Vertex_handle v1(vit);
      if ( v1 != vclosest /*&& !is_infinite(v1)*/ ) {
	Weighted_point_2 p1 = vclosest->point();
	Weighted_point_2 p2 = v1->point();
	if ( side_of_bisector(p1, p2, p) == ON_NEGATIVE_SIDE ) {
	  vclosest = v1;
	}
      }
    }
    return vclosest;
  }

  do {
    vclosest = v;
    Weighted_point_2 p1 = v->point();
    Vertex_circulator vc_start = incident_vertices(v);
    Vertex_circulator vc = vc_start;
    do {
      if ( !is_infinite(vc) ) {
	Vertex_handle v1(vc);
	Weighted_point_2 p2 = v1->point();
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

template< class Gt, bool StoreHidden, class Agds >
inline
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
is_hidden(const Weighted_point_2 &p, const Weighted_point_2 &q) const
{
  return geom_traits().is_hidden_2_object()(p, q);
}

template< class Gt, bool StoreHidden, class Agds >
inline
Oriented_side
Apollonius_graph_2<Gt,StoreHidden,Agds>::
side_of_bisector(const Weighted_point_2 &p1,
		 const Weighted_point_2 &p2,
		 const Point_2 &p) const
{
  return geom_traits().oriented_side_of_bisector_2_object()(p1, p2, p);
}


template< class Gt, bool StoreHidden, class Agds >
inline
Sign
Apollonius_graph_2<Gt,StoreHidden,Agds>::
incircle(const Weighted_point_2 &p1, const Weighted_point_2 &p2,
	 const Weighted_point_2 &p3,	const Weighted_point_2 &q) const
{
  return geom_traits().vertex_conflict_2_object()(p1, p2, p3, q);
}

template< class Gt, bool StoreHidden, class Agds >
inline
Sign
Apollonius_graph_2<Gt,StoreHidden,Agds>::
incircle(const Weighted_point_2 &p1, const Weighted_point_2 &p2,
	 const Weighted_point_2 &q) const
{
  return
    geom_traits().vertex_conflict_2_object()(p1, p2, q);
}


template< class Gt, bool StoreHidden, class Agds >
inline
Sign
Apollonius_graph_2<Gt,StoreHidden,Agds>::
incircle(const Face_handle& f, const Weighted_point_2& q) const
{
  if ( !is_infinite(f) ) {
    return incircle(f->vertex(0)->point(),
		    f->vertex(1)->point(),
		    f->vertex(2)->point(), q);
  }

  int inf_i(-1); // to avoid compiler warning
  for (int i = 0; i < 3; i++) {
    if ( is_infinite(f->vertex(i)) ) {
      inf_i = i;
      break;
    }
  }
  return incircle( f->vertex( ccw(inf_i) )->point(),
		   f->vertex(  cw(inf_i) )->point(), q );
}


template< class Gt, bool StoreHidden, class Agds >
inline
Sign
Apollonius_graph_2<Gt,StoreHidden,Agds>::
incircle(const Vertex_handle& v0, const Vertex_handle& v1,
	 const Vertex_handle& v) const
{
  CGAL_precondition( !is_infinite(v0) && !is_infinite(v1)
		     && !is_infinite(v) );

  return incircle( v0->point(), v1->point(), v->point());
}

template< class Gt, bool StoreHidden, class Agds >
inline
Sign
Apollonius_graph_2<Gt,StoreHidden,Agds>::
incircle(const Vertex_handle& v0, const Vertex_handle& v1,
	 const Vertex_handle& v2, const Vertex_handle& v) const
{
  CGAL_precondition( !is_infinite(v) );

  if ( !is_infinite(v0) && !is_infinite(v1) &&
       !is_infinite(v2) ) {
    return incircle(v0->point(), v1->point(),
		    v2->point(), v->point());
  }

  if ( is_infinite(v0) ) {
    CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) );
    return incircle( v1->point(), v2->point(), v->point());
  }
  if ( is_infinite(v1) ) {
    CGAL_precondition( !is_infinite(v0) && !is_infinite(v2) );
    return incircle( v2->point(), v0->point(), v->point());
  }

  CGAL_assertion( is_infinite(v2) );
  CGAL_precondition( !is_infinite(v0) && !is_infinite(v1) );
  return incircle( v0->point(), v1->point(), v->point());
}


template< class Gt, bool StoreHidden, class Agds >
inline
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
finite_edge_interior(const Weighted_point_2& p1,
		     const Weighted_point_2& p2,
		     const Weighted_point_2& p3,
		     const Weighted_point_2& p4,
		     const Weighted_point_2& q, bool b) const
{
  if ( is_hidden(q, p1) ) { return true; }
  if ( is_hidden(q, p2) ) { return true; }
  return
    geom_traits().finite_edge_interior_conflict_2_object()(p1,p2,p3,p4,q,b);
}

template< class Gt, bool StoreHidden, class Agds >
inline
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
finite_edge_interior(const Face_handle& f, int i,
		     const Weighted_point_2& p, bool b) const
{
  CGAL_precondition( !is_infinite(f) &&
		     !is_infinite(f->neighbor(i)) );
  return finite_edge_interior( f->vertex( ccw(i) )->point(),
			       f->vertex(  cw(i) )->point(),
			       f->vertex(     i  )->point(),
			       f->mirror_vertex(i)->point(), p, b);
}

template< class Gt, bool StoreHidden, class Agds >
inline
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
finite_edge_interior(const Vertex_handle& v1,
		     const Vertex_handle& v2,
		     const Vertex_handle& v3,
		     const Vertex_handle& v4,
		     const Vertex_handle& v, bool b) const
{
  CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) &&
		     !is_infinite(v3) && !is_infinite(v4) &&
		     !is_infinite(v) );
  return finite_edge_interior( v1->point(), v2->point(),
			       v3->point(), v4->point(),
			       v->point(), b);
}

template< class Gt, bool StoreHidden, class Agds >
inline
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
finite_edge_interior_degenerated(const Weighted_point_2& p1,
				 const Weighted_point_2& p2,
				 const Weighted_point_2& p3,
				 const Weighted_point_2& q,
				 bool b) const
{
  if ( is_hidden(q, p1) ) { return true; }
  if ( is_hidden(q, p2) ) { return true; }
  return
    geom_traits().finite_edge_interior_conflict_2_object()(p1,p2,p3,q,b);
}

template< class Gt, bool StoreHidden, class Agds >
inline
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
finite_edge_interior_degenerated(const Weighted_point_2& p1,
				 const Weighted_point_2& p2,
				 const Weighted_point_2& q,
				 bool b) const
{
  if ( is_hidden(q, p1) ) { return true; }
  if ( is_hidden(q, p2) ) { return true; }
  return
    geom_traits().finite_edge_interior_conflict_2_object()(p1, p2, q, b);
}


template< class Gt, bool StoreHidden, class Agds >
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
finite_edge_interior_degenerated(const Face_handle& f, int i,
				 const Weighted_point_2& p, bool b) const
{
  if ( !is_infinite( f->mirror_vertex(i) ) ) {
    CGAL_precondition( is_infinite(f->vertex(i)) );

    Face_handle g = f->neighbor(i);
    int j = f->mirror_index(i);

    return finite_edge_interior_degenerated(g, j, p, b);
  }

  CGAL_precondition( is_infinite( f->mirror_vertex(i) ) );

  Weighted_point_2 p1 = f->vertex( ccw(i) )->point();
  Weighted_point_2 p2 = f->vertex(  cw(i) )->point();

  if ( is_infinite(f->vertex(i)) ) {
    return finite_edge_interior_degenerated(p1, p2, p, b);
  }

  Weighted_point_2 p3 = f->vertex(i)->point();
  return finite_edge_interior_degenerated(p1, p2, p3, p, b);
}

template< class Gt, bool StoreHidden, class Agds >
inline
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
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

  Weighted_point_2 p1 = v1->point();
  Weighted_point_2 p2 = v2->point();
  Weighted_point_2 p = v->point();

  if ( is_infinite(v3) ) {
    return finite_edge_interior_degenerated(p1, p2, p, b);
  }

  Weighted_point_2 p3 = v3->point();
  return finite_edge_interior_degenerated(p1, p2, p3, p, b);
}

template< class Gt, bool StoreHidden, class Agds >
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
infinite_edge_interior(const Weighted_point_2& p2,
		       const Weighted_point_2& p3,
		       const Weighted_point_2& p4,
		       const Weighted_point_2& q,
		       bool b) const
{
  if ( is_hidden(q, p2) ) { return true; }
  return
    geom_traits().infinite_edge_interior_conflict_2_object()(p2,p3,p4,q,b);
}

template< class Gt, bool StoreHidden, class Agds >
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
infinite_edge_interior(const Face_handle& f, int i,
		       const Weighted_point_2& p, bool b) const
{
  if ( !is_infinite( f->vertex(ccw(i)) ) ) {
    CGAL_precondition( is_infinite( f->vertex(cw(i)) ) );
    Face_handle g = f->neighbor(i);
    int j = f->mirror_index(i);

    return infinite_edge_interior(g, j, p, b);
  }

  CGAL_precondition( is_infinite( f->vertex(ccw(i)) ) );

  Weighted_point_2 p2 = f->vertex(  cw(i) )->point();
  Weighted_point_2 p3 = f->vertex(     i  )->point();
  Weighted_point_2 p4 = f->mirror_vertex(i)->point();

  return infinite_edge_interior(p2, p3, p4, p, b);
}


template< class Gt, bool StoreHidden, class Agds >
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
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

  Weighted_point_2 p2 = v2->point();
  Weighted_point_2 p3 = v3->point();
  Weighted_point_2 p4 = v4->point();
  Weighted_point_2 p = v->point();

  return infinite_edge_interior(p2, p3, p4, p, b);
}




template< class Gt, bool StoreHidden, class Agds >
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
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


template< class Gt, bool StoreHidden, class Agds >
bool
Apollonius_graph_2<Gt,StoreHidden,Agds>::
edge_interior(const Face_handle& f, int i,
	      const Weighted_point_2& p, bool b) const
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

template< class Gt, bool StoreHidden, class Agds >
typename Apollonius_graph_2<Gt,StoreHidden,Agds>::Conflict_type
Apollonius_graph_2<Gt,StoreHidden,Agds>::
finite_edge_conflict_type_degenerated(const Weighted_point_2& p1,
				      const Weighted_point_2& p2,
				      const Weighted_point_2& q) const
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


//----------------------------------------------------------------------
// methods for disk removal
//----------------------------------------------------------------------


template< class Gt, bool StoreHidden, class Agds >
inline
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
remove_first(Vertex_handle v)
{
  Delaunay_graph::remove_first(v);
}

template< class Gt, bool StoreHidden, class Agds >
inline
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
remove_second(Vertex_handle v)
{
  Delaunay_graph::remove_second(v);
}

template< class Gt, bool StoreHidden, class Agds >
inline
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
remove_third(Vertex_handle v)
{
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

  this->_tds.remove_dim_down( &(*v) );
}


template< class Gt, bool StoreHidden, class Agds >
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle(NULL) );
  CGAL_triangulation_precondition( !is_infinite(v) );

  // find a neighbor of v to use for point location of hidden sites to
  // be re-inserted
  Vertex_handle vnear(NULL);
  if ( StoreHidden ) {
    if ( number_of_vertices() > 10 ) {
      Vertex_circulator vc_start = v->incident_vertices();
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

  Weighted_point_list wp_list;
  typename Vertex_base::Hidden_weighted_point_iterator wpit;
  for (wpit = v->hidden_weighted_points_begin();
       wpit != v->hidden_weighted_points_end(); ++wpit) {
    wp_list.push_back(*wpit);
  }

  int n = number_of_vertices();
  if ( n == 1 ) {
    remove_first(v);
  } else if ( n == 2 ) {
    remove_second(v);
  } else if ( n == 3 ) {
    remove_third(v);
  } else {
    int degree = v->degree();
    if ( degree == 2 ) {
      remove_degree_2(v);
    } else if ( degree == 3 ) {
      remove_degree_3(v);
    } else {
      remove_degree_d_vertex(v);
    }
  }

  Weighted_point_less_than_comparator less_than(geom_traits());
  std::sort(wp_list.begin(), wp_list.end(), less_than);
  for (unsigned int i = 0; i < wp_list.size(); i++) {
    vnear = insert(wp_list[i], vnear);
  }
}


template< class Gt, bool StoreHidden, class Agds >
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
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
  
  Apollonius_graph_2<Gt,StoreHidden,Agds> ag_small;

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
      if ( &(*vh_small) != NULL ) {
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
  CGAL_triangulation_precondition( v->degree() == 3 );

  this->_tds.remove_degree_3( &(*v), NULL);

  for (unsigned int i = 0; i < num_fe; i++) {
    delete flipped_edges[i];
  }
}


template< class Gt, bool StoreHidden, class Agds >
void
Apollonius_graph_2<Gt,StoreHidden,Agds>::
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

CGAL_END_NAMESPACE
