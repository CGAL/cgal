// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>


// class implementation continued
//=================================

namespace CGAL {

//====================================================================
//====================================================================
//                   CONSTRUCTORS
//====================================================================
//====================================================================

// copy constructor
template<class Gt, class ST, class D_S, class LTag>
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
Segment_Delaunay_graph_2(const Segment_Delaunay_graph_2& other)
  : DG(other.geom_traits())
{
  Segment_Delaunay_graph_2&
    non_const_other = const_cast<Segment_Delaunay_graph_2&>(other);
  copy(non_const_other);
  CGAL_postcondition( is_valid() );
}

// assignment operator
template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Self&
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
operator=(const Self& other)
{
  if ( this != &other ) {
    Segment_Delaunay_graph_2&
      non_const_other = const_cast<Segment_Delaunay_graph_2&>(other);
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

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_first(const Storage_site_2& ss, const Point_2& )
{
  CGAL_precondition( number_of_vertices() == 0 );

  Vertex_handle v = this->_tds.insert_second();
  v->set_site(ss);
  return v;
}

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_second(const Storage_site_2& ss, const Point_2& p)
{
  CGAL_precondition( number_of_vertices() == 1 );
  // p0 is actually a point
  Site_2 p0 = finite_vertices_begin()->site();
  // MK: change the equality test between points by the functor in
  // geometric traits
  Site_2 tp = Site_2::construct_site_2(p);
  if ( same_points(tp,p0) ) {
    // merge info of identical sites
    merge_info(finite_vertices_begin(), ss);
    return finite_vertices_begin();
  }

  return create_vertex_dim_up(ss);
}

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_third(const Storage_site_2& ss, const Point_2& p)
{
  Site_2 t = Site_2::construct_site_2(p);
  return insert_third(t, ss);
}

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_third(const Site_2& t, const Storage_site_2& ss)
{
  CGAL_precondition( number_of_vertices() == 2 );

  // p0 and p1 are actually points
  Vertex_handle v0 = finite_vertices_begin();
  Vertex_handle v1 = ++finite_vertices_begin();
  Site_2 t0 = v0->site();
  Site_2 t1 = v1->site();

  if ( same_points(t, t0) ) {
    // merge info of identical sites
    merge_info(v0, ss);
    return v0;
  }
  if ( same_points(t, t1) ) {
    // merge info of identical sites
    merge_info(v1, ss);
    return v1;
  }

  Vertex_handle v = create_vertex_dim_up(ss);

  Face_handle f(finite_faces_begin());

  Site_2 s1 = f->vertex(0)->site();
  Site_2 s2 = f->vertex(1)->site();
  Site_2 s3 = f->vertex(2)->site();

  Sign s12i3 = geom_traits().vertex_conflict_2_object()(s1, s2, s3);
  Sign s21i3 = geom_traits().vertex_conflict_2_object()(s2, s1, s3);

  CGAL_assertion(s12i3 != ZERO);
  CGAL_assertion(s21i3 != ZERO);

  if ( s12i3 != s21i3 ) {
    if ( s21i3 == NEGATIVE ) {
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


template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_third(const Storage_site_2& ss, Vertex_handle , Vertex_handle )
{
  CGAL_precondition( number_of_vertices() == 2 );

  //  this can only be the case if the first site is a segment
  CGAL_precondition( dimension() == 1 );

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

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_point(const Storage_site_2& ss, const Point_2& p, Vertex_handle vnear)
{
  size_type n = number_of_vertices();
  if ( n == 0 ) {
    return insert_first(ss, p);
  } else if ( n == 1 ) {
    return insert_second(ss, p);
  } else if ( n == 2 ) {
    return insert_third(ss, p);
  }

  Site_2 t = Site_2::construct_site_2(p);
  return insert_point(ss, t, vnear);
}


template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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
      // merge info of identical sites
      merge_info(vnearest, ss);
      return vnearest;
    }
  } else {
    CGAL_assertion( vnearest->is_segment() );
    CGAL_assertion( at_res != AT2::TOUCH_1 );
    CGAL_assertion( at_res != AT2::TOUCH_2 );
    CGAL_assertion( at_res == AT2::DISJOINT || at_res == AT2::INTERIOR );
    if ( at_res == AT2::INTERIOR ) {
      CGAL_assertion( t.is_input() );

      Vertex_triple vt = (this->*insert_exact_point_on_segment_ptr)(
          ss, t, vnearest);
      return vt.first;
    } else {
      // the point to be inserted does not belong to the interior of a
      // segment
      CGAL_assertion( at_res == AT2::DISJOINT );
    }
  }

  return insert_point2(ss, t, vnearest);
}

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_point2(const Storage_site_2& ss, const Site_2& t,
              Vertex_handle vnearest)
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

#ifndef CGAL_SDG_NO_FACE_MAP
  std::map<Face_handle,Sign> sign_map;
#endif

  do {
    Face_handle f(fc);

    s = incircle(f, t);

#ifdef CGAL_SDG_NO_FACE_MAP
    f->tds_data().set_incircle_sign(s);
#else
    sign_map[f] = s;
#endif

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

#ifdef CGAL_SDG_NO_FACE_MAP
      Sign s1 = e.first->tds_data().incircle_sign();
      Sign s2 = e.first->neighbor(e.second)->tds_data().incircle_sign();
#else
      Sign s1 = sign_map[e.first];
      Sign s2 = sign_map[e.first->neighbor(e.second)];
#endif

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

#ifndef CGAL_SDG_NO_FACE_MAP
    sign_map.clear();
#endif

    CGAL_assertion( interior_in_conflict );

    return insert_degree_2(e, ss);
  }


  // we are in conflict with a Voronoi vertex; start from that and
  // find the entire conflict region and then repair the diagram
  List l;
#ifndef CGAL_SDG_NO_FACE_MAP
  Face_map fm;
#endif

  Triple<bool, Vertex_handle, Arrangement_type>
    vcross(false, Vertex_handle(), AT2::DISJOINT);

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT
  // REGION AND EXPANDS THE CONFLICT REGION.
  initialize_conflict_region(start_f, l);
#ifdef CGAL_SDG_NO_FACE_MAP
  expand_conflict_region(start_f, t, l, vcross);
#else
  expand_conflict_region(start_f, t, l, fm, sign_map, vcross);
#endif

  CGAL_assertion( !vcross.first );

  Vertex_handle v = create_vertex(ss);

#ifdef CGAL_SDG_NO_FACE_MAP
  retriangulate_conflict_region(v, l);
#else
  retriangulate_conflict_region(v, l, fm);
#endif

  return v;
}

//--------------------------------------------------------------------
// insertion of a point that lies on a segment
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Face_pair
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
find_faces_to_split(const Vertex_handle& v, const Site_2& t) const
{
  CGAL_precondition( v->is_segment() );

#ifndef CGAL_NO_ASSERTIONS
  {
    // count number of adjacent infinite faces
    Face_circulator fc = incident_faces(v);
    Face_circulator fc_start = fc;
    int n_inf = 0;
    do {
      if ( is_infinite(fc) ) { n_inf++; }
      fc++;
    } while ( fc != fc_start );
    CGAL_assertion( n_inf % 2 == 0 );
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
        os1 = oriented_side(v->site(), sv_ep, sitev_supp, t);
      } else {
        CGAL_assertion(  !is_infinite( ff1->vertex( cw_v) )  );
        CGAL_assertion( ff1->vertex( cw_v)->site().is_point() );
        sv_ep = ff1->vertex( cw_v)->site();
        os1 = oriented_side(sv_ep, v->site(), sitev_supp, t);
      }
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
        os2 = oriented_side(v->site(), sv_ep, sitev_supp, t);
      } else {
        CGAL_assertion(  !is_infinite( ff2->vertex( cw_v) )  );
        CGAL_assertion( ff2->vertex( cw_v)->site().is_point() );
        sv_ep = ff2->vertex( cw_v)->site();
        os2 = oriented_side(sv_ep, v->site(), sitev_supp, t);
      }
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
         os1 != ON_NEGATIVE_SIDE && os2 == ON_NEGATIVE_SIDE ) {
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

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_triple
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_exact_point_on_segment(const Storage_site_2& ss, const Site_2& t,
                              Vertex_handle v)
{
  // splits the segment site v->site() in two and inserts represented by t
  // on return the three vertices are, respectively, the vertex
  // corresponding to t and the two subsegments of v->site()

  CGAL_assertion( t.is_point() );
  CGAL_assertion( t.is_input() );

  Storage_site_2 ssitev = v->storage_site();

  CGAL_assertion( ssitev.is_segment() );

  Face_pair fpair = find_faces_to_split(v, t);

  boost::tuples::tuple<Vertex_handle,Vertex_handle,Face_handle,Face_handle>
    qq = this->_tds.split_vertex(v, fpair.first, fpair.second);

  // now I need to update the sites for vertices v1 and v2
  Vertex_handle v1 = boost::tuples::get<0>(qq); //qq.first;
  Storage_site_2 ssv1 = split_storage_site(ssitev, ss, true);
  v1->set_site( ssv1 );

  Vertex_handle v2 = boost::tuples::get<1>(qq); //qq.second;
  Storage_site_2 ssv2 = split_storage_site(ssitev, ss, false);
  v2->set_site( ssv2 );

  Face_handle qqf = boost::tuples::get<2>(qq); //qq.third;
  Vertex_handle vsx =
    this->_tds.insert_in_edge(qqf, cw(qqf->index(v1)));

  vsx->set_site(ss);
  // merge info of point and segment; the point lies on the segment
  merge_info(vsx, ssitev);

  return Vertex_triple(vsx, v1, v2);
}

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_triple
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_point_on_segment(const Storage_site_2& ss, const Site_2& ,
                        Vertex_handle v, const Tag_true&)
{
  // splits the segment site v->site() in two and inserts the point of
  // intersection of t and v->site()
  // on return the three vertices are, respectively, the point of
  // intersection and the two subsegments of v->site()

  Storage_site_2 ssitev = v->storage_site();
  Storage_site_2 ssx = st_.construct_storage_site_2_object()(ss, ssitev);

  Face_pair fpair = find_faces_to_split(v, ssx.site());

  boost::tuples::tuple<Vertex_handle,Vertex_handle,Face_handle,Face_handle>
    qq = this->_tds.split_vertex(v, fpair.first, fpair.second);

  // now I need to update the sites for vertices v1 and v2
  Vertex_handle v1 = boost::tuples::get<0>(qq); //qq.first;
  Vertex_handle v2 = boost::tuples::get<1>(qq); //qq.second;

  Storage_site_2 ssv1 =
    st_.construct_storage_site_2_object()(ssitev, ss, true);

  Storage_site_2 ssv2 =
    st_.construct_storage_site_2_object()(ssitev, ss, false);

  v1->set_site( ssv1 );
  v2->set_site( ssv2 );

  Face_handle qqf = boost::tuples::get<2>(qq); //qq.third;
  Vertex_handle vsx =
    this->_tds.insert_in_edge(qqf, cw(qqf->index(v1)));

  vsx->set_site(ssx);

  return Vertex_triple(vsx, v1, v2);
}

//--------------------------------------------------------------------
// insertion of a segment
//--------------------------------------------------------------------
template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_segment(const Storage_site_2& ss, const Site_2& t, Vertex_handle vnear)
{
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( t.is_input() );

  if ( is_degenerate_segment(t) ) {
    Storage_site_2 ss_src = ss.source_site();
    convert_info(ss_src, ss, true);
    return insert_point(ss_src, t.source(), vnear);
  }

  Storage_site_2 ss_src = ss.source_site();
  convert_info(ss_src, ss, true);
  Storage_site_2 ss_trg = ss.target_site();
  convert_info(ss_trg, ss, false);

  Vertex_handle v0 = insert_point( ss_src, t.source(), vnear );
  CGAL_assertion( is_valid() );
  Vertex_handle v1 = insert_point( ss_trg, t.target(), v0 );
  CGAL_assertion( is_valid() );

  if ( number_of_vertices() == 2 ) {
    return insert_third(ss, v0, v1);
  }

  return insert_segment_interior(t, ss, v0);
}


template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_segment_interior(const Site_2& t, const Storage_site_2& ss,
                        Vertex_handle vnearest)
{
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( number_of_vertices() > 2 );

  CGAL_assertion( vnearest != Vertex_handle() );

  // find the first conflict

  // first look if there are intersections...
  Vertex_circulator vc = incident_vertices(vnearest);
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
        // merge info of identical items
        merge_info(vv, ss);
        return vv;
      } else if ( at_res == AT2::CROSSING ) {
        Intersections_tag itag;
        return insert_intersecting_segment(ss, t, vv, itag);
      } else if ( at_res == AT2::TOUCH_11_INTERIOR_1 ) {
        Vertex_handle vp = second_endpoint_of_segment(vv);
        Storage_site_2 ssvp = vp->storage_site();
        Storage_site_2 sss = split_storage_site(ss, ssvp, false);

        Storage_site_2 sss1 = split_storage_site(ss, ssvp, true);
        // merge the info of the first (common) subsegment
        merge_info(vv, sss1);
        // merge the info of the (common) splitting endpoint
        merge_info(vp, ss);

        return insert_segment_interior(sss.site(), sss, vp);
      } else if ( at_res == AT2::TOUCH_12_INTERIOR_1 ) {
        Vertex_handle vp = first_endpoint_of_segment(vv);
        Storage_site_2 ssvp = vp->storage_site();
        Storage_site_2 sss = split_storage_site(ss, ssvp, true);

        /*Storage_site_2 sss1 =*/
  split_storage_site(ss, ssvp, false);
        // merge the info of the second (common) subsegment
        //        merge_info(vv, sss);
        // merge the info of the (common) splitting endpoint
        merge_info(vp, ss);

        return insert_segment_interior(sss.site(), sss, vp);
      } else {
        // this should never be reached; the only possible values for
        // at_res are DISJOINT, CROSSING, TOUCH_11_INTERIOR_1
        // and TOUCH_12_INTERIOR_1
        CGAL_error();
      }
    } else {
      CGAL_assertion( vv->is_point() );
      if ( at_res == AT2::INTERIOR ) {
        Storage_site_2 ssvv = vv->storage_site();
        if ( ssvv.is_input() ) {
          Storage_site_2 ss1 = split_storage_site(ss, ssvv, true);
          Storage_site_2 ss2 = split_storage_site(ss, ssvv, false);
          // merge the info of the splitting point and the segment
          merge_info(vv, ss);
          insert_segment_interior(ss1.site(), ss1, vv);
          return insert_segment_interior(ss2.site(), ss2, vv);
        } else {
          // this should never be reached; the only possible values for
          // at_res are DISJOINT and INTERIOR
          CGAL_error();
        }
      }
    }
    ++vc;
  } while ( vc != vc_start );

  // first look for conflict with vertex
  Face_circulator fc_start = incident_faces(vnearest);
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;

#ifndef CGAL_SDG_NO_FACE_MAP
  std::map<Face_handle,Sign> sign_map;
#endif

  do {
    Face_handle f(fc);

    s = incircle(f, t);

#ifdef CGAL_SDG_NO_FACE_MAP
    f->tds_data().set_incircle_sign(s);
#else
    sign_map[f] = s;
#endif

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
#ifndef CGAL_SDG_NO_FACE_MAP
  Face_map fm;
#endif

  Triple<bool, Vertex_handle, Arrangement_type>
    vcross(false, Vertex_handle(), AT2::DISJOINT);

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT
  // REGION AND EXPANDS THE CONFLICT REGION.
  initialize_conflict_region(start_f, l);
#ifdef CGAL_SDG_NO_FACE_MAP
  expand_conflict_region(start_f, t, l, vcross);
#else
  expand_conflict_region(start_f, t, l, fm, sign_map, vcross);
#endif

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
        Storage_site_2 ss1 = split_storage_site(ss, ssvv, true);
        Storage_site_2 ss2 = split_storage_site(ss, ssvv, false);
        // merge the info of the splitting point and the segment
        merge_info(vcross.second, ss);
        insert_segment_interior(ss1.site(), ss1, vcross.second);
        return insert_segment_interior(ss2.site(), ss2, vcross.second);
      } else {
        // this should never be reached; the only possible values for
        // vcross.third are CROSSING, INTERIOR and DISJOINT
        CGAL_error();
      }
    }
  }

  // no intersecting segment has been found; we insert the segment as
  // usual...
  Vertex_handle v = create_vertex(ss);

#ifdef CGAL_SDG_NO_FACE_MAP
  retriangulate_conflict_region(v, l);
#else
  retriangulate_conflict_region(v, l, fm);
#endif

  return v;
}


//--------------------------------------------------------------------
// insertion of an intersecting segment
//--------------------------------------------------------------------
template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& /* ss */,
                                     const Site_2& /* t */, Vertex_handle /* v */,
                                     Tag_false)
{
  print_error_message();
  return Vertex_handle();
}

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
                                     const Site_2& t, Vertex_handle v,
                                     Tag_true tag)
{
  CGAL_precondition( t.is_segment() && v->is_segment() );

  const Storage_site_2& ssitev = v->storage_site();
  Site_2 sitev = ssitev.site();

  if ( same_segments(t, sitev) ) {
    merge_info(v, ss);
    return v;
  }

  Vertex_triple vt = (this->*insert_point_on_segment_ptr)(ss, t, v, tag);

  Vertex_handle vsx = vt.first;

  Storage_site_2 ss3 = st_.construct_storage_site_2_object()(ss, ssitev, true);
  Storage_site_2 ss4 = st_.construct_storage_site_2_object()(ss, ssitev, false);
  Site_2 s3 = ss3.site();
  Site_2 s4 = ss4.site();

  insert_segment_interior(s3, ss3, vsx);
  insert_segment_interior(s4, ss4, vsx);
  return vsx;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// helper methods for insertion (find conflict region)
//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
initialize_conflict_region(const Face_handle& f, List& l)
{


  l.clear();
  for (int i = 0; i < 3; i++) {
    l.push_back(sym_edge(f, i));
  }
}


template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
check_unregistered_face(const Face_handle& n,
                        const Site_2& t,
                        List& l,
#ifndef CGAL_SDG_NO_FACE_MAP
                        Face_map& fm,
#endif
                        Triple<bool, Vertex_handle, Arrangement_type>& vcross)
{
  for (int j = 0; j < 3; ++j)
  {
    Vertex_handle vf = n->vertex(j);

    if ( is_infinite(vf) )
      continue;

    Arrangement_type at_res = arrangement_type(t, vf);

    CGAL_assertion( vcross.third == AT2::DISJOINT ||
                    vcross.third == AT2::CROSSING ||
                    vcross.third == AT2::INTERIOR );

    if ( vf->is_segment() )
    {
      CGAL_assertion( at_res != AT2::IDENTICAL );
      CGAL_assertion( at_res != AT2::TOUCH_11_INTERIOR_1 );
      CGAL_assertion( at_res != AT2::TOUCH_12_INTERIOR_1 );

      if ( at_res == AT2::CROSSING )
      {
        vcross.first = true;
        vcross.second = vf;
        vcross.third = AT2::CROSSING;
        l.clear();
#ifdef CGAL_SDG_NO_FACE_MAP
        fhc_.clear();
#else
        fm.clear();
#endif
        return true;
      }
      else // at_res != AT2::CROSSING
      {
        CGAL_assertion ( at_res == AT2::DISJOINT ||
                         at_res == AT2::TOUCH_1 ||
                         at_res == AT2::TOUCH_2 ||
                         at_res == AT2::TOUCH_11 ||
                         at_res == AT2::TOUCH_12 ||
                         at_res == AT2::TOUCH_21 ||
                         at_res == AT2::TOUCH_22 );
        // we do nothing in these cases
      }
    }
    else // ! vf->is_segment()
    {
      CGAL_assertion( vf->is_point() );
      if ( at_res == AT2::INTERIOR )
      {
        vcross.first = true;
        vcross.second = vf;
        vcross.third = AT2::INTERIOR;
        l.clear();
#ifdef CGAL_SDG_NO_FACE_MAP
        fhc_.clear();
#else
        fm.clear();
#endif
        return true;
      }
    }
  }

  return false;
}

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
expand_conflict_region(const Face_handle& in_f,
                       const Site_2& t,
                       List& l,
#ifndef CGAL_SDG_NO_FACE_MAP
                       Face_map& fm,
                       std::map<Face_handle, Sign>& sign_map,
#endif
                       Triple<bool, Vertex_handle, Arrangement_type>& vcross)
{
  std::stack<Face_handle> face_stack;
  face_stack.push(in_f);

  while(!face_stack.empty())
  {
    // Stop as soon as intersecting segments are found
    if ( vcross.first )
      break;

    const Face_handle curr_f = face_stack.top();
    face_stack.pop();

#ifdef CGAL_SDG_NO_FACE_MAP
    if ( curr_f->tds_data().is_in_conflict() )
      continue;
#else
    if ( fm.find(curr_f) != fm.end() )
      continue;
#endif

    // setting fm[f] to true means that the face has been reached and
    // that the face is available for recycling. If we do not want the
    // face to be available for recycling we must set this flag to false.
#ifdef CGAL_SDG_NO_FACE_MAP
    curr_f->tds_data().mark_in_conflict();
    fhc_.push_back(curr_f);
#else
    fm[curr_f] = true;
#endif

//    CGAL_assertion( fm.find(curr_f) != fm.end() );

    for (int i = 0; i < 3; ++i)
    {
      Face_handle n = curr_f->neighbor(i);

#ifdef CGAL_SDG_NO_FACE_MAP
      bool face_registered = n->tds_data().is_in_conflict();
#else
      bool face_registered = (fm.find(n) != fm.end());
#endif

      if ( !face_registered )
      {
#ifdef CGAL_SDG_NO_FACE_MAP
        if(check_unregistered_face(n, t, l, vcross))
#else
        if(check_unregistered_face(n, t, l, fm, vcross))
#endif
        {
          CGAL_assertion(vcross.first);
          break; // intersecting segments were found
        }
      }

      Sign s = incircle(n, t);

#ifdef CGAL_SDG_NO_FACE_MAP
      n->tds_data().set_incircle_sign(s);

      Sign s_f = curr_f->tds_data().incircle_sign();
#else
      sign_map[n] = s;

      Sign s_f = sign_map[curr_f];
#endif

      if ( s == POSITIVE )
        continue;
      if ( s != s_f )
        continue;
      if ( face_registered )
        continue;

      bool interior_in_conflict = edge_interior(curr_f, i, t, s);
      if ( !interior_in_conflict )
        continue;

      Edge e = sym_edge(curr_f, i);
      CGAL_assertion( l.is_in_list(e) );

      int j = this->_tds.mirror_index(curr_f, i);
      Edge e_before = sym_edge(n, ccw(j));
      Edge e_after = sym_edge(n, cw(j));

      if ( !l.is_in_list(e_before) )
        l.insert_before(e, e_before);

      if ( !l.is_in_list(e_after) )
        l.insert_after(e, e_after);

      l.remove(e);

      face_stack.push(n);
    } // neighbor for-loop
  }
}


//--------------------------------------------------------------------
// retriangulate conflict region
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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


template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_list
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
add_bogus_vertices(List& l)
{
#if defined(USE_INPLACE_LIST) && defined(CGAL_SDG_NO_FACE_MAP)
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

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
remove_bogus_vertices(Vertex_list& vl)
{
  while ( vl.size() > 0 ) {
    Vertex_handle v = vl.front();
    vl.pop_front();
    remove_degree_2(v);
  }
}


template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
#ifdef CGAL_SDG_NO_FACE_MAP
retriangulate_conflict_region(Vertex_handle v, List& l)
#else
retriangulate_conflict_region(Vertex_handle v, List& l,
                              Face_map& fm)
#endif

{
  // 1. add the bogus vetrices
  Vertex_list dummy_vertices = add_bogus_vertices(l);

  // 2. repair the face pointers...
  Edge e_start = l.front();
  Edge eit = e_start;
  do {
    CGAL_assertion_code(Edge esym =) sym_edge(eit);
    Face_handle f = eit.first;
    int k = eit.second;
    CGAL_assertion( !l.is_in_list(esym) );
#ifdef CGAL_SDG_NO_FACE_MAP
    CGAL_assertion( !f->tds_data().is_in_conflict() );
#else
    CGAL_assertion( fm.find(f) == fm.end() );
#endif
    f->vertex(ccw(k))->set_face(f);
    f->vertex( cw(k))->set_face(f);
    eit = l.next(eit);
  } while ( eit != e_start );

  // 3. copy the edge list to a vector of edges and clear the edge list
  // MK:: here I actually need to copy the edges to an std::list<Edge>, or
  // even better add iterators to the list of type List
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

  // 5. remove the bogus vertices
  remove_bogus_vertices(dummy_vertices);

  // 6. remove the unused faces
#ifdef CGAL_SDG_NO_FACE_MAP
  typename std::vector<Face_handle>::iterator it;
  for (it = fhc_.begin(); it != fhc_.end(); ++it) {
    (*it)->tds_data().clear();
    this->_tds.delete_face( *it );
  }

  fhc_.clear();
#else
  typename Face_map::iterator it;
  for (it = fm.begin(); it != fm.end(); ++it) {
    Face_handle fh = (*it).first;
    this->_tds.delete_face(fh);
  }

  fm.clear();
#endif

  // 7. DONE!!!!
}

//====================================================================
//====================================================================
//                   METHODS FOR REMOVAL
//====================================================================
//====================================================================

template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
is_star(const Vertex_handle& v) const
{
  CGAL_precondition( v->storage_site().is_point() );

  Vertex_circulator vc_start = incident_vertices(v);
  Vertex_circulator vc = vc_start;
  Storage_site_2 p = v->storage_site();

  size_type count = 0;
  do {
    Storage_site_2 ss = vc->storage_site();
    if ( ss.is_segment() && is_endpoint_of_segment(p, ss) ) {
      ++count;
      //      if ( count == 3 ) { break; }
      if ( count == 3 ) { return true; }
    }
    ++vc;
  } while ( vc != vc_start );

  return false;
}


template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
is_linear_chain(const Vertex_handle& v0, const Vertex_handle& v1,
                const Vertex_handle& v2) const
{
  Site_2 tt[3] = { v0->site(), v1->site(), v2->site() };

  if ( tt[1].is_point() &&
       tt[0].is_segment() &&
       tt[2].is_segment() &&
       is_endpoint_of_segment(tt[1], tt[0]) &&
       is_endpoint_of_segment(tt[1], tt[2]) ) {
    typename Geom_traits::Equal_2 are_equal = geom_traits().equal_2_object();
    Site_2 s_end[2];
    if (  are_equal( tt[1], tt[0].source_site() )  ) {
      s_end[0] = tt[0].target_site();
    } else {
      s_end[0] = tt[0].source_site();
    }

    if (  are_equal( tt[1], tt[2].source_site() )  ) {
      s_end[1] = tt[2].target_site();
    } else {
      s_end[1] = tt[2].source_site();
    }

    typename Geom_traits::Orientation_2 orientation =
      geom_traits().orientation_2_object();

    return orientation(s_end[0], s_end[1], tt[1]) == COLLINEAR;
  }
  return false;
}


template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
is_flippable(const Face_handle& f, int i) const
{
  CGAL_assertion( !is_infinite(f->vertex( cw(i) )) );

  Vertex_handle v_other = f->vertex( ccw(i) );
  Vertex_handle v0 = f->vertex( i );
  Vertex_handle v1 = this->_tds.mirror_vertex( f, i );

  if ( is_infinite(v_other) || is_infinite(v0) || is_infinite(v1) ) {
    return false;
  }

  Vertex_handle v = f->vertex( cw(i) );

  Storage_site_2 ss = v->storage_site();
  Storage_site_2 ss_other = v_other->storage_site();
  if ( ss_other.is_point() && ss.is_segment() &&
       is_endpoint_of_segment(ss_other,        ss) && is_star(v_other) ) {
    return false;
  }

  if ( is_linear_chain(v0, v_other, v1) ) { return false; }

  return (v0 != v1) && is_degenerate_edge(f, i);
}


template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
minimize_degree(const Vertex_handle& v)
{
  CGAL_precondition ( degree(v) > 3 );

  Face_circulator fc_start = incident_faces(v);
  Face_circulator fc = incident_faces(v);
  bool found(false);
  do {
    Face_handle f = Face_handle(fc);
    int i = ccw( f->index(v) );

    CGAL_assertion( f->vertex( cw(i) ) == v );

    if ( is_flippable(f,i) ) {
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

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
equalize_degrees(const Vertex_handle& v, Self& small_d,
                 std::map<Vertex_handle,Vertex_handle>& vmap,
                 List& l) const
{
  size_type deg = degree(v);
  CGAL_assertion( l.size() <= deg );
  if ( l.size() == deg ) { return; }
#if 0
  std::cerr << "size of l  : " << l.size() << std::endl;
  std::cerr << "degree of v: " << deg << std::endl;
#endif

  //  typedef std::map<Edge,Edge>  Edge_map;
  // maps edges on the boundary of the conflict region from the small
  // diagram to edges of the star of v in the small diagram
  //  Edge_map emap;

  Edge e_small_start = l.front();
  Edge e_small = e_small_start;
  bool found;
  Edge e_small_begin;
  Edge e_large_begin;
  do {
    found = false;
    Vertex_handle v_sml_src = e_small.first->vertex(cw(e_small.second));
    Vertex_handle v_sml_trg = e_small.first->vertex(ccw(e_small.second));

    // first we find a first edge in common
    Face_circulator fc_start = incident_faces(v);
    Face_circulator fc = fc_start;

    do {
      int id = fc->index(v);
      Vertex_handle v_lrg_src = fc->vertex(ccw(id));
      Vertex_handle v_lrg_trg = fc->vertex(cw(id));
      if ( vmap[v_sml_src] == v_lrg_src && vmap[v_sml_trg] == v_lrg_trg ) {
        found = true;
        e_large_begin = Edge(fc, id);
        e_small_begin = e_small;
        break;
      }
    } while ( ++fc != fc_start );
    if ( found ) { break; }
    e_small = l.next(e_small);
  } while ( e_small != e_small_start );

  CGAL_assertion( found );

  Face_circulator fc_start = incident_faces(v, e_large_begin.first);
  Face_circulator fc = fc_start;
  e_small = e_small_begin;
  do {
    int id = fc->index(v);
    Vertex_handle vsml_src = e_small.first->vertex(cw(e_small.second));
    Vertex_handle vsml_trg = e_small.first->vertex(ccw(e_small.second));
    Vertex_handle vlrg_src = fc->vertex(ccw(id));
    Vertex_handle vlrg_trg = fc->vertex(cw(id));
    if ( vmap[vsml_src] != vlrg_src || vmap[vsml_trg] != vlrg_trg ) {
      Edge e_small_prev = l.previous(e_small);
      std::cerr << "size of l: " << l.size() << std::endl;
      l.remove(e_small);

      std::cerr << "size of l: " << l.size() << std::endl;

      Edge e_small_new = small_d.flip(e_small);
      Edge e_small_new_sym = small_d.sym_edge(e_small_new);
      Face_handle f1 = e_small_new.first;
      Face_handle f2 = e_small_new_sym.first;

      if ( f2->vertex(e_small_new_sym.second) == vsml_src ) {
        std::swap(f1, f2);
        std::swap(e_small_new, e_small_new_sym);
        CGAL_assertion( f1->vertex(e_small_new.second) == vsml_src );
        CGAL_assertion( f2->vertex(e_small_new_sym.second) == vsml_trg );
      }

      Edge to_list1(f1, cw(e_small_new.second));
      Edge to_list2(f2, ccw(e_small_new_sym.second));

      e_small = small_d.sym_edge(to_list1);

      l.insert_after(e_small_prev, e_small);
      std::cerr << "size of l: " << l.size() << std::endl;
      l.insert_after(e_small, small_d.sym_edge(to_list2));
      std::cerr << "size of l: " << l.size() << std::endl;
    } else {
      e_small = l.next(e_small);
      ++fc;
    }
    CGAL_assertion( l.size() <= deg );
  } while ( fc != fc_start );

#if 0
  std::cerr << "size of l  : " << l.size() << std::endl;
  std::cerr << "degree of v: " << deg << std::endl;
#endif

#if !defined(CGAL_NO_ASSERTIONS) && !defined(NDEBUG)
  // we go around the boundary of the conflict region verify that all
  // edges are there
  CGAL_assertion( l.size() == degree(v) );
  e_small = e_small_begin;
  fc_start = incident_faces(v, e_large_begin.first);
  fc = fc_start;
  do {
    int id = fc->index(v);
    Vertex_handle vsml_src = e_small.first->vertex(cw(e_small.second));
    Vertex_handle vsml_trg = e_small.first->vertex(ccw(e_small.second));
    Vertex_handle vlrg_src = fc->vertex(ccw(id));
    Vertex_handle vlrg_trg = fc->vertex(cw(id));
    CGAL_assertion(vmap[vsml_src] == vlrg_src && vmap[vsml_trg] == vlrg_trg );

    // go to next edge
    e_small = l.next(e_small);
  } while ( ++fc != fc_start );
#endif
}

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
expand_conflict_region_remove(const Face_handle& in_f,
                              const Site_2& t,
                              List& l,
                              Face_map& fm,
                              Sign_map& sign_map)
{
  std::stack<Face_handle> face_stack;
  face_stack.push(in_f);

  while(!face_stack.empty())
  {
    const Face_handle curr_f = face_stack.top();
    face_stack.pop();

    if ( fm.find(curr_f) != fm.end() )
      continue;

    // setting fm[curr_f] to true means that the face has been reached and
    // that the face is available for recycling. If we do not want the
    // face to be available for recycling we must set this flag to
    // false.
    fm[curr_f] = true;

    //  CGAL_assertion( fm.find(f) != fm.end() );

    for (int i = 0; i < 3; ++i)
    {
      Face_handle n = curr_f->neighbor(i);

      bool face_registered = (fm.find(n) != fm.end());

      Sign s = incircle(n, t);

      sign_map[n] = s;

      Sign s_f = sign_map[curr_f];

      if ( s == POSITIVE )
        continue;
      if ( s != s_f )
        continue;
      if ( face_registered )
        continue;

      bool interior_in_conflict = edge_interior(curr_f, i, t, s);
      if ( !interior_in_conflict )
        continue;

      Edge e = sym_edge(curr_f, i);
      CGAL_assertion( l.is_in_list(e) );

      int j = this->_tds.mirror_index(curr_f, i);
      Edge e_before = sym_edge(n, ccw(j));
      Edge e_after = sym_edge(n, cw(j));

      if ( !l.is_in_list(e_before) )
        l.insert_before(e, e_before);

      if ( !l.is_in_list(e_after) )
        l.insert_after(e, e_after);

      l.remove(e);

      face_stack.push(n);
    } // neighbor for-loop
  }
}


template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
find_conflict_region_remove(const Vertex_handle& v,
                            const Vertex_handle& vnearest,
                            List& l, Face_map& fm, Sign_map&)
{
  CGAL_precondition( vnearest != Vertex_handle() );
  Storage_site_2 ss = v->storage_site();
  Site_2 t = ss.site();

  // find the first conflict

  // first look for conflict with vertex
  Face_circulator fc_start = incident_faces(vnearest);
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;

  Sign_map sign_map;

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

  CGAL_assertion( s == NEGATIVE );

  // we are not in conflict with a Voronoi vertex, so we have to
  // be in conflict with the interior of a Voronoi edge
  if ( s != NEGATIVE ) {
    Edge_circulator ec_start = incident_edges(vnearest);
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

    l.push_back(e);
    l.push_back(sym_edge(e));
    return;
  }

  initialize_conflict_region(start_f, l);
  expand_conflict_region_remove(start_f, t, l, fm, sign_map);
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::size_type
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
count_faces(const List& l) const
{
  std::vector<Face_handle> flist;
  get_faces(l, std::back_inserter(flist));
  return flist.size();
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
fill_hole(const Self& small_d, const Vertex_handle& v, const List& l,
          std::map<Vertex_handle,Vertex_handle>& vmap)
{
#if 0
  std::cerr << "size of l  : " << l.size() << std::endl;
  std::cerr << "degree of v: " << degree(v) << std::endl;
#endif

  typedef std::map<Edge,Edge>  Edge_map;
  // maps edges on the boundary of the conflict region from the small
  // diagram to edges of the star of v in the small diagram
  Edge_map emap;

  Edge e_sml_sym = sym_edge(l.front());

  Vertex_handle v_sml_src = e_sml_sym.first->vertex(ccw(e_sml_sym.second));
  Vertex_handle v_sml_trg = e_sml_sym.first->vertex(cw(e_sml_sym.second));

  // first we find a first edge in common
  Face_circulator fc_start = incident_faces(v);
  Face_circulator fc = fc_start;
  Face_circulator fc_begin;
  CGAL_assertion_code( bool found = false; )
  do {
    int id = fc->index(v);
    Vertex_handle v_lrg_src = fc->vertex(ccw(id));
    Vertex_handle v_lrg_trg = fc->vertex(cw(id));
    if ( vmap[v_sml_src] == v_lrg_src && vmap[v_sml_trg] == v_lrg_trg ) {
      CGAL_assertion_code( found = true; )
      fc_begin = fc;
      break;
    }
  } while ( ++fc != fc_start );
  CGAL_assertion( found );

  // container of faces to delete
  std::list<Face_handle> to_delete;

  // next we go around the boundary of the conflict region and map all edges
  Edge e_small = l.front();
  fc_start = incident_faces(v, fc_begin);
  fc = fc_start;
  do {
    int id = fc->index(v);
#if !defined(CGAL_NO_ASSERTIONS) && !defined(NDEBUG)
    Vertex_handle vsml_src = e_small.first->vertex(cw(e_small.second));
    Vertex_handle vsml_trg = e_small.first->vertex(ccw(e_small.second));
    Vertex_handle vlrg_src = fc->vertex(ccw(id));
    Vertex_handle vlrg_trg = fc->vertex(cw(id));
    CGAL_assertion(vmap[vsml_src] == vlrg_src && vmap[vsml_trg] == vlrg_trg );
#endif
    // set mapping
    emap[e_small] = sym_edge(fc, id);
    // keep face for deletion
    to_delete.push_back(fc);
    // go to next edge
    e_small = l.next(e_small);
  } while ( ++fc != fc_start );


  // map the faces of the small diagram with the new faces of the
  // large diagram
  std::map<Face_handle,Face_handle> fmap;
  std::vector<Face_handle> f_small, f_large;

  small_d.get_faces(l, std::back_inserter(f_small));
  for (unsigned int i = 0; i < f_small.size(); i++) {
    Face_handle f = this->_tds.create_face();
    fmap[ f_small[i] ] = f;
    f_large.push_back(f);
  }

  CGAL_assertion( l.size() == degree(v) );
  CGAL_assertion( f_small.size() == f_large.size() );
  CGAL_assertion( f_small.size() == l.size() - 2 );
  CGAL_assertion( f_small.size() == count_faces(l) );

  // set the vertices for the new faces; also set the face for each
  // vertex
  for (typename std::vector<Face_handle>::iterator fit = f_small.begin();
       fit != f_small.end(); ++fit) {
    Face_handle ff_small = *fit;
    Face_handle f = fmap[ff_small];

    for (int i = 0; i < 3; ++i) {
      CGAL_assertion( vmap.find(ff_small->vertex(i)) != vmap.end() );
      f->set_vertex(i, vmap[ff_small->vertex(i)]);
      // we are setting the face for each vertex a lot of times; we
      // could possibly do it faster if we use the edges on the
      // boundary of the conflict region
      // in fact we may not even need it since the make_hole method
      // updates the faces of the vertices of the boundary of the
      // hole, and we do not introduce any new vertices
      f->vertex(i)->set_face(f);
    }
  }

  // set the neighbors for each face
  for (typename std::vector<Face_handle>::iterator fit = f_small.begin();
       fit != f_small.end(); ++fit) {
    Face_handle ff_small = *fit;

    for (int i = 0; i < 3; i++) {
      Face_handle f = fmap[ff_small];
      Face_handle n_small = ff_small->neighbor(i);
      if ( fmap.find(n_small) != fmap.end() ) {
        // this is one of the new faces
        f->set_neighbor(i, fmap[n_small]);
      } else {
        // otherwise it is one of the old faces outside the conflict
        // region; we have to use the edge map in this case
        Edge e_small_sym = small_d.sym_edge(ff_small, i);
        CGAL_assertion(emap.find(e_small_sym) != emap.end());

        Edge e_large = emap[e_small_sym];
        f->set_neighbor(i, e_large.first);
        e_large.first->set_neighbor(e_large.second, f);
      }
    }
  }

  // delete the unused faces and the vertex
  while ( !to_delete.empty() ) {
    delete_face(to_delete.front());
    to_delete.pop_front();
  }
  delete_vertex(v);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
remove_first(const Vertex_handle& v)
{
  Delaunay_graph::remove_first(v);
  return true;
}

template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
remove_second(const Vertex_handle& v)
{
  Delaunay_graph::remove_second(v);
  return true;
}

template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
remove_third(const Vertex_handle& v)
{
  if ( is_degree_2(v) ) {
    CGAL_assertion( v->storage_site().is_point() );
    Face_handle fh( incident_faces(v) );
    int i = fh->index(v);
    flip(fh, i);
  } else if ( degree(v) == 4 ) {
    Edge_circulator ec = incident_edges(v);
    for (int i = 0; i < 4; i++) {
      Edge e = *ec;
      Edge sym = sym_edge(e);
      if ( e.first->vertex(e.second) !=        sym.first->vertex(sym.second) ) {
        flip(e);
        break;
      }
      ++ec;
    }
  }

  this->_tds.remove_dim_down( v );
  return true;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
compute_small_diagram(const Vertex_handle& v, Self& small_d) const
{
  Vertex_circulator vc_start = incident_vertices(v);
  Vertex_circulator vc = vc_start;

  // insert all neighboring sites
  do {
    if ( !is_infinite(vc) ) {
      Site_2 t = vc->site();
      if ( t.is_input() ) {
        small_d.insert(t);
      } else if ( t.is_point() ) {
        small_d.insert(t.supporting_site(0));
        /*Vertex_handle vnear =*/ small_d.insert(t.supporting_site(1));
        //        vh_small = sdg_small.nearest_neighbor(t, vnear);
      } else {
        CGAL_assertion( t.is_segment() );
        /*Vertex_handle vnear =*/ small_d.insert(t.supporting_site());
        //        vh_small = sdg_small.nearest_neighbor(t, vnear);
      }
      //      CGAL_assertion( vh_small != Vertex_handle() );
      //      vmap[vh_small] = vh_large;
    }
    ++vc;
  } while ( vc != vc_start );
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
compute_vertex_map(const Vertex_handle& v, const Self& small_d,
                   std::map<Vertex_handle,Vertex_handle>& vmap) const
{
  Vertex_circulator vc_start = incident_vertices(v);
  Vertex_circulator vc = vc_start;
  Vertex_handle vh_large, vh_small;

  do {
    vh_large = Vertex_handle(vc);
    if ( is_infinite(vh_large) ) {
      vh_small = small_d.infinite_vertex();
      vmap[vh_small] = vh_large;
    } else {
#if !defined(CGAL_NO_ASSERTIONS) && !defined(NDEBUG)
      vh_small = Vertex_handle();
#endif
      Site_2 t = vc->site();
      if ( t.is_input() ) {
        if ( t.is_segment() ) {
          Vertex_handle vv = small_d.nearest_neighbor( t.source() );
          Vertex_circulator vvc_start = small_d.incident_vertices(vv);
          Vertex_circulator vvc = vvc_start;
          do {
            if ( small_d.same_segments(t, vvc) ) {
              vh_small = vvc;
              break;
            }
          } while ( ++vvc != vvc_start );
          CGAL_assertion( small_d.same_segments(t, vh_small) );
        } else {
          CGAL_assertion( t.is_point() );
          vh_small = small_d.nearest_neighbor( t.point() );
          CGAL_assertion( small_d.same_points(t, vh_small->site()) );
        }
      } else if ( t.is_segment() ) {
        Vertex_handle vv =
          small_d.nearest_neighbor( t.source_site(), Vertex_handle() );
        Vertex_circulator vvc_start = small_d.incident_vertices(vv);
        Vertex_circulator vvc = vvc_start;
        do {
          if ( small_d.same_segments(t, vvc) ) {
            vh_small = vvc;
            break;
          }
        } while ( ++vvc != vvc_start );
        CGAL_assertion( small_d.same_segments(t, vh_small) );
      } else {
        CGAL_assertion( t.is_point() );
        vh_small = small_d.nearest_neighbor( t, Vertex_handle() );
        CGAL_assertion( small_d.same_points(t, vh_small->site()) );
      }
      CGAL_assertion( vh_small != Vertex_handle() );
      vmap[vh_small] = vh_large;
    }
    ++vc;
  } while ( vc != vc_start );
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
remove_degree_d_vertex(const Vertex_handle& v)
{
#if 0
  Self sdg_small;
  compute_small_diagram(v, sdg_small);
  std::map<Vertex_handle,Vertex_handle> vmap;
  compute_vertex_map(v, sdg_small, vmap);

  // find nearest neighbor of v in small diagram
  Site_2 t = v->site();
  Vertex_handle vn;

  CGAL_assertion( t.is_input() );

  if ( t.is_point() ) {
    vn = sdg_small.nearest_neighbor( t.point() );
  } else {
    vn = sdg_small.nearest_neighbor( t.source() );
  }
  CGAL_assertion( vn != Vertex_handle() );

  List l;
  Face_map fm;
  Sign_map sign_map;

  sdg_small.find_conflict_region_remove(v, vn, l, fm, sign_map);

  fm.clear();
  sign_map.clear();

  equalize_degrees(v, sdg_small, vmap, l);

  fill_hole(sdg_small, v, l, vmap);

  l.clear();
  return;

#else
  minimize_degree(v);
  size_type deg = degree(v);
  if ( deg == 3 ) {
    remove_degree_3(v);
    return;
  }
  if ( deg == 2 ) {
    remove_degree_2(v);
    return;
  }

  Self sdg_small;
  compute_small_diagram(v, sdg_small);

  if ( sdg_small.number_of_vertices() <= 2 ) {
    CGAL_assertion( sdg_small.number_of_vertices() == 2 );
    CGAL_assertion( deg == 4 );
    Edge_circulator ec_start = incident_edges(v);
    Edge_circulator ec = ec_start;
    do {
      if ( is_infinite(*ec) ) { break; }
      ++ec;
    } while ( ec != ec_start );
    CGAL_assertion( is_infinite(ec) );
    flip(*ec);
    remove_degree_3(v);
    return;
  }

  std::map<Vertex_handle,Vertex_handle> vmap;
  compute_vertex_map(v, sdg_small, vmap);

  // find nearest neighbor of v in small diagram
  Site_2 t = v->site();
  Vertex_handle vn;

  CGAL_assertion( t.is_input() );

  // here we find a site in the small diagram that serves as a
  // starting point for finding all conflicts.
  // To do that we find the nearest neighbor of t if t is a point;
  // t is guarranteed to have a conflict with its nearest neighbor
  // If t is a segment, then one endpoint of t is enough; t is
  // guarranteed to have a conflict with the Voronoi edges around
  // this endpoint
  if ( t.is_point() ) {
    vn = sdg_small.nearest_neighbor( t.point() );
  } else {
    vn = sdg_small.nearest_neighbor( t.source() );
  }
  CGAL_assertion( vn != Vertex_handle() );

  List l;
  Face_map fm;
  Sign_map sign_map;

  sdg_small.find_conflict_region_remove(v, vn, l, fm, sign_map);

  fill_hole(sdg_small, v, l, vmap);

  l.clear();
  fm.clear();
  sign_map.clear();
#endif
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
remove_base(const Vertex_handle& v)
{
  Storage_site_2 ssv = v->storage_site();
  CGAL_precondition( ssv.is_input() );

  // first consider the case where we have up to 2 points
  size_type n = number_of_vertices();

  if ( n == 1 ) {
    return remove_first(v);
  } else if ( n == 2 ) {
    return remove_second(v);
  }

  // secondly check if the point to be deleted is adjacent to a segment
  if ( ssv.is_point() ) {
    Vertex_circulator vc_start = incident_vertices(v);
    Vertex_circulator vc = vc_start;

    do {
      Storage_site_2 ss = vc->storage_site();
      if ( ss.is_segment() && is_endpoint_of_segment(ssv, ss) ) {
        return false;
      }
      ++vc;
    } while ( vc != vc_start );
  }

  // now do the deletion
  if ( n == 3 ) {
    return remove_third(v);
  }

  size_type deg = degree(v);
  if ( deg == 2 ) {
    remove_degree_2(v);
  } else if ( deg == 3 ) {
    remove_degree_3(v);
  } else {
    remove_degree_d_vertex(v);
  }

  return true;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
remove(const Vertex_handle& v)
{
  CGAL_precondition( !is_infinite(v) );
  const Storage_site_2& ss = v->storage_site();

  if ( !ss.is_input() ) { return false; }

  Point_handle h1, h2;
  bool is_point = ss.is_point();
  if ( is_point ) {
    h1 = ss.point();
  } else {
    CGAL_assertion( ss.is_segment() );
    h1 = ss.source_of_supporting_site();
    h2 = ss.target_of_supporting_site();
  }

  bool success = remove_base(v);

  if ( success ) {
    // unregister the input site
    if ( is_point ) {
      unregister_input_site( h1 );
    } else {
      unregister_input_site( h1, h2 );
    }
  }
  return success;
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
template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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

  //  if ( start_vertex == nullptr ) { return start_vertex; }

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

template<class Gt, class ST, class D_S, class LTag>
Sign
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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


template<class Gt, class ST, class D_S, class LTag>
Sign
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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
template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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

  Site_2 t1 = f->vertex( ccw(i) )->site();
  Site_2 t2 = f->vertex(  cw(i) )->site();

  if ( is_infinite(f->vertex(i)) ) {
    return finite_edge_interior(t1, t2, q, sgn);
  }

  Site_2 t3 = f->vertex(i)->site();
  return finite_edge_interior(t1, t2, t3, q, sgn);
}

template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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

template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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

  Site_2 t2 = f->vertex(  cw(i) )->site();
  Site_2 t3 = f->vertex(     i  )->site();
  Site_2 t4 = this->_tds.mirror_vertex(f, i)->site();

  return infinite_edge_interior(t2, t3, t4, q, sgn);
}


template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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




template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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


template<class Gt, class ST, class D_S, class LTag>
bool
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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


template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Arrangement_type
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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
template<class Gt, class ST, class D_S, class LTag>
Object
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
primal(const Edge e) const
{
  CGAL_precondition( !is_infinite(e) );

  if ( this->dimension() == 1 ) {
    Site_2 p = (e.first)->vertex(cw(e.second))->site();
    Site_2 q = (e.first)->vertex(ccw(e.second))->site();

    return make_object(construct_sdg_bisector_2_object()(p,q));
  }

  // dimension == 2
  // none of the two adjacent faces is infinite
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    Site_2 p = (e.first)->vertex( ccw(e.second) )->site();
    Site_2 q = (e.first)->vertex(  cw(e.second) )->site();
    Site_2 r = (e.first)->vertex(     e.second  )->site();
    Site_2 s = this->_tds.mirror_vertex(e.first, e.second)->site();
    return construct_sdg_bisector_segment_2_object()(p,q,r,s);
  }

  // both of the adjacent faces are infinite
  if ( is_infinite(e.first) &&
       is_infinite(e.first->neighbor(e.second)) )  {
    Site_2 p = (e.first)->vertex(cw(e.second))->site();
    Site_2 q = (e.first)->vertex(ccw(e.second))->site();
    return make_object(construct_sdg_bisector_2_object()(p,q));
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
  Site_2 p = ee.first->vertex( ccw(ee.second) )->site();
  Site_2 q = ee.first->vertex(  cw(ee.second) )->site();
  Site_2 r = ee.first->vertex(     ee.second  )->site();

  return make_object(construct_sdg_bisector_ray_2_object()(p,q,r));
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// validity test method
//--------------------------------------------------------------------
//--------------------------------------------------------------------
template<class Gt, class ST, class D_S, class LTag>
bool Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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


template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
print_error_message() const
{
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(int, once, 0);
  if(once == 0){
    ++once;
    std::cerr << std::endl;
    std::cerr << "ATTENTION:" << std::endl;
    std::cerr << "A segment-segment intersection was found."
              << std::endl;
    std::cerr << "The Segment_Delaunay_graph_2 class is not configured"
              << " to handle this situation." << std::endl;
    std::cerr << "Please look at the documentation on how to handle"
              << " this behavior." << std::endl;
    std::cerr << std::endl;
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// the copy method
//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Storage_site_2
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
copy_storage_site(const Storage_site_2& ss_other, Handle_map& hm,
                  const Tag_false&)
{
  if ( ss_other.is_segment() ) {
    Point_handle p0 = hm[ ss_other.source_of_supporting_site() ];
    Point_handle p1 = hm[ ss_other.target_of_supporting_site() ];
    return st_.construct_storage_site_2_object()(p0, p1);
  } else {
    Point_handle p0 = hm[ ss_other.point() ];
    return st_.construct_storage_site_2_object()(p0);
  }
}

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Storage_site_2
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
copy_storage_site(const Storage_site_2& ss_other, Handle_map& hm,
                  const Tag_true&)
{
  if ( ss_other.is_segment() ) {
    if ( ss_other.is_input() ) {
      Point_handle p0 = hm[ ss_other.source_of_supporting_site() ];
      Point_handle p1 = hm[ ss_other.target_of_supporting_site() ];
      return st_.construct_storage_site_2_object()(p0, p1);
    } else if ( ss_other.is_input(0) ) {
      Point_handle p0 = hm[ ss_other.source_of_supporting_site() ];
      Point_handle p1 = hm[ ss_other.target_of_supporting_site() ];
      Point_handle p4 = hm[ ss_other.source_of_crossing_site(1) ];
      Point_handle p5 = hm[ ss_other.target_of_crossing_site(1) ];
      return st_.construct_storage_site_2_object()(p0, p1, p4, p5, true);
    } else if ( ss_other.is_input(1) ) {
      Point_handle p0 = hm[ ss_other.source_of_supporting_site() ];
      Point_handle p1 = hm[ ss_other.target_of_supporting_site() ];
      Point_handle p2 = hm[ ss_other.source_of_crossing_site(0) ];
      Point_handle p3 = hm[ ss_other.target_of_crossing_site(0) ];
      return st_.construct_storage_site_2_object()(p0, p1, p2, p3, false);
    } else {
      Point_handle p0 = hm[ ss_other.source_of_supporting_site() ];
      Point_handle p1 = hm[ ss_other.target_of_supporting_site() ];
      Point_handle p2 = hm[ ss_other.source_of_crossing_site(0) ];
      Point_handle p3 = hm[ ss_other.target_of_crossing_site(0) ];
      Point_handle p4 = hm[ ss_other.source_of_crossing_site(1) ];
      Point_handle p5 = hm[ ss_other.target_of_crossing_site(1) ];
      return st_.construct_storage_site_2_object()(p0, p1, p2, p3, p4, p5);
    }
  } else {
    if ( ss_other.is_input() ) {
      Point_handle p0 = hm[ ss_other.point() ];
      return st_.construct_storage_site_2_object()(p0);
    } else {
      Point_handle p2 = hm[ ss_other.source_of_supporting_site(0) ];
      Point_handle p3 = hm[ ss_other.target_of_supporting_site(0) ];
      Point_handle p4 = hm[ ss_other.source_of_supporting_site(1) ];
      Point_handle p5 = hm[ ss_other.target_of_supporting_site(1) ];
      return st_.construct_storage_site_2_object()(p2, p3, p4, p5);
    }
  }
}

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
copy(Segment_Delaunay_graph_2& other)
{
  // copy the insert_on method pointers
  insert_point_on_segment_ptr = other.insert_point_on_segment_ptr;
  insert_exact_point_on_segment_ptr =
    other.insert_exact_point_on_segment_ptr;

  // first copy the point container
  pc_ = other.pc_;

  // copy storage traits
  st_ = other.st_;

  // first create a map between the old point handles and the new ones
  Handle_map hm;

  Point_handle it_other = other.pc_.begin();
  Point_handle it_this = pc_.begin();
  for (; it_other != other.pc_.end(); ++it_other, ++it_this) {
    hm.insert( Point_handle_pair(it_other, it_this) );
  }

  copy(other, hm);
}

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
copy(Segment_Delaunay_graph_2& other, Handle_map& hm)
{
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

  // then copy the diagram
  DG::operator=(other);

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
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// getting endpoints of segments
//--------------------------------------------------------------------
//--------------------------------------------------------------------

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
first_endpoint_of_segment(const Vertex_handle& v) const
{
  CGAL_assertion( v->is_segment() );
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

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
file_output(std::ostream& os, const Storage_site_2& t,
            Point_handle_mapper& P) const
{
  CGAL_precondition( t.is_defined() );

  if ( t.is_point() ) {
    // 0 for point
    os << 0;
    if ( is_ascii(os) ) { os << ' '; }
    if ( t.is_input() ) {
      // 0 for input
      os << 0;
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.point()];
    } else {
      // 1 for non-input
      os << 1;
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_supporting_site(0)];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_supporting_site(0)];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_supporting_site(1)];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_supporting_site(1)];
    }
  } else { // t is a segment
    // 1 for segment
    os << 1;
    if ( is_ascii(os) ) { os << ' '; }
    if ( t.is_input() ) {
      // 0 for input
      os << 0;
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_supporting_site()];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_supporting_site()];
    } else if ( t.is_input(0) ) {
      // 1 for input source
      os << 1;
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_supporting_site()];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_supporting_site()];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_crossing_site(1)];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_crossing_site(1)];
    } else if ( t.is_input(1) ) {
      // 2 for input target
      os << 2;
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_supporting_site()];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_supporting_site()];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_crossing_site(0)];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_crossing_site(0)];
    } else {
      // 3 for non-input src & trg
      os << 3;
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_supporting_site()];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_supporting_site()];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_crossing_site(0)];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_crossing_site(0)];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.source_of_crossing_site(1)];
      if ( is_ascii(os) ) { os << ' '; }
      os << P[t.target_of_crossing_site(1)];
    }
  }
}

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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

template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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

  std::map<Vertex_handle,size_type> V;
  std::map<Face_handle,size_type> F;

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



template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>::
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


} //namespace CGAL

// EOF
