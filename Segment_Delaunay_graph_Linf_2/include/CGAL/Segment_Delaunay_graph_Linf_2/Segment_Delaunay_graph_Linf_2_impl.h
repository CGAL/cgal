// Copyright (c) 2015  Università della Svizzera italiana.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

namespace CGAL {

// print face in standard output
template<class Gt, class ST, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag>::
face_output(const char *before, Face_handle f,
            const char *after) const
{
  std::cout << before;
  if (is_infinite(f->vertex(0))) {
    std::cout << " infv";
  } else {
    std::cout << ' ' << f->vertex(0)->site();
  }
  if (is_infinite(f->vertex(1))) {
    std::cout << " infv";
  } else {
    std::cout << ' ' << f->vertex(1)->site();
  }
  if (is_infinite(f->vertex(2))) {
    std::cout << " infv";
  } else {
    std::cout << ' ' << f->vertex(2)->site();
  }
  std::cout << after;
}


//--------------------------------------------------------------------
// insertion of a point that lies on a segment
//--------------------------------------------------------------------

// tiebreak with face oriented side
template<class Gt, class ST, class D_S, class LTag>
Oriented_side
Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag>::
oriented_side_face_tiebreak(
    Face_handle f, const Vertex_handle& v, const Site_2 & sitev,
    const Site_2 & sitev_supp, const Site_2 & t) const
{
  Oriented_side os;
  if ( is_infinite(f) ) {
    int id_v = f->index(v);
    int cw_v = this->cw( id_v );
    int ccw_v = this->ccw( id_v );

    Site_2 sv_ep;
    if ( is_infinite( f->vertex(cw_v) ) ) {
      CGAL_assertion(  !is_infinite( f->vertex(ccw_v) )  );
      CGAL_assertion( f->vertex(ccw_v)->site().is_point() );
      sv_ep = f->vertex(ccw_v)->site();
      CGAL_SDG_DEBUG( std::cout <<
          "debug sv_ep = " << sv_ep << std::endl;
          );
      os = oriented_side(sitev, sv_ep, sitev_supp, t, t.point());
    } else {
      CGAL_assertion(  !is_infinite( f->vertex( cw_v) )  );
      CGAL_assertion( f->vertex( cw_v)->site().is_point() );
      sv_ep = f->vertex( cw_v)->site();
      CGAL_SDG_DEBUG( std::cout <<
          "debug sv_ep = " << sv_ep << std::endl;
          );
      os = oriented_side(sv_ep, sitev, sitev_supp, t, t.point());
    }
  } else {
    os = oriented_side(f->vertex(0)->site(),
        f->vertex(1)->site(),
        f->vertex(2)->site(),
        sitev_supp, t, t.point());
  }
  return os;
}


template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag>::Face_pair
Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag>::
find_faces_to_split(const Vertex_handle& v, const Site_2& t,
    unsigned int & flips_nop,
    unsigned int & flips_pon) const
{
  CGAL_precondition( v->is_segment() );

  CGAL_SDG_DEBUG(std::cout << "debug: impl find_faces_to_split enters "
      << "with v=" << v->site() << " t=" << t << std::endl; );

#ifndef CGAL_NO_ASSERTIONS
  {
    // count number of adjacent infinite faces
    Face_circulator fc = incident_faces(v);
    Face_circulator fc_start = fc;
    int n_inf = 0;
    int n_faces = 0;
    do {
      if ( is_infinite(fc) ) { n_inf++; }
      fc++;
      n_faces++;
    } while ( fc != fc_start );

    CGAL_SDG_DEBUG( std::cout
        << "debug: impl find_faces_to_split: n_faces="
        << n_faces << std::endl; );
    CGAL_SDG_DEBUG( std::cout
        << "debug: impl find_faces_to_split: n_inf="
        << n_inf << std::endl; );

    // philaris: the (n_inf in {0,2,4}) assertion is very
    // specific to L2; it seems that in general, this should
    // be an even number (for Linf, at least)
    //CGAL_assertion( n_inf == 0 || n_inf == 2 || n_inf == 4 );

    // philaris: tocheck
    CGAL_assertion( n_inf % 2 == 0 );
  }
#endif

  Face_circulator fc1 = incident_faces(v);
  Face_circulator fc2 = fc1; ++fc2;
  Face_circulator fc_start = fc1;
  Face_handle f1, f2;
  bool found_f1 = false, found_f2 = false;
  bool os0_fc_start = false;
  bool first_found_f1 = false;
  CGAL_assertion_code( bool first_found_f2 = false; );
  Face_handle f0_no, f0_op, f0_po, f0_on;

#ifndef CGAL_NO_ASSERTIONS
  bool is_nop(false);
  bool is_pon(false);

  bool is_set_f0_on(false);
  bool is_set_f0_op(false);
  bool is_set_f0_no(false);
  bool is_set_f0_po(false);
#endif

  bool connect_all_nop(false);
  bool connect_all_pon(false);

  unsigned int count_nop_zeros(0);
  unsigned int count_pon_zeros(0);
  unsigned int count_zeros(0);

  Site_2 sitev = v->site();
  Site_2 sitev_supp = v->site().supporting_site();

  CGAL_SDG_DEBUG( std::cout
      << "debug: impl find_faces_to_split: vsite="
      << v->site() << " with supporting site="
      << v->site().supporting_site() << std::endl; );

  do {

    CGAL_SDG_DEBUG(std::cout << "debug impl start (found_f1 found_f2)= "
                   << found_f1 << " " << found_f2 << std::endl;);

    Face_handle ff1(fc1), ff2(fc2);
    Oriented_side os1, os2;

#ifdef CGAL_SDG_VERBOSE
    face_output("debug ff1=[", ff1, "]\n");
    face_output("debug ff2=[", ff2, "]\n");
#endif

    if ( is_infinite(ff1) ) {
      int id_v = ff1->index(v);
      int cw_v = this->cw( id_v );
      int ccw_v = this->ccw( id_v );

      Site_2 sv_ep;
      if ( is_infinite( ff1->vertex(cw_v) ) ) {
        CGAL_assertion(  !is_infinite( ff1->vertex(ccw_v) )  );
        CGAL_assertion( ff1->vertex(ccw_v)->site().is_point() );
        sv_ep = ff1->vertex(ccw_v)->site();
        CGAL_SDG_DEBUG( std::cout <<
            "debug sv_ep = " << sv_ep << std::endl;
            );
        os1 = oriented_side(sitev, sv_ep, sitev_supp, t);
      } else {
        CGAL_assertion(  !is_infinite( ff1->vertex( cw_v) )  );
        CGAL_assertion( ff1->vertex( cw_v)->site().is_point() );
        sv_ep = ff1->vertex( cw_v)->site();
        CGAL_SDG_DEBUG( std::cout <<
            "debug sv_ep = " << sv_ep << std::endl;
            );
        os1 = oriented_side(sv_ep, sitev, sitev_supp, t);
      }
    } else {
      os1 = oriented_side(fc1->vertex(0)->site(),
                          fc1->vertex(1)->site(),
                          fc1->vertex(2)->site(),
                          sitev_supp, t);
    }

    CGAL_SDG_DEBUG( std::cout <<
        "debug os1 = " << os1 << std::endl; );

    if ( is_infinite(ff2) ) {
      int id_v = ff2->index(v);
      int cw_v = this->cw( id_v );
      int ccw_v = this->ccw( id_v );

      Site_2 sv_ep;
      if ( is_infinite( ff2->vertex(cw_v) ) ) {
        CGAL_assertion(  !is_infinite( ff2->vertex(ccw_v) )  );
        CGAL_assertion( ff2->vertex(ccw_v)->site().is_point() );
        sv_ep = ff2->vertex(ccw_v)->site();
        os2 = oriented_side(sitev, sv_ep, sitev_supp, t);
      } else {
        CGAL_assertion(  !is_infinite( ff2->vertex( cw_v) )  );
        CGAL_assertion( ff2->vertex( cw_v)->site().is_point() );
        sv_ep = ff2->vertex( cw_v)->site();
        os2 = oriented_side(sv_ep, sitev, sitev_supp, t);
      }
    } else {
      os2 = oriented_side(fc2->vertex(0)->site(),
                          fc2->vertex(1)->site(),
                          fc2->vertex(2)->site(),
                          sitev_supp, t);
    }

    CGAL_SDG_DEBUG( std::cout <<
        "debug os2 = " << os2 << std::endl; );

#ifdef CGAL_SDG_VERBOSE
    face_output("debug ff1=[", ff1, "] has os1=");
    std::cout << os1 << std::endl;
    face_output("debug ff2=[", ff2, "] has os2=");
    std::cout << os2 << std::endl;
#endif

    if (os1 == ON_ORIENTED_BOUNDARY) {
      ++count_zeros;
      if (fc_start == fc1) {
        CGAL_SDG_DEBUG(std::cout << "debug impl start face has os=0"
            << std::endl;);
        os0_fc_start = true;
      }
      if (os2 == ON_NEGATIVE_SIDE) {
        CGAL_assertion(! is_set_f0_on);
        CGAL_assertion_code( is_set_f0_on = true ) ;
        f0_on = ff1;
        CGAL_SDG_DEBUG(std::cout << "debug impl found f0_on "
          << std::endl;);
      } else if (os2 == ON_POSITIVE_SIDE) {
        CGAL_assertion(! is_set_f0_op);
        CGAL_assertion_code( is_set_f0_op = true ) ;
        f0_op = ff1;
        CGAL_SDG_DEBUG(std::cout << "debug impl found f0_op "
          << std::endl;);
      }
    }

    if ( !found_f1 &&
         os1 != ON_POSITIVE_SIDE && os2 == ON_POSITIVE_SIDE ) {
      f1 = ff2;
      found_f1 = true;
      count_nop_zeros = count_zeros;
      count_zeros = 0;
      CGAL_SDG_DEBUG(std::cout << "debug impl found_f1 set to true "
          << "with os1=" << os1 << std::endl;);
#ifndef CGAL_NO_ASSERTIONS
      if (os1 == ON_ORIENTED_BOUNDARY) {
        CGAL_assertion(! is_nop);
        is_nop = true;
      }
#endif
      if ( !found_f2 ) {
        first_found_f1 = true;
      }
    }

    // philaris: change to be more symmetric:
    // before it was:
    //   os1 == ON_POSITIVE_SIDE && os2 != ON_POSITIVE_SIDE
    if ( !found_f2 &&
         os1 != ON_NEGATIVE_SIDE && os2 == ON_NEGATIVE_SIDE ) {
      f2 = ff2;
      found_f2 = true;
      count_pon_zeros = count_zeros;
      count_zeros = 0;
      CGAL_SDG_DEBUG(std::cout << "debug impl found_f2 set to true "
          << "with os1=" << os1 << std::endl;);
#ifndef CGAL_NO_ASSERTIONS
      if (os1 == ON_ORIENTED_BOUNDARY) {
        CGAL_assertion(! is_pon);
        is_pon = true;
      }
#endif
      CGAL_assertion_code(
      if ( !found_f1 ) {
        first_found_f2 = true;
      }
      );
    }

    if (os2 == ON_ORIENTED_BOUNDARY) {
      if (os1 == ON_NEGATIVE_SIDE) {
        CGAL_assertion(! is_set_f0_no);
        CGAL_assertion_code( is_set_f0_no = true );
        f0_no = ff2;
        CGAL_SDG_DEBUG(std::cout << "debug impl found f0_no "
          << std::endl;);
      } else if (os1 == ON_POSITIVE_SIDE) {
        CGAL_assertion(! is_set_f0_po);
        CGAL_assertion_code( is_set_f0_po = true );
        f0_po = ff2;
        CGAL_SDG_DEBUG(std::cout << "debug impl found f0_po "
          << std::endl;);
      }
    }

    CGAL_SDG_DEBUG(std::cout << "debug impl end (found_f1 found_f2)= "
                   << found_f1 << " " << found_f2 << std::endl;);

    if ( found_f1 && found_f2 && (!os0_fc_start) ) { break; }

    ++fc1, ++fc2;
  } while ( fc_start != fc1 );

  CGAL_SDG_DEBUG(std::cout << "debug after loop count_zeros="
      << count_zeros << std::endl;);

#ifndef CGAL_NO_ASSERTIONS
  if (is_nop) {
    CGAL_assertion( f1 != f0_op );
    CGAL_assertion( f1 != f0_no );
    CGAL_assertion( f2 != f0_no );
    CGAL_assertion( f2 != f0_op );
  }
  if (is_pon) {
    CGAL_assertion( f2 != f0_on );
    CGAL_assertion( f2 != f0_po );
    CGAL_assertion( f1 != f0_po );
    CGAL_assertion( f1 != f0_on );
  }
  if (is_nop && is_pon) {
    CGAL_assertion( f0_op != f0_on );
    CGAL_assertion( f0_no != f0_po );
  }
#endif

  // add missing counted zeros
  if (count_zeros > 0) {
    if (first_found_f1) {
      count_nop_zeros = count_nop_zeros + count_zeros;
    } else { // here first_found_f2
      CGAL_assertion(first_found_f2);
      count_pon_zeros = count_pon_zeros + count_zeros;
    }
  }
#ifdef CGAL_SDG_VERBOSE
  std::cout << "debug count_nop_zeros=" << count_nop_zeros << std::endl;
  std::cout << "debug count_pon_zeros=" << count_pon_zeros << std::endl;
#endif

  Oriented_side os_f0_no(ON_ORIENTED_BOUNDARY);
  Oriented_side os_f0_op(ON_ORIENTED_BOUNDARY);
  if (count_nop_zeros > 0) {
    os_f0_no = oriented_side_face_tiebreak(f0_no, v, sitev, sitev_supp, t);
    if (count_nop_zeros > 1) {
      os_f0_op = oriented_side_face_tiebreak(f0_op, v, sitev, sitev_supp, t);
    } else {
      os_f0_op = os_f0_no;
    }
    if (os_f0_op == ON_NEGATIVE_SIDE) {
      // do nothing
    } else if (os_f0_no == ON_POSITIVE_SIDE) {
      CGAL_SDG_DEBUG(std::cerr << "debug f1 updated"
          << std::endl;);
      f1 = f0_no;
    } else {
      connect_all_nop = true;
    }
  }
#ifdef CGAL_SDG_VERBOSE
  std::cout << "debug os_f0_no=" << os_f0_no
    << " os_f0_op=" << os_f0_op
    << " connect_all_nop=" << connect_all_nop << std::endl;
#endif

  Oriented_side os_f0_po(ON_ORIENTED_BOUNDARY);
  Oriented_side os_f0_on(ON_ORIENTED_BOUNDARY);
  if (count_pon_zeros > 0) {
    os_f0_po = oriented_side_face_tiebreak(f0_po, v, sitev, sitev_supp, t);
    if (count_pon_zeros > 1) {
      os_f0_on = oriented_side_face_tiebreak(f0_on, v, sitev, sitev_supp, t);
    } else {
      os_f0_on = os_f0_po;
    }
    if (os_f0_on == ON_POSITIVE_SIDE) {
      // do nothing
    } else if (os_f0_po == ON_NEGATIVE_SIDE) {
      CGAL_SDG_DEBUG(std::cerr << "debug f2 updated"
          << std::endl;);
      f2 = f0_po;
    } else {
      connect_all_pon = true;
    }
  }
#ifdef CGAL_SDG_VERBOSE
  std::cout << "debug os_f0_po=" << os_f0_po
    << " os_f0_on=" << os_f0_on
    << " connect_all_pon=" << connect_all_pon << std::endl;
#endif

#ifdef CGAL_SDG_VERBOSE
  std::cout << "debug find_faces_to_split results" << std::endl;
  face_output("debug f1=[", f1, "]\n");
  face_output("debug f2=[", f2, "]\n");
  if (is_nop) {
    face_output("debug f0_no=[", f0_no, "]\n");
    face_output("debug f0_op=[", f0_op, "]\n");
    std::cout << "debug f0_no, f0_op "
      << ((f0_no == f0_op)? "equal" : "not equal") << std::endl;
  } else {
    std::cout << "debug f0_no=UNSET" << std::endl;
    std::cout << "debug f0_op=UNSET" << std::endl;
    std::cout << "debug f0_no, f0_op UNSET" << std::endl;
  }
  if (is_pon) {
    face_output("debug f0_po=[", f0_po, "]\n");
    face_output("debug f0_on=[", f0_on, "]\n");
    std::cout << "debug f0_po, f0_on "
      << ((f0_po == f0_on)? "equal" : "not equal") << std::endl;
  } else {
    std::cout << "debug f0_po=UNSET" << std::endl;
    std::cout << "debug f0_on=UNSET" << std::endl;
    std::cout << "debug f0_po, f0_on UNSET" << std::endl;
  }
#endif

  CGAL_assertion( found_f1 && found_f2 );
  CGAL_assertion( f1 != f2 );

  CGAL_assertion( is_nop == is_set_f0_no );
  CGAL_assertion( is_nop == is_set_f0_op );
  CGAL_assertion( is_pon == is_set_f0_po );
  CGAL_assertion( is_pon == is_set_f0_on );

  CGAL_assertion( first_found_f1 || first_found_f2);
  CGAL_assertion( ! (first_found_f1 && first_found_f2) );

  CGAL_assertion((! connect_all_nop) || (count_nop_zeros > 0));
  if (connect_all_nop) {
    flips_nop = count_nop_zeros;
  } else {
    flips_nop = 0;
  }

  CGAL_assertion((! connect_all_pon) || (count_pon_zeros > 0));
  if (connect_all_pon) {
    flips_pon = count_pon_zeros;
  } else {
    flips_pon = 0;
  }

  return Face_pair(f1, f2);
}

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag>::Vertex_triple
Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag>::
insert_exact_point_on_segment(const Storage_site_2& ss, const Site_2& t,
                              Vertex_handle v)
{
  // splits the segment site v->site() in two and inserts represented by t
  // on return the three vertices are, respectively, the vertex
  // corresponding to t and the two subsegments of v->site()

  CGAL_assertion( t.is_point() );
  CGAL_assertion( t.is_input() );

  unsigned int flips_nop, flips_pon;

  Storage_site_2 ssitev = v->storage_site();

  CGAL_assertion( ssitev.is_segment() );

  Face_pair fpair = find_faces_to_split(v, t, flips_nop, flips_pon);

#ifdef CGAL_SDG_VERBOSE
  std::cout << "debug ins exact flips_nop=" << flips_nop
    << " flips_pon=" << flips_pon
    << std::endl;
#endif

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

  Face_handle otherf;
  int f_i = -1;
  if (flips_nop > 0) {
    otherf = qqf->neighbor(qqf->index(v1));
    f_i = this->_tds.mirror_index(qqf, qqf->index(v1));
  }

  Face_handle otherg;
  int g_i = -1;
  if (flips_pon > 0) {
    Face_handle qqg = boost::tuples::get<3>(qq); //qq.fourth;
    otherg = qqg->neighbor(qqg->index(v2));
    g_i = this->_tds.mirror_index(qqg, qqg->index(v2));
  }

  Vertex_handle vsx =
    this->_tds.insert_in_edge(qqf, cw(qqf->index(v1)));

  if (flips_nop > 0) {
    CGAL_assertion( f_i != -1 );
#ifdef CGAL_SDG_VERBOSE
    std::cout << "debug flip_nop>0 otherf=";
    face_output("[", otherf, "]");
    std::cout << " at f_i=" << f_i << std::endl;
#endif
    unsigned int remaining_nop = flips_nop;
    Face_handle next_face = otherf;
    int next_i = f_i;
    while (remaining_nop > 0) {
      const Face_handle current_face = next_face;
      const int current_i = next_i;
      if (remaining_nop > 1) {
        next_face = current_face->neighbor(ccw(current_i));
        next_i = this->_tds.mirror_index(current_face, ccw(current_i));
      }
#ifdef CGAL_SDG_VERBOSE
      std::cout << "debug remaining_nop=" << remaining_nop;
      face_output(" flip curf=[", current_face, "]");
      std::cout << " at f_i=" << current_i << std::endl;
#endif
      this->_tds.flip(current_face, current_i);
      --remaining_nop;
    }
  }

  if (flips_pon > 0) {
    CGAL_assertion( g_i != -1 );
#ifdef CGAL_SDG_VERBOSE
    std::cout << "debug flip otherg=";
    face_output("[", otherg, "]");
    std::cout << " at g_i=" << g_i << std::endl;
#endif
    unsigned int remaining_pon = flips_pon;
    Face_handle next_face = otherg;
    int next_i = g_i;
    while (remaining_pon > 0) {
      const Face_handle current_face = next_face;
      const int current_i = next_i;
      if (remaining_pon > 1) {
        next_face = current_face->neighbor(ccw(current_i));
        next_i = this->_tds.mirror_index(current_face, ccw(current_i));
      }
#ifdef CGAL_SDG_VERBOSE
      std::cout << "debug remaining_pon=" << remaining_pon;
      face_output(" flip curg=[", current_face, "]");
      std::cout << " at g_i=" << current_i << std::endl;
#endif
      this->_tds.flip(current_face, current_i);
      --remaining_pon;
    }
  }

  vsx->set_site(ss);
  // merge info of point and segment; the point lies on the segment
  merge_info(vsx, ssitev);

#ifndef CGAL_NO_ASSERTIONS
  // check if vsx has exactly: 4 + flips_nop + flips_pon neighbors
  {
    // count number of adjacent faces of vsx
    Face_circulator fc = incident_faces(vsx);
    Face_circulator fc_start = fc;
    unsigned int n_faces = 0;
    do {
      fc++;
      n_faces++;
    } while ( fc != fc_start );

    CGAL_SDG_DEBUG( std::cout
        << "debug: impl insert_exact_point_on_segment: n_faces="
        << n_faces << std::endl; );
    CGAL_assertion( n_faces == 4 + flips_nop + flips_pon );
  }
#endif

  return Vertex_triple(vsx, v1, v2);
}

template<class Gt, class ST, class D_S, class LTag>
typename Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag>::Vertex_triple
Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag>::
insert_point_on_segment(const Storage_site_2& ss, const Site_2& ,
                        Vertex_handle v, const Tag_true&)
{
  // splits the segment site v->site() in two and inserts the point of
  // intersection of t and v->site()
  // on return the three vertices are, respectively, the point of
  // intersection and the two subsegments of v->site()

  unsigned int flips_nop(0), flips_pon(0);

  Storage_site_2 ssitev = v->storage_site();
  Storage_site_2 ssx =
    storage_traits().construct_storage_site_2_object()(ss, ssitev);

  #ifdef CGAL_SDG_VERBOSE
  Site_2 sitev = ssitev.site();
  std::cout << "debug insert_point_on_segment intsec="
    << ssx.site() << " v=" << sitev << std::endl ;
  #endif
  Face_pair fpair = find_faces_to_split(v, ssx.site(), flips_nop, flips_pon);

#ifdef CGAL_SDG_VERBOSE
  std::cout << "debug ins flips_nop=" << flips_nop
    << " flips_pon=" << flips_pon
    << std::endl;
#endif

  boost::tuples::tuple<Vertex_handle,Vertex_handle,Face_handle,Face_handle>
    qq = this->_tds.split_vertex(v, fpair.first, fpair.second);

  // now I need to update the sites for vertices v1 and v2
  Vertex_handle v1 = boost::tuples::get<0>(qq); //qq.first;
  Vertex_handle v2 = boost::tuples::get<1>(qq); //qq.second;

  Storage_site_2 ssv1 =
    storage_traits().construct_storage_site_2_object()(ssitev, ss, true);

  Storage_site_2 ssv2 =
    storage_traits().construct_storage_site_2_object()(ssitev, ss, false);

  v1->set_site( ssv1 );
  v2->set_site( ssv2 );

  Face_handle qqf = boost::tuples::get<2>(qq); //qq.third;

  Face_handle otherf;
  int f_i = -1;
  if (flips_nop > 0) {
    otherf = qqf->neighbor(qqf->index(v1));
    f_i = this->_tds.mirror_index(qqf, qqf->index(v1));
  }

  Face_handle otherg;
  int g_i = -1;
  if (flips_pon > 0) {
    Face_handle qqg = boost::tuples::get<3>(qq); //qq.fourth;
    otherg = qqg->neighbor(qqg->index(v2));
    g_i = this->_tds.mirror_index(qqg, qqg->index(v2));
  }

  Vertex_handle vsx =
    this->_tds.insert_in_edge(qqf, cw(qqf->index(v1)));

  if (flips_nop > 0) {
    CGAL_assertion( f_i != -1 );
#ifdef CGAL_SDG_VERBOSE
    std::cout << "debug flip_nop>0 otherf=";
    face_output("[", otherf, "]");
    std::cout << " at f_i=" << f_i << std::endl;
#endif
    unsigned int remaining_nop = flips_nop;
    Face_handle next_face = otherf;
    int next_i = f_i;
    while (remaining_nop > 0) {
      const Face_handle current_face = next_face;
      const int current_i = next_i;
      if (remaining_nop > 1) {
        next_face = current_face->neighbor(ccw(current_i));
        next_i = this->_tds.mirror_index(current_face, ccw(current_i));
      }
#ifdef CGAL_SDG_VERBOSE
      std::cout << "debug remaining_nop=" << remaining_nop;
      face_output(" flip curf=[", current_face, "]");
      std::cout << " at f_i=" << current_i << std::endl;
#endif
      this->_tds.flip(current_face, current_i);
      --remaining_nop;
    }
  }

  if (flips_pon > 0) {
    CGAL_assertion( g_i != -1 );
#ifdef CGAL_SDG_VERBOSE
    std::cout << "debug flip otherg=";
    face_output("[", otherg, "]");
    std::cout << " at g_i=" << g_i << std::endl;
#endif
    unsigned int remaining_pon = flips_pon;
    Face_handle next_face = otherg;
    int next_i = g_i;
    while (remaining_pon > 0) {
      const Face_handle current_face = next_face;
      const int current_i = next_i;
      if (remaining_pon > 1) {
        next_face = current_face->neighbor(ccw(current_i));
        next_i = this->_tds.mirror_index(current_face, ccw(current_i));
      }
#ifdef CGAL_SDG_VERBOSE
      std::cout << "debug remaining_pon=" << remaining_pon;
      face_output(" flip curg=[", current_face, "]");
      std::cout << " at g_i=" << current_i << std::endl;
#endif
      this->_tds.flip(current_face, current_i);
      --remaining_pon;
    }
  }

  vsx->set_site(ssx);

#ifndef CGAL_NO_ASSERTIONS
  // check if vsx has exactly: 4 + flips_nop + flips_pon neighbors
  {
    // count number of adjacent faces of vsx
    Face_circulator fc = incident_faces(vsx);
    Face_circulator fc_start = fc;
    unsigned int n_faces = 0;
    do {
      fc++;
      n_faces++;
    } while ( fc != fc_start );

    CGAL_SDG_DEBUG( std::cout
        << "debug: impl insert_point_on_segment: n_faces="
        << n_faces << std::endl; );
    CGAL_assertion( n_faces == 4 + flips_nop + flips_pon );
  }
#endif

  return Vertex_triple(vsx, v1, v2);
}

} //namespace CGAL

// EOF
