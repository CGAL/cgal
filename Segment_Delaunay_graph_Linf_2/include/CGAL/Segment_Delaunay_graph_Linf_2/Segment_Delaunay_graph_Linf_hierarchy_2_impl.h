// Copyright (c) 2015  Università della Svizzera italiana.
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

// class implementation continued
//=================================

namespace CGAL {


//===========================================================================

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  constructors
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
init_hierarchy(const Geom_traits& gt)
{
  hierarchy[0] = this; 
  for(unsigned int i = 1; i < sdg_hierarchy_2__maxlevel; ++i) {
    hierarchy[i] = new Base(gt);
  }
}

template<class Gt, class ST, class STag, class D_S, class LTag>
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
Segment_Delaunay_graph_Linf_hierarchy_2(const Gt& gt)
  : Base(gt)
{ 
  init_hierarchy(gt);
}


// copy constructor duplicates vertices and faces
template<class Gt, class ST, class STag, class D_S, class LTag>
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
Segment_Delaunay_graph_Linf_hierarchy_2
(const Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag> &sdg)
    : Base(sdg.geom_traits())
{ 
  // create an empty triangulation to be able to delete it !
  init_hierarchy(sdg.geom_traits());
  copy(sdg);
} 

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  destructor
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class ST, class STag, class D_S, class LTag>
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>:: 
~Segment_Delaunay_graph_Linf_hierarchy_2()
{
  clear();
  for(unsigned int i = 1; i < sdg_hierarchy_2__maxlevel; ++i){ 
    delete hierarchy[i];
  }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// assignment operator
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class ST, class STag, class D_S, class LTag>
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag> &
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
operator=(const Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag> &sdg)
{
  if ( this != &sdg ) {
    copy(sdg);
  }
  return *this;
}

//====================================================================
//====================================================================
//                   METHODS FOR INSERTION
//====================================================================
//====================================================================

//--------------------------------------------------------------------
// insertion of a point
//--------------------------------------------------------------------

template<class Gt, class ST, class STag, class D_S, class LTag>
std::pair<bool,int>
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
insert_point(const Point_2& p, const Storage_site_2& ss, int level,
	     Vertex_handle* vertices)
{
  CGAL_SDG_DEBUG(std::cout << "debug hier insert_point " << p << " at level "
                 << level << std::endl;);

  CGAL_precondition( level != UNDEFINED_LEVEL );

  int new_level = level;

  bool lies_on_seg = false;

  Vertex_handle vertex;
  Vertex_handle vs1, vs2;
  Vertex_handle vnear[sdg_hierarchy_2__maxlevel];

  Site_2 t = Site_2::construct_site_2(p);

  nearest_neighbor(t, vnear, false);

  size_type n = hierarchy[0]->number_of_vertices();
  if ( n > 2 ) {
    Arrangement_type at_res = this->arrangement_type(t, vnear[0]);

    CGAL_assertion( at_res == AT2::DISJOINT ||
		    at_res == AT2::INTERIOR ||
		    at_res == AT2::IDENTICAL );

    if ( vnear[0]->is_point() ) {
      CGAL_SDG_DEBUG(std::cout << "debug hier insert_point nearest is point"
                     << std::endl;);
      if ( at_res == AT2::IDENTICAL ) {
	vertex = vnear[0];
	merge_info(vertex, ss);
      } else {
	vertex = hierarchy[0]->insert_point2(ss, t, vnear[0]);
      }
    } else { // nearest neighbor is a segment
      CGAL_assertion( vnear[0]->is_segment() );
      CGAL_assertion( at_res == AT2::DISJOINT ||
		      at_res == AT2::INTERIOR );

      if ( at_res == AT2::INTERIOR ) {
	CGAL_assertion( t.is_input() );
	lies_on_seg = true;

	int vnear_level = find_level(vnear[0]);

	// I need to find the level of the nearest neighbor that t
	// lies on and update the level of t
	if ( new_level < vnear_level ) {
	  new_level = vnear_level;
	}

	Vertex_triple vt =
	  hierarchy[0]->insert_exact_point_on_segment(ss, t, vnear[0]);
	vertex = vt.first;
	vs1 = vt.second;
	vs2 = vt.third;
      } else {
	vertex = hierarchy[0]->insert_point2(ss, t, vnear[0]);
      }
    }
  } else if ( n == 0 ) {
    vertex = hierarchy[0]->insert_first(ss, p);
  } else if ( n == 1 ) {
    vertex = hierarchy[0]->insert_second(ss, p);
  } else if ( n == 2 ) {
    vertex = hierarchy[0]->insert_third(ss, p);
  }

  CGAL_assertion( vertex != Vertex_handle() );

  if ( vertices != NULL ) { vertices[0] = vertex; }

  // insert at other levels
  Vertex_handle previous = vertex;
  Vertex_handle vs1_prev = vs1;
  Vertex_handle vs2_prev = vs2;

  //  Storage_site_2 ss = vertex->storage_site();

  int k = 1;
  while ( k <= new_level ) {
    size_type nv = hierarchy[k]->number_of_vertices();
    if ( nv > 2 ) {
      if ( nv < sdg_hierarchy_2__minsize ) {
	CGAL_assertion( vnear[k] == Vertex_handle() );
	vnear[k] = hierarchy[k]->nearest_neighbor(p);
      }
      Arrangement_type at_res = this->arrangement_type(t, vnear[k]->site());
      if ( at_res == AT2::INTERIOR ) {
	Vertex_triple vt =
	  hierarchy[k]->insert_exact_point_on_segment(ss, t, vnear[k]);
	vertex = vt.first;
	vs1 = vt.second;
	vs2 = vt.third;

	CGAL_assertion( (same_segments(vs1_prev->site(), vs1->site()) &&
			 same_segments(vs2_prev->site(), vs2->site())) ||
			(same_segments(vs1_prev->site(), vs2->site()) &&
			 same_segments(vs2_prev->site(), vs1->site())) );
	if ( same_segments(vs1_prev->site(), vs1->site()) ) {
	  vs1->set_down(vs1_prev);
	  vs1_prev->set_up(vs1);
	  vs1_prev = vs1;

	  vs2->set_down(vs2_prev);
	  vs2_prev->set_up(vs2);
	  vs2_prev = vs2;
	} else {
	  vs1->set_down(vs2_prev);
	  vs2_prev->set_up(vs1);
	  vs2_prev = vs1;

	  vs2->set_down(vs1_prev);
	  vs1_prev->set_up(vs2);
	  vs1_prev = vs2;
	}
      } else {
	vertex = hierarchy[k]->insert_point(ss, t, vnear[k]);	
      }
    } else if ( nv == 2 ) {
      vertex = hierarchy[k]->insert_third(t, ss);
    } else { // nv == 0 || nv == 1
      vertex = hierarchy[k]->insert_no_register(ss, p, vnear[k]);
    }

    CGAL_assertion( vertex != Vertex_handle() );

    if ( vertices != NULL ) { vertices[k] = vertex; }

    vertex->set_down(previous); // link with other levels
    previous->set_up(vertex);
    previous = vertex;
    k++;
  }
  return std::make_pair(lies_on_seg, new_level);
}


template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
insert_point(const Site_2& t, const Storage_site_2& ss,
	     int low, int high, Vertex_handle vbelow,
	     Vertex_handle* vertices)
{
  CGAL_precondition( low >= 1 && low <= high );
  CGAL_precondition( vbelow != Vertex_handle() );

  Vertex_handle vertex;
  Vertex_handle vnear[sdg_hierarchy_2__maxlevel];

  nearest_neighbor(t, vnear, false);

  Vertex_handle previous = vbelow;

  // insert at all levels
  int k = low;
  while ( k <= high ) {
    // MK::ERROR: this is a hack. I need to change the code in the
    // segment Delaunay graph class, so that I can insert sites as
    // the first, second, or third site...; actually this is
    // problematic only if the number of vertices is exactly 2; it
    // cannot be smaller since we already have added the endpoints
    // of the segment at level k.
    size_type n = hierarchy[k]->number_of_vertices();
    if ( n > 2 ) {
      vertex = hierarchy[k]->insert_point(ss, t, vnear[k]);
    } else if ( n == 2 ) {
      vertex = hierarchy[k]->insert_third(t, ss);
    } else {
      if ( ss.is_input() ) {
	vertex = hierarchy[k]->insert_no_register(ss, t.point(), vnear[k]);
      } else {
	break;
      }
      // ideally what should be instead of the break-statement above is the
      // following statement(s) in the #if-#endif block
#if 0
    } else if ( n == 1 ) {
      vertex = hierarchy[k]->insert_second(t, ss);
    } else {
      CGAL_assertion( n == 0 );
      vertex = hierarchy[k]->insert_first(t, ss);
#endif
    }

    CGAL_assertion( vertex != Vertex_handle() );

    if ( vertices != NULL ) { vertices[k] = vertex; }

    vertex->set_down(previous); // link with other levels
    previous->set_up(vertex);
    previous = vertex;
    k++;
  }
}


//--------------------------------------------------------------------
// insertion of a segment
//--------------------------------------------------------------------
template<class Gt, class ST, class STag, class D_S, class LTag>
typename
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
insert_segment(const Point_2& p0, const Point_2& p1,
	       const Storage_site_2& ss, int level)
{

  CGAL_SDG_DEBUG(std::cout << "debug hier insert_segment " 
                 << p0 << ' ' << p1 << std::endl;);

  // the tag is true so we DO insert segments in hierarchy
  if ( level == UNDEFINED_LEVEL ) {
    level = random_level();
  }

  Site_2 t = Site_2::construct_site_2(p0, p1);

  if ( is_degenerate_segment(t) ) {
    Storage_site_2 ss_src = ss.source_site();
    convert_info(ss_src, ss, true);
    return insert_point(p0, ss_src, level);
  }

  Vertex_handle vertices0[sdg_hierarchy_2__maxlevel];
  Vertex_handle vertices1[sdg_hierarchy_2__maxlevel];

  Storage_site_2 ss_src = ss.source_site();
  convert_info(ss_src, ss, true);
#if 1
  insert_point(p0, ss_src, level, vertices0);
#else
  std::pair<bool,int> res = insert_point(p0, ss_src, level, vertices0);
  if ( res.first ) { level = res.second; }
#endif

  Storage_site_2 ss_trg = ss.target_site();
  convert_info(ss_trg, ss, false);
  insert_point(p1, ss_trg, level, vertices1);

  CGAL_assertion( vertices0[0] != Vertex_handle() );
  CGAL_assertion( vertices1[0] != Vertex_handle() );

  Vertex_handle vertex;

  if ( hierarchy[0]->number_of_vertices() == 2 ) {
    Segments_in_hierarchy_tag stag;

    vertex = hierarchy[0]->insert_third(ss, vertices0[0], vertices1[0]);
    insert_segment_in_upper_levels(t, vertex->storage_site(),
				   vertex, vertices0, level, stag);
  } else {
    vertex = insert_segment_interior(t, ss, vertices0, level);
  }

  return vertex;
}

template<class Gt, class ST, class STag, class D_S, class LTag>
typename
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
insert_segment_interior(const Site_2& t, const Storage_site_2& ss,
			const Vertex_handle* vertices, int level)
{
  // insert the interior of a segment, and DO insert segments in
  // upper levels of the hierarchy
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( this->number_of_vertices() >= 2 );

  CGAL_assertion( vertices[0] != Vertex_handle() );
  // MK: add here code that checks if the inserted segment has already
  // been inserted; MAYBE THIS IS NOT NEEDED; I ALREADY DO IT IN
  // arrangement_type

  // the tags
  Intersections_tag          itag;
  Segments_in_hierarchy_tag  stag;

  CGAL_SDG_DEBUG(std::cout << "debug insert_segment_interior t=" << t 
                 << " at level " << level << std::endl;);

  // find the first conflict

  // first look if there are intersections...
  Vertex_circulator vc = hierarchy[0]->incident_vertices(vertices[0]);
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    if ( is_infinite(vv) ) {
      vc++;
      continue;
    }

    Arrangement_type at_res = this->arrangement_type(t, vv);
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
        CGAL_SDG_DEBUG(std::cout << "debug insert_segment_interior crossing "
                       << t << " and " << *vv << std::endl;);
	return insert_intersecting_segment_with_tag(ss, t, vv, level,
						    itag, stag);
      } else if ( at_res == AT2::TOUCH_11_INTERIOR_1 ) {
	Vertex_handle vp = second_endpoint_of_segment(vv);

	Storage_site_2 sss = split_storage_site(ss, vp->storage_site(), true);
	// merge the info of the first (common) subsegment
	merge_info(vv, sss);

	return insert_segment_on_point(ss, vp, level, 1);
      } else if ( at_res == AT2::TOUCH_12_INTERIOR_1 ) {
	Vertex_handle vp = first_endpoint_of_segment(vv);
	return insert_segment_on_point(ss, vp, level, 0);
      } else {
	// this should never be reached; the only possible values for
	// at_res are DISJOINT, CROSSING, TOUCH_11_INTERIOR_1
	// and TOUCH_12_INTERIOR_1
	CGAL_error();
      }
    } else {
      CGAL_assertion( vv->is_point() );
      if ( at_res == AT2::INTERIOR ) {
	Storage_site_2 svv = vv->storage_site();
	if ( svv.is_input() ) {
	  return insert_segment_on_point(ss, vv, level, 2);
	} else {
	  // MK::ERROR:: not ready yet
	  CGAL_error();
	}
      }
    }
    ++vc;
  } while ( vc != vc_start );



  CGAL_SDG_DEBUG(std::cout << "debug insert_segment_interior t=" << t 
                 << " start looking for vert cf " << std::endl;);

  // first look for conflict with vertex
  Face_circulator fc_start = hierarchy[0]->incident_faces(vertices[0]);
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

#ifdef CGAL_SDG_VERBOSE
  // debug
  std::cout << "debug hier insert_segment_interior t=" 
    << t << " with " 
    << (is_infinite(start_f)? "infinite" : "finite") 
    << " face [" ;
  if (is_infinite(start_f->vertex(0))) {
    std::cout << " infv";
  } else {
    std::cout << ' ' << start_f->vertex(0)->site();
  }
  if (is_infinite(start_f->vertex(1))) {
    std::cout << " infv";
  } else {
    std::cout << ' ' << start_f->vertex(1)->site();
  }
  if (is_infinite(start_f->vertex(2))) {
    std::cout << " infv";
  } else {
    std::cout << ' ' << start_f->vertex(2)->site();
  }
  std::cout << "] has sign s=" << s << std::endl; 

#endif

  // we are in conflict with a Voronoi vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;
  Face_map fm;

  Triple<bool, Vertex_handle, Arrangement_type>
    vcross(false, Vertex_handle(), AT2::DISJOINT);

  hierarchy[0]->initialize_conflict_region(start_f, l);
  hierarchy[0]->expand_conflict_region(start_f, t, ss, l, fm,
				       sign_map, vcross);

  CGAL_assertion( vcross.third == AT2::DISJOINT ||
		  vcross.third == AT2::CROSSING ||
		  vcross.third == AT2::INTERIOR );

  // the following condition becomes true only if intersecting
  // segments are found
  if ( vcross.first ) {
    if ( t.is_segment() ) {
      if ( vcross.third == AT2::CROSSING ) {
	Intersections_tag itag;
	return insert_intersecting_segment_with_tag(ss, t, vcross.second,
						    level, itag, stag);
      } else if ( vcross.third == AT2::INTERIOR ) {
	return insert_segment_on_point(ss, vcross.second, level, 2);
      } else {
	// this should never be reached; the only possible values for
	// vcross.third are CROSSING, INTERIOR and DISJOINT
	CGAL_error();
      }
    }
  }

  // no intersecting segment has been found; we insert the segment as
  // usual...
  Vertex_handle v = hierarchy[0]->create_vertex(ss);

  hierarchy[0]->retriangulate_conflict_region(v, l, fm);

  insert_segment_in_upper_levels(t, ss, v, vertices, level, stag);

  return v;
}


template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
insert_segment_in_upper_levels(const Site_2& t, const Storage_site_2& ss,
			       Vertex_handle vbelow,
			       const Vertex_handle* vertices,
			       int level, Tag_true /* stag */)
{
  CGAL_precondition( vertices != NULL );
  CGAL_precondition( vbelow != Vertex_handle() );

  // insert at all upper levels
  Vertex_handle previous = vbelow;
  Vertex_handle vertex = vbelow;

  int k = 1;
  while ( k <= level ) {
    if ( hierarchy[k]->number_of_vertices() == 2 ) {
      Vertex_handle v0(hierarchy[k]->finite_vertices_begin());
      Vertex_handle v1(++(hierarchy[k]->finite_vertices_begin()));
      CGAL_precondition( v0 != Vertex_handle() &&
			 v1 != Vertex_handle() );
      vertex = hierarchy[k]->insert_third(ss, v0, v1);
    } else {
      vertex = hierarchy[k]->insert_segment_interior(t, ss, vertices[k]);
    }

    CGAL_assertion( vertex != Vertex_handle() );

    vertex->set_down(previous); // link with level above
    previous->set_up(vertex);
    previous = vertex;
    k++;
  }
}

//--------------------------------------------------------------------
// insertion of a segment that goes through a point
//--------------------------------------------------------------------

template<class Gt, class ST, class STag, class D_S, class LTag>
typename
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
insert_segment_on_point(const Storage_site_2& ss,
			const Vertex_handle& v,
			int level, int which)
{  
  // inserts the segment represented by ss in the case where this
  // segment goes through a point which has already been inserted and
  // corresponds to the vertex handle v
  CGAL_precondition( ss.is_segment() );
  CGAL_precondition( v->is_point() );

  Storage_site_2 ssv = v->storage_site();
  Site_2 sv = ssv.site();

  Vertex_handle vnear[sdg_hierarchy_2__maxlevel];

  nearest_neighbor(sv, vnear, false);

  int v_level = find_level(v);

  if ( v_level < level ) {
    Vertex_handle vbelow = vnear[v_level];

    insert_point(sv, ssv, v_level + 1, level, vbelow, vnear);
  }

  // merge the info of the (common) splitting endpoint
  merge_info(v, ss);

  if ( which == 2 ) {
    Storage_site_2 ss1 = this->split_storage_site(ss, ssv, true);
    Storage_site_2 ss2 = this->split_storage_site(ss, ssv, false);

    insert_segment_interior(ss1.site(), ss1, vnear, level);
    return insert_segment_interior(ss2.site(), ss2, vnear, level);
  }

  Storage_site_2 ss1 = this->split_storage_site(ss, ssv, which == 0);
  return insert_segment_interior(ss1.site(), ss1, vnear, level);
}

//--------------------------------------------------------------------
// insertion of an intersecting segment
//--------------------------------------------------------------------
template<class Gt, class ST, class STag, class D_S, class LTag>
typename
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     int level,
				     Tag_true itag, Tag_false /* stag */)
{
  CGAL_precondition( t.is_segment() && v->is_segment() );

  const Storage_site_2& ssitev = v->storage_site();
  Site_2 sitev = ssitev.site();

  if ( same_segments(t, sitev) ) {
    merge_info(v, ss);
    return v;
  }

  Vertex_triple vt = hierarchy[0]->insert_point_on_segment(ss, t, v, itag);

  Vertex_handle verticesx[sdg_hierarchy_2__maxlevel];

  Vertex_handle vsx = vt.first;
  verticesx[0] = vsx;

  bool compute_new_level = true;
  int new_level = compute_new_level ? random_level() : level;
  if ( new_level > 0 ) {
    Storage_site_2 ssx = vsx->storage_site();
    Site_2 sx = ssx.site();
    insert_point(sx, ssx, 1, new_level, vsx, verticesx);
  }

  Storage_site_2 ss3 =
    this->st_.construct_storage_site_2_object()(ss, ssitev, true);
  Storage_site_2 ss4 =
    this->st_.construct_storage_site_2_object()(ss, ssitev, false);

  Site_2 s3 = ss3.site();
  Site_2 s4 = ss4.site();

  insert_segment_interior(s3, ss3, verticesx, level);
  insert_segment_interior(s4, ss4, verticesx, level);
  return vsx;
}

template<class Gt, class ST, class STag, class D_S, class LTag>
typename
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::Vertex_handle
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     int level,
				     Tag_true itag, Tag_true /* stag */)
{
  CGAL_precondition( t.is_segment() && v->is_segment() );

  CGAL_expensive_precondition( arrangement_type(t, v->site()) );

  const Storage_site_2& ssitev = v->storage_site();
  Site_2 sitev = ssitev.site();

  if ( same_segments(t, sitev) ) {
    // MK::ERROR: I may need to insert it to levels higher than its
    // previous level...
    merge_info(v, ss);
    return v;
  }

  Vertex_handle verticesx[sdg_hierarchy_2__maxlevel];

  int levelv = find_level(v);

  Vertex_handle vcross = v;
  Vertex_handle v1_old, v2_old, vsx_old;

  int k = 0;
  while ( k <= levelv ) {
    Vertex_handle vcross_up = vcross->up();

    Vertex_triple vt =
      hierarchy[k]->insert_point_on_segment(ss, t, vcross, itag);

    // now I need to update the sites for vertices v1 and v2
    Vertex_handle vsx = vt.first;
    Vertex_handle v1 = vt.second;
    Vertex_handle v2 = vt.third;

    CGAL_assertion( v1->is_segment() && v2->is_segment() );

    if ( k > 0 ) {
      CGAL_assertion( (same_segments(v1->site(), v1_old->site()) &&
		       same_segments(v2->site(), v2_old->site())) ||
		      (same_segments(v1->site(), v2_old->site()) &&
		       same_segments(v2->site(), v1_old->site()))
		      );
      if ( same_segments(v1->site(), v1_old->site()) ) {
	v1->set_down(v1_old);
	v2->set_down(v2_old);
	v1_old->set_up(v1);
	v2_old->set_up(v2);
      } else {
	v1->set_down(v2_old);
	v2->set_down(v1_old);
	v1_old->set_up(v2);
	v2_old->set_up(v1);
      }
      vsx_old->set_up(vsx);
      vsx->set_down(vsx_old);
    }

    v1_old = v1;
    v2_old = v2;
    vsx_old = vsx;

    verticesx[k] = vsx;

    vcross = vcross_up;
    k++;
  }

  if ( levelv < level ) {
    Storage_site_2 ssx = verticesx[0]->storage_site();
    Site_2 sx = ssx.site();

    insert_point(sx, ssx, levelv + 1, level, verticesx[levelv], verticesx);
  }

  Storage_site_2 ss3 =
    this->st_.construct_storage_site_2_object()(ss, ssitev, true);

  Storage_site_2 ss4 =
    this->st_.construct_storage_site_2_object()(ss, ssitev, false);

  Site_2 s3 = ss3.site();
  Site_2 s4 = ss4.site();

  insert_segment_interior(s3, ss3, verticesx, level);
  insert_segment_interior(s4, ss4, verticesx, level);
  return verticesx[0];
}


//====================================================================
//====================================================================
//                   METHODS FOR REMOVAL
//====================================================================
//====================================================================

template<class Gt, class ST, class STag, class D_S, class LTag>
bool
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
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

  // remove at level 0
  Vertex_handle u = v->up();
  bool success = hierarchy[0]->remove_base(v);
  if ( !success ) { return false; }

  // remove at higher levels
  unsigned int l = 1;
  Vertex_handle vv = u;
  while ( u != Vertex_handle() ) {
    if ( l >= sdg_hierarchy_2__maxlevel) { break; }
    vv = u;
    u = vv->up();
    success = hierarchy[l++]->remove_base(vv);
    CGAL_assertion( success );
  }

  // unregister the input site
  if ( is_point ) {
    this->unregister_input_site( h1 );
  } else {
    this->unregister_input_site( h1, h2 );
  }

  return true;
}

//===========================================================================
//===========================================================================

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  nearest neighbor location
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class ST, class STag, class D_S, class LTag>
typename
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::Vertex_handle 
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
nearest_neighbor(const Point_2& p, bool force_point) const
{
  Vertex_handle vnear[sdg_hierarchy_2__maxlevel];
  //  nearest_neighbor(Site_2(p), vnear, force_point);

  Site_2 t = Site_2::construct_site_2(p);

  nearest_neighbor(t, vnear, force_point);
  return vnear[0];
}

template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
nearest_neighbor(const Site_2& t,
		 Vertex_handle vnear[sdg_hierarchy_2__maxlevel],
		 bool /* force_point */) const
{
  CGAL_precondition( t.is_point() );

  Vertex_handle nearest;
  int level  = sdg_hierarchy_2__maxlevel;

  // find the highest level with enough vertices
  while ( hierarchy[--level]->number_of_vertices() 
	  < sdg_hierarchy_2__minsize ) {
    if ( !level ) break;  // do not go below 0
  }
  for (unsigned int i = level + 1; i < sdg_hierarchy_2__maxlevel; i++) {
    vnear[i] = Vertex_handle();
  }

  while ( level > 0 ) {
    vnear[level] = nearest =
      hierarchy[level]->nearest_neighbor(t, nearest);  

    CGAL_assertion( !hierarchy[level]->is_infinite(vnear[level]) );
    CGAL_assertion( vnear[level] != Vertex_handle() );
    // go at the same vertex on level below
    nearest = nearest->down();
    --level;
  }
  vnear[0] = hierarchy[0]->nearest_neighbor(t, nearest);
  // at level 0
}



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  miscellaneous methods
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::   
copy(const Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag> &sdg)
{
#ifndef CGAL_NO_ASSERTIONS
  for (unsigned int i = 1; i < sdg_hierarchy_2__maxlevel; ++i) {
    CGAL_assertion( hierarchy[i]->pc_.size() == 0 );
    CGAL_assertion( hierarchy[i]->isc_.size() == 0 );

    CGAL_assertion( sdg.hierarchy[i]->pc_.size() == 0 );
    CGAL_assertion( sdg.hierarchy[i]->isc_.size() == 0 );
  }
#endif

  // first copy the point container and input point container
  this->pc_ = sdg.pc_;

  // first create a map between the old point handles and the new ones
  Handle_map hm;

  Point_handle it_other = sdg.hierarchy[0]->pc_.begin();
  Point_handle it_this = hierarchy[0]->pc_.begin();
  for (; it_other != sdg.hierarchy[0]->pc_.end(); ++it_other, ++it_this) {
    hm.insert( Point_handle_pair(it_other, it_this) );
  }

  {
    for(unsigned int i = 0; i < sdg_hierarchy_2__maxlevel; ++i) {
      hierarchy[i]->copy(*sdg.hierarchy[i], hm);
    }
  }

  //up and down have been copied in straightforward way
  // compute a map at lower level
  std::map< Vertex_handle, Vertex_handle > V;
  {
    for(Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(); 
	it != hierarchy[0]->finite_vertices_end(); ++it) {
      if ( it->up() != Vertex_handle() ) {
	V[ it->up()->down() ] = it;
      }
    }
  }
  {
    for(unsigned int i = 1; i < sdg_hierarchy_2__maxlevel; ++i) {
      for(Finite_vertices_iterator it = hierarchy[i]->finite_vertices_begin(); 
	  it != hierarchy[i]->finite_vertices_end(); ++it) {
	// down pointer goes in original instead in copied triangulation
	it->set_down(V[it->down()]);
	// make reverse link
	it->down()->set_up( it );
	// make map for next level
	if ( it->up() != Vertex_handle() ) {
	  V[ it->up()->down() ] = it;
	}
      }
    }
  }

  // the point container and the input sites container are copied by
  // the operator= of the one-level classes
  //  hierarchy[0]->pc_ = sdg.hierarchy[0]->pc_;
  //  hierarchy[0]->isc_ = sdg.hierarchy[0]->isc_;
}

template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>:: 
clear()
{
  for(unsigned int i = 0; i < sdg_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->clear();
  }
}

template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>:: 
swap(Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>& other)
{
  Base* temp;
  Base::swap(other);
  for(unsigned int i = 1; i < sdg_hierarchy_2__maxlevel; ++i) {
    temp = hierarchy[i];
    hierarchy[i] = other.hierarchy[i];
    other.hierarchy[i]= temp;
  }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  validity check
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class ST, class STag, class D_S, class LTag>
bool
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>:: 
is_valid(bool verbose, int level) const
{
  // philaris: the level argument has to do with debugging level
  // philaris: the level argument is not related to the hierarchy levels

  bool result(true);

#define DEBUGVALIDHIER true

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>

//#undef DEBUGVALIDHIER

#ifdef DEBUGVALIDHIER  
  std::cout << "debug hierarchy is_valid entering" 
    << " level=" << level << std::endl;
#endif

  //verify correctness of triangulation at all levels
  for(unsigned int i = 0; i < sdg_hierarchy_2__maxlevel; ++i) {
#ifdef DEBUGVALIDHIER  
    std::cout << "debug hierarchy is_valid level " << i << std::endl;
#endif
    if ( verbose ) {
      std::cerr << "Level " << i << ": " << std::flush;
    }
    bool is_valid_level = hierarchy[i]->is_valid(verbose, level);

#ifdef DEBUGVALIDHIER  
    std::cout << "debug level=" << i 
      << " validity=" << is_valid_level << std::endl;
#endif

    result = result && is_valid_level;
    if ( verbose ) {
      std::cerr << std::endl;
    }
  }
  //verify that lower level has no down pointers
  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->finite_vertices_end(); ++it) {
    result = result && ( it->down() == Vertex_handle() );
  }

  //verify that other levels has down pointer and reciprocal link is fine
  for(unsigned int i = 1; i < sdg_hierarchy_2__maxlevel; ++i) {
    for( Finite_vertices_iterator it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->finite_vertices_end(); ++it) {
      Vertex_handle vit(it);
      result = result && ( it->down()->up() == vit );
    }
  }
#ifdef DEBUGVALIDHIER  
  std::cout << "debug hierarchy is_valid returns " << result << std::endl;
#endif
  return result;
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  local helper methods
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
print_error_message() const
{
  std::cerr << std::endl;
  std::cerr << "WARNING:" << std::endl;
  std::cerr << "A segment-segment intersection was found."
	    << std::endl;
  std::cerr << "The Segment_Delaunay_graph_Linf_hierarchy_2 class is not"
	    << " configured to handle this situation." << std::endl;
  std::cerr << "Please look at the documentation on how to handle"
	    << " this behavior." << std::endl;
  std::cerr << std::endl;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  file I/O
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
file_output(std::ostream& os) const
{
  typedef std::map<Vertex_handle,int>        Vertex_map;
  typedef typename Base::const_Point_handle  const_Point_handle;

  // compute the point handle mapper
  typename Base::Point_handle_mapper P;
  size_type inum = 0;
  for (const_Point_handle ph = this->pc_.begin(); ph != this->pc_.end(); ++ph) {
    P[ph] = inum++;
  }

  // write each level of the hierarchy
  hierarchy[0]->file_output(os, P, true);
  if ( is_ascii(os) ) { os << std::endl << std::endl; }
  for (unsigned int i = 1; i < sdg_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->file_output(os, P, false);
    if ( is_ascii(os) ) { os << std::endl << std::endl; }
  }

  Vertex_map* V = new Vertex_map[sdg_hierarchy_2__maxlevel];

  // create the map of vertex indices
  for (unsigned int i = 0; i < sdg_hierarchy_2__maxlevel; ++i) {
    int inum = 0;
    for (Finite_vertices_iterator vit = hierarchy[i]->finite_vertices_begin();
	 vit != hierarchy[i]->finite_vertices_end(); ++vit) {
      V[i][vit] = inum++;
    }
  }

  Vertex_map* V_up   = new Vertex_map[sdg_hierarchy_2__maxlevel];
  Vertex_map* V_down = new Vertex_map[sdg_hierarchy_2__maxlevel];

  // create the map of up and down pointers
  for (unsigned int i = 0; i < sdg_hierarchy_2__maxlevel; ++i) {
    for (Finite_vertices_iterator vit = hierarchy[i]->finite_vertices_begin();
	 vit != hierarchy[i]->finite_vertices_end(); ++vit) {
      if ( vit->up() != Vertex_handle() ) {
	V_up[i][vit] = V[i+1][vit->up()];
      } else {
	V_up[i][vit] = -1;
      }

      if ( vit->down() != Vertex_handle() ) {
	V_down[i][vit] = V[i-1][vit->down()];
      } else {
	V_down[i][vit] = -1;
      }
    }
  }

  // write up and down pointer info
  if ( is_ascii(os) ) { os << std::endl << std::endl; }
  for (unsigned int i = 0; i < sdg_hierarchy_2__maxlevel; ++i) {
    os << i;
    if ( is_ascii(os) ) { os << " "; }
    os << hierarchy[i]->number_of_vertices();
    if ( is_ascii(os) ) { os << std::endl; }
    for (Finite_vertices_iterator vit = hierarchy[i]->finite_vertices_begin();
	 vit != hierarchy[i]->finite_vertices_end(); ++vit) {
      os << V[i][vit];
      if ( is_ascii(os) ) { os << " "; }
      os << V_down[i][vit];
      if ( is_ascii(os) ) { os << " "; }
      os << V_up[i][vit];
      if ( is_ascii(os) ) { os << std::endl; }
    }
    if ( is_ascii(os) ) { os << std::endl << std::endl; }
  }

  delete[] V;
  delete[] V_up;
  delete[] V_down;
}

template<class Gt, class ST, class STag, class D_S, class LTag>
void
Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST,STag,D_S,LTag>::
file_input(std::istream& is)
{
  typedef std::vector<Vertex_handle>  Vertex_vector;

  // firstly, read the segment Delaunay graph at each level
  clear();
  typename Base::Point_handle_vector P;
  hierarchy[0]->file_input(is, true, P);
  for (unsigned int i = 1; i < sdg_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->file_input(is, false, P);
  }

  Vertex_vector* V = new Vertex_vector[sdg_hierarchy_2__maxlevel];

  // secondly, create the map of vertex indices
  for (unsigned int i = 0; i < sdg_hierarchy_2__maxlevel; ++i) {
    V[i].resize(hierarchy[i]->number_of_vertices());
    int j = 0;
    for (Finite_vertices_iterator vit = hierarchy[i]->finite_vertices_begin();
	 vit != hierarchy[i]->finite_vertices_end(); ++vit, ++j) {
      V[i][j] = vit;
    }
  }

  // read the correspondences between up and down pointers and set
  // them appropriately
  for (unsigned int i = 0; i < sdg_hierarchy_2__maxlevel; ++i) {
    unsigned int level;
    int vnum;
    is >> level >> vnum;
    for (int k = 0; k < vnum; k++) {
      int ithis, iup, idown;
      is >> ithis >> idown >> iup;
      if ( idown != -1 ) { V[i][ithis]->set_down(V[i-1][idown]); }
      if ( iup != -1 )   { V[i][ithis]->set_up(V[i+1][iup]); }
    }
  }

  delete[] V;
}


//--------------------------------------------------------------------


} //namespace CGAL


// EOF
