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


//===========================================================================

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  constructors
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template<class Gt, class STag, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
init_hierarchy(const Geom_traits& gt)
{
  hierarchy[0] = this; 
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i) {
    hierarchy[i] = new Base(gt);
  }
}

template<class Gt, class STag, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
Segment_Voronoi_diagram_hierarchy_2(const Gt& gt)
  : Base(gt), random((long)0)
{ 
  init_hierarchy(gt);
}


// copy constructor duplicates vertices and faces
template<class Gt, class STag, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
Segment_Voronoi_diagram_hierarchy_2
(const Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> &svd)
    : Base(svd.geom_traits()), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  init_hierarchy(svd.geom_traits());
  copy(svd);
} 

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  destructor
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class STag, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>:: 
~Segment_Voronoi_diagram_hierarchy_2()
{
  clear();
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i){ 
    delete hierarchy[i];
  }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// assignment operator
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class STag, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> &
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
operator=(const Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> &svd)
{
  if ( this != &svd ) {
    copy(svd);
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

template<class Gt, class STag, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
insert_point(const Point_2& p, int level, Vertex_handle* vertices)
{
  CGAL_precondition( level != UNDEFINED_LEVEL );

  int new_level = level;

  Vertex_handle vertex;
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];

  Site_2 t = Site_2::construct_site_2(p);

  nearest_neighbor(t, vnear, false);

  size_type n = hierarchy[0]->number_of_vertices();
  if ( n > 2 ) {
    Arrangement_type at_res = this->arrangement_type(t, vnear[0]);

    CGAL_assertion( at_res == AT2::DISJOINT ||
		    at_res == AT2::INTERIOR ||
		    at_res == AT2::IDENTICAL );

    if ( vnear[0]->is_point() ) {
      if ( at_res == AT2::IDENTICAL ) {
	vertex = vnear[0];
      } else {
	Storage_site_2 ss = create_storage_site(p);
	vertex = hierarchy[0]->insert_point2(ss, t, vnear[0]);
      }
    } else { // nearest neighbor is a segment
      CGAL_assertion( vnear[0]->is_segment() );
      CGAL_assertion( at_res == AT2::DISJOINT ||
		      at_res == AT2::INTERIOR );

      Storage_site_2 ss = create_storage_site(p);

      if ( at_res == AT2::INTERIOR ) {
	CGAL_assertion( t.is_input() );

	int vnear_level = find_level(vnear[0]);

	// I need to find the level of the nearest neighbor that t
	// lies on and update the level of t
	if ( new_level < vnear_level ) {
	  new_level = vnear_level;
	}

	Vertex_triple vt =
	  hierarchy[0]->insert_exact_point_on_segment(ss, t, vnear[0]);
	vertex = vt.first;
      } else {
	vertex = hierarchy[0]->insert_point2(ss, t, vnear[0]);
      }
    }
  } else if ( n == 0 ) {
    vertex = hierarchy[0]->insert_first(p);
  } else if ( n == 1 ) {
    vertex = hierarchy[0]->insert_second(p);
  } else if ( n == 2 ) {
    vertex = hierarchy[0]->insert_third(p);
  }

  CGAL_assertion( vertex != Vertex_handle() );

  if ( vertices != NULL ) { vertices[0] = vertex; }

  // insert at other levels
  Vertex_handle previous = vertex;

  Storage_site_2 ss = vertex->storage_site();

  int k = 1;
  while ( k <= new_level ) {
    int nv = hierarchy[k]->number_of_vertices();
    if ( nv > 2 ) {
      vertex = hierarchy[k]->insert_point(ss, t, vnear[k]);
    } else if ( nv == 2 ) {
      vertex = hierarchy[k]->insert_third(t, ss);
    } else { // nv == 0 || nv == 1
      vertex = hierarchy[k]->insert_no_register(p, vnear[k]);
    }

    CGAL_assertion( vertex != Vertex_handle() );

    if ( vertices != NULL ) { vertices[k] = vertex; }

    vertex->set_down(previous); // link with other levels
    previous->set_up(vertex);
    previous = vertex;
    k++;
  }
}


template<class Gt, class STag, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
insert_point(const Site_2& t, const Storage_site_2& ss,
	     int low, int high, Vertex_handle vbelow,
	     Vertex_handle* vertices)
{
  CGAL_precondition( low >= 0 && low <= high );
  CGAL_precondition( (low == 0 && vbelow == Vertex_handle()) ||
		     (low > 0  && vbelow != Vertex_handle()) );

  Vertex_handle vertex;
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];

  nearest_neighbor(t, vnear, false);

  Vertex_handle previous = vbelow;

  // insert at all levels
  int k = low;
  while ( k <= high ) {
    // MK::ERROR: this is a hack. I need to change the code in the
    // segment Voronoi diagram class, so that I can insert sites as
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
	vertex = hierarchy[k]->insert_no_register(t.point(), vnear[k]);
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
template<class Gt, class STag, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
insert_segment(const Point_2& p0, const Point_2& p1,
	       int level/*, Tag_true stag*/)
{
  // the tag is true so we DO insert segments in hierarchy
  if ( level == UNDEFINED_LEVEL ) {
    level = random_level();
  }

  Site_2 t = Site_2::construct_site_2(p0, p1);

  if ( is_degenerate_segment(t) ) {
    return insert_point(p0, level);
  }

  Vertex_handle vertices0[svd_hierarchy_2__maxlevel];
  Vertex_handle vertices1[svd_hierarchy_2__maxlevel];

  insert_point(p0, level, vertices0);

#if 0
  insert_point(p1, level, vertices1);
#else // this way may be faster...
  vertices1[0] = hierarchy[0]->insert_no_register(p1, vertices0[0]);
  if ( level >= 1 ) {
    Storage_site_2 ss1 = vertices1[0]->storage_site();
    insert_point(ss1.site(), ss1, 1, level, vertices1[0], vertices1);
  }
#endif

  CGAL_assertion( vertices0[0] != Vertex_handle() );
  CGAL_assertion( vertices1[0] != Vertex_handle() );

  Vertex_handle vertex;

  if ( hierarchy[0]->number_of_vertices() == 2 ) {
    static Segments_in_hierarchy_tag stag;

    vertex = hierarchy[0]->insert_third(vertices0[0], vertices1[0]);
    insert_segment_in_upper_levels(t, vertex->storage_site(),
				   vertex, vertices0, level, stag);
  } else {
    Storage_site_2 ss = create_storage_site(vertices0[0], vertices1[0]);
    vertex = insert_segment_interior(t, ss, vertices0, level);
  }

  return vertex;
}

template<class Gt, class STag, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
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
  static Intersections_tag          itag;
  static Segments_in_hierarchy_tag  stag;

  // find the first conflict

  // first look if there are intersections...
  Vertex_circulator vc = vertices[0]->incident_vertices();
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
	return vv;
      } else if ( at_res == AT2::CROSSING ) {
	return insert_intersecting_segment_with_tag(ss, t, vv, level,
						    itag, stag);
      } else if ( at_res == AT2::TOUCH_11_INTERIOR_1 ) {
#if 1
	Vertex_handle vp = second_endpoint_of_segment(vv);
	return insert_segment_on_point(ss, vp, level, stag, 1);
#else
	Vertex_handle vp = second_endpoint_of_segment(vv);
	Storage_site_2 ssvp = vp->storage_site();
	Site_2 svp = ssvp.site();

	Vertex_handle vp_vnear[svd_hierarchy_2__maxlevel];

	nearest_neighbor(svp, vp_near, false);

	int vp_level = find_level(vp);

	if ( vp_level < level ) {
	  Vertex_handle vbelow = vp_vertices[vp_level];

	  insert_point(svp, ssvp, vp_level + 1,	level, vbelow, vp_vnear);
	}

	Storage_site_2 sss = split_storage_site(ss, svp, 1, itag);
	return insert_segment_interior(sss.site(), sss, vp_vnear, level);
#endif
      } else if ( at_res == AT2::TOUCH_12_INTERIOR_1 ) {
#if 1
	Vertex_handle vp = first_endpoint_of_segment(vv);
	return insert_segment_on_point(ss, vp, level, stag, 0);
#else
	Vertex_handle vp = first_endpoint_of_segment(vv);
	Storage_site_2 ssvp = vp->storage_site();
	Site_2 svp = ssvp.site();

	Vertex_handle vp_vnear[svd_hierarchy_2__maxlevel];

	nearest_neighbor(svp, vp_near, false);

	int vp_level = find_level(vp);

	if ( vp_level < level ) {
	  Vertex_handle vbelow = vp_vertices[vp_level];

	  insert_point(svp, ssvp, vp_level + 1, level, vbelow, vp_vnear);
	}

	Storage_site_2 sss = split_storage_site(ss, ssvp, 0, itag);
	return insert_segment_interior(sss.site(), sss, vp_vnear, level);
#endif
      } else {
	// this should never be reached; the only possible values for
	// at_res are DISJOINT, CROSSING, TOUCH_11_INTERIOR_1
	// and TOUCH_12_INTERIOR_1
	CGAL_assertion( false );
      }
    } else {
      CGAL_assertion( vv->is_point() );
      if ( at_res == AT2::INTERIOR ) {
	Storage_site_2 svv = vv->storage_site();
	if ( svv.is_input() ) {
#if 1
	  return insert_segment_on_point(ss, vv, level, stag, 2);
#else
	  //************************************************************
	  //************************************************************
	  //************************************************************
	  //************************************************************
	  // here I need a new method: insert_segment_on_point
	  // which will do the splitting and arrange the hierarchy
	  // pointers
	  //************************************************************
	  //************************************************************
	  //************************************************************
	  //************************************************************
	  Storage_site_2 ss1 = this->split_storage_site(ss, svv, 0, itag);
	  Storage_site_2 ss2 = this->split_storage_site(ss, svv, 1, itag);
	  CGAL_assertion( false );
	  return Vertex_handle();
	  //	  insert_segment_interior(ss1.site(), ss1, vv);
	  //	  return insert_segment_interior(ss2.site(), ss2, vv);
#endif
	} else {
	  // MK::ERROR:: not ready yet
	  CGAL_assertion( false );
	}
      }
    }
    ++vc;
  } while ( vc != vc_start );

  // first look for conflict with vertex
  Face_circulator fc_start = vertices[0]->incident_faces();
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
#if 1
	return insert_segment_on_point(ss, vcross.second, level, stag, 2);
#else
	Storage_site_2 ssvv = vcross.second->storage_site();
	Intersections_tag itag;
	//************************************************************
	//************************************************************
	//************************************************************
	//************************************************************
	// here I need a new method: insert_segment_on_point
	// which will do the splitting and arrange the hierarchy
	// pointers
	//************************************************************
	//************************************************************
	//************************************************************
	//************************************************************
	Storage_site_2 ss1 = this->split_storage_site(ss, ssvv, 0, itag);
	Storage_site_2 ss2 = this->split_storage_site(ss, ssvv, 1, itag);
	CGAL_assertion( false );
	return Vertex_handle();
	//	insert_segment_interior(ss1.site(), ss1, vcross.second);
	//	return insert_segment_interior(ss2.site(), ss2, vcross.second);
#endif
      } else {
	// this should never be reached; the only possible values for
	// vcross.third are CROSSING, INTERIOR and DISJOINT
	CGAL_assertion( false );
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


template<class Gt, class STag, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
insert_segment_in_upper_levels(const Site_2& t, const Storage_site_2& ss,
			       Vertex_handle vbelow,
			       const Vertex_handle* vertices,
			       int level, Tag_true stag)
{
  CGAL_precondition( vertices != NULL );

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
      vertex = hierarchy[k]->insert_third(v0, v1);
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
template<class Gt, class STag, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
insert_segment_on_point(const Storage_site_2& ss,
			const Vertex_handle& v,
			int level, Tag_true stag, int which)
{  
  // inserts the segment represented by ss in the case where this
  // segment goes through a point which has already been inserted and
  // corresponds to the vertex handle v
  CGAL_precondition( ss.is_segment() );
  CGAL_precondition( v->is_point() );

  Storage_site_2 ssv = v->storage_site();
  Site_2 sv = ssv.site();

  Vertex_handle vnear[svd_hierarchy_2__maxlevel];

  nearest_neighbor(sv, vnear, false);

  int v_level = find_level(v);

  if ( v_level < level ) {
    Vertex_handle vbelow = vnear[v_level];

    insert_point(sv, ssv, v_level + 1, level, vbelow, vnear);
  }

  Intersections_tag itag;
  if ( which == 2 ) {
    Storage_site_2 ss1 = this->split_storage_site(ss, ssv, 0, itag);
    Storage_site_2 ss2 = this->split_storage_site(ss, ssv, 1, itag);

    insert_segment_interior(ss1.site(), ss1, vnear, level);
    return insert_segment_interior(ss2.site(), ss2, vnear, level);
  } else {
    Storage_site_2 ss1 = this->split_storage_site(ss, ssv, which, itag);

    return insert_segment_interior(ss1.site(), ss1, vnear, level);
  }
}

//--------------------------------------------------------------------
// insertion of an intersecting segment
//--------------------------------------------------------------------
template<class Gt, class STag, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     int level,
				     Tag_true itag, Tag_false stag)
{
  CGAL_precondition( t.is_segment() && v->is_segment() );

  const Storage_site_2& ssitev = v->storage_site();
  Site_2 sitev = ssitev.site();

  if ( same_segments(t, sitev) ) {
    return v;
  }

  Vertex_triple vt = hierarchy[0]->insert_point_on_segment(ss, t, v, itag);

  Vertex_handle verticesx[svd_hierarchy_2__maxlevel];

  Vertex_handle vsx = vt.first;
  verticesx[0] = vsx;

  bool compute_new_level = true;
  int new_level = compute_new_level ? random_level() : level;
  if ( new_level > 0 ) {
    Storage_site_2 ssx = vsx->storage_site();
    Site_2 sx = ssx.site();
    insert_point(sx, ssx, 1, new_level, vsx, verticesx);
  }

  Storage_site_2 ss3, ss4;
  Site_2 s3, s4;
  if ( t.is_input(0) ) {
    ss3 = create_storage_site(ss, ssitev, true);
  } else {
    ss3 = create_storage_site_type1(ss, ss, ssitev);
  }
  s3 = ss3.site();

  if ( t.is_input(1) ) {
    ss4 = create_storage_site(ss, ssitev, false);
  } else {
    ss4 = create_storage_site_type2(ss, ssitev, ss);
  }
  s4 = ss4.site();

  insert_segment_interior(s3, ss3, verticesx, level);
  insert_segment_interior(s4, ss4, verticesx, level);
  return vsx;
}

template<class Gt, class STag, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     int level,
				     Tag_true itag, Tag_true stag)
{
  CGAL_precondition( t.is_segment() && v->is_segment() );

  CGAL_expensive_precondition( arrangement_type(t, v->site()) );

  const Storage_site_2& ssitev = v->storage_site();
  Site_2 sitev = ssitev.site();

  if ( same_segments(t, sitev) ) {
    // MK::ERROR: I may need to insert it to levels higher than its
    // previous level...
    return v;
  }

  Vertex_handle verticesx[svd_hierarchy_2__maxlevel];

  int levelv = find_level(v);

  Vertex_handle vcross = v;
  Vertex_handle v1_old, v2_old, vsx_old;

  int k = 0;
  while ( k <= levelv ) {
    // MK::ERROR: I have to remove this; too expensive...
    CGAL_expensive_precondition( arrangement_type(t, vertex->site()) );

    Vertex_handle vcross_up = vcross->up();

    Vertex_triple vt =
      hierarchy[k]->insert_point_on_segment(ss, t, vcross, itag);

    // now I need to update the sites for vertices v1 and v2
    Vertex_handle vsx = vt.first;
    Vertex_handle v1 = vt.second;
    Vertex_handle v2 = vt.third;

    CGAL_assertion( v1->is_segment() && v2->is_segment() );

    if ( k > 0 ) {
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

  Storage_site_2 ss3, ss4;
  Site_2 s3, s4;
  if ( t.is_input(0) ) {
    ss3 = create_storage_site(ss, ssitev, true);
  } else {
    ss3 = create_storage_site_type1(ss, ss, ssitev);
  }
  s3 = ss3.site();

  if ( t.is_input(1) ) {
    ss4 = create_storage_site(ss, ssitev, false);
  } else {
    ss4 = create_storage_site_type2(ss, ssitev, ss);
  }
  s4 = ss4.site();

  insert_segment_interior(s3, ss3, verticesx, level);
  insert_segment_interior(s4, ss4, verticesx, level);
  return verticesx[0];
}


//===========================================================================
//===========================================================================

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  nearest neighbor location
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class STag, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::Vertex_handle 
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
nearest_neighbor(const Point_2& p, bool force_point) const
{
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];
  //  nearest_neighbor(Site_2(p), vnear, force_point);

  Site_2 t = Site_2::construct_site_2(p);

  nearest_neighbor(t, vnear, force_point);
  return vnear[0];
}

template<class Gt, class STag, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
nearest_neighbor(const Site_2& t,
		 Vertex_handle vnear[svd_hierarchy_2__maxlevel],
		 bool force_point) const
{
  CGAL_precondition( t.is_point() );

  Vertex_handle nearest;
  int level  = svd_hierarchy_2__maxlevel;

  // find the highest level with enough vertices
  while ( hierarchy[--level]->number_of_vertices() 
	  < svd_hierarchy_2__minsize ) {
    if ( !level ) break;  // do not go below 0
  }
  for (unsigned int i = level + 1; i < svd_hierarchy_2__maxlevel; i++) {
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
template<class Gt, class STag, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::   
copy(const Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> &svd)
{
  std::map< Vertex_handle, Vertex_handle > V;
  {
    for(unsigned int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
      *(hierarchy[i]) = *svd.hierarchy[i];
    }
  }

  //up and down have been copied in straightforward way
  // compute a map at lower level
  {
    for(Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(); 
	it != hierarchy[0]->finite_vertices_end(); ++it) {
      if ( it->up() != Vertex_handle() ) {
	V[ it->up()->down() ] = it;
      }
    }
  }
  {
    for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i) {
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
  //  hierarchy[0]->pc_ = svd.hierarchy[0]->pc_;
  //  hierarchy[0]->isc_ = svd.hierarchy[0]->isc_;
}

template<class Gt, class STag, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>:: 
clear()
{
  for(unsigned int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->clear();
  }
}

template<class Gt, class STag, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>:: 
swap(Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>& other)
{
  Base* temp;
  Base::swap(other);
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i) {
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
template<class Gt, class STag, class DS, class LTag>
bool
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>:: 
is_valid(bool verbose, int level) const
{
  bool result(true);

  //verify correctness of triangulation at all levels
  for(unsigned int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
    if ( verbose ) {
      std::cerr << "Level " << i << ": " << std::flush;
    }
    result = result && hierarchy[i]->is_valid(verbose, level);
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
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i) {
    for( Finite_vertices_iterator it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->finite_vertices_end(); ++it) {
      Vertex_handle vit(it);
      result = result && ( it->down()->up() == vit );
    }
  }
  return result;
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  local helper methods
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template<class Gt, class STag, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>::
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


CGAL_END_NAMESPACE


// EOF
