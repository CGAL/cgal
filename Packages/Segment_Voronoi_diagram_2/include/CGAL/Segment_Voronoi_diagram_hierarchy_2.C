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
template<class Gt, class STag, class PC, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
Segment_Voronoi_diagram_hierarchy_2(const Geom_traits& traits)
  : Base(traits), random((long)0)
{ 
  hierarchy[0] = this; 
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Base(traits);
}


// copy constructor duplicates vertices and faces
template<class Gt, class STag, class PC, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
Segment_Voronoi_diagram_hierarchy_2
(const Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag> &svd)
    : Base(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i = 1; i < svd_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Base(svd.geom_traits());
  copy_triangulation(svd);
} 

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  destructor
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class STag, class PC, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>:: 
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
template<class Gt, class STag, class PC, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag> &
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
operator=(const Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag> &svd)
{
  copy_triangulation(svd);
  return *this;
}

//====================================================================
//====================================================================
//                   METHODS FOR INSERTION
//====================================================================
//====================================================================

template<class Gt, class STag, class PC, class DS, class LTag>
inline typename 
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert(const Site_2& t) {
  // the intended use is to unify the calls to insert(...);
  // thus the site must be an exact one; 
  CGAL_precondition( t.is_exact() );
  if ( t.is_segment() ) {
    return insert_segment(t.source(), t.target(), UNDEFINED_LEVEL);
  } else if ( t.is_point() ) {
    return insert_point(t.point(), UNDEFINED_LEVEL);
  } else {
    CGAL_precondition ( t.is_defined() );
    return Vertex_handle(); // to avoid compiler error
  }
}

//--------------------------------------------------------------------
// insertion of a point
//--------------------------------------------------------------------
template<class Gt, class STag, class PC, class DS, class LTag>
inline typename 
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_point(const Point_2& p, int level)

{
  if ( level == UNDEFINED_LEVEL ) {
    level = random_level();
  }

  Vertex_handle vertices[svd_hierarchy_2__maxlevel];
  
  insert_point(p, level, vertices);

  return vertices[0];
}

template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_point(const Point_2& p, int level, Vertex_handle* vertices)
{
  CGAL_precondition( level != UNDEFINED_LEVEL );

  Vertex_handle vertex;
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];

  nearest_neighbor(p, vnear, false);

  vertex = hierarchy[0]->insert(p, vnear[0]);

  if ( vertices != NULL ) { vertices[0] = vertex; }

  CGAL_assertion( vertex != Vertex_handle() );

  // insert at other levels
  Vertex_handle previous = vertex;
      
  int k = 1;
  while ( k <= level ) {
    vertex = hierarchy[k]->insert(p, vnear[k]);

    CGAL_assertion( vertex != Vertex_handle() );

    if ( vertices != NULL ) { vertices[k] = vertex; }

    vertex->set_down(previous); // link with other levels
    previous->set_up(vertex);
    previous = vertex;
    k++;
  }
}


template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
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
      if ( ss.is_exact() ) {
	vertex = hierarchy[k]->insert(t.point(), vnear[k]);
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
template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_segment(const Point_2& p0, const Point_2& p1,
	       int level/*, Tag_true stag*/)
{
  // the tag is true so we DO insert segments in hierarchy
  if ( level == UNDEFINED_LEVEL ) {
    level = random_level();
  }

  Site_2 t(p0, p1);

  if ( is_degenerate_segment(t) ) {
    return insert_point(p0, level);
  }

  Vertex_handle vertices0[svd_hierarchy_2__maxlevel];
  Vertex_handle vertices1[svd_hierarchy_2__maxlevel];

  insert_point(p0, level, vertices0);
  insert_point(p1, level, vertices1);

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

template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_segment_interior(const Site_2& t, const Storage_site_2& ss,
			const Vertex_handle* vertices, int level)
{
  // insert the interior of a segment, and DO insert segments in
  // upper levels of the hierarchy
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( number_of_vertices() >= 2 );

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
    if ( same_segments(t, vv) ) {
      return vv;
    }
    if ( arrangement_type(t, vv) ) {
      if ( t.is_segment() ) {
	return insert_intersecting_segment_with_tag(ss, t, vv, level,
						    itag, stag);
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

  std::pair<bool, Vertex_handle> vcross(false, Vertex_handle());

  hierarchy[0]->initialize_conflict_region(start_f, l);
  hierarchy[0]->expand_conflict_region(start_f, t, ss, l, fm,
				       sign_map, vcross);

  // the following condition becomes true only if intersecting
  // segments are found
  if ( vcross.first ) {
    if ( t.is_segment() ) {
      Intersections_tag itag;
      return insert_intersecting_segment_with_tag(ss, t, vcross.second,
						  level, itag, stag);
    }
  }

  // no intersecting segment has been found; we insert the segment as
  // usual...
  Vertex_handle v = hierarchy[0]->create_vertex(ss);

  hierarchy[0]->retriangulate_conflict_region(v, l, fm);

  insert_segment_in_upper_levels(t, ss, v, vertices, level, stag);
  return v;
}


template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
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
// insertion of an intersecting segment
//--------------------------------------------------------------------
template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
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
  if ( t.is_exact(0) ) {
    s3.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), true);
    ss3 = create_storage_site(ss, ssitev, true);
  } else {
    s3.set_segment(t.point(0), t.point(1),
		   t.point(2), t.point(3),
		   sitev.point(0), sitev.point(1));
    ss3 = create_storage_site_type1(ss, ss, ssitev);
  }

  if ( t.is_exact(1) ) {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), false);
    ss4 = create_storage_site(ss, ssitev, false);
  } else {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1),
		   t.point(4), t.point(5));
    ss4 = create_storage_site_type2(ss, ssitev, ss);
  }

  insert_segment_interior(s3, ss3, verticesx, level);
  insert_segment_interior(s4, ss4, verticesx, level);
  return vsx;
}

template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
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
  if ( t.is_exact(0) ) {
    s3.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), true);
    ss3 = create_storage_site(ss, ssitev, true);
  } else {
    s3.set_segment(t.point(0), t.point(1),
		   t.point(2), t.point(3),
		   sitev.point(0), sitev.point(1));
    ss3 = create_storage_site_type1(ss, ss, ssitev);
  }

  if ( t.is_exact(1) ) {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), false);
    ss4 = create_storage_site(ss, ssitev, false);
  } else {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1),
		   t.point(4), t.point(5));
    ss4 = create_storage_site_type2(ss, ssitev, ss);
  }

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
template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle 
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
nearest_neighbor(const Point_2& p, bool force_point) const
{
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];
  nearest_neighbor(Site_2(p), vnear, force_point);
  return vnear[0];
}

template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
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
template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::   
copy_triangulation
(const Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag> &svd)
{
  std::map< Vertex_handle, Vertex_handle > V;
  {
    for(int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
      *(hierarchy[i]) = *svd.hierarchy[i];
    }
  }

  //up and down have been copied in straightforward way
  // compute a map at lower level
  {
    for(All_vertices_iterator it = hierarchy[0]->all_vertices_begin(); 
	it != hierarchy[0]->all_vertices_end(); ++it) {
      if ( it->up() != Vertex_handle() ) {
	V[ it->up()->down() ] = it;
      }
    }
  }
  {
    for(int i = 1; i < svd_hierarchy_2__maxlevel; ++i) {
      for(All_vertices_iterator it = hierarchy[i]->all_vertices_begin(); 
	  it != hierarchy[i]->all_vertices_end(); ++it) {
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

  // copy the point container
  hierarchy[0]->pc_ = svd.hierarchy[0]->pc_;
}

template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>:: 
clear()
{
  for(unsigned int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->clear();
  }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//  validity check
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template<class Gt, class STag, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>:: 
is_valid(bool verbose, int level) const
{
  bool result(true);

  //verify correctness of triangulation at all levels
  for(unsigned int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
    if ( verbose ) {
      std::cout << "Level " << i << ": " << std::flush;
    }
    result = result && hierarchy[i]->is_valid(verbose, level);
    if ( verbose ) {
      std::cout << std::endl;
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
template<class Gt, class STag, class PC, class DS, class LTag>
int
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
random_level()
{
  unsigned int l = 0;
  while ( true ) {
    if ( random(svd_hierarchy_2__ratio) ) break;
    ++l;
  }
  if (l >= svd_hierarchy_2__maxlevel)
    l = svd_hierarchy_2__maxlevel -1;
  return l;
}

template<class Gt, class STag, class PC, class DS, class LTag>
inline typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::size_type
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
find_level(Vertex_handle v) const
{
  CGAL_precondition( v != Vertex_handle() );
  size_type level = 0;
  Vertex_handle vertex = v;
  while ( vertex->up() != Vertex_handle() ) {
    vertex = vertex->up();
    level++;
  }

  return level;
}

template<class Gt, class STag, class PC, class DS, class LTag>
inline
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
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
