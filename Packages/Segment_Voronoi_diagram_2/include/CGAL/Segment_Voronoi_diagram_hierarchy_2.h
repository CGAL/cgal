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



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_HIERARCHY_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_HIERARCHY_2_H

#include <CGAL/Random.h>
#include <map>
// the following include should be removed
#include <CGAL/Triangulation_hierarchy_2.h>

#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_data_structure_2.h>
#include <CGAL/Segment_Voronoi_diagram_vertex_base_2.h>
#include <CGAL/Segment_Voronoi_diagram_hierarchy_vertex_base_2.h>
#include <CGAL/Segment_Voronoi_diagram_face_base_2.h>


CGAL_BEGIN_NAMESPACE

//--------------------------------------------------------------------
//--------------------------------------------------------------------

// parameterization of the hierarchy
const unsigned int svd_hierarchy_2__ratio    = 30;
const unsigned int svd_hierarchy_2__minsize  = 20;
const unsigned int svd_hierarchy_2__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template < class Gt,
  class PContainer = Point_container<typename Gt::Point_2>,
  class Tds = Segment_Voronoi_diagram_data_structure_2<
    Segment_Voronoi_diagram_hierarchy_vertex_base_2<
       Segment_Voronoi_diagram_vertex_base_2<Gt,
                                 typename PContainer::Point_handle> >,
    Segment_Voronoi_diagram_face_base_2<Gt> > >
class Segment_Voronoi_diagram_hierarchy_2
  : public Segment_Voronoi_diagram_2<Gt,PContainer,Tds>
{
public:
  typedef Segment_Voronoi_diagram_2<Gt,PContainer,Tds>  Base;

  typedef typename Base::Geom_traits        Geom_traits;

  typedef typename Geom_traits::Point_2     Point;
  typedef typename Geom_traits::Segment_2   Segment;
  typedef typename Geom_traits::Site_2      Site;

  typedef typename Base::Vertex_handle      Vertex_handle;
  typedef typename Base::Face_handle        Face_handle;
  typedef typename Base::Edge               Edge;

  typedef typename Base::Vertex_circulator         Vertex_circulator;
  typedef typename Base::Edge_circulator           Edge_circulator;
  typedef typename Base::Face_circulator           Face_circulator;

  typedef typename Base::All_faces_iterator        All_faces_iterator;
  typedef typename Base::Finite_faces_iterator     Finite_faces_iterator;
  typedef typename Base::All_vertices_iterator     All_vertices_iterator;
  typedef typename Base::Finite_vertices_iterator  Finite_vertices_iterator;
  typedef typename Base::All_edges_iterator        All_edges_iterator;
  typedef typename Base::Finite_edges_iterator     Finite_edges_iterator;

  struct Vertex_iterator {};

private:
  static const int UNDEFINED_LEVEL;

  typedef typename Base::Storage_site_2            Storage_site_2;

private:
  // here is the stack of triangulations which form the hierarchy
  Base*   hierarchy[svd_hierarchy_2__maxlevel];
  Random random; // random generator

public:
  Segment_Voronoi_diagram_hierarchy_2
  (const Geom_traits& traits = Geom_traits());
  Segment_Voronoi_diagram_hierarchy_2
  (const Segment_Voronoi_diagram_hierarchy_2& svd);

  Segment_Voronoi_diagram_hierarchy_2 &operator=
  (const  Segment_Voronoi_diagram_hierarchy_2& svd);

  ~Segment_Voronoi_diagram_hierarchy_2();

  //Helping
  void copy_triangulation
  (const Segment_Voronoi_diagram_hierarchy_2 &svd);

  void clear();

  // CHECKING
  bool is_valid(bool verbose = true, int level = 1) const;

public:
   // insertion of a point/segment

  Vertex_handle  insert(const Point& p);
  Vertex_handle  insert(const Point& p0, const Point& p1);

  template< class Input_iterator >
  void insert(Input_iterator first, Input_iterator beyond,
	      bool do_shuffle = true)
  {
    if ( do_shuffle ) {
      std::random_shuffle(first, beyond);
    }

    for (Input_iterator it = first; it != beyond; ++it) {
      insert(*it, UNDEFINED_LEVEL);
    }
  }

  Vertex_handle insert(const Point& p, Vertex_handle)
  {
    return insert(p);
  }

  Vertex_handle insert(const Point& p1, const Point& p2,
		       Vertex_handle)
  {
    return insert(p1, p2);
  }

private:
  Vertex_handle insert_point(const Point& p, int level);
  void          insert_point(const Point& p, int level,
			     Vertex_handle* vertices);

  Vertex_handle insert_segment(const Point& p0, const Point& p1,
			       int level); 
  void          insert(const Site& p, int level,
		       Vertex_handle* vertices);

#if 0
public:
  // removal
  void remove(Vertex_handle  v, bool remove_endpoints = true);
#endif

public:
  // nearest neighbor
  Vertex_handle  nearest_neighbor(const Point& p,
				  bool force_point = false) const;

  Vertex_handle  nearest_neighbor(const Point& p, Vertex_handle)
  {
    return nearest_neighbor(p);
  }

private:
  void nearest_neighbor(const Site& p,
			Vertex_handle vnear[svd_hierarchy_2__maxlevel],
			bool force_point) const; 
  int random_level();
};

template<class Gt, class P, class Tds>
const int Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::UNDEFINED_LEVEL = -1;

//**************************************************************************
//**************************************************************************

template<class Gt, class P, class Tds>
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
Segment_Voronoi_diagram_hierarchy_2(const Geom_traits& traits)
  : Base(traits), random((long)0)
{ 
  hierarchy[0] = this; 
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Base(traits);
}


// copy constructor duplicates vertices and faces
template<class Gt, class P, class Tds>
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
Segment_Voronoi_diagram_hierarchy_2
(const Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds> &svd)
    : Base(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i = 1; i < svd_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Base(svd.geom_traits());
  copy_triangulation(svd);
} 
 

//Assignement
template<class Gt, class P, class Tds>
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds> &
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
operator=(const Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds> &svd)
{
  copy_triangulation(svd);
  return *this;
}

template<class Gt, class P, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::   
copy_triangulation
(const Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds> &svd)
{
  std::map< Vertex_handle, Vertex_handle > V;
  {
    for(int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
      //      hierarchy[i]->copy_triangulation(*svd.hierarchy[i]);
      *(hierarchy[i]) = *svd.hierarchy[i];
    }
  }
  //up and down have been copied in straightforward way
  // compute a map at lower level
  {
    for(All_vertices_iterator it = hierarchy[0]->all_vertices_begin(); 
	it != hierarchy[0]->all_vertices_end(); ++it) {
      if ( it->up() != NULL ) V[ it->up()->down() ] = it;
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
	if ( it->up() != NULL ) V[ it->up()->down() ] = it;
      }
    }
  }

  // copy the point container
  hierarchy[0]->pc_ = svd.hierarchy[0]->pc_;
}

template<class Gt, class P, class Tds>
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>:: 
~Segment_Voronoi_diagram_hierarchy_2()
{
  clear();
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i){ 
    delete hierarchy[i];
  }
}

template<class Gt, class P, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>:: 
clear()
{
  for(unsigned int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->clear();
  }
}

template<class Gt, class P, class Tds>
inline
typename Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
insert(const Point& p)
{
  return insert_point(p, UNDEFINED_LEVEL);
}

template<class Gt, class P, class Tds>
inline
typename Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
insert_point(const Point& p, int level)

{
  if ( level == UNDEFINED_LEVEL ) {
    level = random_level();
  }

  Vertex_handle vertices[svd_hierarchy_2__maxlevel];
  
  insert_point(p, level, vertices);

  return vertices[0];
}

template<class Gt, class P, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
insert_point(const Point& p, int level,	Vertex_handle* vertices)
{
  CGAL_precondition( level != UNDEFINED_LEVEL );

  Vertex_handle vertex;
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];

  nearest_neighbor(p, vnear, false);

  vertex = hierarchy[0]->insert(p, vnear[0]);

  if ( vertices != NULL ) { vertices[0] = vertex; }

  CGAL_assertion( vertex != Vertex_handle(NULL) );

  // insert at other levels
  Vertex_handle previous = vertex;
      
  int k = 1;
  while (k <= level ) {
    vertex = hierarchy[k]->insert(p, vnear[k]);

    CGAL_assertion( vertex != Vertex_handle(NULL) );

    if ( vertices != NULL ) { vertices[k] = vertex; }

    vertex->set_down(previous); // link with other levels
    previous->set_up(vertex);
    previous = vertex;
    k++;
  }
}


template<class Gt, class P, class Tds>
inline
typename Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
insert(const Point& p0, const Point& p1)
{
  return insert_segment(p0, p1, UNDEFINED_LEVEL);
}


template<class Gt, class P, class Tds>
typename Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
insert_segment(const Point& p0, const Point& p1, int level)
{
  if ( level == UNDEFINED_LEVEL ) {
    level = random_level();
  }

  Site t(Segment(p0, p1));

  if ( is_degenerate_segment(t) ) {
    return insert_point(p0, level);
  }

  Vertex_handle vertices0[svd_hierarchy_2__maxlevel];
  Vertex_handle vertices1[svd_hierarchy_2__maxlevel];

  insert_point(p0, level, vertices0);
  insert_point(p1, level, vertices1);

  CGAL_assertion( vertices0[0] != Vertex_handle(NULL) );
  CGAL_assertion( vertices1[0] != Vertex_handle(NULL) );

  Storage_site_2 ss = create_storage_site(vertices0[0], vertices1[0]);

  Vertex_handle vertex;

  if ( hierarchy[0]->number_of_vertices() == 2 ) {
    vertex = hierarchy[0]->insert_third(vertices0[0], vertices1[0]);
  } else {
    vertex = hierarchy[0]->insert_segment2(t, ss, vertices0[0], false);
  }

  CGAL_assertion( vertex != Vertex_handle(NULL) );

#if 0
  // this is the case when the new site is a segment and it intersects
  // existing segments
  if ( vertex == Vertex_handle(NULL) ) {
    return vertex;
  }
#endif

  // MK:: by doing this the hierarchy stores segments only at the
  //      bottom-most level
  if ( this->intersection_flag ) {
    return vertex;
  }

  // insert at other levels
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;
      
  int k = 1;
  while (k <= level ){
    if ( hierarchy[k]->number_of_vertices() == 2 ) {
      vertex = hierarchy[k]->insert_third(vertices0[k], vertices1[k]);
    } else {
      vertex = hierarchy[k]->insert_segment2(t, ss, vertices0[k], false);
    }

    CGAL_assertion( vertex != Vertex_handle(NULL) );

    vertex->set_down(previous); // link with level above
    previous->set_up(vertex);
    previous = vertex;
    k++;
  }
  return first;
}

#if 0
template<class Gt, class P, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
remove(Vertex_handle v, bool remove_endpoints)
{
  void* u = v->up();
  int l = 0;
  while ( true ) {
    hierarchy[l++]->remove(v, remove_endpoints);
    if ( !u )  { break; }
    if( l > svd_hierarchy_2__maxlevel )   { break; }
    v = u;
    u = v->up();
  }
}
#endif

template<class Gt, class P, class Tds>
typename Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::Vertex_handle 
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
nearest_neighbor(const Point& p, bool force_point) const
{
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];
  nearest_neighbor(Site(p), vnear, force_point);
  return vnear[0];
}

template<class Gt, class P, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
nearest_neighbor(const Site& p,
		 Vertex_handle vnear[svd_hierarchy_2__maxlevel],
		 bool force_point) const
{
  CGAL_precondition( p.is_point() );

  Vertex_handle nearest(NULL);
  int level  = svd_hierarchy_2__maxlevel;

  // find the highest level with enough vertices
  while ( hierarchy[--level]->number_of_vertices() 
	  < svd_hierarchy_2__minsize ) {
    if ( !level ) break;  // do not go below 0
  }
  for (unsigned int i = level + 1; i < svd_hierarchy_2__maxlevel; i++) {
    vnear[i] = Vertex_handle(NULL);
  }

  while ( level > 0 ) {
    vnear[level] = nearest =
      hierarchy[level]->nearest_neighbor(p, nearest);  

    CGAL_assertion( !hierarchy[level]->is_infinite(vnear[level]) );
    CGAL_assertion( vnear[level] != Vertex_handle(NULL) );
    // go at the same vertex on level below
    nearest = nearest->down();
    --level;
  }
  vnear[0] = hierarchy[level]->nearest_neighbor(p, nearest);
  // at level 0
}

template<class Gt, class P, class Tds>
int
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>::
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


template<class Gt, class P, class Tds>
bool
Segment_Voronoi_diagram_hierarchy_2<Gt,P,Tds>:: 
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
    result = result && ( it->down() == 0 );
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



CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_HIERARCHY_2_H

