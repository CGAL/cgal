// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/Segment_Voronoi_diagram_hierarchy_2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_HIERARCHY_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_HIERARCHY_2_H

#include <CGAL/Random.h>
#include <map>
#include <CGAL/Triangulation_hierarchy_2.h>

#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_data_structure_2.h>
#include <CGAL/Segment_Voronoi_diagram_vertex_base_2.h>
#include <CGAL/Segment_Voronoi_diagram_face_base_2.h>

CGAL_BEGIN_NAMESPACE

template < class Vbb>
class Segment_Voronoi_diagram_hierarchy_vertex_base_2
 : public Vbb
{
public:
  typedef Vbb V_Base;
  typedef typename V_Base::Segment_Voronoi_diagram_data_structure_2 Svdds;

  typedef typename V_Base::Point_2            Point_2;
  typedef typename V_Base::Segment_2          Segment_2;
  typedef typename V_Base::Site_2             Site_2;

  typedef Svdds         Segment_Voronoi_diagram_data_structure_2;
  typedef typename Svdds::Vertex_handle       Vertex_handle;
  typedef typename Svdds::Face_handle         Face_handle;

  template < typename SVDDS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<SVDDS2>::Other      Vb2;
    typedef Segment_Voronoi_diagram_hierarchy_vertex_base_2<Vb2>  Other;
  };

  Segment_Voronoi_diagram_hierarchy_vertex_base_2()
    : V_Base(), _up(0), _down(0) {}

  Segment_Voronoi_diagram_hierarchy_vertex_base_2(const Site_2& t,
						  Face_handle f)
    : V_Base(t,f), _up(0), _down(0) {}

#if 0
  Segment_Voronoi_diagram_hierarchy_vertex_base_2(const Point& p,
						  void* f = NULL)
    : V_Base(p,f), _up(0), _down(0) {}

  Segment_Voronoi_diagram_hierarchy_vertex_base_2(const Segment& s,
						  void* f = NULL)
    : V_Base(s,f), _up(0), _down(0) {}
#endif

public:  // for use in hierarchy only
  Vertex_handle up() {return _up;}
  Vertex_handle down() {return _down;}
  void set_up(Vertex_handle u) {_up=u;}
  void set_down(Vertex_handle d) {if (this) _down=d;}

private:
  Vertex_handle _up;    // same vertex one level above
  Vertex_handle _down;  // same vertex one level below
};

//--------------------------------------------------------------------
//--------------------------------------------------------------------

// parameterization of the  hierarchy
const int svd_hierarchy_2__ratio    = 30;
const int svd_hierarchy_2__minsize  = 20;
const int svd_hierarchy_2__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template < class Gt,
  class Tds = Segment_Voronoi_diagram_data_structure_2<
    Segment_Voronoi_diagram_hierarchy_vertex_base_2<
       Segment_Voronoi_diagram_vertex_base_2<Gt> >,
    Segment_Voronoi_diagram_face_base_2<Gt> > >
class Segment_Voronoi_diagram_hierarchy_2
  : public Segment_Voronoi_diagram_2< Gt, Tds >
{
public:
  typedef Segment_Voronoi_diagram_2<Gt, Tds>  Base;

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
  Vertex_handle  insert(const Segment& p);
  Vertex_handle  insert(const Site& t);

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

private:
  Vertex_handle insert(const Site& t, int hierarchy_level);
  void          insert(const Site& p, int hierarchy_level,
		       Vertex_handle* vertices);

public:
  // removal
  void remove(Vertex_handle  v, bool remove_endpoints = true);

public:
  // find nearest neighbor
  Vertex_handle  find_nearest_neighbor(const Point& p,
				       bool force_point = false) const;

private:
  void 
  find_nearest_neighbor(const Site& p,
			Vertex_handle vnear[svd_hierarchy_2__maxlevel],
			bool force_point) const; 
  int random_level();
};

template<class Gt, class Tds>
const int Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::UNDEFINED_LEVEL = -1;

//**************************************************************************
//**************************************************************************

template < class Gt, class Tds>
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
Segment_Voronoi_diagram_hierarchy_2(const Geom_traits& traits)
  : Base(traits), random((long)0)
{ 
  hierarchy[0] = this; 
  for(int i = 1; i < svd_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Base(traits);
}


// copy constructor duplicates vertices and faces
template <class Gt, class Tds>
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
Segment_Voronoi_diagram_hierarchy_2
(const Segment_Voronoi_diagram_hierarchy_2<Gt,Tds> &svd)
    : Base(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i = 1; i < svd_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Base(svd.geom_traits());
  copy_triangulation(svd);
} 
 

//Assignement
template <class Gt, class Tds>
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds> &
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
operator=(const Segment_Voronoi_diagram_hierarchy_2<Gt,Tds> &svd)
{
  copy_triangulation(svd);
  return *this;
}

template <class Gt, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::   
copy_triangulation
(const Segment_Voronoi_diagram_hierarchy_2<Gt,Tds> &svd)
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
}

template <class Gt, class Tds>
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>:: 
~Segment_Voronoi_diagram_hierarchy_2()
{
  clear();
  for(int i = 1; i < svd_hierarchy_2__maxlevel; ++i){ 
    delete hierarchy[i];
  }
}

template <class Gt, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>:: 
clear()
{
  for(int i = 0; i < svd_hierarchy_2__maxlevel; ++i)
    hierarchy[i]->clear();
}


template <class Gt, class Tds>
inline
typename Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
insert(const Point& p)
{
  Site t(p);
  return insert(t, UNDEFINED_LEVEL);
}

template <class Gt, class Tds>
inline
typename Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
insert(const Segment& s)
{
  Site t(s);
  return insert(t, UNDEFINED_LEVEL);
}

template <class Gt, class Tds>
inline
typename Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
insert(const Site& t)
{
  return insert(t, UNDEFINED_LEVEL);
}

template <class Gt, class Tds>
typename Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
insert(const Site& t, int hierarchy_level)
{
  int vertex_level(hierarchy_level);
  if ( hierarchy_level == UNDEFINED_LEVEL ) {
    vertex_level = random_level();
  }

  Vertex_handle vertices[svd_hierarchy_2__maxlevel];

  if ( t.is_point() ) {
    insert(t, vertex_level, vertices);
    return vertices[0];
  }

  if ( t.segment().is_degenerate() ) {
    insert(t.source_site(), vertex_level, vertices);
    return vertices[0];
  }

  CGAL_assertion( t.is_segment() );

  insert(t.source_site(), vertex_level, vertices);
  insert(t.target_site(), vertex_level, NULL);

  Vertex_handle vertex = hierarchy[0]->insert(t, vertices[0], false);

  // this is the case when the new site is a segment and it intersects
  // existing segments
  if ( vertex == Vertex_handle(NULL) ) {
    return vertex;
  }

  // MK:: by doing this the hierarchy stores segments only at the
  //      bottom-most level
  if ( this->intersection_flag && t.is_segment() ) {
    return vertex;
  }

  // insert at other levels
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;
      
  int level = 1;
  while (level <= vertex_level ){
    vertex = hierarchy[level]->insert(t, vertices[level], false);
    vertex->set_down(previous); // link with level above
    previous->set_up(vertex);
    previous = vertex;
    level++;
  }
  return first;
}


template <class Gt, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
insert(const Site& p, int hierarchy_level,
       Vertex_handle* vertices)
{
  CGAL_precondition( p.is_point() );

  int vertex_level(hierarchy_level);
  if ( hierarchy_level == UNDEFINED_LEVEL ) {
    vertex_level = random_level();
  }

  Vertex_handle vertex;
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];

  find_nearest_neighbor(p, vnear, false);

  vertex = hierarchy[0]->insert(p, vnear[0]);
  if ( vertices != NULL ) { vertices[0] = vertex; }

  // insert at other levels
  Vertex_handle previous = vertex;
      
  int level = 1;
  while (level <= vertex_level ){
    vertex = hierarchy[level]->insert(p, vnear[level]);
    if ( vertices != NULL ) { vertices[level] = vertex; }
    vertex->set_down(previous); // link with level above
    previous->set_up(vertex);
    previous = vertex;
    level++;
  }
}

template <class Gt, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
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

template <class Gt, class Tds>
typename Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::Vertex_handle 
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
find_nearest_neighbor(const Point& p, bool force_point) const
{
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];
  find_nearest_neighbor(p, vnear, force_point);
  return vnear[0];
}

template <class Gt, class Tds>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
find_nearest_neighbor(const Site& p,
		      Vertex_handle vnear[svd_hierarchy_2__maxlevel],
		      bool force_point)
  const
{
  CGAL_precondition( p.is_point() );

  Vertex_handle nearest = 0;
  int level  = svd_hierarchy_2__maxlevel;

  // find the highest level with enough vertices
  while ( hierarchy[--level]->number_of_vertices() 
	  < svd_hierarchy_2__minsize ) {
    if ( !level ) break;  // do not go below 0
  }
  for (int i = level+1; i < svd_hierarchy_2__maxlevel; ++i){
    vnear[i] = NULL;
  }

  while ( level > 0 ) {
    vnear[level] = nearest =
      hierarchy[level]->nearest_neighbor(p, nearest);  

    CGAL_assertion( !hierarchy[level]->is_infinite(vnear[level]) );
    // go at the same vertex on level below
    nearest = nearest->down();
    --level;
  }
  vnear[0] =
    hierarchy[level]->nearest_neighbor(p, nearest);  
  // at level 0
}

template <class Gt, class Tds>
int
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>::
random_level()
{
  int l = 0;
  while (1) {
    if ( random(svd_hierarchy_2__ratio) ) break;
    ++l;
  }
  if (l >= svd_hierarchy_2__maxlevel)
    l = svd_hierarchy_2__maxlevel -1;
  return l;
}


template<class Gt, class Tds>
bool
Segment_Voronoi_diagram_hierarchy_2<Gt,Tds>:: 
is_valid(bool verbose, int level) const
{
  bool result(true);

  //verify correctness of triangulation at all levels
  for(int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
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
  for(int i = 1; i < svd_hierarchy_2__maxlevel; ++i) {
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

