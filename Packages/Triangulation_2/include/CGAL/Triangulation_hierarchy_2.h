// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_hierarchy_2.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Triangulation_2
// author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//                 Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_HIERARCHY_2_H
#define CGAL_TRIANGULATION_HIERARCHY_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>
#include <map>

CGAL_BEGIN_NAMESPACE

// parameterization of the  hierarchy
//const float Triangulation_hierarchy_2__ratio    = 30.0;
const int Triangulation_hierarchy_2__ratio      = 30;
const int   Triangulation_hierarchy_2__minsize  = 20;
const int   Triangulation_hierarchy_2__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

template < class Tr>
class Triangulation_hierarchy_2
: public Tr
{
 public:
  typedef Tr                                   Tr_Base;
  typedef typename Tr_Base::Geom_traits        Geom_traits;
  typedef typename Tr_Base::Point              Point;
  typedef typename Tr_Base::Vertex_handle      Vertex_handle;
  typedef typename Tr_Base::Face_handle        Face_handle;
  typedef typename Tr_Base::Vertex             Vertex;
  typedef typename Tr_Base::Locate_type        Locate_type;
  typedef typename Tr_Base::Finite_vertices_iterator  Finite_vertices_iterator;
  //typedef typename Tr_Base::Finite_faces_iterator     Finite_faces_iterator;

 private:
  // here is the stack of triangulations which form the hierarchy
  Tr_Base*   hierarchy[Triangulation_hierarchy_2__maxlevel];
  Random random; // random generator

public:
  Triangulation_hierarchy_2(const Geom_traits& traits = Geom_traits());
  Triangulation_hierarchy_2(const Triangulation_hierarchy_2& tr);

  Triangulation_hierarchy_2 &operator=(const  Triangulation_hierarchy_2& tr);
  ~Triangulation_hierarchy_2();

  //Helping
  void copy_triangulation(const Triangulation_hierarchy_2 &tr);
  void swap(Triangulation_hierarchy_2 &tr);
  void clear();

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;

  // INSERT REMOVE
  Vertex_handle insert(const Point &p, Face_handle start = Face_handle(NULL) );
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, int li );
  Vertex_handle push_back(const Point &p);
 
  template < class InputIterator >
  int insert(InputIterator first, InputIterator last)
    {
      int n = number_of_vertices();
      while(first != last){
	insert(*first);
	++first;
      }
      return number_of_vertices() - n;
    }

  void remove_degree_3(Vertex_handle  v);
  void remove_first(Vertex_handle  v);
  void remove_second(Vertex_handle v);
  void remove(Vertex_handle  v);

  //LOCATE
  Face_handle
  locate(const Point& p,
	 Locate_type& lt,
	 int& li,
	 Face_handle start = Face_handle(NULL)) const;

  Face_handle
  locate(const Point &p,
	 Face_handle start = Face_handle(NULL)) const;

private:
  void  locate_in_all(const Point& p,
		      Locate_type& lt,
		      int& li,
		      Face_handle
		      pos[Triangulation_hierarchy_2__maxlevel]) const;
  int random_level();
};



template <class Tr >
Triangulation_hierarchy_2<Tr>::
Triangulation_hierarchy_2(const Geom_traits& traits)
  : Tr_Base(traits), random((long)0)
{ 
  hierarchy[0] = this; 
  for(int i=1;i<Triangulation_hierarchy_2__maxlevel;++i)
    hierarchy[i] = new Tr_Base(traits);
}


// copy constructor duplicates vertices and faces
template <class Tr>
Triangulation_hierarchy_2<Tr>::
Triangulation_hierarchy_2(const Triangulation_hierarchy_2<Tr> &tr)
    : Tr_Base(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i=1;i<Triangulation_hierarchy_2__maxlevel;++i)
    hierarchy[i] = new Tr_Base(tr.geom_traits());
  copy_triangulation(tr);
} 
 

//Assignement
template <class Tr>
Triangulation_hierarchy_2<Tr> &
Triangulation_hierarchy_2<Tr>::
operator=(const Triangulation_hierarchy_2<Tr> &tr)
{
  copy_triangulation(tr);
  return *this;
}


template <class Tr>
void
Triangulation_hierarchy_2<Tr>::   
copy_triangulation(const Triangulation_hierarchy_2<Tr> &tr)
{
  std::map<Vertex_handle, Vertex_handle > V;
  {
    for(int i=0;i<Triangulation_hierarchy_2__maxlevel;++i)
    hierarchy[i]->copy_triangulation(*tr.hierarchy[i]);
  }
  //up and down have been copied in straightforward way
  // compute a map at lower level
  {
    for( Finite_vertices_iterator it=hierarchy[0]->finite_vertices_begin(); 
	 it != hierarchy[0]->finite_vertices_end(); ++it) {
      if (it->up() != NULL) V[ it->up()->down() ] = it;
    }
  }
  {
    for(int i=1;i<Triangulation_hierarchy_2__maxlevel;++i) {
      for( Finite_vertices_iterator it=hierarchy[i]->finite_vertices_begin(); 
	   it != hierarchy[i]->finite_vertices_end(); ++it) {
	// down pointer goes in original instead in copied triangulation
	it->set_down(V[it->down()]);
	// make reverse link
	it->down()->set_up(it);
	// make map for next level
	if (it->up()!= NULL ) V[ it->up()->down() ] = it;
      }
    }
  }
}

template <class Tr>
void
Triangulation_hierarchy_2<Tr>:: 
swap(Triangulation_hierarchy_2<Tr> &tr)
{
  Tr_Base* temp;
  Tr_Base::swap(tr);
  for(int i= 1; i<Triangulation_hierarchy_2__maxlevel; ++i){
    temp = hierarchy[i];
    hierarchy[i] = tr.hierarchy[i];
    tr.hierarchy[i]= temp;
  }
}

template <class Tr>
Triangulation_hierarchy_2<Tr>:: 
~Triangulation_hierarchy_2()
{
  clear();
  for(int i= 1; i<Triangulation_hierarchy_2__maxlevel; ++i){ 
    delete hierarchy[i];
  }
}

template <class Tr>
void
Triangulation_hierarchy_2<Tr>:: 
clear()
{
        for(int i=0;i<Triangulation_hierarchy_2__maxlevel;++i)
	hierarchy[i]->clear();
}


template <class Tr>
bool
Triangulation_hierarchy_2<Tr>:: 
is_valid(bool verbose, int level) const
{
  bool result = true;
  int i;
  Finite_vertices_iterator it;
  //verify correctness of triangulation at all levels
  for(i=0;i<Triangulation_hierarchy_2__maxlevel;++i)
	result = result && hierarchy[i]->is_valid(verbose,level);
  //verify that lower level has no down pointers
  for( it = hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->finite_vertices_end(); ++it) 
    result = result && ( it->down() == NULL );
  //verify that other levels have down pointer and reciprocal link is fine
  for(i=1;i<Triangulation_hierarchy_2__maxlevel;++i)
    for( it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->finite_vertices_end(); ++it) 
      result = result && 
	       ( &*(it->down()->up())  ==  &*(it) );
  //verify that levels have up pointer and reciprocal link is fine
  for(i=0;i<Triangulation_hierarchy_2__maxlevel-1;++i)
    for( it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->finite_vertices_end(); ++it) 
      result = result && ( it->up() == NULL ||
	        &*it == &*(it->up())->down() );
  return result;
}

  
template <class Tr>
typename Triangulation_hierarchy_2<Tr>::Vertex_handle
Triangulation_hierarchy_2<Tr>::
insert(const Point &p, Face_handle)
{
  int vertex_level = random_level();
  Locate_type lt;
  int i;
  // locate using hierarchy
  Face_handle positions[Triangulation_hierarchy_2__maxlevel];
  locate_in_all(p,lt,i,positions);
  //insert at level 0
  Vertex_handle vertex=hierarchy[0]->Tr_Base::insert(p,lt,positions[0],i);
  Vertex_handle previous=vertex;
  Vertex_handle first = vertex;
      
  int level  = 1;
  while (level <= vertex_level ){
    vertex=hierarchy[level]->Tr_Base::insert(p,positions[level]);
    vertex->set_down(previous);// link with level above
    previous->set_up(vertex);
    previous=vertex;
    level++;
  }
  return first;
}

template <class Tr>
typename Triangulation_hierarchy_2<Tr>::Vertex_handle
Triangulation_hierarchy_2<Tr>::
insert(const Point& p,
       Locate_type lt,
       Face_handle loc, 
       int li )
{
  int vertex_level = random_level();
  //insert at level 0
  Vertex_handle vertex=hierarchy[0]->Tr_Base::insert(p,lt,loc,li);
  Vertex_handle previous=vertex;
  Vertex_handle first = vertex;

  if (vertex_level > 0) {
    // locate using hierarchy
    Locate_type ltt;
    int lii;
    Face_handle positions[Triangulation_hierarchy_2__maxlevel];
    locate_in_all(p,ltt,lii,positions);
    //insert in higher levels
    int level  = 1;
    while (level <= vertex_level ){
      vertex=hierarchy[level]->Tr_Base::insert(p,positions[level]);
      vertex->set_down(previous);// link with level above
      previous->set_up(vertex);
      previous=vertex;
      level++;
    }
  }
  return first;
}


template <class Tr>
inline
typename Triangulation_hierarchy_2<Tr>::Vertex_handle
Triangulation_hierarchy_2<Tr>::
push_back(const Point &p)
{
  return insert(p);
}

template <class Tr>
void 
Triangulation_hierarchy_2<Tr>::
remove(Vertex_handle v )
{
  Vertex_handle u=v->up();
  int l = 0 ;
  while(1){
    hierarchy[l++]->remove(v);
    if (u == NULL) break; 
    if(l>Triangulation_hierarchy_2__maxlevel) break;
    v=u; u=v->up();
  }
}

template <class Tr>
inline void 
Triangulation_hierarchy_2<Tr>::
remove_degree_3(Vertex_handle v )
{
  remove(v);
}

template <class Tr>
inline void 
Triangulation_hierarchy_2<Tr>::
remove_first(Vertex_handle v )
{
  remove(v);
}

template <class Tr>
inline void 
Triangulation_hierarchy_2<Tr>::
remove_second(Vertex_handle v )
{
  remove(v);
}


template <class Tr>
typename Triangulation_hierarchy_2<Tr>::Face_handle 
Triangulation_hierarchy_2<Tr>::
locate(const Point& p, Locate_type& lt, int& li, Face_handle) const
{
  Face_handle positions[Triangulation_hierarchy_2__maxlevel];
  locate_in_all(p,lt,li,positions);
  return positions[0];
}

template <class Tr>
typename Triangulation_hierarchy_2<Tr>::Face_handle 
Triangulation_hierarchy_2<Tr>::
locate(const Point& p, Face_handle ) const
{
  Locate_type lt;
  int li;
  return locate(p, lt, li);
}

template <class Tr>
void
Triangulation_hierarchy_2<Tr>::
locate_in_all(const Point& p,
    Locate_type& lt,
    int& li,
    Face_handle pos[Triangulation_hierarchy_2__maxlevel]) const
{
  Face_handle position(NULL);
  Vertex_handle nearest(NULL);
  int level  = Triangulation_hierarchy_2__maxlevel;
  typename Geom_traits::Compare_distance_2 
    closer = geom_traits().compare_distance_2_object();

  // find the highest level with enough vertices
  while (hierarchy[--level]->number_of_vertices() 
	 < Triangulation_hierarchy_2__minsize){
    if ( ! level) break;  // do not go below 0
  }
  for (int i=level+1; i<Triangulation_hierarchy_2__maxlevel;++i) pos[i]=0;
  while(level > 0) {
    pos[level]=position=hierarchy[level]->locate(p,position);  
    // locate at that level from "position"
    // result is stored in "position" for the next level
    // find the nearest between vertices 0 and 1
    if (hierarchy[level]->is_infinite(position->vertex(0)))
      nearest = position->vertex(1);
    else if (hierarchy[level]->is_infinite(position->vertex(1)))
      nearest = position->vertex(0);
     else if ( closer(p,
		      position->vertex(0)->point(),
		      position->vertex(1)->point()) == SMALLER)
      nearest = position->vertex(0);
    else
      nearest = position->vertex(1);
    // compare to vertex 2
    if ( !  hierarchy[level]->is_infinite(position->vertex(2)))
      if ( closer( p, 
		   position->vertex(2)->point(),
		   nearest->point()) == SMALLER )
	nearest = position->vertex(2);
    // go at the same vertex on level below
    nearest  = nearest->down();
    position = nearest->face();                // incident face
    --level;
  }
  pos[0]=hierarchy[level]->locate(p,lt,li,position);  // at level 0
}


template <class Tr>
int
Triangulation_hierarchy_2<Tr>::
random_level()
{
  int l = 0;
  while (1) {
    if ( random(Triangulation_hierarchy_2__ratio) ) break;
    ++l;
  }
  if (l >= Triangulation_hierarchy_2__maxlevel)
    l = Triangulation_hierarchy_2__maxlevel -1;
  return l;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_HIERARCHY_2_H
