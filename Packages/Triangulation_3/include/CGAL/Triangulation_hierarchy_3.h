// ============================================================================
//
// Copyright (c) 1998, 2001 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_hierarchy_3.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Triangulation3
// author(s)     : Olivier Devillers <Olivier.Devillers@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_HIERARCHY_3_H
#define CGAL_TRIANGULATION_HIERARCHY_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_3.h>

CGAL_BEGIN_NAMESPACE

// parameterization of the  hierarchy
//const float Triangulation_hierarchy_3__ratio    = 30.0;
const int Triangulation_hierarchy_3__ratio    = 30;
const int Triangulation_hierarchy_3__minsize  = 20;
const int Triangulation_hierarchy_3__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

template < class Tr >
class Triangulation_hierarchy_3
  : public Tr
{
public:
  typedef Tr                                   Tr_Base;
  typedef typename Tr_Base::Geom_traits        Geom_traits;
  typedef typename Tr_Base::Point              Point;
  typedef typename Tr_Base::Vertex_handle      Vertex_handle;
  typedef typename Tr_Base::Cell_handle        Cell_handle;
  typedef typename Tr_Base::Vertex_iterator    Vertex_iterator;
  typedef typename Tr_Base::Vertex             Vertex;
  typedef typename Tr_Base::Locate_type        Locate_type;
  typedef typename Tr_Base::Finite_vertices_iterator  Finite_vertices_iterator;
  typedef typename Tr_Base::Finite_cells_iterator     Finite_cells_iterator;
  typedef typename Tr_Base::Finite_facets_iterator    Finite_facets_iterator;
  typedef typename Tr_Base::Finite_edges_iterator     Finite_edges_iterator;

private:
  // here is the stack of triangulations which form the hierarchy
  Tr_Base*   hierarchy[Triangulation_hierarchy_3__maxlevel];
  Random     random; // random generator

public:
  Triangulation_hierarchy_3(const Geom_traits& traits = Geom_traits());
  Triangulation_hierarchy_3(const Triangulation_hierarchy_3& tr);

  Triangulation_hierarchy_3 &operator=(const  Triangulation_hierarchy_3& tr);
  ~Triangulation_hierarchy_3();

  //Helping
  void copy_triangulation(const Triangulation_hierarchy_3 &tr);
  void swap(Triangulation_hierarchy_3 &tr);
  void clear();

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;

  // INSERT REMOVE
  Vertex_handle insert(const Point &p);

  Vertex_handle push_back(const Point &p)
  {
      return insert(p);
  }
 
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

  bool remove(Vertex_handle v);

  //LOCATE
  Cell_handle locate(const Point& p, Locate_type& lt, int& li, int& lj) const;
  Cell_handle locate(const Point& p) const;

private:

  struct locs {
      Cell_handle pos;
      int li, lj;
      Locate_type lt;
  };

  void locate(const Point& p, Locate_type& lt, int& li, int& lj,
	      locs pos[Triangulation_hierarchy_3__maxlevel]) const;
  int random_level();

  // added to make the test program of usual triangulations work
  // undocumented
public:

  Vertex_handle insert(const Point& p, Locate_type lt, Cell_handle loc,
	               int li, int lj)
  {
    return Tr_Base::insert(p, lt, loc, li, lj);
  }
  
  Vertex_handle insert(const Point &p, Cell_handle start)
  {
    return Tr_Base::insert(p, start);
  }

  Cell_handle locate(const Point& p, Locate_type& lt, int& li, int& lj,
		     Cell_handle start) const
  {
    return Tr_Base::locate(p, lt, li, lj, start);
  }
};


template <class Tr >
Triangulation_hierarchy_3<Tr>::
Triangulation_hierarchy_3(const Geom_traits& traits)
  : Tr_Base(traits), random((long)0)
{ 
  hierarchy[0] = this; 
  for(int i=1;i<Triangulation_hierarchy_3__maxlevel;++i)
    hierarchy[i] = new Tr_Base(traits);
}

// copy constructor duplicates vertices and cells
template <class Tr>
Triangulation_hierarchy_3<Tr>::
Triangulation_hierarchy_3(const Triangulation_hierarchy_3<Tr> &tr)
    : Tr_Base(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i=1;i<Triangulation_hierarchy_3__maxlevel;++i)
    hierarchy[i] = new Tr_Base(tr.geom_traits());
  copy_triangulation(tr);
} 
 

//Assignement
template <class Tr>
Triangulation_hierarchy_3<Tr> &
Triangulation_hierarchy_3<Tr>::
operator=(const Triangulation_hierarchy_3<Tr> &tr)
{
  copy_triangulation(tr);
  return *this;
}


template <class Tr>
void
Triangulation_hierarchy_3<Tr>::   
copy_triangulation(const Triangulation_hierarchy_3<Tr> &tr)
{
  std::map< const void*, void* > V;

  for(int i=0; i<Triangulation_hierarchy_3__maxlevel; ++i)
    hierarchy[i]->copy_triangulation(*tr.hierarchy[i]);

  // up and down have been copied in straightforward way
  // compute a map at lower level

  for( Finite_vertices_iterator it=hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->finite_vertices_end(); ++it)
    if (it->up())
      V[ ((Vertex*)(it->up()))->down() ] = &(*it);

  for(int j=1; j<Triangulation_hierarchy_3__maxlevel; ++j) {
    for( Finite_vertices_iterator it=hierarchy[j]->finite_vertices_begin();
	 it != hierarchy[j]->finite_vertices_end(); ++it) {
	// down pointer goes in original instead in copied triangulation
	it->set_down(V[it->down()]);
	// make reverse link
	((Vertex*)(it->down()))->set_up( &(*it) );
	// make map for next level
	if (it->up())
	    V[ ((Vertex*)(it->up()))->down() ] = &(*it);
    }
  }
}

template <class Tr>
void
Triangulation_hierarchy_3<Tr>:: 
swap(Triangulation_hierarchy_3<Tr> &tr)
{
  Tr_Base::swap(tr);
  for(int i=1; i<Triangulation_hierarchy_3__maxlevel; ++i)
      std::swap(hierarchy[i], tr.hierarchy[i]);
}

template <class Tr>
Triangulation_hierarchy_3<Tr>:: 
~Triangulation_hierarchy_3()
{
  clear();
  for(int i=1; i<Triangulation_hierarchy_3__maxlevel; ++i)
    delete hierarchy[i];
}

template <class Tr>
void
Triangulation_hierarchy_3<Tr>:: 
clear()
{
        for(int i=0;i<Triangulation_hierarchy_3__maxlevel;++i)
	hierarchy[i]->clear();
}

template <class Tr>
bool
Triangulation_hierarchy_3<Tr>:: 
is_valid(bool verbose, int level) const
{
  bool result = true;
  
  // verify correctness of triangulation at all levels
  for(int i=0; i<Triangulation_hierarchy_3__maxlevel; ++i)
	result = result && hierarchy[i]->is_valid(verbose, level);

  // verify that lower level has no down pointers
  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->finite_vertices_end(); ++it) 
    result = result && (it->down() == NULL);

  // verify that other levels has down pointer and reciprocal link is fine
  for(int j=1; j<Triangulation_hierarchy_3__maxlevel; ++j)
    for( Finite_vertices_iterator it = hierarchy[j]->finite_vertices_begin(); 
	 it != hierarchy[j]->finite_vertices_end(); ++it) 
      result = result && (&*it == (Vertex*)((Vertex*)it->down())->up());

  // verify that other levels has down pointer and reciprocal link is fine
  for(int k=0; k<Triangulation_hierarchy_3__maxlevel-1; ++k)
    for( Finite_vertices_iterator it = hierarchy[k]->finite_vertices_begin(); 
	 it != hierarchy[k]->finite_vertices_end(); ++it) 
      result = result && ( ((Vertex*)it->up() == NULL) ||
	       ( &*it == (Vertex*)((Vertex*)it->up())->down() ));

  return result;
}
  
template <class Tr>
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
insert(const Point &p)
{
  int vertex_level = random_level();
  Locate_type lt;
  int i, j;
  // locate using hierarchy
  locs positions[Triangulation_hierarchy_3__maxlevel];
  locate(p, lt, i, j, positions);
  // insert at level 0
  Vertex_handle vertex = hierarchy[0]->insert(p,
	                                      positions[0].lt,
	                                      positions[0].pos,
	                                      positions[0].li,
	                                      positions[0].lj);
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;

  int level = 1;
  while (level <= vertex_level ){
      if (positions[level].pos == NULL)
          vertex = hierarchy[level]->insert(p);
      else
          vertex = hierarchy[level]->insert(p,
	                                    positions[level].lt,
	                                    positions[level].pos,
	                                    positions[level].li,
	                                    positions[level].lj);
    vertex->set_down((void *) &*previous);// link with level above
    previous->set_up((void *) &*vertex);
    previous=vertex;
    level++;
  }
  return first;
}

template <class Tr>
bool
Triangulation_hierarchy_3<Tr>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition(v != NULL);
  void * u = v->up();
  int l = 0;
  bool result = true;
  while (1) {
    if (! hierarchy[l++]->remove(v))
	result = false;
    if (!u || l>Triangulation_hierarchy_3__maxlevel)
	break;
    v = (Vertex*) u;
    u = v->up();
  }
  return result;
}

template <class Tr>
inline
typename Triangulation_hierarchy_3<Tr>::Cell_handle
Triangulation_hierarchy_3<Tr>::
locate(const Point& p, Locate_type& lt, int& li, int& lj) const
{
  locs positions[Triangulation_hierarchy_3__maxlevel];
  locate(p, lt, li, lj, positions);
  return positions[0].pos;
}

template <class Tr>
inline
typename Triangulation_hierarchy_3<Tr>::Cell_handle 
Triangulation_hierarchy_3<Tr>::
locate(const Point& p) const
{
  Locate_type lt;
  int li, lj;
  return locate(p, lt, li, lj);
}

template <class Tr>
void
Triangulation_hierarchy_3<Tr>::
locate(const Point& p, Locate_type& lt, int& li, int& lj,
       locs pos[Triangulation_hierarchy_3__maxlevel]) const
{
  int level = Triangulation_hierarchy_3__maxlevel;

  // find the highest level with enough vertices
  while (hierarchy[--level]->number_of_vertices() 
	 < Triangulation_hierarchy_3__minsize) {
    if ( ! level)
	break;  // do not go below 0
  }

  for (int i=level+1; i<Triangulation_hierarchy_3__maxlevel; ++i)
      pos[i].pos=0;

  Cell_handle position(NULL);
  while(level > 0) {
    // locate at that level from "position"
    // result is stored in "position" for the next level
    pos[level].pos = position = hierarchy[level]->locate(p,
	                                                 pos[level].lt,
	                                                 pos[level].li,
	                                                 pos[level].lj,
	                                                 position);
    // find the nearest vertex.
    Vertex_handle nearest =
	hierarchy[level]->nearest_vertex_in_cell(p, position);

    // go at the same vertex on level below
    nearest = (Vertex*) nearest->down();
    position = nearest->cell();                // incident cell
    --level;
  }
  pos[0].pos = hierarchy[level]->locate(p, lt, li, lj, position); // at level 0
  pos[0].lt = lt;
  pos[0].li = li;
  pos[0].lj = lj;
}


template <class Tr>
int
Triangulation_hierarchy_3<Tr>::
random_level()
{
  int l = 0;
  while ( ! random(Triangulation_hierarchy_3__ratio) )
    ++l;

  if (l >= Triangulation_hierarchy_3__maxlevel)
    l = Triangulation_hierarchy_3__maxlevel - 1;

  return l;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_HIERARCHY_3_H
