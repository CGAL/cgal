// Copyright (c) 1998, 2001, 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL$
// $Id$
//
// Author(s)     : Olivier Devillers <Olivier.Devillers@sophia.inria.fr>
//                 Sylvain Pion

#ifndef CGAL_TRIANGULATION_HIERARCHY_3_H
#define CGAL_TRIANGULATION_HIERARCHY_3_H

#include <CGAL/basic.h>
#include <CGAL/Random.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_3.h>

CGAL_BEGIN_NAMESPACE

template < class Tr >
class Triangulation_hierarchy_3
  : public Tr
{
  // parameterization of the hierarchy
  // maximal number of points is 30^5 = 24 millions !
  enum { ratio = 30 };
  enum { minsize = 20};
  enum { maxlevel = 5};

public:
  typedef Tr                                   Tr_Base;
  typedef typename Tr_Base::Geom_traits        Geom_traits;
  typedef typename Tr_Base::Point              Point;
  typedef typename Tr_Base::size_type          size_type;
  typedef typename Tr_Base::Vertex_handle      Vertex_handle;
  typedef typename Tr_Base::Cell_handle        Cell_handle;
  typedef typename Tr_Base::Vertex_iterator    Vertex_iterator;
  typedef typename Tr_Base::Vertex             Vertex;
  typedef typename Tr_Base::Locate_type        Locate_type;
  typedef typename Tr_Base::Finite_vertices_iterator  Finite_vertices_iterator;
  typedef typename Tr_Base::Finite_cells_iterator     Finite_cells_iterator;
  typedef typename Tr_Base::Finite_facets_iterator    Finite_facets_iterator;
  typedef typename Tr_Base::Finite_edges_iterator     Finite_edges_iterator;

  using Tr_Base::number_of_vertices;
  using Tr_Base::geom_traits;

private:
  // here is the stack of triangulations which form the hierarchy
  Tr_Base*   hierarchy[maxlevel];
  Random     random; // random generator

public:
  Triangulation_hierarchy_3(const Geom_traits& traits = Geom_traits());
  Triangulation_hierarchy_3(const Triangulation_hierarchy_3& tr);

  template < typename InputIterator >
  Triangulation_hierarchy_3(InputIterator first, InputIterator last,
                            const Geom_traits& traits = Geom_traits())
    : Tr_Base(traits), random((long)0)
  {
      hierarchy[0] = this;
      for(int i=1; i<maxlevel; ++i)
          hierarchy[i] = new Tr_Base(traits);
      insert(first, last);
  }

  Triangulation_hierarchy_3 & operator=(const Triangulation_hierarchy_3& tr)
  {
    Triangulation_hierarchy_3 tmp(tr);
    swap(tmp);
    return *this;
  }

  ~Triangulation_hierarchy_3();

  void swap(Triangulation_hierarchy_3 &tr);
  void clear();

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;

  // INSERT REMOVE
  Vertex_handle insert(const Point &p, Vertex_handle hint)
  {
    return insert(p, hint == Vertex_handle() ? this->infinite_cell() : hint->cell());
  }
  Vertex_handle insert(const Point &p, Cell_handle start = Cell_handle ());
  Vertex_handle insert(const Point &p, Locate_type lt, Cell_handle loc,
      int li, int lj);

  template < class InputIterator >
  int insert(InputIterator first, InputIterator last)
  {
      int n = number_of_vertices();

      std::vector<Point> points (first, last);
      std::random_shuffle (points.begin(), points.end());
      spatial_sort (points.begin(), points.end(), geom_traits());

      // hints[i] is the vertex of the previously inserted point in level i.
      // Thanks to spatial sort, they are better hints than what the hierarchy
      // would give us.
      Vertex_handle hints[maxlevel];
      for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
              p != end; ++p)
      {
          int vertex_level = random_level();

          Vertex_handle v = hints[0] = hierarchy[0]->insert (*p, hints[0]);
          Vertex_handle prev = v;

          for (int level = 1; level <= vertex_level; ++level) {
              v = hints[level] = hierarchy[level]->insert (*p, hints[level]);
              v->set_down (prev);
              prev->set_up (v);
              prev = v;
          }
      }

      return number_of_vertices() - n;
  }

  // bool only for backward compatibility, we document void.
  bool remove(Vertex_handle v);

  template < typename InputIterator >
  int remove(InputIterator first, InputIterator beyond)
  {
    int n = number_of_vertices();
    while (first != beyond) {
      remove(*first);
      ++first;
    }
    return n - number_of_vertices();
  }

  Vertex_handle move_point(Vertex_handle v, const Point & p);

  //LOCATE
  Cell_handle locate(const Point& p, Locate_type& lt, int& li, int& lj,
                     Vertex_handle hint) const
  {
    return locate(p, lt, li, lj, hint == Vertex_handle() ? this->infinite_cell() : hint->cell());
  }
  Cell_handle locate(const Point& p, Vertex_handle hint) const
  {
    return locate(p, hint == Vertex_handle() ? this->infinite_cell() : hint->cell());
  }

  Cell_handle locate(const Point& p, Locate_type& lt, int& li, int& lj,
          Cell_handle start = Cell_handle ()) const;
  Cell_handle locate(const Point& p, Cell_handle start = Cell_handle ()) const;

  Vertex_handle
  nearest_vertex(const Point& p, Cell_handle start = Cell_handle()) const;

private:

  struct locs {
      Cell_handle pos;
      int li, lj;
      Locate_type lt;
  };

  void locate(const Point& p, Locate_type& lt, int& li, int& lj,
	      locs pos[maxlevel], Cell_handle start = Cell_handle ()) const;
  int random_level();
};


template <class Tr >
Triangulation_hierarchy_3<Tr>::
Triangulation_hierarchy_3(const Geom_traits& traits)
  : Tr_Base(traits), random((long)0)
{
  hierarchy[0] = this;
  for(int i=1;i<maxlevel;++i)
    hierarchy[i] = new Tr_Base(traits);
}

// copy constructor duplicates vertices and cells
template <class Tr>
Triangulation_hierarchy_3<Tr>::
Triangulation_hierarchy_3(const Triangulation_hierarchy_3<Tr> &tr)
    : Tr_Base(tr), random((long)0)
{
  hierarchy[0] = this;
  for(int i=1; i<maxlevel; ++i)
    hierarchy[i] = new Tr_Base(*tr.hierarchy[i]);

  // up and down have been copied in straightforward way
  // compute a map at lower level

  std::map< Vertex_handle, Vertex_handle > V;

  for( Finite_vertices_iterator it=hierarchy[0]->finite_vertices_begin();
       it != hierarchy[0]->finite_vertices_end(); ++it)
    if (it->up() != Vertex_handle())
      V[ it->up()->down() ] = it;

  for(int j=1; j<maxlevel; ++j) {
    for( Finite_vertices_iterator it=hierarchy[j]->finite_vertices_begin();
	 it != hierarchy[j]->finite_vertices_end(); ++it) {
	// down pointer goes in original instead in copied triangulation
	it->set_down(V[it->down()]);
	// make reverse link
	it->down()->set_up( it );
	// make map for next level
	if (it->up() != Vertex_handle())
	    V[ it->up()->down() ] = it;
    }
  }
}

template <class Tr>
void
Triangulation_hierarchy_3<Tr>::
swap(Triangulation_hierarchy_3<Tr> &tr)
{
  Tr_Base::swap(tr);
  for(int i=1; i<maxlevel; ++i)
      std::swap(hierarchy[i], tr.hierarchy[i]);
}

template <class Tr>
Triangulation_hierarchy_3<Tr>::
~Triangulation_hierarchy_3()
{
  clear();
  for(int i=1; i<maxlevel; ++i)
    delete hierarchy[i];
}

template <class Tr>
void
Triangulation_hierarchy_3<Tr>::
clear()
{
        for(int i=0;i<maxlevel;++i)
	hierarchy[i]->clear();
}

template <class Tr>
bool
Triangulation_hierarchy_3<Tr>::
is_valid(bool verbose, int level) const
{
  bool result = true;

  // verify correctness of triangulation at all levels
  for(int i=0; i<maxlevel; ++i)
	result = result && hierarchy[i]->is_valid(verbose, level);

  // verify that lower level has no down pointers
  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin();
       it != hierarchy[0]->finite_vertices_end(); ++it)
    result = result && (it->down() == Vertex_handle());

  // verify that other levels has down pointer and reciprocal link is fine
  for(int j=1; j<maxlevel; ++j)
    for( Finite_vertices_iterator it = hierarchy[j]->finite_vertices_begin();
	 it != hierarchy[j]->finite_vertices_end(); ++it)
      result = result && &*(it) == &*(it->down()->up());

  // verify that other levels has down pointer and reciprocal link is fine
  for(int k=0; k<maxlevel-1; ++k)
    for( Finite_vertices_iterator it = hierarchy[k]->finite_vertices_begin();
	 it != hierarchy[k]->finite_vertices_end(); ++it)
      result = result && ( it->up() == Vertex_handle() ||
	        &*it == &*(it->up())->down() );

  return result;
}

template <class Tr>
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
insert(const Point &p, Cell_handle start)
{
  int vertex_level = random_level();
  Locate_type lt;
  int i, j;
  // locate using hierarchy
  locs positions[maxlevel];
  locate(p, lt, i, j, positions, start);
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
      if (positions[level].pos == Cell_handle())
          vertex = hierarchy[level]->insert(p);
      else
          vertex = hierarchy[level]->insert(p,
	                                    positions[level].lt,
	                                    positions[level].pos,
	                                    positions[level].li,
	                                    positions[level].lj);
    vertex->set_down(previous);// link with level above
    previous->set_up(vertex);
    previous=vertex;
    level++;
  }
  return first;
}

template <class Tr>
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
insert(const Point &p, Locate_type lt, Cell_handle loc, int li, int lj)
{
  int vertex_level = random_level();
  // insert at level 0
  Vertex_handle vertex = hierarchy[0]->insert(p,lt,loc,li,lj);
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;

  if (vertex_level > 0) {
    Locate_type lt;
    int i, j;
    // locate using hierarchy
    locs positions[maxlevel];
    locate(p, lt, i, j, positions, loc);

    int level = 1;
    while (level <= vertex_level ){
      if (positions[level].pos == Cell_handle())
	vertex = hierarchy[level]->insert(p);
      else
	vertex = hierarchy[level]->insert(p,
	    positions[level].lt,
	    positions[level].pos,
	    positions[level].li,
	    positions[level].lj);
      vertex->set_down(previous);// link with level above
      previous->set_up(vertex);
      previous=vertex;
      level++;
    }
  }
  return first;
}

template <class Tr>
bool
Triangulation_hierarchy_3<Tr>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition(v != Vertex_handle());
  for (int l = 0; l < maxlevel; ++l) {
    Vertex_handle u = v->up();
    hierarchy[l]->remove(v);
    if (u == Vertex_handle())
	break;
    v = u;
  }
  return true;
}

template < class Tr >
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
move_point(Vertex_handle v, const Point & p)
{
  CGAL_triangulation_precondition(v != Vertex_handle());
  Vertex_handle old, ret;

  for (int l = 0; l < maxlevel; ++l) {
    Vertex_handle u = v->up();
    Vertex_handle w = hierarchy[l]->move_point(v, p);
    if (l == 0) {
	ret = w;
    }
    else {
	old->set_up(w);
	w->set_down(old);
    }
    if (u == Vertex_handle())
	break;
    old = w;
    v = u;
  }

  return ret;
}

template <class Tr>
inline
typename Triangulation_hierarchy_3<Tr>::Cell_handle
Triangulation_hierarchy_3<Tr>::
locate(const Point& p, Locate_type& lt, int& li, int& lj, Cell_handle start) const
{
  if (start != Cell_handle ()) return Tr_Base::locate (p, lt, li, lj, start);
  locs positions[maxlevel];
  locate(p, lt, li, lj, positions);
  return positions[0].pos;
}

template <class Tr>
inline
typename Triangulation_hierarchy_3<Tr>::Cell_handle
Triangulation_hierarchy_3<Tr>::
locate(const Point& p, Cell_handle start) const
{
  if (start != Cell_handle ()) return Tr_Base::locate (p, start);
  Locate_type lt;
  int li, lj;
  return locate(p, lt, li, lj);
}

template <class Tr>
void
Triangulation_hierarchy_3<Tr>::
locate(const Point& p, Locate_type& lt, int& li, int& lj,
       locs pos[maxlevel], Cell_handle start) const
{
  int level = maxlevel;

  // find the highest level with enough vertices
  while (hierarchy[--level]->number_of_vertices() < (size_type) minsize) {
    if ( ! level)
	break;  // do not go below 0
  }

  for (int i=level+1; i<maxlevel; ++i)
      pos[i].pos = Cell_handle();

  Cell_handle position = Cell_handle();
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
    nearest = nearest->down();
    position = nearest->cell();                // incident cell
    --level;
  }
  if (start != Cell_handle ()) position = start;
  pos[0].pos = hierarchy[0]->locate(p, lt, li, lj, position); // at level 0
  pos[0].lt = lt;
  pos[0].li = li;
  pos[0].lj = lj;
}

template <class Tr>
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
nearest_vertex(const Point& p, Cell_handle start) const
{
    return Tr_Base::nearest_vertex(p, start != Cell_handle() ? start
	                                                     : locate(p));
}

template <class Tr>
int
Triangulation_hierarchy_3<Tr>::
random_level()
{
  int l = 0;
  while ( ! random(ratio) && l < maxlevel-1 )
    ++l;

  return l;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_HIERARCHY_3_H
