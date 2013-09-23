// Copyright (c) 1998, 2001, 2003, 2009  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Olivier Devillers <Olivier.Devillers@sophia.inria.fr>
//                 Sylvain Pion
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_TRIANGULATION_HIERARCHY_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_HIERARCHY_3_H

#include <CGAL/basic.h>
#include <CGAL/array.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_3.h>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace CGAL {

template < class PTr >
class Periodic_3_triangulation_hierarchy_3
  : public PTr
{
  // parameterization of the hierarchy
  // maximal number of points is 30^5 = 24 millions !
  enum { ratio = 30 };
  enum { minsize = 20};
  enum { maxlevel = 5};

public:
  typedef PTr                                  PTr_Base;
  typedef typename PTr_Base::Geom_traits       Geom_traits;
  typedef typename PTr_Base::Point             Point;
  typedef typename PTr_Base::Iso_cuboid        Iso_cuboid;
  typedef typename PTr_Base::size_type         size_type;
  typedef typename PTr_Base::Vertex_handle     Vertex_handle;
  typedef typename PTr_Base::Cell_handle       Cell_handle;
  typedef typename PTr_Base::Vertex_iterator   Vertex_iterator;
  typedef typename PTr_Base::Vertex            Vertex;
  typedef typename PTr_Base::Locate_type       Locate_type;
  typedef typename PTr_Base::Cell_iterator     Cell_iterator;
  typedef typename PTr_Base::Facet_iterator    Facet_iterator;
  typedef typename PTr_Base::Edge_iterator     Edge_iterator;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using PTr_Base::number_of_vertices;
  using PTr_Base::geom_traits;
  using PTr_Base::is_virtual;
#endif

private:
  // here is the stack of triangulations which form the hierarchy
  PTr_Base*  hierarchy[maxlevel];
  boost::rand48 random;
  int level_mult_cover;

public:
  Periodic_3_triangulation_hierarchy_3(
      const Iso_cuboid& domain = Iso_cuboid(0,0,0,1,1,1),
      const Geom_traits& traits = Geom_traits());

  Periodic_3_triangulation_hierarchy_3(
      const Periodic_3_triangulation_hierarchy_3& tr);

  template < typename InputIterator >
  Periodic_3_triangulation_hierarchy_3(InputIterator first, InputIterator last,
      const Iso_cuboid& domain = Iso_cuboid(0,0,0,1,1,1),
      const Geom_traits& traits = Geom_traits())
    : PTr_Base(domain,traits), level_mult_cover(0)
  {
      hierarchy[0] = this; 
      for(int i=1; i<maxlevel; ++i)
	hierarchy[i] = new PTr_Base(domain,traits);
      insert(first, last);
  }
 
  Periodic_3_triangulation_hierarchy_3 & operator=(
      const Periodic_3_triangulation_hierarchy_3& tr)
  {
    Periodic_3_triangulation_hierarchy_3 tmp(tr);
    swap(tmp);
    return *this;
  }

  ~Periodic_3_triangulation_hierarchy_3();

  void swap(Periodic_3_triangulation_hierarchy_3 &tr);
  void clear();

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;

  // INSERT REMOVE
  Vertex_handle insert(const Point &p, Cell_handle start = Cell_handle ());
  Vertex_handle insert(const Point &p, Locate_type lt, Cell_handle loc,
      int li, int lj);

  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last, bool = false)
  {
    size_type n = number_of_vertices();

    std::vector<Point> points (first, last);
    std::random_shuffle (points.begin(), points.end());
    spatial_sort (points.begin(), points.end(), geom_traits());

    // hints[i] is the cell of the previously inserted point in level i.
    // Thanks to spatial sort, they are better hints than what the hierarchy
    // would give us.
    Cell_handle hints[maxlevel];
    for (typename std::vector<Point>::const_iterator p = points.begin(),
	   end = points.end(); p != end; ++p) {
      int vertex_level = random_level();
	
      Vertex_handle v = hierarchy[0]->insert (*p, hints[0]);
      hints[0] = v->cell();
	
      Vertex_handle prev = v;
	
      for (int level = 1; level <= vertex_level; ++level) {
	v = hierarchy[level]->insert (*p, hints[level]);
	hints[level] = v->cell();
	  
	v->set_down (prev);
	if (hierarchy[level]->number_of_sheets()[0] != 1) {
	  std::vector<Vertex_handle> vtc 
	    = hierarchy[level]->periodic_copies(v);
	  for (unsigned int i=0 ; i<vtc.size() ; i++) vtc[i]->set_down(prev);
	}

	prev->set_up (v);
	prev = v;
      }
    }
    return number_of_vertices() - n;
  }

  void remove(Vertex_handle v);

  template < typename InputIterator >
  std::ptrdiff_t remove(InputIterator first, InputIterator beyond) {
    size_type n = number_of_vertices();
    while (first != beyond) {
      remove(*first);
      ++first;
    }
    return n-number_of_vertices();
  }

  Vertex_handle move_point(Vertex_handle v, const Point & p);

  //LOCATE
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

  // added to make the test program of usual triangulations work
  // undocumented
public:

};


template <class PTr >
Periodic_3_triangulation_hierarchy_3<PTr>::
Periodic_3_triangulation_hierarchy_3(
    const Iso_cuboid& domain, const Geom_traits& traits)
  : PTr_Base(domain, traits), level_mult_cover(0)
{ 
  hierarchy[0] = this; 
  for(int i=1;i<maxlevel;++i)
    hierarchy[i] = new PTr_Base(domain,traits);
}

// copy constructor duplicates vertices and cells
template <class PTr>
Periodic_3_triangulation_hierarchy_3<PTr>::
Periodic_3_triangulation_hierarchy_3(
    const Periodic_3_triangulation_hierarchy_3<PTr> &tr)
  : PTr_Base(tr), level_mult_cover(tr.level_mult_cover)
{ 
  hierarchy[0] = this;
  for(int i=1; i<maxlevel; ++i)
    hierarchy[i] = new PTr_Base(*tr.hierarchy[i]);

  // up and down have been copied in straightforward way
  // compute a map at lower level

  std::map< Vertex_handle, Vertex_handle > V;

  for( Vertex_iterator it=hierarchy[0]->vertices_begin(); 
       it != hierarchy[0]->vertices_end(); ++it) {
    if (hierarchy[0]->is_virtual(it)) continue;
    if (it->up() != Vertex_handle())
      V[ it->up()->down() ] = it;
  }

  for(int j=1; j<maxlevel; ++j) {
    for( Vertex_iterator it=hierarchy[j]->vertices_begin();
	 it != hierarchy[j]->vertices_end(); ++it) {
      if (hierarchy[j]->is_virtual(it)) {
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
}

template <class PTr>
void
Periodic_3_triangulation_hierarchy_3<PTr>:: 
swap(Periodic_3_triangulation_hierarchy_3<PTr> &tr)
{
  PTr_Base::swap(tr);
  for(int i=1; i<maxlevel; ++i)
      std::swap(hierarchy[i], tr.hierarchy[i]);
}

template <class PTr>
Periodic_3_triangulation_hierarchy_3<PTr>:: 
~Periodic_3_triangulation_hierarchy_3()
{
  clear();
  for(int i=1; i<maxlevel; ++i)
    delete hierarchy[i];
}

template <class PTr>
void
Periodic_3_triangulation_hierarchy_3<PTr>:: 
clear()
{
        for(int i=0;i<maxlevel;++i)
	hierarchy[i]->clear();
}

template <class PTr>
bool
Periodic_3_triangulation_hierarchy_3<PTr>:: 
is_valid(bool verbose, int level) const
{
  bool result = true;
  
  // verify correctness of triangulation at all levels
  for(int i=0; i<maxlevel; ++i)
	result = result && hierarchy[i]->is_valid(verbose, level);

  // verify that lower level has no down pointers
  for( Vertex_iterator it = hierarchy[0]->vertices_begin(); 
       it != hierarchy[0]->vertices_end(); ++it) 
    if (!hierarchy[0]->is_virtual(it))
      result = result && (it->down() == Vertex_handle());

  // verify that other levels has down pointer and reciprocal link is fine
  for(int j=1; j<maxlevel; ++j)
    for( Vertex_iterator it = hierarchy[j]->vertices_begin(); 
	 it != hierarchy[j]->vertices_end(); ++it) 
      if (!hierarchy[j]->is_virtual(it))
	result = result && &*(it) == &*(it->down()->up());

  // verify that other levels has down pointer and reciprocal link is fine
  for(int k=0; k<maxlevel-1; ++k)
    for( Vertex_iterator it = hierarchy[k]->vertices_begin(); 
	 it != hierarchy[k]->vertices_end(); ++it) 
      if (!hierarchy[k]->is_virtual(it))
	result = result && ( it->up() == Vertex_handle() ||
	    &*it == &*(it->up())->down() );

  return result;
}
  
template <class PTr>
typename Periodic_3_triangulation_hierarchy_3<PTr>::Vertex_handle
Periodic_3_triangulation_hierarchy_3<PTr>::
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
    if (hierarchy[level]->number_of_sheets()[0] != 1) {
      std::vector<Vertex_handle> vtc 
	= hierarchy[level]->periodic_copies(vertex);
      for (unsigned int i=0 ; i<vtc.size() ; i++) vtc[i]->set_down(previous);
    }
    previous->set_up(vertex);
    previous=vertex;
    level++;
  }
  return first;
}

template <class PTr>
typename Periodic_3_triangulation_hierarchy_3<PTr>::Vertex_handle
Periodic_3_triangulation_hierarchy_3<PTr>::
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
      if (hierarchy[level]->number_of_sheets()[0] != 1) {
	std::vector<Vertex_handle> vtc 
	  = hierarchy[level]->periodic_copies(vertex);
	for (unsigned int i=0 ; i<vtc.size() ; i++) vtc[i]->set_down(previous);
      }
      previous->set_up(vertex);
      previous=vertex;
      level++;
    }
  }
  return first;
}

template <class PTr>
void
Periodic_3_triangulation_hierarchy_3<PTr>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_expensive_precondition(is_vertex(v));
  for (int l = 0; l < maxlevel; ++l) {
    Vertex_handle u = v->up();
    hierarchy[l]->remove(v);
    if (u == Vertex_handle())
	break;
    v = u;
  }
}

template < class PTr >
typename Periodic_3_triangulation_hierarchy_3<PTr>::Vertex_handle
Periodic_3_triangulation_hierarchy_3<PTr>::
move_point(Vertex_handle v, const Point & p)
{
  CGAL_triangulation_precondition(v != Vertex_handle());
  Vertex_handle old, ret;

  for (int l = 0; l < maxlevel; ++l) {
    Vertex_handle u = v->up();
    CGAL_triangulation_assertion(hierarchy[l]->is_valid());
    Vertex_handle w = hierarchy[l]->move_point(v, p);
    if (l == 0) {
	ret = w;
    }
    else {
	old->set_up(w);
	w->set_down(old);
	if (hierarchy[l]->number_of_sheets()[0] != 1) {
	  std::vector<Vertex_handle> vtc = hierarchy[l]->periodic_copies(w);
	  for (unsigned int i=0 ; i<vtc.size() ; i++) vtc[i]->set_down(old);
	}
    }
    if (u == Vertex_handle())
	break;
    old = w;
    v = u;
  }

  return ret;
}
 
template <class PTr>
inline
typename Periodic_3_triangulation_hierarchy_3<PTr>::Cell_handle
Periodic_3_triangulation_hierarchy_3<PTr>::
locate(const Point& p, Locate_type& lt, int& li, int& lj, Cell_handle start) const
{
  if (start != Cell_handle ()) return PTr_Base::locate (p, lt, li, lj, start);
  locs positions[maxlevel];
  locate(p, lt, li, lj, positions);
  return positions[0].pos;
}

template <class PTr>
inline
typename Periodic_3_triangulation_hierarchy_3<PTr>::Cell_handle 
Periodic_3_triangulation_hierarchy_3<PTr>::
locate(const Point& p, Cell_handle start) const
{
  if (start != Cell_handle ()) return PTr_Base::locate (p, start);
  Locate_type lt;
  int li, lj;
  return locate(p, lt, li, lj);
}

template <class PTr>
void
Periodic_3_triangulation_hierarchy_3<PTr>::
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
      hierarchy[level]->nearest_vertex_in_cell(position,p);

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

template <class PTr>
typename Periodic_3_triangulation_hierarchy_3<PTr>::Vertex_handle 
Periodic_3_triangulation_hierarchy_3<PTr>::
nearest_vertex(const Point& p, Cell_handle start) const
{
    return PTr_Base::nearest_vertex(p, start != Cell_handle() ? start
	                                                     : locate(p));
}

template <class PTr>
int
Periodic_3_triangulation_hierarchy_3<PTr>::
random_level()
{
  if ( level_mult_cover < maxlevel
      && hierarchy[level_mult_cover]->number_of_sheets() == make_array(1,1,1) )
    ++level_mult_cover;

   boost::geometric_distribution<> proba(1.0/ratio);
   boost::variate_generator<boost::rand48&, boost::geometric_distribution<> >
     die(random, proba);
   return (std::min)(die()-1, level_mult_cover);
}

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_HIERARCHY_3_H
