// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//                 Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//                 Nico Kruithof  <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_HIERARCHY_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_HIERARCHY_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>


#include <CGAL/basic.h>
#include <CGAL/array.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>
#include <map>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace CGAL
{

template < class PTr>
class Periodic_2_triangulation_hierarchy_2
  : public PTr
{
  // parameterization of the  hierarchy
  static const int m_ratio    = 30;
  static const int m_minsize  = 20;
  static const int m_maxlevel = 5;
  // maximal number of points is 30^5 = 24 millions !


public:
  typedef PTr                              PTr_Base;
  typedef typename PTr::Geom_traits        Geom_traits;
  typedef typename PTr::Point              Point;
  typedef typename PTr::Iso_rectangle      Iso_rectangle;
  typedef typename PTr::size_type          size_type;
  typedef typename PTr::Vertex_handle      Vertex_handle;
  typedef typename PTr::Face_handle        Face_handle;
  typedef typename PTr::Vertex             Vertex;
  typedef typename PTr::Locate_type        Locate_type;
  typedef typename PTr::Finite_vertices_iterator  Finite_vertices_iterator;
  //typedef typename PTr::Finite_faces_iterator     Finite_faces_iterator;

  typedef typename PTr::Weighted_tag       Weighted_tag;
  typedef typename PTr::Periodic_tag       Periodic_tag;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using PTr_Base::geom_traits;
#endif

private:
  // here is the stack of triangulations which form the hierarchy
  PTr_Base*   hierarchy[m_maxlevel];
  boost::rand48  random;
  int level_mult_cover;

public:
  Periodic_2_triangulation_hierarchy_2(
    const Iso_rectangle& domain = Iso_rectangle(0, 0, 1, 1),
    const Geom_traits& traits = Geom_traits());

  Periodic_2_triangulation_hierarchy_2(
    const Periodic_2_triangulation_hierarchy_2& tr);

  template < typename InputIterator >
  Periodic_2_triangulation_hierarchy_2(InputIterator first, InputIterator last,
                                       const Iso_rectangle& domain = Iso_rectangle(0, 0, 1, 1),
                                       const Geom_traits& traits = Geom_traits())
    : PTr_Base(domain, traits), level_mult_cover(0)
  {
    hierarchy[0] = this;
    for(int i = 1; i < m_maxlevel; ++i)
      hierarchy[i] = new PTr_Base(domain, traits);
    insert(first, last);
  }

  Periodic_2_triangulation_hierarchy_2 &operator=(const Periodic_2_triangulation_hierarchy_2& tr);
  ~Periodic_2_triangulation_hierarchy_2();

  //Helping
  void copy_triangulation(const Periodic_2_triangulation_hierarchy_2 &tr);
  void swap(Periodic_2_triangulation_hierarchy_2 &tr);
  void clear();

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;

  // INSERT REMOVE MOVE
  Vertex_handle insert(const Point &p,
                       Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle loc, int li );
  Vertex_handle push_back(const Point &p);

  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last, bool /* is_large_point_set */ = false)
  {
    std::ptrdiff_t n = this->number_of_vertices();

    std::vector<Point> points (first, last);
    CGAL::spatial_sort (points.begin(), points.end(), geom_traits());

    // hints[i] is the face of the previously inserted point in level i.
    // Thanks to spatial sort, they are better hints than what the hierarchy
    // would give us.
    Face_handle hints[m_maxlevel];
    for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
         p != end; ++p)
      {
        int vertex_level = random_level();

        Vertex_handle v = hierarchy[0]->insert (*p, hints[0]);
        hints[0] = v->face();

        Vertex_handle prev = v;

        for (int level = 1; level <= vertex_level; ++level)
          {
            v = hierarchy[level]->insert (*p, hints[level]);
            hints[level] = v->face();

            v->set_down (prev);
            if (hierarchy[level]->number_of_sheets()[0] != 1)
              {
                std::vector<Vertex_handle> vtc = hierarchy[level]->periodic_copies(v);
                for (unsigned int i = 0 ; i < vtc.size() ; i++) vtc[i]->set_down(prev);
              }

            prev->set_up (v);
            prev = v;
          }
      }
    std::ptrdiff_t m = this->number_of_vertices();
    return m - n;
  }

  void remove(Vertex_handle  v);
  void remove_first(Vertex_handle  v);
  void remove_degree_3(Vertex_handle  v);

  Vertex_handle move_point(Vertex_handle v, const Point &p);
  Vertex_handle move_if_no_collision(Vertex_handle v, const Point &p);

  //LOCATE
  Face_handle
  locate(const Point& p,
         Locate_type& lt,
         int& li,
         Face_handle start = Face_handle()) const;

  Face_handle
  locate(const Point &p,
         Face_handle start = Face_handle()) const;

  Vertex_handle
  nearest_vertex(const Point& p, Face_handle start = Face_handle()) const
  {
    return PTr_Base::nearest_vertex(p, start != Face_handle() ? start : locate(p));
  }

private:
  void  locate_in_all(const Point& p,
                      Locate_type& lt,
                      int& li,
                      Face_handle loc,
                      Face_handle pos[m_maxlevel]) const;
  int random_level();
};



template <class PTr >
Periodic_2_triangulation_hierarchy_2<PTr>::
Periodic_2_triangulation_hierarchy_2(const Iso_rectangle& domain, const Geom_traits& traits)
  : PTr_Base(domain, traits)
{
  level_mult_cover = 0;
  hierarchy[0] = this;
  for(int i = 1; i < m_maxlevel; ++i)
    hierarchy[i] = new PTr_Base(domain, traits);
}


// copy constructor duplicates vertices and faces
template <class PTr>
Periodic_2_triangulation_hierarchy_2<PTr>::
Periodic_2_triangulation_hierarchy_2(const Periodic_2_triangulation_hierarchy_2<PTr> &tr)
  : PTr_Base()
{
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this;
  for(int i = 1; i < m_maxlevel; ++i)
    hierarchy[i] = new PTr_Base(tr.domain(), tr.geom_traits());
  copy_triangulation(tr);
}


//Assignement
template <class PTr>
Periodic_2_triangulation_hierarchy_2<PTr> &
Periodic_2_triangulation_hierarchy_2<PTr>::
operator=(const Periodic_2_triangulation_hierarchy_2<PTr> &tr)
{
  Periodic_2_triangulation_hierarchy_2 tmp(tr);
  swap(tmp);
  return *this;
}


template <class PTr>
void
Periodic_2_triangulation_hierarchy_2<PTr>::
copy_triangulation(const Periodic_2_triangulation_hierarchy_2<PTr> &tr)
{
  {
    for(int i = 0; i < m_maxlevel; ++i)
      hierarchy[i]->copy_triangulation(*tr.hierarchy[i]);
  }


  //up and down have been copied in straightforward way
  // compute a map at lower level
  std::map<Vertex_handle, Vertex_handle > V;
  {
    for(Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin();
        it != hierarchy[0]->finite_vertices_end(); ++it)
      {
        if (it->up() != Vertex_handle()) V[ it->up()->down() ] = it;
      }
  }

  {
    for(int i = 1; i < m_maxlevel; ++i)
      {
        for( Finite_vertices_iterator it = hierarchy[i]->finite_vertices_begin();
             it != hierarchy[i]->finite_vertices_end(); ++it)
          {
            if (hierarchy[i]->is_virtual(it))
              {
                // down pointer goes in original instead in copied triangulation
                it->set_down(V[it->down()]);
                // make reverse link
                it->down()->set_up(it);
                // I think the next line is unnecessary (my)
                // make map for next level
                if (it->up() !=  Vertex_handle() ) V[ it->up()->down() ] = it;
              }
          }
      }
  }
}

template <class PTr>
void
Periodic_2_triangulation_hierarchy_2<PTr>::
swap(Periodic_2_triangulation_hierarchy_2<PTr> &tr)
{
  PTr_Base::swap(tr);
  for(int i = 1; i < m_maxlevel; ++i)
    std::swap(hierarchy[i], tr.hierarchy[i]);
}

template <class PTr>
Periodic_2_triangulation_hierarchy_2<PTr>::
~Periodic_2_triangulation_hierarchy_2()
{
  clear();
  for(int i = 1; i < m_maxlevel; ++i)
      delete hierarchy[i];
}

template <class PTr>
void
Periodic_2_triangulation_hierarchy_2<PTr>::
clear()
{
  for(int i = 0; i < m_maxlevel; ++i) 
  {
      CGAL_assertion(hierarchy[i] != NULL);
      hierarchy[i]->clear();
  }
}


template <class PTr>
bool
Periodic_2_triangulation_hierarchy_2<PTr>::
is_valid(bool verbose, int level) const
{
  bool result = true;
  int i;
  Finite_vertices_iterator it;
  //verify correctness of triangulation at all levels
  for(i = 0; i < m_maxlevel; ++i)
    {
      if(verbose) // print  number of vertices at each level
        std::cout << "number_of_vertices "
                  <<  hierarchy[i]->number_of_vertices() << std::endl;
      result = result && hierarchy[i]->is_valid(verbose, level);
    }
  //verify that lower level has no down pointers
  for( it = hierarchy[0]->finite_vertices_begin();
       it != hierarchy[0]->finite_vertices_end(); ++it)
    if (!hierarchy[0]->is_virtual(it))
      result = result && (it->down() == Vertex_handle());

  //verify that other levels have down pointer and reciprocal link is fine
  for(i = 1; i < m_maxlevel; ++i)
    for( it = hierarchy[i]->finite_vertices_begin();
         it != hierarchy[i]->finite_vertices_end(); ++it)
      if (!hierarchy[i]->is_virtual(it))
        result = result && (&*(it->down()->up()) == &*(it));

  //verify that levels have up pointer and reciprocal link is fine
  for(i = 0; i < m_maxlevel - 1; ++i)
    for( it = hierarchy[i]->finite_vertices_begin();
         it != hierarchy[i]->finite_vertices_end(); ++it)
      if (!hierarchy[i]->is_virtual(it))
        result = result && ( it->up() == Vertex_handle() || &*it == &*(it->up())->down() );

  return result;
}


template <class PTr>
typename Periodic_2_triangulation_hierarchy_2<PTr>::Vertex_handle
Periodic_2_triangulation_hierarchy_2<PTr>::
insert(const Point &p, Face_handle loc)
{
  int vertex_level = random_level();
  Locate_type lt;
  int i;
  // locate using hierarchy
  Face_handle positions[m_maxlevel];
  locate_in_all(p, lt, i, loc, positions);
  Vertex_handle vertex = hierarchy[0]->PTr_Base::insert(p, lt, positions[0], i);
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;

  int level  = 1;
  while (level <= vertex_level )
    {
      vertex = hierarchy[level]->PTr_Base::insert(p, positions[level]);
      vertex->set_down(previous);// link with level above
      if (hierarchy[level]->number_of_sheets()[0] != 1)
        {
          std::vector<Vertex_handle> vtc = hierarchy[level]->periodic_copies(vertex);
          for (unsigned int i = 0 ; i < vtc.size() ; i++) vtc[i]->set_down(previous);
        }
      previous->set_up(vertex);
      previous = vertex;
      level++;
    }
  return first;
}

template <class PTr>
typename Periodic_2_triangulation_hierarchy_2<PTr>::Vertex_handle
Periodic_2_triangulation_hierarchy_2<PTr>::
insert(const Point& p,
       Locate_type lt,
       Face_handle loc,
       int li )
{
  int vertex_level = random_level();
  //insert at level 0
  Vertex_handle vertex = hierarchy[0]->PTr_Base::insert(p, lt, loc, li);
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;

  if (vertex_level > 0)
    {
      // locate using hierarchy
      Locate_type ltt;
      int lii;
      Face_handle positions[m_maxlevel];
      locate_in_all(p, ltt, lii, loc, positions);
      //insert in higher levels
      int level  = 1;
      while (level <= vertex_level )
        {
          vertex = hierarchy[level]->PTr_Base::insert(p, positions[level]);
          vertex->set_down(previous);// link with level above
          if (hierarchy[level]->number_of_sheets()[0] != 1)
            {
              std::vector<Vertex_handle> vtc = hierarchy[level]->periodic_copies(vertex);
              for (unsigned int i = 0 ; i < vtc.size() ; i++) vtc[i]->set_down(previous);
            }
          previous->set_up(vertex);
          previous = vertex;
          level++;
        }
    }
  return first;
}

template <class PTr>
inline
typename Periodic_2_triangulation_hierarchy_2<PTr>::Vertex_handle
Periodic_2_triangulation_hierarchy_2<PTr>::
push_back(const Point &p)
{
  return insert(p);
}

template <class PTr>
void
Periodic_2_triangulation_hierarchy_2<PTr>::
remove(Vertex_handle v )
{
  Vertex_handle u = v->up();
  int l = 0 ;
  while(1)
    {
      hierarchy[l++]->remove(v);
      if (u == Vertex_handle()) break;
      if (l >= m_maxlevel) break;
      v = u;
      u = v->up();
    }
}

template <class PTr>
inline void
Periodic_2_triangulation_hierarchy_2<PTr>::
remove_degree_3(Vertex_handle v )
{
  remove(v);
}

template <class PTr>
inline void
Periodic_2_triangulation_hierarchy_2<PTr>::
remove_first(Vertex_handle v )
{
  remove(v);
}

template <class PTr>
typename Periodic_2_triangulation_hierarchy_2<PTr>::Vertex_handle
Periodic_2_triangulation_hierarchy_2<PTr>::
move_if_no_collision(Vertex_handle v, const Point &p)
{
  CGAL_triangulation_precondition(v != Vertex_handle());
  Vertex_handle old, ret;

  for (int l = 0; l < m_maxlevel; ++l)
    {
      Vertex_handle u = v->up();
      CGAL_triangulation_assertion(hierarchy[l]->is_valid());
      Vertex_handle w = hierarchy[l]->move_if_no_collision(v, p);
      if (l == 0)
        {
          ret = w;
        }
      else
        {
          old->set_up(w);
          w->set_down(old);
          if (hierarchy[l]->number_of_sheets()[0] != 1)
            {
              std::vector<Vertex_handle> vtc = hierarchy[l]->periodic_copies(w);
              for (unsigned int i = 0 ; i < vtc.size() ; i++) vtc[i]->set_down(old);
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
typename Periodic_2_triangulation_hierarchy_2<PTr>::Vertex_handle
Periodic_2_triangulation_hierarchy_2<PTr>::
move_point(Vertex_handle v, const Point &p)
{
  CGAL_triangulation_precondition(v != Vertex_handle());
  Vertex_handle old, ret;

  for (int l = 0; l < m_maxlevel; ++l)
    {
      Vertex_handle u = v->up();
      CGAL_triangulation_assertion(hierarchy[l]->is_valid());
      Vertex_handle w = hierarchy[l]->move_point(v, p);
      if (l == 0)
        {
          ret = w;
        }
      else
        {
          old->set_up(w);
          w->set_down(old);
          if (hierarchy[l]->number_of_sheets()[0] != 1)
            {
              std::vector<Vertex_handle> vtc = hierarchy[l]->periodic_copies(w);
              for (unsigned int i = 0 ; i < vtc.size() ; i++) vtc[i]->set_down(old);
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
typename Periodic_2_triangulation_hierarchy_2<PTr>::Face_handle
Periodic_2_triangulation_hierarchy_2<PTr>::
locate(const Point& p, Locate_type& lt, int& li, Face_handle loc) const
{
  Face_handle positions[m_maxlevel];
  locate_in_all(p, lt, li, loc, positions);
  return positions[0];
}

template <class PTr>
typename Periodic_2_triangulation_hierarchy_2<PTr>::Face_handle
Periodic_2_triangulation_hierarchy_2<PTr>::
locate(const Point& p, Face_handle loc ) const
{
  if (loc != Face_handle ()) return PTr_Base::locate (p, loc);
  Locate_type lt;
  int li;
  return locate(p, lt, li, loc);
}

template <class PTr>
void
Periodic_2_triangulation_hierarchy_2<PTr>::
locate_in_all(const Point& p,
              Locate_type& lt,
              int& li,
              Face_handle loc,
              Face_handle pos[m_maxlevel]) const
{
  Face_handle position;
  Vertex_handle nearest;
  int level  = m_maxlevel;
  typename Geom_traits::Compare_distance_2
  closer = this->geom_traits().compare_distance_2_object();

  // find the highest level with enough vertices that is at the same time 2D
  while ( (hierarchy[--level]->number_of_vertices()
           < static_cast<size_type> (m_minsize ))
          || (hierarchy[level]->dimension() < 2) )
    {
      if ( ! level) break;  // do not go below 0
    }
  if((level > 0) && (hierarchy[level]->dimension() < 2))
    {
      level--;
    }

  for (int i = level + 1; i < m_maxlevel; ++i) pos[i] = 0;
  while(level > 0)
    {
      pos[level] = position = hierarchy[level]->locate(p, position);
      // locate at that level from "position"
      // result is stored in "position" for the next level
      // find the nearest between vertices 0 and 1
      if (hierarchy[level]->is_infinite(position->vertex(0)))
        {

          nearest = position->vertex(1);
        }
      else if (hierarchy[level]->is_infinite(position->vertex(1)))
        {
          nearest = position->vertex(0);
        }
      else if ( closer(p,
                       position->vertex(0)->point(),
                       position->vertex(1)->point()) == SMALLER)
        {
          nearest = position->vertex(0);
        }
      else
        {
          nearest = position->vertex(1);
        }
      // compare to vertex 2, but only if the triangulation is 2D, because otherwise vertex(2) is  NULL
      if ( (hierarchy[level]->dimension() == 2) && (!  hierarchy[level]->is_infinite(position->vertex(2))))
        {
          if ( closer( p,
                       position->vertex(2)->point(),
                       nearest->point()) == SMALLER )
            {
              nearest = position->vertex(2);
            }
        }
      // go at the same vertex on level below
      nearest  = nearest->down();
      position = nearest->face();                // incident face
      --level;
    }
  pos[0] = hierarchy[0]->locate(p, lt, li, loc == Face_handle() ? position : loc); // at level 0
}

template <class PTr>
int
Periodic_2_triangulation_hierarchy_2<PTr>::
random_level()
{
  if ( level_mult_cover < m_maxlevel
       && hierarchy[level_mult_cover]->number_of_sheets() == make_array(1, 1) )
    ++level_mult_cover;

  boost::geometric_distribution<> proba(1.0 / m_ratio);
  boost::variate_generator<boost::rand48&, boost::geometric_distribution<> >
  die(random, proba);
  return (std::min)(die() - 1, level_mult_cover);
}

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_HIERARCHY_2_H
