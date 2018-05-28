// Copyright (c) 1998, 2001, 2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Olivier Devillers <Olivier.Devillers@sophia.inria.fr>
//                 Sylvain Pion

#ifndef CGAL_TRIANGULATION_HIERARCHY_3_H
#define CGAL_TRIANGULATION_HIERARCHY_3_H

#include <CGAL/license/Triangulation_3.h>

// Commented because the class is actually used by Delaunay_triangulation_hierarchy_3.h
// #define CGAL_DEPRECATED_HEADER "<CGAL/Triangulation_hierarchy_3.h>"
// #include <CGAL/internal/deprecation_warning.h>

// This class is deprecated, but must be kept for backward compatibility.
//
// It would be better to move its content to the Delaunay_triangulation_3
// specializations for Fast_location and make Triangulation_hierarchy_3 the
// empty nutshell instead.
//
// Then, later, maybe merge the Compact/Fast codes in a cleaner factorized way.

#include <CGAL/basic.h>
#include <CGAL/internal/Has_nested_type_Bare_point.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_3.h>
#include <CGAL/Location_policy.h>

#include <CGAL/internal/boost/function_property_map.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/utility/result_of.hpp>

#ifndef CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/internal/info_check.h>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/if.hpp>

#endif //CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO

namespace CGAL {

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
  typedef Fast_location                        Location_policy;
  typedef typename Tr_Base::Geom_traits        Geom_traits;
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

  // this may be weighted or not
  typedef typename Tr_Base::Point              Point;

  typedef typename Tr_Base::Weighted_tag       Weighted_tag;
  typedef typename Tr_Base::Periodic_tag       Periodic_tag;

  using Tr_Base::number_of_vertices;
  using Tr_Base::geom_traits;

private:

  // here is the stack of triangulations which form the hierarchy
  Tr_Base*       hierarchy[maxlevel];
  boost::rand48  random;

  void set_up_down(Vertex_handle up, Vertex_handle down)
  {
    up->set_down(down);
    down->set_up(up);
  }

public:

  Triangulation_hierarchy_3(const Geom_traits& traits = Geom_traits());

  Triangulation_hierarchy_3(const Triangulation_hierarchy_3& tr);

  template < typename InputIterator >
  Triangulation_hierarchy_3(InputIterator first, InputIterator last,
                            const Geom_traits& traits = Geom_traits())
    : Tr_Base(traits)
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

#ifndef CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template < class InputIterator >
  std::ptrdiff_t
  insert( InputIterator first, InputIterator last,
          typename boost::enable_if<
            boost::is_convertible<
                typename std::iterator_traits<InputIterator>::value_type,
                Point
            >
          >::type* = NULL
  )
#else
  template < class InputIterator >
  std::ptrdiff_t
  insert( InputIterator first, InputIterator last)
#endif //CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  {
    size_type n = number_of_vertices();

      std::vector<Point> points (first, last);

      // Spatial sort can only be used with Geom_traits::Point_3: we need an adapter
      typedef typename Geom_traits::Construct_point_3 Construct_point_3;
      typedef typename boost::result_of<const Construct_point_3(const Point&)>::type Ret;
      typedef CGAL::internal::boost_::function_property_map<Construct_point_3, Point, Ret> fpmap;
      typedef CGAL::Spatial_sort_traits_adapter_3<Geom_traits, fpmap> Search_traits_3;

      spatial_sort(points.begin(), points.end(),
                   Search_traits_3(
                     CGAL::internal::boost_::make_function_property_map<Point, Ret, Construct_point_3>(
                       geom_traits().construct_point_3_object()), geom_traits()));

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
	      set_up_down(v, prev);
              prev = v;
          }
      }
      return number_of_vertices() - n;
  }

#ifndef CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
private:
  //top stands for tuple-or-pair
  template <class Info>
  const Point& top_get_first(const std::pair<Point,Info>& pair) const { return pair.first; }
  template <class Info>
  const Info& top_get_second(const std::pair<Point,Info>& pair) const { return pair.second; }
  template <class Info>
  const Point& top_get_first(const boost::tuple<Point,Info>& tuple) const { return boost::get<0>(tuple); }
  template <class Info>
  const Info& top_get_second(const boost::tuple<Point,Info>& tuple) const { return boost::get<1>(tuple); }

  template<class Construct_bare_point, class Container>
  struct Index_to_Bare_point
  {
    const typename Geom_traits::Point_3& operator()(const std::size_t& i) const
    {
      return cp(c[i]);
    }

    Index_to_Bare_point(const Container& c, const Construct_bare_point& cp)
      : c(c), cp(cp) { }

    const Container& c;
    const Construct_bare_point cp;
  };

  template <class Tuple_or_pair,class InputIterator>
  std::ptrdiff_t insert_with_info(InputIterator first,InputIterator last)
  {
    size_type n = number_of_vertices();
    std::vector<std::size_t> indices;
    std::vector<Point> points;
    std::vector<typename Vertex::Info> infos;
    std::size_t index=0;
    for (InputIterator it=first;it!=last;++it){
      Tuple_or_pair value=*it;
      points.push_back( top_get_first(value)  );
      infos.push_back ( top_get_second(value) );
      indices.push_back(index++);
    }

    // We need to sort the points and their info at the same time through
    // the `indices` vector AND spatial sort can only handle Geom_traits::Point_3.
    typedef typename Geom_traits::Construct_point_3 Construct_point_3;
    typedef Index_to_Bare_point<Construct_point_3,
                                std::vector<Point> > Access_bare_point;
    typedef typename boost::result_of<const Construct_point_3(const Point&)>::type Ret;
    typedef CGAL::internal::boost_::function_property_map<Access_bare_point, std::size_t, Ret> fpmap;
    typedef CGAL::Spatial_sort_traits_adapter_3<Geom_traits, fpmap> Search_traits_3;

    Access_bare_point accessor(points, geom_traits().construct_point_3_object());
    spatial_sort(indices.begin(), indices.end(),
                 Search_traits_3(
                   CGAL::internal::boost_::make_function_property_map<
                     std::size_t, Ret, Access_bare_point>(accessor),
                   geom_traits()));


    // hints[i] is the vertex of the previously inserted point in level i.
    // Thanks to spatial sort, they are better hints than what the hierarchy
    // would give us.
    Vertex_handle hints[maxlevel];
    for (typename std::vector<std::size_t>::const_iterator
      it = indices.begin(), end = indices.end();
      it != end; ++it)
    {
        int vertex_level = random_level();

        Vertex_handle v = hints[0] = hierarchy[0]->insert (points[*it], hints[0]);
        v->info()=infos[*it];
        Vertex_handle prev = v;

        for (int level = 1; level <= vertex_level; ++level) {
            v = hints[level] = hierarchy[level]->insert (points[*it], hints[level]);
            set_up_down(v, prev);
            prev = v;
        }
    }
    return number_of_vertices() - n;
  }

public:

  template < class InputIterator >
  std::ptrdiff_t
  insert( InputIterator first,
          InputIterator last,
          typename boost::enable_if<
            boost::is_convertible<
              typename std::iterator_traits<InputIterator>::value_type,
              std::pair<Point,typename internal::Info_check<Vertex>::type>
            > >::type* =NULL
  )
  {
    return insert_with_info< std::pair<Point,typename internal::Info_check<Vertex>::type> >(first,last);
  }

  template <class  InputIterator_1,class InputIterator_2>
  std::ptrdiff_t
  insert( boost::zip_iterator< boost::tuple<InputIterator_1,InputIterator_2> > first,
          boost::zip_iterator< boost::tuple<InputIterator_1,InputIterator_2> > last,
          typename boost::enable_if<
            boost::mpl::and_<
              boost::is_convertible< typename std::iterator_traits<InputIterator_1>::value_type, Point >,
              boost::is_convertible< typename std::iterator_traits<InputIterator_2>::value_type, typename internal::Info_check<Vertex>::type >
            >
          >::type* =NULL
  )
  {
    return insert_with_info< boost::tuple<Point,typename internal::Info_check<Vertex>::type> >(first,last);
  }
#endif //CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO

  void remove(Vertex_handle v);

  template < typename InputIterator >
  size_type remove(InputIterator first, InputIterator beyond)
  {
    size_type n = number_of_vertices();
    while (first != beyond) {
      remove(*first);
      ++first;
    }
    return n - number_of_vertices();
  }

  template < typename InputIterator >
  size_type remove_cluster(InputIterator first, InputIterator beyond)
  {
    CGAL_triangulation_precondition(!this->does_repeat_in_range(first, beyond));
    CGAL_triangulation_precondition(!this->infinite_vertex_in_range(first, beyond));
    size_type n = this->number_of_vertices();
    std::vector<Vertex_handle> vo(first, beyond), vc;
    int l=0;
    while(1) {
      size_type n = vo.size();
      if(n == 0) break;
      for(size_type i=0; i<n; i++) {
        if(vo[i]->up() != Vertex_handle()) vc.push_back(vo[i]->up());
      }
      hierarchy[l++]->remove_cluster(vo.begin(), vo.end());
      std::swap(vo,vc);
      vc.clear();
    }
    return n - this->number_of_vertices();
  }

  Vertex_handle move_if_no_collision(Vertex_handle v, const Point &p);
  Vertex_handle move(Vertex_handle v, const Point &p);

public: // some internal methods

  // INSERT REMOVE DISPLACEMENT
  // GIVING NEW FACES

  template <class OutputItCells>
  Vertex_handle insert_and_give_new_cells(const Point  &p, 
                                          OutputItCells fit,
                                          Cell_handle start = Cell_handle() );
		
  template <class OutputItCells>
  Vertex_handle insert_and_give_new_cells(const Point& p,
                                          OutputItCells /* fit */,
                                          Vertex_handle hint)
  {
    return insert_and_give_new_cells(p, hint == Vertex_handle() ? 
                                     this->infinite_cell() : hint->cell());			
  }

  template <class OutputItCells>
  Vertex_handle insert_and_give_new_cells(const Point& p,
                                          Locate_type lt,
                                          Cell_handle c, int li, int lj, 
                                          OutputItCells fit);

  template <class OutputItCells>
  void remove_and_give_new_cells(Vertex_handle v, 
                                 OutputItCells fit);

  template <class OutputItCells>
  Vertex_handle move_if_no_collision_and_give_new_cells(Vertex_handle v, 
                                                        const Point &p, OutputItCells fit);
	
public:	


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

protected:

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
  : Tr_Base(traits)
{
  hierarchy[0] = this;
  for(int i=1;i<maxlevel;++i)
    hierarchy[i] = new Tr_Base(traits);
}

// copy constructor duplicates vertices and cells
template <class Tr>
Triangulation_hierarchy_3<Tr>::
Triangulation_hierarchy_3(const Triangulation_hierarchy_3<Tr> &tr)
    : Tr_Base(tr)
{
  hierarchy[0] = this;
  for(int i=1; i<maxlevel; ++i)
    hierarchy[i] = new Tr_Base(*tr.hierarchy[i]);

  // up and down have been copied in straightforward way
  // compute a map at lower level

  std::map< Vertex_handle, Vertex_handle > V;

  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(),
       end = hierarchy[0]->finite_vertices_end(); it != end; ++it)
    if (it->up() != Vertex_handle())
      V[ it->up()->down() ] = it;

  for(int j=1; j<maxlevel; ++j) {
    for( Finite_vertices_iterator it = hierarchy[j]->finite_vertices_begin(),
	 end = hierarchy[j]->finite_vertices_end(); it != end; ++it) {
	// current it->down() pointer goes in original instead in copied triangulation
	set_up_down(it, V[it->down()]);
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
  for(int i=1; i<maxlevel; ++i) {
    delete hierarchy[i];
  }
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
  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(),
       end = hierarchy[0]->finite_vertices_end(); it != end; ++it)
    result = result && (it->down() == Vertex_handle());

  // verify that other levels has down pointer and reciprocal link is fine
  for(int j=1; j<maxlevel; ++j)
    for( Finite_vertices_iterator it = hierarchy[j]->finite_vertices_begin(),
	 end = hierarchy[j]->finite_vertices_end(); it != end; ++it)
      result = result && &*(it) == &*(it->down()->up());

  // verify that other levels has down pointer and reciprocal link is fine
  for(int k=0; k<maxlevel-1; ++k)
    for( Finite_vertices_iterator it = hierarchy[k]->finite_vertices_begin(),
	 end = hierarchy[k]->finite_vertices_end(); it != end; ++it)
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
    set_up_down(vertex, previous);
    previous=vertex;
    level++;
  }
  return first;
}

template <class Tr>
template <class OutputItCells>
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
insert_and_give_new_cells(const Point &p, OutputItCells fit, Cell_handle start)
{
  int vertex_level = random_level();
  Locate_type lt;
  int i, j;
  // locate using hierarchy
  locs positions[maxlevel];
  locate(p, lt, i, j, positions, start);
  // insert at level 0
  Vertex_handle vertex = hierarchy[0]->insert_and_give_new_cells(p,
                                                                 positions[0].lt,
                                                                 positions[0].pos,
                                                                 positions[0].li,
                                                                 positions[0].lj,fit);
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
    set_up_down(vertex, previous);
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
    locate(p, lt, i, j, positions, vertex->cell());

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
      set_up_down(vertex, previous);
      previous=vertex;
      level++;
    }
  }
  return first;
}

template <class Tr>
template <class OutputItCells>
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
insert_and_give_new_cells(const Point &p, Locate_type lt, Cell_handle loc, 
  int li, int lj, OutputItCells fit)
{
  int vertex_level = random_level();
  // insert at level 0
  Vertex_handle vertex = 
    hierarchy[0]->insert_and_give_new_cells(p,lt,loc,li,lj,fit);
  Vertex_handle previous = vertex;
  Vertex_handle first = vertex;

  if (vertex_level > 0) {
    Locate_type lt;
    int i, j;
    // locate using hierarchy
    locs positions[maxlevel];
    locate(p, lt, i, j, positions, vertex->cell());

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
      set_up_down(vertex, previous);
      previous=vertex;
      level++;
    }
  }
  return first;
}

template <class Tr>
void
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
}

template <class Tr>
template <class OutputItCells>
void
Triangulation_hierarchy_3<Tr>::
remove_and_give_new_cells(Vertex_handle v, OutputItCells fit)
{
  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_triangulation_precondition(!is_infinite(v));
  for (int l = 0; l < maxlevel; ++l) {
    Vertex_handle u = v->up();
    if(l) hierarchy[l]->remove(v);
    else hierarchy[l]->remove_and_give_new_cells(v, fit);
    if (u == Vertex_handle())
	break;
    v = u;
  }
}

template <class Tr>
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
move_if_no_collision(Vertex_handle v, const Point & p)
{
  CGAL_triangulation_precondition(!this->is_infinite(v));	
  if(v->point() == p) return v;
  Vertex_handle ans;
  for (int l = 0; l < maxlevel; ++l) {
    Vertex_handle u = v->up();
    if(l) hierarchy[l]->move_if_no_collision(v, p);
    else ans = hierarchy[l]->move_if_no_collision(v, p);
    if(ans != v) return ans;
    if (u == Vertex_handle())
      break;
    v = u;
  }
  return ans;
}

template <class Tr>
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
move(Vertex_handle v, const Point & p)
{
  CGAL_triangulation_precondition(!this->is_infinite(v));
  if(v->point() == p) return v;
  Vertex_handle w = move_if_no_collision(v,p);
  if(w != v) {
    remove(v);
    return w;
  }
  return v;
}

template <class Tr>
template <class OutputItCells>
typename Triangulation_hierarchy_3<Tr>::Vertex_handle
Triangulation_hierarchy_3<Tr>::
move_if_no_collision_and_give_new_cells(
  Vertex_handle v, const Point & p, OutputItCells fit)
{
  CGAL_triangulation_precondition(!is_infinite(v));	
  if(v->point() == p) return v;
  Vertex_handle ans;
  for (int l = 0; l < maxlevel; ++l) {
    Vertex_handle u = v->up();
    if(l) hierarchy[l]->move_if_no_collision(v, p);
    else ans = 
           hierarchy[l]->move_if_no_collision_and_give_new_cells(v, p, fit);
    if(ans != v) return ans;
    if (u == Vertex_handle())
      break;
    v = u;
  }
  return ans;
}

template <class Tr>
inline
typename Triangulation_hierarchy_3<Tr>::Cell_handle
Triangulation_hierarchy_3<Tr>::
locate(const Point& p, Locate_type& lt, int& li, int& lj, Cell_handle start) const
{
  if (start != Cell_handle ())
    return Tr_Base::locate (p, lt, li, lj, start);
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
  if (start != Cell_handle ())
    return Tr_Base::locate (p, start);
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
    Vertex_handle nearest = hierarchy[level]->nearest_vertex_in_cell(p, position);

    // go at the same vertex on level below
    nearest = nearest->down();
    position = nearest->cell();                // incident cell
    --level;
  }

  if (start != Cell_handle())
    position = start;

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
    return Tr_Base::nearest_vertex(p, start != Cell_handle() ? start : locate(p));
}

template <class Tr>
int
Triangulation_hierarchy_3<Tr>::
random_level()
{
  boost::geometric_distribution<> proba(1.0/ratio);
  boost::variate_generator<boost::rand48&, boost::geometric_distribution<> > die(random, proba);

  return (std::min)(die(), (int)maxlevel)-1;
}

} //namespace CGAL

#endif // CGAL_TRIANGULATION_HIERARCHY_3_H
