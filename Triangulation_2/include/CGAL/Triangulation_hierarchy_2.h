// Copyright (c) 1998-2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//                 Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_TRIANGULATION_HIERARCHY_2_H
#define CGAL_TRIANGULATION_HIERARCHY_2_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/STL_Extension/internal/Has_nested_type_Bare_point.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>
#include <CGAL/assertions.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>

#include <boost/mpl/identity.hpp>
#include <boost/mpl/if.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/utility/result_of.hpp>

#include <iostream>
#include <map>
#include <vector>
#include <array>
#include <CGAL/array.h>

namespace CGAL {

// parameterization of the  hierarchy
//const float Triangulation_hierarchy_2__ratio    = 30.0;
const int Triangulation_hierarchy_2__ratio      = 30;
const int   Triangulation_hierarchy_2__minsize  = 20;
const int   Triangulation_hierarchy_2__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

template <class Tr_>
class Triangulation_hierarchy_2
  : public Tr_
{
public:
  typedef Tr_                                  Tr_Base;
  typedef typename Tr_Base::Geom_traits        Geom_traits;
  typedef typename Tr_Base::size_type          size_type;
  typedef typename Tr_Base::Vertex_handle      Vertex_handle;
  typedef typename Tr_Base::Face_handle        Face_handle;
  typedef typename Tr_Base::Vertex             Vertex;
  typedef typename Tr_Base::Locate_type        Locate_type;
  typedef typename Tr_Base::Finite_vertices_iterator  Finite_vertices_iterator;

  // this one may be weighted or not
  typedef typename Tr_Base::Point                  Point;

  // If the triangulation has defined the `Bare_point` typename, use it.
  typedef typename boost::mpl::eval_if_c<
    internal::Has_nested_type_Bare_point<Tr_Base>::value,
    typename internal::Bare_point_type<Tr_Base>,
    boost::mpl::identity<typename Tr_Base::Point>
  >::type                                          Bare_point;
  typedef typename Geom_traits::Weighted_point_2   Weighted_point;

  typedef typename Tr_Base::Weighted_tag           Weighted_tag;
  typedef typename Tr_Base::Periodic_tag           Periodic_tag;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Tr_Base::geom_traits;
  using Tr_Base::number_of_vertices;
  using Tr_Base::is_infinite;
#endif

 private:
  void init_hierarchy() {
    hierarchy[0] = this;
    for(int i=1; i<Triangulation_hierarchy_2__maxlevel; ++i)
      hierarchy[i] = &hierarchy_triangulations[i-1];
  }

 // here is the stack of triangulations which form the hierarchy
  std::array<Tr_Base,Triangulation_hierarchy_2__maxlevel-1> hierarchy_triangulations;
  std::array<Tr_Base*,Triangulation_hierarchy_2__maxlevel> hierarchy;
  boost::rand48  random;

public:
  Triangulation_hierarchy_2(const Geom_traits& traits = Geom_traits());
  Triangulation_hierarchy_2(const Triangulation_hierarchy_2& tr);

  Triangulation_hierarchy_2(Triangulation_hierarchy_2&& other)
    noexcept( noexcept(Tr_Base(std::move(other))) )
    : Tr_Base(std::move(other))
    , hierarchy_triangulations(std::move(other.hierarchy_triangulations))
    , random(std::move(other.random))
  {
    init_hierarchy();
  }

  template<class InputIterator>
  Triangulation_hierarchy_2(InputIterator first, InputIterator beyond,
                            const Geom_traits& traits = Geom_traits())
    : Tr_Base(traits)
  {
    init_hierarchy();
    insert (first, beyond);
  }

  Triangulation_hierarchy_2 &operator=(const  Triangulation_hierarchy_2& tr);

  Triangulation_hierarchy_2 & operator=(Triangulation_hierarchy_2&& other)
    noexcept( noexcept(Triangulation_hierarchy_2(std::move(other))) )
  {
    static_cast<Tr_Base&>(*this) = std::move(other);
    hierarchy_triangulations = std::move(other.hierarchy_triangulations);
    return *this;
  }

  ~Triangulation_hierarchy_2() = default;

  //Helping
  void copy_triangulation(const Triangulation_hierarchy_2 &tr);
  void swap(Triangulation_hierarchy_2 &tr);
  void clear();

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;

  // INSERT REMOVE MOVE
  Vertex_handle insert(const Point &p, Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle loc, int li );
  Vertex_handle push_back(const Point &p);

  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last)
  {
    std::ptrdiff_t n = number_of_vertices();

      std::vector<Point> points (first, last);

      // Spatial sort can only be used with Gt::Point_2: we need an adapter
      typedef typename Geom_traits::Construct_point_2 Construct_point_2;
      typedef typename boost::result_of<const Construct_point_2(const Point&)>::type Ret;
      typedef boost::function_property_map<Construct_point_2, Point, Ret> fpmap;
      typedef CGAL::Spatial_sort_traits_adapter_2<Geom_traits, fpmap> Search_traits_2;

      spatial_sort(points.begin(), points.end(),
                   Search_traits_2(
                     boost::make_function_property_map<Point, Ret, Construct_point_2>(
                       geom_traits().construct_point_2_object()), geom_traits()));

      // hints[i] is the face of the previously inserted point in level i.
      // Thanks to spatial sort, they are better hints than what the hierarchy
      // would give us.
      Face_handle hints[Triangulation_hierarchy_2__maxlevel];
      for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
              p != end; ++p)
      {
          int vertex_level = random_level();

          Vertex_handle v = hierarchy[0]->insert (*p, hints[0]);
          hints[0] = v->face();

          Vertex_handle prev = v;

          for (int level = 1; level <= vertex_level; ++level) {
              v = hierarchy[level]->insert (*p, hints[level]);
              hints[level] = v->face();

              v->set_down (prev);
              prev->set_up (v);
              prev = v;
          }
      }
      std::ptrdiff_t m = number_of_vertices();
      return m - n;
  }

  void remove_degree_3(Vertex_handle  v);
  void remove_first(Vertex_handle  v);
  void remove_second(Vertex_handle v);
  void remove(Vertex_handle  v);

  Vertex_handle move_if_no_collision(Vertex_handle v, const Point &p);
  Vertex_handle move(Vertex_handle v, const Point &p);

protected: // some internal methods

  // INSERT REMOVE DISPLACEMENT
  // GIVING NEW FACES

  template <class OutputItFaces>
  Vertex_handle insert_and_give_new_faces(const Point  &p,
                                          OutputItFaces fit,
                                          Face_handle start = Face_handle() );

  template <class OutputItFaces>
  Vertex_handle insert_and_give_new_faces(const Point& p,
                                          Locate_type lt,
                                          Face_handle loc, int li,
                                          OutputItFaces fit);

  template <class OutputItFaces>
  void remove_and_give_new_faces(Vertex_handle v,
                                 OutputItFaces fit);

  template <class OutputItFaces>
  Vertex_handle move_if_no_collision_and_give_new_faces(Vertex_handle v,
                                                        const Point &p, OutputItFaces fit);

public:

  //LOCATE
  Face_handle
  locate(const Point& p,
         Locate_type& lt,
         int& li,
         Face_handle start = Face_handle()) const;

  Face_handle
  locate(const Point&p, Face_handle start = Face_handle()) const;

  Vertex_handle
  nearest_vertex(const Point& p, Face_handle start = Face_handle()) const
  {
    return nearest_vertex_dispatch<Tr_Base>(p, start, Weighted_tag());
  }

private:

  template < typename T >
  Vertex_handle
  nearest_vertex_dispatch(const Point& p, Face_handle start, Tag_false) const
  {
    return T::nearest_vertex(p, start != Face_handle() ? start : locate(p));
  }

  // There is no nearest_vertex() for Regular.
  template < typename T >
  Vertex_handle
  nearest_vertex_dispatch(const Point&, Face_handle, Tag_true) const
  {
    CGAL_assertion(false);
    return Vertex_handle();
  }

  void  locate_in_all(const Point& p,
                      Locate_type& lt,
                      int& li,
                      Face_handle loc,
                      Face_handle
                      pos[Triangulation_hierarchy_2__maxlevel]) const;

  int random_level();

  // helping function to copy_triangulation
  // the version to be used with Tag_true is templated to avoid
  // systematic instantiation
  template <class Tag>
  void add_hidden_vertices_into_map(Tag,
                                    std::map<Vertex_handle,Vertex_handle >& V)
  {
    for (typename Tr_Base::Hidden_vertices_iterator
           it=hierarchy[0]->hidden_vertices_begin();
         it != hierarchy[0]->hidden_vertices_end(); ++it) {
      if (it->up() != Vertex_handle()) V[ it->up()->down() ] = it;
    }
  }


  void add_hidden_vertices_into_map(Tag_false ,
                                    std::map<Vertex_handle,Vertex_handle >& )
  {return;}
};



template <class Tr_>
Triangulation_hierarchy_2<Tr_>::
Triangulation_hierarchy_2(const Geom_traits& traits)
  : Tr_Base(traits)
  , hierarchy_triangulations(
           make_filled_array<Triangulation_hierarchy_2__maxlevel-1,
                             Tr_Base>(traits))
{
  init_hierarchy();
}


// copy constructor duplicates vertices and faces
template <class Tr_>
Triangulation_hierarchy_2<Tr_>::
Triangulation_hierarchy_2(const Triangulation_hierarchy_2<Tr_> &tr)
  : Triangulation_hierarchy_2(tr.geom_traits())
{
  copy_triangulation(tr);
}


//Assignment
template <class Tr_>
Triangulation_hierarchy_2<Tr_> &
Triangulation_hierarchy_2<Tr_>::
operator=(const Triangulation_hierarchy_2<Tr_> &tr)
{
  copy_triangulation(tr);
  return *this;
}


template <class Tr_>
void
Triangulation_hierarchy_2<Tr_>::
copy_triangulation(const Triangulation_hierarchy_2<Tr_> &tr)
{
  {
    for(int i=0;i<Triangulation_hierarchy_2__maxlevel;++i)
    hierarchy[i]->copy_triangulation(*tr.hierarchy[i]);
  }


  //up and down have been copied in straightforward way
  // compute a map at lower level
  std::map<Vertex_handle, Vertex_handle > V;
  {
    for( Finite_vertices_iterator it=hierarchy[0]->finite_vertices_begin();
         it != hierarchy[0]->finite_vertices_end(); ++it) {
      if (it->up() != Vertex_handle()) V[ it->up()->down() ] = it;
    }
  }

  add_hidden_vertices_into_map(Weighted_tag(), V);

  {
    for(int i=1;i<Triangulation_hierarchy_2__maxlevel;++i) {
      for( Finite_vertices_iterator it=hierarchy[i]->finite_vertices_begin();
           it != hierarchy[i]->finite_vertices_end(); ++it) {
        // down pointer goes in original instead in copied triangulation
        it->set_down(V[it->down()]);
        // make reverse link
        it->down()->set_up(it);
        // I think the next line is unnecessary (my)
        // make map for next level
        if (it->up()!=  Vertex_handle() ) V[ it->up()->down() ] = it;
      }
    }
  }
}

/* template <class Tr_> */
/* void */
/* Triangulation_hierarchy_2<Tr_>::  */
/* add_hidden_vertices_into_map(Tag_false, */
/*                              std::map<Vertex_handle,Vertex_handle >& V) { */
/*   return; */
/* } */


/* template <class Tr_> */
/* void */
/* Triangulation_hierarchy_2<Tr_>::  */
/* add_hidden_vertices_into_map(Tag_true, */
/*                              std::map<Vertex_handle,Vertex_handle >& V)  */
/* { */
/*   for (typename Tr_Base::Hidden_vertices_iterator  */
/*          it=hierarchy[0]->hidden_vertices_begin();  */
/*        it != hierarchy[0]->hidden_vertices_end(); ++it) { */
/*     if (it->up() != Vertex_handle()) V[ it->up()->down() ] = it; */
/*   } */
/* } */


template <class Tr_>
void
Triangulation_hierarchy_2<Tr_>::
swap(Triangulation_hierarchy_2<Tr_> &tr)
{
  Tr_Base::swap(tr);
  using std::swap;
  swap(hierarchy_triangulations, tr.hierarchy_triangulations);
}

template <class Tr_>
void
Triangulation_hierarchy_2<Tr_>::
clear()
{
  for(int i=0;i<Triangulation_hierarchy_2__maxlevel;++i)
    if(hierarchy[i]) hierarchy[i]->clear();
}


template <class Tr_>
bool
Triangulation_hierarchy_2<Tr_>::
is_valid(bool verbose, int level) const
{
  bool result = true;
  int i;
  Finite_vertices_iterator it;
  //verify correctness of triangulation at all levels
  for(i=0;i<Triangulation_hierarchy_2__maxlevel;++i) {
    if(verbose) // pirnt  number of vertices at each level
      std::cout << "number_of_vertices "
                <<  hierarchy[i]->number_of_vertices() << std::endl;
    result = result && hierarchy[i]->is_valid(verbose,level);
  }
    //verify that lower level has no down pointers
  for( it = hierarchy[0]->finite_vertices_begin();
       it != hierarchy[0]->finite_vertices_end(); ++it)
    result = result && ( it->down() ==   Vertex_handle());
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
      result = result && ( it->up() ==  Vertex_handle() ||
                &*it == &*(it->up())->down() );
  return result;
}


template <class Tr_>
typename Triangulation_hierarchy_2<Tr_>::Vertex_handle
Triangulation_hierarchy_2<Tr_>::
insert(const Point &p, Face_handle loc)
{
  int vertex_level = random_level();
  Locate_type lt;
  int i;
  // locate using hierarchy
  Face_handle positions[Triangulation_hierarchy_2__maxlevel];
  locate_in_all(p,lt,i,loc,positions);
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

template <class Tr_>
typename Triangulation_hierarchy_2<Tr_>::Vertex_handle
Triangulation_hierarchy_2<Tr_>::
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
    locate_in_all(p,ltt,lii,loc,positions);
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

template <class Tr_>
inline
typename Triangulation_hierarchy_2<Tr_>::Vertex_handle
Triangulation_hierarchy_2<Tr_>::
push_back(const Point &p)
{
  return insert(p);
}

template <class Tr_>
void
Triangulation_hierarchy_2<Tr_>::
remove(Vertex_handle v )
{
  Vertex_handle u=v->up();
  int l = 0 ;
  while(1){
    hierarchy[l++]->remove(v);
    if (u == Vertex_handle()) break;
    if (l >= Triangulation_hierarchy_2__maxlevel) break;
    v=u; u=v->up();
  }
}

template <class Tr_>
template <class OutputItFaces>
void
Triangulation_hierarchy_2<Tr_>::
remove_and_give_new_faces(Vertex_handle v, OutputItFaces fit)
{
  Vertex_handle u=v->up();
  int l = 0 ;
  while(1){
    if(l==0) hierarchy[l++]->remove_and_give_new_faces(v, fit);
    else hierarchy[l++]->remove(v);
    if (u == Vertex_handle()) break;
    if (l >= Triangulation_hierarchy_2__maxlevel) break;
    v=u; u=v->up();
  }
}


template <class Tr_>
inline void
Triangulation_hierarchy_2<Tr_>::
remove_degree_3(Vertex_handle v )
{
  remove(v);
}

template <class Tr_>
inline void
Triangulation_hierarchy_2<Tr_>::
remove_first(Vertex_handle v )
{
  remove(v);
}

template <class Tr_>
inline void
Triangulation_hierarchy_2<Tr_>::
remove_second(Vertex_handle v )
{
  remove(v);
}

template <class Tr_>
typename Triangulation_hierarchy_2<Tr_>::Vertex_handle
Triangulation_hierarchy_2<Tr_>::
move_if_no_collision(Vertex_handle v, const Point &p) {
  Vertex_handle u=v->up(), norm = v;
  int l = 0 ;
  while(1) {
    Vertex_handle w = hierarchy[l++]->move_if_no_collision(v, p);
    if(w != v) return w;
    if (u == Vertex_handle()) break;
    if (l >= Triangulation_hierarchy_2__maxlevel) break;
    v=u; u=v->up();
  }
  return norm;
}

template <class Tr_>
typename Triangulation_hierarchy_2<Tr_>::Vertex_handle
Triangulation_hierarchy_2<Tr_>::
move(Vertex_handle v, const Point &p) {
  CGAL_precondition(!is_infinite(v));
  Vertex_handle w = move_if_no_collision(v,p);
  if(w != v) {
    remove(v);
    return w;
  }
  return v;
}

template <class Tr_>
template <class OutputItFaces>
typename Triangulation_hierarchy_2<Tr_>::Vertex_handle
Triangulation_hierarchy_2<Tr_>::
move_if_no_collision_and_give_new_faces(Vertex_handle v, const Point &p,
                                        OutputItFaces oif)
{
  Vertex_handle u=v->up(), norm = v;
  int l = 0 ;
  while(1){
    Vertex_handle w;
    if(l == 0)
      w =
        hierarchy[l++]->move_if_no_collision_and_give_new_faces(v, p, oif);
    else w = hierarchy[l++]->move_if_no_collision(v, p);

    if(w != v) return w;

    if (u == Vertex_handle()) break;
    if (l >= Triangulation_hierarchy_2__maxlevel) break;
    v=u; u=v->up();
  }
  return norm;
}

template <class Tr_>
template < class OutputItFaces >
inline
typename Triangulation_hierarchy_2<Tr_>::Vertex_handle
Triangulation_hierarchy_2<Tr_>::insert_and_give_new_faces(const Point  &p,
                                                         OutputItFaces oif,
                                                         Face_handle loc)
{
  int vertex_level = random_level();
  Locate_type lt;
  int i;
  // locate using hierarchy
  Face_handle positions[Triangulation_hierarchy_2__maxlevel];
  locate_in_all(p,lt,i,loc,positions);
  Vertex_handle vertex=
    hierarchy[0]->Tr_Base::insert_and_give_new_faces(p,lt,positions[0],i,oif);
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

template <class Tr_>
template < class OutputItFaces >
inline
typename Triangulation_hierarchy_2<Tr_>::Vertex_handle
Triangulation_hierarchy_2<Tr_>::
insert_and_give_new_faces(const Point  &p,
                          Locate_type lt,
                          Face_handle loc,
                          int li,
                          OutputItFaces oif)
{
  int vertex_level = random_level();
  //insert at level 0
  Vertex_handle vertex=hierarchy[0]->Tr_Base::insert(p,lt,loc,li,oif);
  Vertex_handle previous=vertex;
  Vertex_handle first = vertex;

  if (vertex_level > 0) {
    // locate using hierarchy
    Locate_type ltt;
    int lii;
    Face_handle positions[Triangulation_hierarchy_2__maxlevel];
    locate_in_all(p,ltt,lii,loc,positions);
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

template <class Tr_>
typename Triangulation_hierarchy_2<Tr_>::Face_handle
Triangulation_hierarchy_2<Tr_>::
locate(const Point& p, Locate_type& lt, int& li, Face_handle loc) const
{
  Face_handle positions[Triangulation_hierarchy_2__maxlevel];
  locate_in_all(p,lt,li,loc,positions);
  return positions[0];
}

template <class Tr_>
typename Triangulation_hierarchy_2<Tr_>::Face_handle
Triangulation_hierarchy_2<Tr_>::
locate(const Point& p, Face_handle loc ) const
{
  Locate_type lt;
  int li;
  return locate(p, lt, li, loc);
}

template <class Tr_>
void
Triangulation_hierarchy_2<Tr_>::
locate_in_all(const Point& p,
              Locate_type& lt,
              int& li,
              Face_handle loc,
              Face_handle pos[Triangulation_hierarchy_2__maxlevel]) const
{
  Face_handle position;
  Vertex_handle nearest;
  int level  = Triangulation_hierarchy_2__maxlevel;
  typename Geom_traits::Compare_distance_2
    closer = geom_traits().compare_distance_2_object();

  typename Geom_traits::Construct_point_2
    construct_point = geom_traits().construct_point_2_object();

  // find the highest level with enough vertices that is at the same time 2D
  while ( (hierarchy[--level]->number_of_vertices()
           < static_cast<size_type> (Triangulation_hierarchy_2__minsize ))
          || (hierarchy[level]->dimension()<2) ){
    if ( ! level) break;  // do not go below 0
  }
  if((level>0) && (hierarchy[level]->dimension()<2)){
    level--;
  }

  for (int i=level+1; i<Triangulation_hierarchy_2__maxlevel;++i) pos[i]=nullptr;
  while(level > 0) {
    pos[level]=position=hierarchy[level]->locate(p, position);
    // locate at that level from "position"
    // result is stored in "position" for the next level
    // find the nearest between vertices 0 and 1
    if (hierarchy[level]->is_infinite(position->vertex(0))){

      nearest = position->vertex(1);
    }
    else if (hierarchy[level]->is_infinite(position->vertex(1))){
      nearest = position->vertex(0);
    }     else if ( closer(construct_point(p),
                           construct_point(position->vertex(0)->point()),
                           construct_point(position->vertex(1)->point())) == SMALLER){
      nearest = position->vertex(0);
}
    else{
      nearest = position->vertex(1);
}
    // compare to vertex 2, but only if the triangulation is 2D, because otherwise vertex(2) is  nullptr
    if ( (hierarchy[level]->dimension()==2) && (!  hierarchy[level]->is_infinite(position->vertex(2)))){
      if ( closer( construct_point(p),
                   construct_point(position->vertex(2)->point()),
                   construct_point(nearest->point())) == SMALLER ){
        nearest = position->vertex(2);
      }
    }
    // go at the same vertex on level below
    nearest  = nearest->down();
    position = nearest->face();                // incident face
    --level;
  }
  pos[0]=hierarchy[0]->locate(p,lt,li,loc == Face_handle() ? position : loc);  // at level 0
}

template <class Tr_>
int
Triangulation_hierarchy_2<Tr_>::
random_level()
{
  boost::geometric_distribution<> proba(1.0/Triangulation_hierarchy_2__ratio);
  boost::variate_generator<boost::rand48&, boost::geometric_distribution<> > die(random, proba);

  return (std::min)(die(), Triangulation_hierarchy_2__maxlevel)-1;

}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_TRIANGULATION_HIERARCHY_2_H
