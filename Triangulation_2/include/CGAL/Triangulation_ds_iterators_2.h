// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Mariette Yvinec

#ifndef CGAL_TRIANGULATION_DS_ITERATORS_2_H
#define CGAL_TRIANGULATION_DS_ITERATORS_2_H

#include <iterator>
#include <CGAL/triangulation_assertions.h>
//#include <CGAL/Triangulation_iterator_adaptator.h>

namespace CGAL {

template <class Tds>
class Triangulation_ds_edge_iterator_2
{
public:
  typedef typename Tds::Edge           Edge;
  typedef typename Tds::Face_iterator  Face_iterator;
  typedef typename Tds::Face_handle    Face_handle;
    
  typedef Edge            value_type;
  typedef Edge*           pointer;
  typedef Edge&           reference;
  typedef std::size_t                            size_type;
  typedef std::ptrdiff_t                         difference_type;
  typedef std::bidirectional_iterator_tag        iterator_category;

  typedef Triangulation_ds_edge_iterator_2<Tds> Edge_iterator;
  
private:
const Tds* _tds;
Face_iterator pos;
mutable Edge edge;

public:
  Triangulation_ds_edge_iterator_2()  {}
  Triangulation_ds_edge_iterator_2(const Tds * tds);
  Triangulation_ds_edge_iterator_2(const Tds* tds, int );

  bool  operator==(const Edge_iterator& fi) const ;
  bool  operator!=(const Edge_iterator& fi) const {return !(*this== fi);}
  Edge_iterator&    operator++();
  Edge_iterator&    operator--();
  Edge_iterator    operator++(int);
  Edge_iterator    operator--(int);
  Edge*     operator->() const;
  Edge&     operator*() const ;

private: 
  void increment();
  void decrement();
  bool associated_edge();
};


// Edge iterator implementation

template<class Tds>
Triangulation_ds_edge_iterator_2<Tds> ::
Triangulation_ds_edge_iterator_2(const Tds * tds)
 :  _tds(tds) 
{
  edge.second = 0;
  if (_tds->dimension()<= 0) {
    pos = _tds->faces().end();       // there is no edge
    return;
  }
  pos = _tds->faces().begin();
  if (_tds->dimension() == 1) edge.second = 2;
    while ( pos != _tds->faces().end()  
	  && !associated_edge() ) increment();
}

template<class Tds>
Triangulation_ds_edge_iterator_2<Tds> ::
Triangulation_ds_edge_iterator_2(const Tds * tds, int )
  : _tds(tds) 
{
  pos = tds->faces().end();
  edge.second = 0;
  if (_tds->dimension() == 1) {edge.second = 2;}
}


template<class Tds>
inline
bool
Triangulation_ds_edge_iterator_2<Tds> ::
operator==(const Edge_iterator& fi) const
{
  return _tds == fi._tds  && pos == fi.pos  && edge.second == fi.edge.second;
}

template<class Tds>
inline
void
Triangulation_ds_edge_iterator_2<Tds> ::
increment()
{
  CGAL_triangulation_precondition(_tds->dimension() >= 1);
  if (_tds->dimension() == 1) ++pos;
  else  if (edge.second == 2) {edge.second =  0; ++pos;}
  else       edge.second += 1;
  return;
}

template<class Tds>
inline
void
Triangulation_ds_edge_iterator_2<Tds> ::
decrement()
{
  CGAL_triangulation_precondition(_tds->dimension() >= 1);
  if (_tds->dimension() == 1) --pos;
  else    if (edge.second == 0) { edge.second = 2; --pos;}
  else  edge.second -= 1;
  return;
}

template<class Tds>
inline
bool
Triangulation_ds_edge_iterator_2<Tds> ::
associated_edge()
{
  if (_tds->dimension() == 1) {return true;}
  return Face_handle(pos) < pos->neighbor(edge.second);
}

template<class Tds>
inline
Triangulation_ds_edge_iterator_2<Tds>&
Triangulation_ds_edge_iterator_2<Tds> ::
operator++()
{
  //CGAL_triangulation_precondition(pos != Iterator_base() && 
  //			       pos != _tds->faces().end());
  do     increment();
  while( pos != _tds->faces().end() && !associated_edge());
  return *this;
}
    

template<class Tds>
inline
Triangulation_ds_edge_iterator_2<Tds>&
Triangulation_ds_edge_iterator_2<Tds> ::
operator--()
{
  // CGAL_triangulation_precondition(pos != Iterator_base() 
  //                          && *this != Edge_iterator(_tds));
  do      decrement();
  while ( !associated_edge() && *this != Edge_iterator(_tds) ); 
  return *this;
}

    
template<class Tds>
inline
Triangulation_ds_edge_iterator_2<Tds>
Triangulation_ds_edge_iterator_2<Tds> ::    
operator++(int)
{
  Edge_iterator tmp(*this);
  ++(*this);
  return tmp;
}
    
template<class Tds>
inline
Triangulation_ds_edge_iterator_2<Tds>
Triangulation_ds_edge_iterator_2<Tds> ::     
operator--(int)
{
  Edge_iterator tmp(*this);
  --(*this);
  return tmp;
}
    
template<class Tds>
inline
typename Triangulation_ds_edge_iterator_2<Tds>::Edge*
Triangulation_ds_edge_iterator_2<Tds> ::    
operator->() const
{
  edge.first = pos;
  return &edge;
}

template<class Tds>
inline
typename Triangulation_ds_edge_iterator_2<Tds>::Edge&
Triangulation_ds_edge_iterator_2<Tds> ::    
operator*() const 
{
  edge.first = pos;
  return edge;
}
} //namespace CGAL

#endif //CGAL_TRIANGULATION_DS_ITERATORS_2_H
