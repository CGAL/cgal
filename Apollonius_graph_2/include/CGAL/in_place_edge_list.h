// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_IN_PLACE_EDGE_LIST_H
#define CGAL_IN_PLACE_EDGE_LIST_H

namespace CGAL {

template<class Edge>
class In_place_edge_list {
private:
  typedef typename Edge::first_type          Face_handle;
  typedef typename Face_handle::value_type   Face;

private:
  Edge _front;
  unsigned int _size;

private:
  inline void increase_size() {
    _size++;
  }

  inline void decrease_size() {
    _size--;
  }

public:
  bool is_valid() const { return true; }

  inline unsigned int size() const {
    return _size;
  }

  inline void push_first(const Edge& e) {
    _front = e;
    set_next(e, e);
    set_previous(e, e);
    increase_size();
  }

  inline Edge next(const Edge& e) const {
    CGAL_precondition( is_in_list(e) );
    return e.first->next(e.second);
  }

  inline Edge previous(const Edge& e) const {
    CGAL_precondition( is_in_list(e) );
    return e.first->previous(e.second);
  }

  inline void set_next(const Edge& e, const Edge& next) {
    Edge _next(next.first, next.second);
    e.first->set_next(e.second, _next);
  }

  inline void set_previous(const Edge& e, const Edge& prev) {
    Edge _prev(prev.first, prev.second);
    e.first->set_previous(e.second, _prev);
  }

  inline bool is_first(const Edge& e) const {
    return ( (e.first == _front.first &&
	      e.second == _front.second) );
  }

public:
  inline bool is_singleton() const {
    CGAL_precondition( !is_empty() );
    return (size() == 1);
  }

public:
  In_place_edge_list(const Edge& e = Edge(Face_handle(),-1) )
    : _size(0) {
    _front = e;
  }

  inline Edge front() const {
    CGAL_precondition( !is_empty() );
    return _front;
  }

  inline Edge back() const {
    CGAL_precondition( !is_empty() );
    return previous(_front);
  }

  inline bool is_empty() const {
    return ( _front.first == Face_handle() );
  }

  inline void pop() {
    CGAL_precondition( !is_empty() );
    remove(front()); // it is important here that I do not pass the
    // variable _front but rather a copy of it...
  }

  inline void push_front(const Edge& e) {
    CGAL_precondition( !is_in_list(e) );
    push(e);
    _front = e;
  }

  inline void push_back(const Edge& e) {
    push(e);
  }

  void push(const Edge& e) {
    CGAL_precondition( !is_in_list(e) );

    if ( is_empty() ) {
      push_first(e);
      return;
    }
    Edge last_edge = back();
    set_next(last_edge, e);
    set_next(e, _front);
    set_previous(e, last_edge);
    set_previous(_front, e);

    increase_size();
  }

  inline void insert_after(const Edge& e, const Edge& new_e) {
    CGAL_precondition( is_in_list(e) );
    Edge old_front = _front;
    _front = next(e);
    push_front(new_e);
    _front = old_front;
  }

  inline void insert_before(const Edge& e, const Edge& new_e) {
    CGAL_precondition( is_in_list(e) );
    Edge old_front = _front;
    _front = e;
    push(new_e);
    _front = old_front;
  }

  inline void replace(const Edge& e, const Edge& new_e) {
    insert_before(e, new_e);
    remove(e);
  }

  void remove(const Edge& e) {
    CGAL_precondition( is_in_list(e) );
    static Edge SENTINEL_QUEUE_EDGE = Edge(Face_handle(), -1);

    if ( is_singleton() ) {
      _front = SENTINEL_QUEUE_EDGE;
      set_next(e, SENTINEL_QUEUE_EDGE);
      set_previous(e, SENTINEL_QUEUE_EDGE);
      decrease_size();
      return;
    }

    Edge _next = next(e);
    Edge _prev = previous(e);

    if ( is_first(e) ) {
      _front = _next;
    }

    set_next(e, SENTINEL_QUEUE_EDGE);
    set_previous(e, SENTINEL_QUEUE_EDGE);

    set_next(_prev, _next);
    set_previous(_next, _prev);
    decrease_size();
  }

  bool is_in_list(const Edge& e) const {
    return e.first->is_in_list(e.second);
  }

  void clear() {
    while ( !is_empty() ) {
      pop();
    }
  }
};

} //namespace CGAL


#endif // CGAL_IN_EDGE_PLACE_LIST_H
