// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_EDGE_LIST_H
#define CGAL_EDGE_LIST_H


#include <CGAL/Unique_hash_map.h>


CGAL_BEGIN_NAMESPACE

#define USE_BF 0

namespace CGALi {

class Edge_hash_function
  : public Handle_hash_function
{
private:
  typedef Handle_hash_function     Base;

public:
  typedef Base::result_type        result_type;

  template<class Edge>
  result_type operator()(const Edge& e) const
  {
    return Base::operator()(e.first) * (e.second + 1);
  }
};


template<class Edge>
class Edge_list_item
{
private:
  typedef typename Edge::first_type   Face_handle;
  typedef Edge_list_item<Edge>        Self;

private:
  Edge next_;
  Edge prev_;

public:
  // remove the following method and make SENTINEL_EDGE a static const
  // member of the class.
  static Edge sentinel_edge() {
    return Edge(Face_handle(), sentinel_index());
  }
private:
  static int sentinel_index() { return -1; }

private:
  void init() {
    Edge SENTINEL_EDGE = sentinel_edge();
    init( SENTINEL_EDGE, SENTINEL_EDGE );
  }

  void init(const Edge& prev, const Edge& next) {
    next_ = next;
    prev_ = prev;
  }

public:
  Edge_list_item() { init(); }
  Edge_list_item(const Edge& prev, const Edge& next)
  {
    init(prev, next);
  }

  bool is_in_list() const
  {
    return ( next_.second != sentinel_index() ||
	     prev_.second != sentinel_index() );
  }

  void set_next(const Edge& next)
  {
    next_ = next;
  }

  void set_previous(const Edge& prev)
  {
    prev_ = prev;
  }

  const Edge& next()     const { return next_; }
  const Edge& previous() const { return prev_; }

  void reset() {
    Edge SENTINEL_EDGE = sentinel_edge();
    next_ = prev_ = SENTINEL_EDGE;
  }
};

} // namespace CGALi

template<class Edge>
class Edge_list
{
private:
  typedef CGALi::Edge_list_item<Edge>     List_item;
  typedef
  Unique_hash_map<Edge,List_item,CGALi::Edge_hash_function> Edge_map;

public:
  // TYPES
  //======
  typedef std::size_t      size_type;

private:
  // PRIVATE DATA MEMBERS
  //=====================
  Edge_map             emap;
  Edge                 front_;
  size_type            size_;

private:
  // PRIVATE METHODS
  //================
  void increase_size() {
    size_++;
  }

  void decrease_size() {
    size_--;
  }

  void insert_before_nocheck(const Edge& e, const Edge& new_e) {
#if USE_BF
    Edge old_front = front_;
    front_ = e;
    push_back(new_e);
    front_ = old_front;
#else
    List_item& li_e = emap[e];

    const Edge& prev_e = li_e.previous();
    List_item& li_prev_e = emap[prev_e];

    emap[new_e] = List_item(prev_e, e);
    li_e.set_previous(new_e);
    li_prev_e.set_next(new_e);
    increase_size();
#endif
  }

public:
  // CONSTRUCTOR
  //============
  Edge_list(const Edge& e = List_item::sentinel_edge() )
    : emap(), front_(e), size_(0) {}

public:
  // PREDICATES
  //===========
  bool is_valid() const { return true; }

  bool is_in_list(const Edge& e) const {
    if ( !emap.is_defined(e) ) { return false; }
    return emap[e].is_in_list();
  }

public:
  // ACCESS METHODS
  //===============
  size_type size() const {
    return size_;
  }

  const Edge& next(const Edge& e) const {
    CGAL_precondition( is_in_list(e) );
    return emap[e].next();
  }

  const Edge& previous(const Edge& e) const {
    CGAL_precondition( is_in_list(e) );
    return emap[e].previous();
  }

  const Edge& front() const {
    CGAL_precondition( size() > 0 );
    return front_;
  }

  const Edge& back() const {
    CGAL_precondition( size() > 0 );
    return previous(front_);
  }

public:
  // INSERTION
  //==========
  void push_front(const Edge& e) {
    push_back(e);
    front_ = e;
  }

  void push_back(const Edge& e) {
    CGAL_precondition( !is_in_list(e) );

    if ( size() == 0 ) {
      emap[e] = List_item(e,e);
      front_ = e;
      increase_size();
      return;
    }

#if USE_BF
    Edge last_edge = back();
    emap[e] = List_item(last_edge, front_);
    emap[last_edge].set_next(e);
    emap[front_].set_previous(e);
    increase_size();
#else
    insert_before_nocheck(front_, e);
#endif
  }

  void insert_after(const Edge& e, const Edge& new_e) {
    CGAL_precondition( is_in_list(e) );
    CGAL_precondition( !is_in_list(new_e) );
#if USE_BF
    Edge old_front = front_;
    front_ = emap[e].next();
    push_front(new_e);
    front_ = old_front;
#else
    List_item& li_e = emap[e];

    const Edge& next_e = li_e.next();
    List_item& li_next_e = emap[next_e];

    emap[new_e] = List_item(e, next_e);
    li_e.set_next(new_e);
    li_next_e.set_previous(new_e);
    increase_size();
#endif
  }

  void insert_before(const Edge& e, const Edge& new_e) {
    CGAL_precondition( is_in_list(e) );
    CGAL_precondition( !is_in_list(new_e) );
    insert_before_nocheck(e, new_e);
  }

public:
  // REPLACEMENT
  //============
  void replace(const Edge& e, const Edge& new_e) {
#if USE_BF
    insert_before(e, new_e);
    remove(e);
#else
    CGAL_precondition( is_in_list(e) );
    CGAL_precondition( !is_in_list(new_e) );

    List_item& li_e = emap[e];

    if ( size() == 1 ) {
      emap[new_e] = List_item(new_e, new_e);
      front_ = new_e;
      li_e.reset();
    }

#if USE_BF
    Edge next_e = li_e.next();
    Edge prev_e = li_e.previous();
#else
    const Edge& next_e = li_e.next();
    const Edge& prev_e = li_e.previous();
#endif

    emap[prev_e].set_next(new_e);
    emap[next_e].set_previous(new_e);

    emap[new_e] = List_item(prev_e, next_e);

    li_e.reset();

    if ( e == front_ ) {
      front_ = new_e;
    }
#endif
  }

  // REMOVAL
  //========

  void remove(const Edge& e) {
    CGAL_precondition( is_in_list(e) );

    if ( size() == 1 ) {
      front_ = List_item::sentinel_edge();
      emap[e].reset();
      decrease_size();
      return;
    }

    List_item& li_e = emap[e];
#if USE_BF
    Edge next_e = li_e.next();
    Edge prev_e = li_e.previous();
#else
    const Edge& next_e = li_e.next();
    const Edge& prev_e = li_e.previous();
#endif

    //    Edge ne = li_e.next();
    //    Edge pe = li_e.previous();

    //    emap[pe].set_next(ne);
    //    emap[ne].set_previous(pe);

    emap[prev_e].set_next(next_e);
    emap[next_e].set_previous(prev_e);

    if ( e == front_ ) {
      //      front_ = next_e;
      front_ = next_e;
    }

    li_e.reset();

    decrease_size();
  }

  // MISCELLANEOUS
  //==============
  void clear() {
    emap.clear();
  }
};

CGAL_END_NAMESPACE


#endif // CGAL_EDGE_LIST_H
