// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#include <CGAL/boost/graph/internal/graph_traits_2D_TDS_helper.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/iterator.h>
#include <CGAL/use.h>

#include <boost/graph/graph_traits.hpp>

#include <utility>

#ifndef CGAL_GRAPH_TRAITS_2D_TRIANGULATION_HELPERS
#define CGAL_GRAPH_TRAITS_2D_TRIANGULATION_HELPERS

namespace CGAL {
namespace internal {

// just for clarity
template <typename T>
using T2_halfedge_descriptor = TDS2_halfedge_descriptor<T>;

template <typename T, typename EI>
using T2_halfedge_iterator = TDS2_halfedge_iterator<T, EI>;

template <typename T>
using T2_edge_descriptor = TDS2_edge_descriptor<T>;

template <typename T, typename EI>
using T2_edge_iterator = TDS2_edge_iterator<T, EI>;


template <typename Tr>
struct T2_edge_circulator
  : public Tr::Edge_circulator
{
private:
  typedef T2_edge_circulator<Tr>                                Self;
  typedef typename Tr::Edge                                     Edge;
  typedef typename Tr::Edge_circulator                          Base;

public:
  typedef T2_edge_descriptor<Tr>                                value_type;
  typedef value_type*                                           pointer;
  typedef value_type&                                           reference;

  T2_edge_circulator() : Base() { }
  T2_edge_circulator(const Base& c, const Tr& tr) : Base(c), tr(&tr), e() { }

  // Note that the inf check is on the edge in the circulator, not on 'e', which isn't built yet
  Self& operator++() {
    do { this->Base::operator++(); } while(tr->is_infinite(this->Base::operator*()));
    return *this;
  }
  Self& operator--() {
    do { this->Base::operator--(); } while(tr->is_infinite(this->Base::operator*()));
    return *this;
  }
  Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
  Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

protected:
  const Tr* tr;
  mutable value_type e;
};

template <typename Tr>
struct In_edge_circulator
  : public T2_edge_circulator<Tr>
{
private:
  typedef T2_edge_circulator<Tr>                                Base;
  typedef typename Tr::Edge                                     Edge;
  typedef typename Tr::Edge_circulator                          Edge_circulator;

public:
  typedef T2_edge_descriptor<Tr>                                value_type;
  typedef value_type*                                           pointer;
  typedef value_type&                                           reference;

  In_edge_circulator() : Base() { }
  In_edge_circulator(const Edge_circulator& c, const Tr& tr) : Base(c, tr) { }

  const value_type& operator*() const
  {
    this->e = value_type(this->Base::operator*());
    return this->e;
  }
};

template <typename Tr>
struct Out_edge_circulator
  : public T2_edge_circulator<Tr>
{
private:
  typedef T2_edge_circulator<Tr>                                Base;
  typedef typename Tr::Edge                                     Edge;
  typedef typename Tr::Edge_circulator                          Edge_circulator;

public:
  typedef T2_edge_descriptor<Tr>                                value_type;
  typedef value_type*                                           pointer;
  typedef value_type&                                           reference;

  Out_edge_circulator() : Base() { }
  Out_edge_circulator(const Edge_circulator& c, const Tr& tr) : Base(c, tr) { }

  const value_type& operator*() const
  {
    Edge ed(this->Base::operator*());
    this->e = value_type(ed.first->neighbor(ed.second),
                         ed.first->neighbor(ed.second)->index(ed.first));
    return this->e;
  }
};

template <typename Tr>
struct T2_vertex_circulator
  : public In_edge_circulator<Tr>
{
private:
  typedef In_edge_circulator<Tr>                                Base;
  typedef T2_edge_descriptor<Tr>                                edge_descriptor;
  typedef typename Tr::Edge_circulator                          Edge_circulator;
  typedef typename Tr::Vertex_handle                            Vertex_handle;

public:
  typedef Vertex_handle                                         value_type;
  typedef value_type&                                           reference;

  T2_vertex_circulator() : Base() { }
  T2_vertex_circulator(const Edge_circulator& c, const Tr& tr) : Base(c, tr) { }

  const value_type& operator*() const
  {
    const edge_descriptor& edge = this->Base::operator*();
    v = edge.first->vertex(this->tr->ccw(edge.second));
    return v;
  }

private:
  // Because we wrap the iterator with a Counting_iterator, which returns a ref in its operator*()
  mutable Vertex_handle v;
};

template <typename Tr, typename Iterator, typename Handle>
struct Dereference_to_handle_enforcer
  : public boost::iterator_adaptor<
      Dereference_to_handle_enforcer<Tr, Iterator, Handle>,
      Iterator /*base*/,
      Handle /*value*/,
      typename std::iterator_traits<Iterator>::iterator_category,
      Handle /*reference*/
    >
{
public:
  typedef Handle                                                                       value_type;

private:
  typedef Dereference_to_handle_enforcer<Tr, Iterator, Handle>                         Self;
  typedef Iterator                                                                     I;
  typedef typename std::iterator_traits<I>::iterator_category                          Category;
  typedef boost::iterator_adaptor<Self, I, value_type, Category, value_type>           Base;

public:
  Dereference_to_handle_enforcer() { }
  explicit Dereference_to_handle_enforcer(const I& i) : Base(i) { }

private:
  friend class boost::iterator_core_access;
  value_type dereference() const { return value_type(this->base()); }
};

} // namespace internal
} // namespace CGAL

#endif // CGAL_GRAPH_TRAITS_2D_TRIANGULATION_HELPERS
