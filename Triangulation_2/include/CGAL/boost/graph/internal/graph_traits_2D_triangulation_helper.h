// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#include <CGAL/Iterator_range.h>
#include <CGAL/iterator.h>
#include <CGAL/use.h>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>

#include <utility>

#ifndef CGAL_GRAPH_TRAITS_2D_TRIANGULATION_HELPERS
#define CGAL_GRAPH_TRAITS_2D_TRIANGULATION_HELPERS

namespace CGAL {
namespace internal {

// A triangulation edge is a face handle + an int, and is thus actually a halfedge...
template <class Tr>
struct T2_halfedge_descriptor
  : public Tr::Edge
{
  typedef typename Tr::Edge                                   Base;
  typedef typename Tr::Face_handle                            Face_handle;

  T2_halfedge_descriptor() {}
  T2_halfedge_descriptor(Face_handle fh, int i) : Base(fh, i) { }
  explicit T2_halfedge_descriptor(const Base& e) : Base(e) { }
  T2_halfedge_descriptor(const T2_halfedge_descriptor& h) : Base(h) { }

  const Base& base() const { return static_cast<const Base&>(*this); }

  T2_halfedge_descriptor& operator=(const T2_halfedge_descriptor& h)
  {
    this->first = h.first;
    this->second = h.second;
    return *this;
  }

  friend std::size_t hash_value(const T2_halfedge_descriptor& e) {
    return hash_value(e.first);
  }
};

// An edge is just a halfedge, but we give it a complete structure to distinguish it from Tr::Edge
template <typename Tr>
struct T2_edge_descriptor
{
  typedef typename Tr::Face_handle                              Face_handle;

  T2_edge_descriptor() : first(), second(0) { }
  explicit T2_edge_descriptor(const typename Tr::Edge& e) : first(e.first), second(e.second) { }
  T2_edge_descriptor(Face_handle fd, int i) : first(fd), second(i) { }

  // so that we can still do stuff like tr.is_finite(edge_descriptor) without any hassle
  operator std::pair<Face_handle, int>() const { return std::make_pair(first, second); }

  friend std::size_t hash_value(const T2_edge_descriptor& h)
  {
    if(h.first == Face_handle())
      return 0;

    return hash_value(h.first < h.first->neighbor(h.second) ? h.first
                                                            : h.first->neighbor(h.second));
  }

  bool operator==(const T2_edge_descriptor& other) const
  {
    if((first == other.first) && (second == other.second))
      return true;

    Face_handle fh = first->neighbor(second);
    if(other.first != fh)
      return false;

    int i = fh->index(first);
    return (other.second == i);
  }
  bool operator!=(T2_edge_descriptor& other) const { return ! (*this == other); }

  void get_canonical_edge_representation(Face_handle& fh, int& i) const
  {
    Face_handle neigh_fh = fh->neighbor(i);
    Face_handle canonical_fh = (fh < neigh_fh) ? fh : neigh_fh;

    int canonical_i = (fh < neigh_fh) ? i : neigh_fh->index(fh);

    fh = canonical_fh;
    i = canonical_i;
  }

  bool operator<(const T2_edge_descriptor& other) const
  {
    if(*this == other)
      return false;

    Face_handle tfh = first;
    int ti = second;
    get_canonical_edge_representation(tfh, ti);

    Face_handle ofh = other.first;
    int oi = other.second;
    get_canonical_edge_representation(ofh, oi);

    if(tfh < ofh) return true;
    if(tfh > ofh) return false;
    return ti < oi;
  }

  Face_handle first;
  int second;
};

// A halfedge iterator is just an edge iterator that duplicates everything twice,
// to see the edge from either side.
// Could probably be factorized with T2_edge_iterator, but it's clearer this way.
template <typename Tr, typename EdgeIterator>
struct T2_halfedge_iterator
{
private:
  typedef T2_halfedge_iterator<Tr, EdgeIterator>                Self;
  typedef EdgeIterator                                          Edge_iterator;
  typedef T2_halfedge_descriptor<Tr>                            Descriptor;
  typedef typename Tr::Face_handle                              Face_handle;

public:
  typedef Descriptor                                            value_type;
  typedef value_type*                                           pointer;
  typedef value_type&                                           reference;
  typedef std::size_t                                           size_type;
  typedef std::ptrdiff_t                                        difference_type;
  typedef std::bidirectional_iterator_tag                       iterator_category;

  T2_halfedge_iterator() { }
  T2_halfedge_iterator(const Edge_iterator& feit) : it(feit), on_adjacent_face(false) { }

  Self& operator++()
  {
    // If we are on the first face, move to the opposite face. If we are already on the opposite face,
    // then it's time to move on the next edge
    if(on_adjacent_face) {
      ++it;
      on_adjacent_face = false;
    } else {
      on_adjacent_face = true;
    }

    return *this;
  }

  Self& operator--()
  {
    // Note that while decreasing, we start from the opposite face
    if(on_adjacent_face) {
      on_adjacent_face = false;
    } else {
      --it;
      on_adjacent_face = true;
    }

    return *this;
  }

  Self operator++(int) { Self tmp = *this; operator++(); return tmp; }
  Self operator--(int) { Self tmp = *this; operator--(); return tmp; }

  bool operator==(const Self& other) const { return it == other.it; }
  bool operator!=(const Self& other) const { return !(*this == other); }

  reference operator*() const
  {
    if(on_adjacent_face)
    {
      Face_handle neigh_f = it->first->neighbor(it->second);
      hd = Descriptor(neigh_f, neigh_f->index(it->first));
      return hd;
    } else {
      hd = Descriptor(it->first, it->second);
      return hd;
    }
  }

private:
  Edge_iterator it;
  bool on_adjacent_face;
  mutable Descriptor hd;
};

template <typename Tr, typename EdgeIterator>
struct T2_edge_iterator
{
private:
  typedef T2_edge_iterator<Tr, EdgeIterator>                    Self;
  typedef EdgeIterator                                          Edge_iterator;
  typedef T2_edge_descriptor<Tr>                                Descriptor;

public:
  typedef Descriptor                                            value_type;
  typedef value_type*                                           pointer;
  typedef value_type&                                           reference;
  typedef std::size_t                                           size_type;
  typedef std::ptrdiff_t                                        difference_type;
  typedef std::bidirectional_iterator_tag                       iterator_category;

  T2_edge_iterator() { }
  T2_edge_iterator(const Edge_iterator& feit) : it(feit) { }

  bool operator==(const Self& other) const { return it == other.it; }
  bool operator!=(const Self& other) const { return !(*this == other); }
  Self& operator++() { ++it; return *this; }
  Self& operator--() { --it; return *this; }
  Self operator++(int) { Self tmp = *this; operator++(); return tmp; }
  Self operator--(int) { Self tmp = *this; operator--(); return tmp; }

  reference operator*() const
  {
    ed = Descriptor(*it);
    return ed;
  }

private:
  Edge_iterator it;
  mutable Descriptor ed;
};

// Must distinguish TDS and triangulations circulators (later are filtered)
template <class Circ, class E>
class TDS2_Out_edge_circulator
  : public Circ
{
private:
  mutable E e;

public:
  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  TDS2_Out_edge_circulator() : Circ() {}
  TDS2_Out_edge_circulator(Circ c) : Circ(c) {}

  const E& operator*() const
  {
    E ed = static_cast<const Circ*>(this)->operator*();
    e = E(ed.first->neighbor(ed.second), ed.first->neighbor(ed.second)->index(ed.first));
    return e;
  }
};

template <class Circ, class E>
class TDS2_In_edge_circulator
  : public Circ
{
private:
  mutable E e;

public:
  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  TDS2_In_edge_circulator() : Circ() {}
  TDS2_In_edge_circulator(Circ c) : Circ(c) {}

  const E& operator*() const
  {
    typename Circ::value_type ed = static_cast<const Circ*>(this)->operator*();
    e = E(ed);
    return e;
  }
};

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
      boost::use_default,
      Handle /*reference*/
    >
{
public:
  typedef Handle                                                                       value_type;

private:
  typedef Dereference_to_handle_enforcer<Tr, Iterator, Handle>                         Self;
  typedef Iterator                                                                     I;
  typedef boost::iterator_adaptor<Self, I, value_type, boost::use_default, value_type> Base;

public:
  Dereference_to_handle_enforcer() { }
  explicit Dereference_to_handle_enforcer(const I& i) : Base(i) { }

private:
  friend class boost::iterator_core_access;
  value_type dereference() const { return value_type(this->base()); }
};

} // namespace internal
} // namespace CGAL

namespace std {

// workaround a bug detected on at least g++ 4.4 where boost::next(Iterator)
// is picked as a candidate for next(h,g)
template <typename Tr>
struct iterator_traits< CGAL::internal::T2_halfedge_descriptor<Tr> >
{
  typedef void* iterator_category;
  typedef void* difference_type;
  typedef void* value_type;
  typedef void* reference;
};

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash
#endif

#ifndef CGAL_CFG_NO_STD_HASH

template < class Tr>
struct hash<CGAL::internal::T2_halfedge_descriptor<Tr> >
{
  std::size_t operator()(const CGAL::internal::T2_halfedge_descriptor<Tr>& e) const {
    return hash_value(e);
  }
};

#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std

#endif // CGAL_GRAPH_TRAITS_2D_TRIANGULATION_HELPERS
