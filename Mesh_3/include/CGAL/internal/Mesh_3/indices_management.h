// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef CGAL_INTERNAL_MESH_3_INDICES_MANAGEMENT_H
#define CGAL_INTERNAL_MESH_3_INDICES_MANAGEMENT_H

#include <CGAL/license/Triangulation_3.h>


#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <CGAL/Mesh_3/Has_features.h>
#include <CGAL/IO/io.h>

#include <tuple>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

// -----------------------------------
// Index_generator
// Don't use boost::variant if types are the same type
// -----------------------------------
template < typename Subdomain_index, typename Surface_patch_index >
struct Index_generator
{
  typedef boost::variant<Subdomain_index,Surface_patch_index> Index;
  typedef Index                                               type;
};

template < typename T >
struct Index_generator<T, T>
{
  typedef T       Index;
  typedef Index   type;
};

template <typename MD, bool has_feature = Has_features<MD>::value>
struct Indices_tuple_generator
{
  using type = std::tuple<typename MD::Subdomain_index,
                          typename MD::Surface_patch_index,
                          typename MD::Curve_index,
                          typename MD::Corner_index
                          >;
};

template <typename MD>
struct Indices_tuple_generator<MD, false>
{
  using type = std::tuple<typename MD::Subdomain_index,
                          typename MD::Surface_patch_index>;
};

template <typename MD>
using Indices_tuple_t = typename Indices_tuple_generator<MD>::type;

// Nasty meta-programming to get a boost::variant of four types that
// may not be all different.
template <typename T0> struct seq1 {
  typedef T0 type;
};
template <typename T0, typename T1> struct seq2 {
  typedef boost::variant<T0, T1> type;
};
template <typename T0, typename T1, typename T2> struct seq3 {
  typedef boost::variant<T0, T1, T2> type;
};
template <typename T0, typename T1, typename T2, typename T3> struct seq4 {
  typedef boost::variant<T0, T1, T2, T3> type;
};

template <typename T, typename U> struct insert;
template <typename T, typename U> struct insert<seq1<T>, U> {
  typedef seq2<T, U> type;
};
template <typename T, typename U, typename V> struct insert<seq2<T, U>, V> {
  typedef seq3<T, U, V> type;
};
template <typename T, typename U, typename V, typename W>
struct insert<seq3<T, U, V>, W> {
  typedef seq4<T, U, V, W> type;
};
template <typename T> struct insert<seq1<T>, T> {
  typedef seq1<T> type;
};
template <typename T, typename U> struct insert<seq2<T, U>, T> {
  typedef seq2<T, U> type;
};
template <typename T, typename U> struct insert<seq2<T, U>, U> {
  typedef seq2<T, U> type;
};
template <typename T, typename U, typename V> struct insert<seq3<T, U, V>, T> {
  typedef seq3<T, U, V> type;
};
template <typename T, typename U, typename V> struct insert<seq3<T, U, V>, U> {
  typedef seq3<T, U, V> type;
};
template <typename T, typename U, typename V> struct insert<seq3<T, U, V>, V> {
  typedef seq3<T, U, V> type;
};

template < typename Subdomain_index,
           typename Surface_patch_index,
           typename Curves_index,
           typename Corner_index>
struct Index_generator_with_features
{
  typedef typename insert<
    typename insert<
      typename insert<seq1<Subdomain_index>,
                      Surface_patch_index
                      >::type,
      Curves_index
      >::type,
    Corner_index>::type seq;
  typedef typename seq::type Index;
  typedef Index                                        type;
};

template < typename T>
struct Index_generator_with_features<T, T, T, T>
{
  typedef T Index;
  typedef Index                                         type;
};

template <typename T, typename Boost_variant>
const T& get_index(const Boost_variant& x,
                   typename boost::disable_if<boost::is_same<T, Boost_variant> >::type * = 0)
{ return boost::get<T>(x); }

template <typename T>
const T& get_index(const T& x) { return x; }

template <typename Mesh_domain,
          bool has_feature = Has_features<Mesh_domain>::value>
struct Read_mesh_domain_index {
  // here we have has_feature==true

  typedef Mesh_domain MT; // was named "mesh traits" previously

  typename Mesh_domain::Index
  operator()(int dimension, std::istream& is) const {
    switch(dimension) {
    case 0:
      typename MT::Corner_index ci;
      if(is_ascii(is)) is >> ci;
      else CGAL::read(is, ci);
      return  ci;
      break;
    case 1:
      typename MT::Curve_index si;
      if(is_ascii(is)) is >> si;
      else CGAL::read(is, si);
      return  si;
      break;
    default:
      return Read_mesh_domain_index<Mesh_domain, false>()(dimension, is);
    }
  }
}; // end template partial specialization
   // Read_mesh_domain_index<Mesh_domain, true>

template <typename Mesh_domain,
          bool has_feature = Has_features<Mesh_domain>::value>
struct Write_mesh_domain_index {
  // here we have has_feature==true

  typedef Mesh_domain MT; // was named "mesh traits" previously
  typedef typename MT::Corner_index Ci;
  typedef typename MT::Curve_index  Si;

  void
  operator()(std::ostream& os, int dimension,
             const typename Mesh_domain::Index& index) const {
    switch(dimension) {
    case 0: {
      const Ci& ci = get_index<Ci>(index);
      if(is_ascii(os)) os << oformat(ci);
      else CGAL::write(os, ci);
    }
      break;
    case 1: {
      const Si& si = get_index<Si>(index);
      if(is_ascii(os)) os << oformat(si);
      else CGAL::write(os, si);
    }
      break;
    default:
      Write_mesh_domain_index<Mesh_domain, false>()(os, dimension, index);
    }
  }
}; // end template partial specialization
   // Write_mesh_domain_index<Mesh_domain, true>

template <typename Mesh_domain>
struct Read_mesh_domain_index<Mesh_domain, false> {
  // here we have has_feature==false

  typedef Mesh_domain MT; // was named "mesh traits" previously

  typename Mesh_domain::Index
  operator()(int dimension, std::istream& is) const {
    switch(dimension) {
    case 2: {
      typename MT::Surface_patch_index spi;
      if(is_ascii(is)) is >> iformat(spi);
      else CGAL::read(is, spi);
      return  spi;
    }
      break;
    default: {// 3
      typename MT::Subdomain_index di;
      if(is_ascii(is)) is >> iformat(di);
      else CGAL::read(is, di);
      return  di;
    }
      break;
    }
  }
}; // end template partial specialization
   // Read_mesh_domain_index<Mesh_domain, false>

template <typename Mesh_domain>
struct Write_mesh_domain_index<Mesh_domain, false> {
  // here we have has_feature==false

  typedef Mesh_domain MT; // was named "mesh traits" previously
  typedef typename MT::Surface_patch_index Spi;
  typedef typename MT::Subdomain_index Di;

  void
  operator()(std::ostream& os, int dimension,
             const typename Mesh_domain::Index& index) const {
    switch(dimension) {
    case 2: {
      const Spi& spi = get_index<Spi>(index);
      if(is_ascii(os)) os << oformat(spi);
      else CGAL::write(os, spi);
    }
      break;
    default: {// 3
      const Di& di = get_index<Di>(index);
      if(is_ascii(os)) os << oformat(di);
      else CGAL::write(os, di);
    }
      break;
    }
  }
}; // end template partial specialization
   // Write_mesh_domain_index<Mesh_domain, false>

template <typename, typename Index>
struct Read_write_index {
  void operator()(std::ostream& os, int, Index index) const {
    if(is_ascii(os)) os << oformat(index);
    else CGAL::write(os, index);
  }
  Index operator()(std::istream& is, int) const {
    Index index;
    if(is_ascii(is)) is >> iformat(index);
    else CGAL::read(is, index);
    return index;
  }
};

struct Variant_write_visitor {
  std::ostream& os;
  template <typename T>
  void operator()(T v) const {
    if(is_ascii(os)) os << CGAL::oformat(v);
    else CGAL::write(os, v);
  }
};

template <typename Index>
struct Variant_read_visitor {
  std::istream& is;
  Index& variant;
  template <typename T>
  void operator()(T) const {
    T v;
    if(is_ascii(is)) is >> CGAL::iformat(v);
    else CGAL::read(is, v);
    variant = v;
  }
};

template <typename Indices_types, typename... Args>
struct Read_write_index<Indices_types, boost::variant<Args...>> {
  using Index = boost::variant<Args...>;
  using index_seq = std::make_index_sequence<std::tuple_size<Indices_types>::value>;

  template <std::size_t... Is>
  Index get_index(int dimension, std::index_sequence<Is...>) const{
    static const Index variants[] = { std::tuple_element_t<Is, Indices_types>{}... };
    return variants[dimension < 0 ? 0 : 3-dimension];
  }

  void operator()(std::ostream& os, int, Index index) const {
    Variant_write_visitor visitor{os};
    apply_visitor(visitor, index);
  }
  Index operator()(std::istream& is, int dimension) const {
    Index index = get_index(dimension, index_seq{});
    Variant_read_visitor<Index> visitor{is, index};
    apply_visitor(visitor, index);
    return index;
  }
};

} // end namespace internal
} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_INTERNAL_MESH_3_INDICES_MANAGEMENT_H
