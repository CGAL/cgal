// Copyright (c) 2019-2023 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_ALPHA_WRAP_TRIANGULATION_CELL_BASE_3_H
#define CGAL_ALPHA_WRAP_TRIANGULATION_CELL_BASE_3_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

enum class Cell_label
{
  // Cells that have been carved
  OUTSIDE,
  // Cells that have not yet been carved
  INSIDE,
  // OUTSIDE cells that have been labeled "inside" again as to make the result manifold
  MANIFOLD
};

template < typename GT,
           typename Cb = CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<GT> >
class Alpha_wrap_triangulation_cell_base_3
  : public Cb
{
public:
  typedef typename Cb::Vertex_handle                   Vertex_handle;
  typedef typename Cb::Cell_handle                     Cell_handle;

public:
  template < typename TDS2 >
  struct Rebind_TDS
  {
    using Cb2 = typename Cb::template Rebind_TDS<TDS2>::Other;
    using Other = Alpha_wrap_triangulation_cell_base_3<GT, Cb2>;
  };

private:
  Cell_label m_label = Cell_label::INSIDE;

#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
  unsigned int m_erase_counter;
#endif

public:
  Alpha_wrap_triangulation_cell_base_3()
    : Cb()
  {}

  Alpha_wrap_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                                       Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3)
  {}

  Alpha_wrap_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                                       Vertex_handle v2, Vertex_handle v3,
                                       Cell_handle   n0, Cell_handle   n1,
                                       Cell_handle   n2, Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3)
  {}

public:
  Cell_label label() const { return m_label; }
  void set_label(const Cell_label label) { m_label = label; }
  bool is_inside() const { return m_label == Cell_label::INSIDE; }
  bool is_outside() const { return m_label == Cell_label::OUTSIDE; }

#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
  unsigned int erase_counter() const
  {
    return m_erase_counter;
  }
  void set_erase_counter(unsigned int c)
  {
    m_erase_counter = c;
  }
  void increment_erase_counter()
  {
    ++m_erase_counter;
  }
#endif
};

template <typename Cb>
class Cell_base_with_timestamp
  : public Cb
{
  std::size_t time_stamp_;

public:
  using Has_timestamp = CGAL::Tag_true;

  template <class TDS>
  struct Rebind_TDS
  {
    using Cb2 = typename Cb::template Rebind_TDS<TDS>::Other;
    using Other = Cell_base_with_timestamp<Cb2>;
  };

public:
  template <typename... Args>
  Cell_base_with_timestamp(const Args&... args)
    : Cb(args...), time_stamp_(-1)
  { }

  Cell_base_with_timestamp(const Cell_base_with_timestamp& other)
    : Cb(other), time_stamp_(other.time_stamp_)
  { }

public:
  std::size_t time_stamp() const { return time_stamp_; }
  void set_time_stamp(const std::size_t& ts) { time_stamp_ = ts; }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_TRIANGULATION_CELL_BASE_3_H
