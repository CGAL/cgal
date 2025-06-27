// Copyright (c) 2019-2023 Google LLC (USA).
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

#ifndef CGAL_ALPHA_WRAP_TRIANGULATION_FACE_BASE_2_H
#define CGAL_ALPHA_WRAP_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/Delaunay_triangulation_face_base_with_circumcenter_2.h>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

enum class Face_label
{
  // Faces that have been carved
  OUTSIDE = 0,
  // Faces that have not yet been carved
  INSIDE,
  // OUTSIDE faces that have been labeled "inside" again as to make the result manifold
  MANIFOLD
};

std::ostream& operator <<(std::ostream& os, const Face_label& label)
{
   os << static_cast<std::underlying_type<Face_label>::type>(label);
   return os;
}

template < typename GT,
           typename Fb = CGAL::Delaunay_triangulation_face_base_with_circumcenter_2<GT> >
class Alpha_wrap_triangulation_face_base_2
  : public Fb
{
public:
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;

public:
  template < typename TDS2 >
  struct Rebind_TDS
  {
    using Cb2 = typename Fb::template Rebind_TDS<TDS2>::Other;
    using Other = Alpha_wrap_triangulation_face_base_2<GT, Cb2>;
  };

private:
  Face_label m_label = Face_label::INSIDE;

#ifndef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
  unsigned int m_erase_counter;
#endif

public:
  Alpha_wrap_triangulation_face_base_2()
    : Fb()
  {}

  Alpha_wrap_triangulation_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
    : Fb(v0, v1, v2)
  {}

  Alpha_wrap_triangulation_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
                                       Face_handle   n0, Face_handle   n1, Face_handle   n2)
    : Fb(v0, v1, v2, n0, n1, n2)
  {}

public:
  Face_label label() const { return m_label; }
  void set_label(const Face_label label) { m_label = label; }
  bool is_inside() const { return m_label == Face_label::INSIDE; }
  bool is_outside() const { return m_label == Face_label::OUTSIDE; }

#ifndef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
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

template <typename Fb>
class Face_base_with_timestamp
  : public Fb
{
  std::size_t time_stamp_ = std::size_t(-2);

public:
  using Has_timestamp = CGAL::Tag_true;

  template <class TDS>
  struct Rebind_TDS
  {
    using Fb2 = typename Fb::template Rebind_TDS<TDS>::Other;
    using Other = Face_base_with_timestamp<Fb2>;
  };

public:
  template <typename... Args>
  Face_base_with_timestamp(const Args&... args)
    : Fb(args...)
  { }

  Face_base_with_timestamp(const Face_base_with_timestamp& other)
    : Fb(other), time_stamp_(other.time_stamp_)
  { }

public:
  std::size_t time_stamp() const { return time_stamp_; }
  void set_time_stamp(const std::size_t& ts) { time_stamp_ = ts; }
};

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_TRIANGULATION_FACE_BASE_2_H
