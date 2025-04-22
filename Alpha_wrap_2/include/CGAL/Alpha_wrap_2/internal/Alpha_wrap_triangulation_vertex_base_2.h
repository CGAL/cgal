// Copyright (c) 2019-2023 Google LLC (USA).
// Copyright (c) 2025 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_ALPHA_WRAP_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_ALPHA_WRAP_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

enum class Vertex_type
{
  DEFAULT = 0,
  BBOX_VERTEX,
  SEED_VERTEX
};

template <typename GT,
          typename Vb = Triangulation_vertex_base_2<GT> >
class Alpha_wrap_triangulation_vertex_base_2
  : public Vb
{
private:
  Vertex_type vertex_type = Vertex_type::DEFAULT;

public:
  using Face_handle = typename Vb::Face_handle;
  using Point = typename Vb::Point;

  template <typename TDS2>
  struct Rebind_TDS
  {
    using Vb2 = typename Vb::template Rebind_TDS<TDS2>::Other;
    using Other = Alpha_wrap_triangulation_vertex_base_2<GT, Vb2>;
  };

public:
  Alpha_wrap_triangulation_vertex_base_2()
    : Vb() {}

  Alpha_wrap_triangulation_vertex_base_2(const Point& p)
    : Vb(p) {}

  Alpha_wrap_triangulation_vertex_base_2(const Point& p, Face_handle f)
    : Vb(p, f) {}

  Alpha_wrap_triangulation_vertex_base_2(Face_handle f)
    : Vb(f) {}

public:
  const Vertex_type& type() const { return vertex_type; }
  Vertex_type& type() { return vertex_type; }
};

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_TRIANGULATION_VERTEX_BASE_2_H
