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

#ifndef CGAL_ALPHA_WRAP_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_ALPHA_WRAP_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/TDS_3/Indexed_storage.h>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

enum class Vertex_type
{
  DEFAULT = 0,
  BBOX_VERTEX,
  SEED_VERTEX
};

template <typename GT,
          typename Vb = Triangulation_vertex_base_3<GT> >
class Alpha_wrap_triangulation_vertex_base_3
  : public Vb
{
private:
  Vertex_type vertex_type = Vertex_type::DEFAULT;

public:
  using Cell_handle = typename Vb::Cell_handle;
  using Point = typename Vb::Point;

  template <typename TDS2>
  struct Rebind_TDS
  {
    using Vb2 = typename Vb::template Rebind_TDS<TDS2>::Other;
    using Other = Alpha_wrap_triangulation_vertex_base_3<GT, Vb2>;
  };

public:
  Alpha_wrap_triangulation_vertex_base_3()
    : Vb() {}

  Alpha_wrap_triangulation_vertex_base_3(const Point& p)
    : Vb(p) {}

  Alpha_wrap_triangulation_vertex_base_3(const Point& p, Cell_handle c)
    : Vb(p, c) {}

  Alpha_wrap_triangulation_vertex_base_3(Cell_handle c)
    : Vb(c) {}

public:
  const Vertex_type& type() const { return vertex_type; }
  Vertex_type& type() { return vertex_type; }
};

template <typename GT,
            typename Vb = VertexWithPoint<GT>>
class Vertex4Alpha_wrap_3
  : public Vb
  {
    public:
    using Vb::Vb; // inherit constructors
    using Point = typename GT::Point_3;
    using TDS = typename Vb::Triangulation_data_structure;
    using Vertex_handle = typename TDS::Vertex_handle;
    using Cell_handle = typename TDS::Cell_handle;

    struct Storage : public Vb::Storage {
      Vertex_type vertex_type = Vertex_type::DEFAULT;
    };

    template < typename TDS2 >
    struct Rebind_TDS {
      using Vb2 = typename Vb::template Rebind_TDS<TDS2>::Other;
      using Other = Vertex4Alpha_wrap_3<GT,Vb2>;
    };
    auto&& storage() {
      return this->tds()->vertex_storage()[this->index()];
    }

    auto&& storage() const { return this->tds()->vertex_storage()[this->index()]; }

  public:

  const Vertex_type& type() const { return storage().vertex_type; }
  Vertex_type& type() { return storage().vertex_type; }
  };

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_TRIANGULATION_VERTEX_BASE_3_H
