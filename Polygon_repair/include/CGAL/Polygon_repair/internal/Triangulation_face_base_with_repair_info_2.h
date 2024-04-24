// Copyright (c) 2023 GeometryFactory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ken Arroyo Ohori

#ifndef CGAL_TRIANGULATION_WITH_REPAIR_INFO_2_H
#define CGAL_TRIANGULATION_WITH_REPAIR_INFO_2_H

#include <CGAL/license/Polygon_repair.h>

#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL {
namespace Polygon_repair {
namespace internal {

template <typename Kernel, typename FaceBase = Triangulation_face_base_2<Kernel>>
class Triangulation_face_base_with_repair_info_2 : public FaceBase {
  int _label;
  bool _processed;
public:
  using Vertex_handle = typename FaceBase::Vertex_handle;
  using Face_handle = typename FaceBase::Face_handle;

  template <typename TDS2>
  struct Rebind_TDS {
    using FaceBase2 = typename FaceBase::template Rebind_TDS<TDS2>::Other;
    using Other = Triangulation_face_base_with_repair_info_2<Kernel, FaceBase2>;
  };

  Triangulation_face_base_with_repair_info_2() : FaceBase() {}

  Triangulation_face_base_with_repair_info_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
    : FaceBase(v0, v1, v2) {}

  Triangulation_face_base_with_repair_info_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
                                             Face_handle   n0, Face_handle   n1, Face_handle   n2)
    : FaceBase(v0, v1, v2, n0, n1, n2) {}

  const bool& processed() const { return _processed; }
  bool& processed() { return _processed; }
  const int& label() const { return _label; }
  int& label() { return _label; }
};

} // namespace internal
} // namespace Polygon_repair
} //namespace CGAL

#endif // CGAL_TRIANGULATION_WITH_REPAIR_INFO_2_H
