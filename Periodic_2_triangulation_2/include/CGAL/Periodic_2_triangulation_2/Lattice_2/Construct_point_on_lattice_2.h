// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_P2T2_CONSTRUCT_POINT_ON_LATTICE_2_H
#define CGAL_P2T2_CONSTRUCT_POINT_ON_LATTICE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_offset_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Lattice_2.h>

namespace CGAL {
namespace P2T2 {
namespace internal {

template <typename Kernel_,
          typename Offset_ = CGAL::Periodic_2_offset_2,
          typename Construct_point_2_base_ = typename Kernel_::Construct_point_2>
class Construct_point_on_lattice_2
  : public Construct_point_2_base_
{
  typedef Construct_point_2_base_            Base;

  typedef Kernel_                            Kernel;
  typedef Offset_                            Offset;

  typedef typename Kernel::Point_2           Point;

  typedef Lattice_2<Kernel_>                 Lattice;

public:
  Construct_point_on_lattice_2(const Lattice* lattice, const Base& cp)
    : Base(cp), lattice_(lattice)
  { }

  using Base::operator();

  Point operator()(const Point& p, const Offset& o) const
  {
    CGAL_assertion(lattice_ != nullptr);
    return lattice_->translate_by_offset(p, o);
  }

private:
  const Lattice* lattice_;
};

} // namespace internal
} // namespace P2T2
} // namespace CGAL

#endif // CGAL_P2T2_CONSTRUCT_POINT_ON_LATTICE_2_H
