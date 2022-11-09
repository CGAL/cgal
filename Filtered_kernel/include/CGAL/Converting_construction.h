// Copyright (c) 2022  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CONVERTING_CONSTRUCTION_H
#define CGAL_CONVERTING_CONSTRUCTION_H

#include <CGAL/config.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Cartesian/Is_trivial_construction.h>
#include <type_traits>

namespace CGAL {

template <class Source_construction, class Target_construction,
          class Converter, class Backward_converter>
struct Converting_construction
{
  CGAL_NO_UNIQUE_ADDRESS Source_construction source_construction;
  CGAL_NO_UNIQUE_ADDRESS Target_construction construct;
  CGAL_NO_UNIQUE_ADDRESS Converter convert;
  CGAL_NO_UNIQUE_ADDRESS Backward_converter backward_convert;

  template <typename... Args,
            std::enable_if_t<CartesianFunctors::Is_trivial_construction<Source_construction, Args...>::value>* = nullptr>
  decltype(auto) operator()(Args&&... args) const {
    return source_construction(std::forward<Args>(args)...);
  }

  template <typename... Args,
            std::enable_if_t<!CartesianFunctors::Is_trivial_construction<Source_construction, Args...>::value>* = nullptr>
  auto operator()(Args&&... args) const {
    return backward_convert(construct(convert(args)...));
  }
};

template <class Kernel_A, class Kernel_B>
struct Converting_constructions_kernel_adaptor : public Kernel_A
{
  using A_to_B = Cartesian_converter<Kernel_A, Kernel_B>;
  using B_to_A = Cartesian_converter<Kernel_B, Kernel_A>;

#define CGAL_Kernel_cons(C,Cf) \
  using C = Converting_construction<typename Kernel_A::C, \
                                    typename Kernel_B::C, \
                                    A_to_B,               \
                                    B_to_A>;              \
  C Cf() const { return C(); }

#include <CGAL/Kernel/interface_macros.h>

};

} // end namespace CGAL

#endif // CGAL_CONVERTING_CONSTRUCTION_H
