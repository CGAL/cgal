// Copyright (c) 2017  GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri, Laurent Rineau

#ifndef CGAL_EPIC_PREDICATE_IF_CONVERTIBLE_H
#define CGAL_EPIC_PREDICATE_IF_CONVERTIBLE_H

#include <CGAL/Epic_converter.h>

namespace CGAL {

template <typename AK, typename FP, typename EpicP>
class EPIC_predicate_if_convertible {
public:
  FP fp;
  EpicP epicp;

  template <typename... Args>
  auto
  operator()(const Args&... args) const
    -> decltype(fp(args...))
  {
    CGAL::Epic_converter<AK> converter;

    std::tuple<decltype(converter(approx(args)))...> converted_args;

    auto convert = [&](auto index, const auto& arg)
    {
      auto converted = converter(approx(arg));
      if(converted.second)
        std::get<index>(converted_args) = converted;
      return converted.second;
    };

    bool success = [&]<std::size_t... I>(std::index_sequence<I...>) {
      return (... && convert(std::integral_constant<std::size_t, I>{}, args));
    }(std::index_sequence_for<Args...>{});

    if(!success) // failed to convert all arguments, call the base predicate
      return fp(args...);

    return [&]<std::size_t... I>(std::index_sequence<I...>) {
      return epicp(std::get<I>(converted_args).first...);
    }(std::index_sequence_for<Args...>{});
  }
};

} // CGAL

#endif // CGAL_EPIC_PREDICATE_IF_CONVERTIBLE_H
