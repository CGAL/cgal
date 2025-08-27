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
    const CGAL::Epic_converter<AK> converter;

    std::tuple<decltype(converter(approx(args)))...> converted_args;

    // When C++20 is available, check the blame and clean this all up with lambdas
    bool success = convert_all_impl(std::index_sequence_for<Args...>{}, converter, converted_args, args...);
    if(!success) // failed to convert all arguments, call the base predicate
      return fp(args...);

    return call_epicp_impl(std::index_sequence_for<Args...>{}, converted_args);
  }

private:
  template <std::size_t... I, typename Converter, typename Tuple, typename... Args>
  bool convert_all_impl(std::index_sequence<I...>,
                        const Converter& converter,
                        Tuple& converted_args,
                        const Args&... args) const
  {
    auto convert = [&](auto index, const auto& arg) {
      auto converted = converter(approx(arg));
      if(converted.second)
        std::get<index>(converted_args) = converted;
      return converted.second;
    };

    return (... && convert(std::integral_constant<std::size_t, I>{}, args));
  }

  template <std::size_t... I, typename Tuple>
  auto call_epicp_impl(std::index_sequence<I...>, const Tuple& converted_args) const
    -> decltype(epicp(std::get<I>(converted_args).first...))
  {
    return epicp(std::get<I>(converted_args).first...);
  }
};

} // CGAL

#endif // CGAL_EPIC_PREDICATE_IF_CONVERTIBLE_H
