// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#include "ArrangementIO.h"
#include "ArrangementDemoTab.h"
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "Conic_reader.h"

#include <CGAL/IO/Arr_text_formatter.h>
#include <CGAL/IO/Arr_with_history_iostream.h>
#include <CGAL/IO/Arr_with_history_text_formatter.h>

template <
  typename Arrangement,
  typename Traits = typename Arrangement::Geometry_traits_2>
struct ArrReader
{
  Arrangement* operator()(std::ifstream& ifs)
  {
    using Text_formatter = CGAL::Arr_text_formatter<Arrangement>;
    using ArrFormatter = CGAL::Arr_with_history_text_formatter<Text_formatter>;

    ArrFormatter arrFormatter;
    auto arr = new Arrangement();
    CGAL::read(*arr, ifs, arrFormatter);
    return arr;
  }
};

#ifdef CGAL_USE_CORE
template <
  typename Arrangement, typename Rat_kernel_, typename Alg_kernel_,
  typename Nt_traits_>
struct ArrReader<
  Arrangement, CGAL::Arr_conic_traits_2<Rat_kernel_, Alg_kernel_, Nt_traits_>>
{
  using Traits = typename Arrangement::Geometry_traits_2;
  using Curve_2 = typename Arrangement::Curve_2;

  Arrangement* operator()(std::ifstream& ifs)
  {
    Conic_reader<Traits> conicReader;
    std::vector<Curve_2> curve_list;
    CGAL::Bbox_2 bbox;
    conicReader.read_data(ifs, std::back_inserter(curve_list), bbox);
    auto arr = new Arrangement();
    CGAL::insert(*arr, curve_list.begin(), curve_list.end());
    return arr;
  }
};

template <
  typename Arrangement, typename Rat_kernel_, typename Alg_kernel_,
  typename Nt_traits_, typename Bounding_traits_>
struct ArrReader<
  Arrangement, CGAL::Arr_Bezier_curve_traits_2<
                 Rat_kernel_, Alg_kernel_, Nt_traits_, Bounding_traits_>>
{
  Arrangement* operator()(std::ifstream&) { return nullptr; }
};

template <typename Arrangement, typename AlgebraicKernel_d_1_>
struct ArrReader<
  Arrangement, CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1_>>
{
  Arrangement* operator()(std::ifstream&) { return nullptr; }
};
#endif

std::pair<CGAL::Object, demo_types::TraitsType>
ArrangementIO::read(std::ifstream& ifs)
{
  // read type info
  while (ifs.peek() == '#' || std::isspace(ifs.peek())) ifs.get();

  int tt_int;
  ifs >> tt_int;
  auto tt = static_cast<demo_types::TraitsType>(tt_int);

  std::pair<CGAL::Object, demo_types::TraitsType> res;
  demo_types::visitArrangementType(tt, [&](auto type_holder) {
    using Arrangement = typename decltype(type_holder)::type;
    auto arr = ArrReader<Arrangement>{}(ifs);
    res = {CGAL::make_object(arr), demo_types::enumFromArrType<Arrangement>()};
  });
  return res;
}

template <
  typename Arrangement,
  typename Traits = typename Arrangement::Geometry_traits_2>
struct ArrWriter
{
  void operator()(Arrangement* arr, std::ofstream& ofs)
  {
    using TextFormatter = CGAL::Arr_text_formatter<Arrangement>;
    using ArrFormatter = CGAL::Arr_with_history_text_formatter<TextFormatter>;

    ArrFormatter arrFormatter;
    CGAL::write(*arr, ofs, arrFormatter);
  }
};

#ifdef CGAL_USE_CORE
template <
  typename Arrangement, typename Rat_kernel_, typename Alg_kernel_,
  typename Nt_traits_>
struct ArrWriter<
  Arrangement, CGAL::Arr_conic_traits_2<Rat_kernel_, Alg_kernel_, Nt_traits_>>
{
  using Traits = typename Arrangement::Geometry_traits_2;
  using Curve_2 = typename Arrangement::Curve_2;

  void operator()(Arrangement* arr, std::ofstream& ofs)
  {
    Conic_reader<Traits> conicReader;
    conicReader.write_data(ofs, arr->curves_begin(), arr->curves_end());
  }
};

template <
  typename Arrangement, typename Rat_kernel_, typename Alg_kernel_,
  typename Nt_traits_, typename Bounding_traits_>
struct ArrWriter<
  Arrangement, CGAL::Arr_Bezier_curve_traits_2<
                 Rat_kernel_, Alg_kernel_, Nt_traits_, Bounding_traits_>>
{
  void operator()(Arrangement*, std::ofstream&) { }
};

template <typename Arrangement, typename AlgebraicKernel_d_1_>
struct ArrWriter<
  Arrangement, CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1_>>
{
  void operator()(Arrangement*, std::ofstream&) { }
};
#endif

bool ArrangementIO::write(
  const std::pair<CGAL::Object, demo_types::TraitsType>& arr_pair,
  std::ofstream& ofs)
{
  auto tt = arr_pair.second;
  auto arr_obj = arr_pair.first;

  // write type info
  ofs << "# " << static_cast<int>(tt) << std::endl;

  bool result = false;
  demo_types::visitArrangementType(tt, [&](auto type_holder) {
    using Arrangement = typename decltype(type_holder)::type;
    Arrangement* arr;
    if (CGAL::assign(arr, arr_obj))
    {
      ArrWriter<Arrangement>{}(arr, ofs);
      result = true;
    }
  });
  return result;
}
