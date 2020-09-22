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
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "ArrangementDemoTab.h"
#include "Conic_reader.h"

#include <CGAL/IO/Arr_text_formatter.h>
#include <CGAL/IO/Arr_with_history_iostream.h>
#include <CGAL/IO/Arr_with_history_text_formatter.h>

#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Arr_overlay_2.h>

using TraitsType = demo_types::TraitsType;

struct ArrReader
{
  template <typename Arrangement>
  auto operator()(demo_types::TypeHolder<Arrangement>)
  {
    using Text_formatter = CGAL::Arr_text_formatter<Arrangement>;
    using ArrFormatter = CGAL::Arr_with_history_text_formatter<Text_formatter>;

    ArrFormatter arrFormatter;
    auto arr = std::make_unique<Arrangement>();
    CGAL::read(*arr, ifs, arrFormatter);
    return arr;
  }

#ifdef CGAL_USE_Core
  auto operator()(demo_types::TypeHolder<demo_types::Conic_arr>)
  {
    using namespace demo_types;

    Conic_reader<Conic_arr::Geometry_traits_2> conicReader;
    std::vector<Conic_arr::Curve_2> curve_list;
    CGAL::Bbox_2 bbox;
    conicReader.read_data(ifs, std::back_inserter(curve_list), bbox);
    auto arr = std::make_unique<Conic_arr>();
    CGAL::insert(*arr, curve_list.begin(), curve_list.end());
    return arr;
  }

  auto operator()(demo_types::TypeHolder<demo_types::Bezier_arr>)
    -> std::unique_ptr<demo_types::Bezier_arr>
  {
    return nullptr;
  }

  auto operator()(demo_types::TypeHolder<demo_types::Rational_arr>)
    -> std::unique_ptr<demo_types::Rational_arr>
  {
    return nullptr;
  }
#endif

  std::ifstream& ifs;
};

ArrangementDemoTabBase* ArrangementIO::read(std::ifstream& ifs)
{
  // read type info
  while (ifs.peek() == '#' || std::isspace(ifs.peek()))
    ifs.get();

  int tt_int;
  ifs >> tt_int;
  auto tt = static_cast<TraitsType>(tt_int);

  ArrangementDemoTabBase* tab = nullptr;
  demo_types::visitArrangementType(tt, [&](auto type_holder) {
    using Arrangement = typename decltype(type_holder)::type;
    auto arr = ArrReader{ifs}(type_holder);
    tab = new ArrangementDemoTab<Arrangement>(nullptr, std::move(arr));
  });
  return tab;
}

struct ArrWriter
{
  template <typename Arrangement>
  void operator()(Arrangement* arr)
  {
    using TextFormatter = CGAL::Arr_text_formatter<Arrangement>;
    using ArrFormatter = CGAL::Arr_with_history_text_formatter<TextFormatter>;

    ArrFormatter arrFormatter;
    CGAL::write(*arr, ofs, arrFormatter);
  }

#ifdef CGAL_USE_Core
  void operator()(demo_types::Conic_arr* arr)
  {
    Conic_reader<demo_types::Conic_arr::Geometry_traits_2> conicReader;
    conicReader.write_data(ofs, arr->curves_begin(), arr->curves_end());
  }

  void operator()(demo_types::Bezier_arr*) { }
  void operator()(demo_types::Rational_arr*) { }
#endif

  std::ofstream& ofs;
};

bool ArrangementIO::write(ArrangementDemoTabBase* tab, std::ofstream& ofs)
{
  auto tt = tab->traitsType();

  // write type info
  ofs << "# " << static_cast<int>(tt) << std::endl;

  bool result = false;
  demo_types::visitArrangementType(tt, [&](auto type_holder) {
    using Arrangement = typename decltype(type_holder)::type;
    Arrangement* arr;
    if (CGAL::assign(arr, tab->getArrangement()))
    {
      ArrWriter{ofs}(arr);
      result = true;
    }
  });
  return result;
}

ArrangementDemoTabBase* makeOverlayTab(const std::vector<CGAL::Object>& arrs)
{
  ArrangementDemoTabBase* tab = nullptr;
  if (arrs.size() == 2)
  {
    demo_types::forEachArrangementType([&](auto type_holder) {
      using Arrangement = typename decltype(type_holder)::type;

      Arrangement* arr1;
      Arrangement* arr2;
      if (!tab && CGAL::assign(arr1, arrs[0]) && CGAL::assign(arr2, arrs[1]))
      {
        auto overlayArr = std::make_unique<Arrangement>();
        CGAL::Arr_default_overlay_traits<Arrangement> defaultTraits;

        CGAL::overlay(*arr1, *arr2, *overlayArr, defaultTraits);
        tab = new ArrangementDemoTab<Arrangement>(nullptr, std::move(overlayArr));
      }
    });
  }
  return tab;
}
