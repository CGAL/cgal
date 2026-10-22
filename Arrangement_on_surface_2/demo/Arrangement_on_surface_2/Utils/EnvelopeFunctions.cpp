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
//            Efi Fogel <efifogel@gmail.com>

#include "EnvelopeFunctions.h"
#include "ArrangementTypes.h"

#include <CGAL/Envelope_2/envelope_2.h>
#include <CGAL/Envelope_2/Envelope_diagram_1.h>
#include <vector>

template <typename Arr_>
auto EnvelopeFunctions<Arr_>::getXMonotoneCurves(Arrangement* arr)
  -> std::vector<X_monotone_curve_2>
{
  std::vector<X_monotone_curve_2> curves;
  for (auto it = arr->edges_begin(); it != arr->edges_end(); ++it)
    curves.push_back(it->curve());
  return curves;
}

template <typename Arr_>
void EnvelopeFunctions<Arr_>::lowerEnvelope(
  Arrangement* arr, Diagram_1& diagram)
{
  auto curves = getXMonotoneCurves(arr);
  CGAL::Envelope_2::lower_envelope_x_monotone_2(
    curves.begin(), curves.end(), diagram);
}

template <typename Arr_>
void EnvelopeFunctions<Arr_>::upperEnvelope(
  Arrangement* arr, Diagram_1& diagram)
{
  auto curves = getXMonotoneCurves(arr);
  CGAL::Envelope_2::upper_envelope_x_monotone_2(
    curves.begin(), curves.end(), diagram);
}

ARRANGEMENT_DEMO_SPECIALIZE_ARR(EnvelopeFunctions)
