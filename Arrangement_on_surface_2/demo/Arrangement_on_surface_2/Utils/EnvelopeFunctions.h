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

#ifndef ARRANGEMENT_DEMO_ENVELOPE_FUNCTIONS
#define ARRANGEMENT_DEMO_ENVELOPE_FUNCTIONS

#include <vector>
#include <CGAL/Envelope_diagram_1.h>

template <typename Arr_>
class EnvelopeFunctions
{
public:
  using Arrangement = Arr_;
  using Traits = typename Arrangement::Geometry_traits_2;
  using X_monotone_curve_2 = typename Arrangement::X_monotone_curve_2;
  using Diagram_1 = CGAL::Envelope_diagram_1<Traits>;

  void lowerEnvelope(Arrangement* arr, Diagram_1& diagram);
  void upperEnvelope(Arrangement* arr, Diagram_1& diagram);

private:
  std::vector<X_monotone_curve_2> getXMonotoneCurves(Arrangement* arr);
};

#endif
