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

#include "IntersectCurves.h"
#include "ArrangementTypes.h"

template <typename Traits_>
Intersect_curves<Traits_>::Intersect_curves(const Traits* traits) :
    intersect{traits->intersect_2_object()}
{
}

template <typename Traits_>
void Intersect_curves<Traits_>::operator()(
  const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2,
  std::vector<CGAL::Object>& output)
{
  this->intersect(cv1, cv2, std::back_inserter(output));
}

ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(Intersect_curves)
