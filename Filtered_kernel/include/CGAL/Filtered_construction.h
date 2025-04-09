// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_FILTERED_CONSTRUCTION_H
#define CGAL_FILTERED_CONSTRUCTION_H

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>

#include <type_traits>

namespace CGAL {

template <class AC, class EC, class FC, class C2E, class C2F,
          class E2C, class F2C,        bool Protection = true>
class Filtered_construction
{
private:
  EC Exact_construction;
  FC Filter_construction;
  C2E To_Exact;
  C2F To_Filtered;
  E2C From_Exact;
  F2C From_Filtered;

public:
  Filtered_construction() {}

  template <class A1>
  auto operator()(const A1 &a1) const
    -> typename std::invoke_result<AC, A1>::type
  {
    {
      Protect_FPU_rounding<Protection> P1;
      try
      {
        return From_Filtered(Filter_construction(To_Filtered(a1)));
      }
      catch (Uncertain_conversion_exception&)
      {}
    }
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
    return From_Exact(Exact_construction(To_Exact(a1)));
  }

  template <class A1, class A2>
  auto operator()(const A1 &a1, const A2 &a2) const
    -> typename std::invoke_result<AC, A1, A2>::type
  {
    {
      Protect_FPU_rounding<Protection> P1;
      try
      {
        return From_Filtered(Filter_construction(To_Filtered(a1),
                                                 To_Filtered(a2)));
      }
      catch (Uncertain_conversion_exception&)
      {}
    }
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
    return From_Exact(Exact_construction(To_Exact(a1),
                                         To_Exact(a2)));
  }

  template <class A1, class A2, class A3>
  auto operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
    -> typename std::invoke_result<AC, A1, A2, A3>::type
  {
    {
      Protect_FPU_rounding<Protection> P1;
      try
      {
        return From_Filtered(Filter_construction(To_Filtered(a1),
                                                 To_Filtered(a2),
                                                 To_Filtered(a3)));
      }
      catch (Uncertain_conversion_exception&)
      {}
    }
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
    return From_Exact(Exact_construction(To_Exact(a1),
                                         To_Exact(a2),
                                         To_Exact(a3)));
  }
};

} // namespace CGAL

#endif // CGAL_FILTERED_CONSTRUCTION_H
