// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CMAP_3CLOSE_CC_H
#define CMAP_3CLOSE_CC_H

#include <cstddef>

/** 3-close a connected component
 *  @param dh a dart
 *  @return the number of new darts.
 */
template<typename LCC>
std::size_t close_cc_for_beta3(LCC& lcc, typename LCC::Dart_handle dh)
{
  std::size_t res=0;
  typename LCC::Dart_handle d, d2;

  for (auto it=lcc.template darts_of_cell<3>(dh).begin(),
       itend=lcc.template darts_of_cell<3>(dh).end(); it!=itend; ++it)
  {
    if(lcc.template is_free<3>(it))
    {
      d=lcc.create_dart();
      ++res;
      lcc.template link_beta_for_involution<3>(it, d);

      // Special cases for 0 and 1
      if (!lcc.template is_free<1>(it) &&
          !lcc.template is_free<3>(lcc.template beta<1>(it)))
      { lcc.template link_beta<1>(lcc.template beta<1,3>(it),d); }
      if (!lcc.template is_free<0>(it) &&
          !lcc.template is_free<3>(lcc.template beta<0>(it)))
      { lcc.template link_beta<0>(lcc.template beta<0,3>(it),d); }

      d2=lcc.template beta<2>(it);
      while (d2!=lcc.null_dart_handle &&
             !lcc.template is_free<2>(lcc.template beta<3>(d2)))
      { d2=lcc.template beta<3, 2>(d2); }
      if (d2!=lcc.null_dart_handle && !lcc.template is_free<3>(d2))
      {
        lcc.template basic_link_beta_for_involution<2>
            (lcc.template beta<3>(d2), d);
      }
    }
  }
  return res;
}

#endif // CMAP_3CLOSE_CC_H
