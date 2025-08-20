// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
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
