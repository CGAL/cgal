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
#ifndef LCC_GEOMETRY_TRANSFORMATION_H
#define LCC_GEOMETRY_TRANSFORMATION_H

#include <CGAL/Aff_transformation_3.h>
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void lcc_translate_cc(LCC& lcc, typename LCC::Dart_handle dh,
                      const typename LCC::Vector& v)
{
  typename LCC::Traits::Aff_transformation_3 translate(CGAL::TRANSLATION, v);
  for(auto it=lcc.template one_dart_per_incident_cell<0,4>(dh).begin(),
      itend=lcc.template one_dart_per_incident_cell<0,4>(dh).end(); it!=itend; ++it)
  { lcc.point(it)=translate(lcc.point(it)); }
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void lcc_translate_all(LCC& lcc, const typename LCC::Vector& v)
{
  typename LCC::Traits::Aff_transformation_3 translate(CGAL::TRANSLATION, v);
  for(auto itv=lcc.vertex_attribute_range().begin(),
      itvend=lcc.vertex_attribute_range().end(); itv!=itvend; ++itv)
  {
    lcc.point_of_vertex_attribute(itv)=
        translate(lcc.point_of_vertex_attribute(itv));
  }
}
///////////////////////////////////////////////////////////////////////////////
#endif // LCC_GEOMETRY_TRANSFORMATION_H
