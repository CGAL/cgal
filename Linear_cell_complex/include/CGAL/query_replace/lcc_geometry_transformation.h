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
