// Copyright (c) 2011 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_MESH_3_POLYLINE_WITH_CONTEXT_H
#define CGAL_MESH_3_POLYLINE_WITH_CONTEXT_H

#include <CGAL/license/Mesh_3.h>


#include <set>
#include <functional>

namespace CGAL { namespace Mesh_3 {

template <typename Surface_patch_id,
          typename Curve_id_>
struct Polyline_context {
  typedef std::set<Surface_patch_id> Patches_ids;
  typedef Curve_id_ Curve_id;

  Curve_id id;
  Patches_ids adjacent_patches_ids;
};


template <typename Surface_patch_id,
          typename Curve_id,
          typename Bare_polyline_>
struct Polyline_with_context
{
  typedef Bare_polyline_ Bare_polyline;
  typedef typename Bare_polyline::value_type Point_3;
  typedef Polyline_context<Surface_patch_id,Curve_id> Context;

  Context context;
  Bare_polyline polyline_content;
};

/**
 * @class Extract_bare_polyline
 * A functor to get bare polyline from Polyline_with_context (Pwc_).
 */
template <typename Pwc_>
struct Extract_bare_polyline :
  public CGAL::cpp98::unary_function<Pwc_, const typename Pwc_::Bare_polyline&>
{
  typedef CGAL::cpp98::unary_function<Pwc_, const typename Pwc_::Bare_polyline&> Base;
  typedef typename Base::result_type                                     result_type;
  typedef typename Base::argument_type                                   argument_type;

  result_type operator()(const argument_type& p) const { return p.polyline_content; }
};

}} // end namespaces

#endif // CGAL_MESH_3_POLYLINE_WITH_CONTEXT_H
