// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_COST_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_COST_H 1

#include <CGAL/license/Surface_mesh_simplification.h>


#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core.h>

namespace CGAL {

namespace Surface_mesh_simplification  
{

  template<class TM_>
class LindstromTurk_cost
{
public:
    
  typedef TM_ TM ;
  /*

  typedef Edge_profile<TM> Profile ;
  typedef typename Traits::Point_3 Point;
  typedef typename Traits::FT FT ;
  
  typedef optional<FT> result_type ;
  */
public:

  LindstromTurk_cost( LindstromTurk_params const& aParams = LindstromTurk_params() ) : mParams(aParams) {}
     
  template <typename Profile>
  optional<typename Profile::FT>
  operator()( Profile const& aProfile, optional<typename Profile::Point> const& aPlacement ) const
  {
    return LindstromTurkCore<TM,Profile>(mParams,aProfile).compute_cost(aPlacement) ;
  }

private:

  LindstromTurk_params mParams ;    
};

} // namespace Surface_mesh_simplification

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_COST_H //
// EOF //
 
