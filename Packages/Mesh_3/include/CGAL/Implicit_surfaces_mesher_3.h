// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_IMPLICIT_SURFACES_MESHER_3_H
#define CGAL_IMPLICIT_SURFACES_MESHER_3_H

#include <CGAL/Mesh_3/Refine_tets.h>
#include <CGAL/Chew_4_surfaces/Chew_4_surfaces.h>

namespace CGAL {

  namespace Mesh_3 
  {
    

template <
  typename Tr,
  typename Oracle,
  typename Facets_criteria,
  typename Tets_criteria
  >
class Implicit_surfaces_mesher_3
{
public:
  typedef typename
     Chew_4_surfaces::Chew_4_surfaces<Tr,
                                      Oracle,
                                      Facets_criteria> Facets_level;

  typedef typename Mesh_3::Refine_tets<Tr,
                                       Tets_criteria,
                                       Mesh_3::Refine_tets_base<Tr,
                                                                Tets_criteria>,
                                       Facets_level> Tets_level;

private:
  Null_mesher_level null_mesher_level;
  Null_mesh_visitor null_visitor;

  Facets_level facets;
  Tets_level tets;

  bool initialized;

public:
  Implicit_surfaces_mesher_3(Tr& t, Oracle& s,
                             Facets_criteria& c,
                             Tets_criteria tets_crit)
    : facets(t, s, c), tets(t, tets_crit, facets), initialized(false)
  {}

  void init()
  {
    facets.scan_triangulation();
    tets.scan_triangulation();
    initialized = true;
  }

  void refine_mesh()
  {
    if(!initialized)
      init();
    tets.refine(null_visitor);
  }
}; // end Implicit_surfaces_mesher_3

} // end namespace CGAL

#endif // CGAL_IMPLICIT_SURFACES_MESHER_3_H
