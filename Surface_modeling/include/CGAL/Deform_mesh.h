// Copyright (c) 2014 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Yin Xu, Andreas Fabri and Ilker O. Yaz

#ifndef CGAL_DEFORM_MESH_H
#define CGAL_DEFORM_MESH_H

#ifdef DOXYGEN_RUNNING
template <
  class HG,
  class VIM=Default,
  class HIM=Default,
  Deformation_algorithm_tag TAG = SPOKES_AND_RIMS,
  class WC = Default,
  class ST = Default,
  class CR = Default,
  class VPM = Default
  >
class Surface_mesh_deformation;
#endif

#ifndef CGAL_NO_DEPRECATED_CODE

#define CGAL_DEPRECATED_HEADER "<CGAL/Deform_mesh.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Surface_mesh_deformation.h>"
#include <CGAL/internal/deprecation_warning.h>

#include <CGAL/Surface_mesh_deformation.h>

namespace CGAL {

 ///
 /// \ingroup PkgSurfaceModeling
 /// Class renamed to `Surface_mesh_deformation`.
 /// \deprecated This class name is deprecated and has been renamed to `Surface_mesh_deformation`.
template <
  class HG,
  class VIM=Default,
  class HIM=Default,
  Deformation_algorithm_tag TAG = SPOKES_AND_RIMS,
  class WC = Default,
  class ST = Default,
  class CR = Default,
  class VPM = Default
  >
class Deform_mesh :  public Surface_mesh_deformation<HG, VIM, HIM, TAG, WC, ST, CR, VPM>
{
  typedef Deform_mesh<HG, VIM, HIM, TAG, WC, ST, CR, VPM> Self;
  typedef Surface_mesh_deformation<HG, VIM, HIM, TAG, WC, ST, CR, VPM> Base;
#ifndef CGAL_CFG_NO_CPP0X_DELETED_AND_DEFAULT_FUNCTIONS
public:
  Deform_mesh(const Self&) = delete; // no copy
#else
private:
  Deform_mesh(const Self&); // no copy
#endif

public:
  typedef typename Base::Halfedge_graph Halfedge_graph;
  typedef typename Base::Vertex_index_map Vertex_index_map;
  typedef typename Base::Hedge_index_map Hedge_index_map;
  typedef typename Base::Weight_calculator Weight_calculator;
  typedef typename Base::Vertex_point_map Vertex_point_map;

  //vertex_point_map set by default
  Deform_mesh(Halfedge_graph& halfedge_graph,
              Vertex_index_map vertex_index_map,
              Hedge_index_map hedge_index_map
             )
    : Base(halfedge_graph, vertex_index_map, hedge_index_map)
  {}

  //vertex_point_map and hedge_index_map set by default
  Deform_mesh(Halfedge_graph& halfedge_graph,
              Vertex_index_map vertex_index_map)
    : Base(halfedge_graph, vertex_index_map)
  {}
  //vertex_point_map, hedge_index_map and vertex_index_map set by default
  Deform_mesh(Halfedge_graph& halfedge_graph)
    : Base(halfedge_graph)
  {}

  // Constructor with all the parameters provided
  Deform_mesh(Halfedge_graph& halfedge_graph,
              Vertex_index_map vertex_index_map,
              Hedge_index_map hedge_index_map,
              Vertex_point_map vertex_point_map,
              Weight_calculator weight_calculator = Weight_calculator()
             )
    : Base(halfedge_graph, vertex_index_map, hedge_index_map, vertex_point_map, weight_calculator)
  {}
};
} //namespace CGAL

#endif //CGAL_NO_DEPRECATED_CODE

#endif  // CGAL_DEFORM_MESH_H
