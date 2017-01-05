// ======================================================================
//
// Copyright (c) 2005-2017 GeometryFactory (France).  All Rights Reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s): Le-Jeng Shiue <Andy.Shiue@gmail.com>
//
// ======================================================================

#ifndef CGAL_POLYHEDRON_SUBDIVISION_H_01292002
#define CGAL_POLYHEDRON_SUBDIVISION_H_01292002

#include <CGAL/basic.h>

#include <vector>

#include <CGAL/circulator.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Subdivision_method_impl_3.h>
#include <CGAL/Subdivision_mask_3.h>

namespace CGAL {

// ======================================================================
namespace Subdivision_method_3 {

  namespace parameters = Polygon_mesh_processing::parameters;
  // 
  template <class PolygonMesh, class Mask>
  void PQQ(PolygonMesh& pmesh, Mask mask, int step = 1) {
    // todo:  static assert that PolygonMesh == Mask::PolygonMesh
    for (int i = 0; i < step; i++) Private::PQQ_1step(pmesh, get(vertex_point,pmesh), mask);
  }  


  template <class PolygonMesh, class Mask, class NamedParameters>
  void PQQ(PolygonMesh& pmesh, Mask mask, const NamedParameters& np) {
    // todo:  static assert that PolygonMesh == Mask::PolygonMesh  
    using boost::choose_param;
    using boost::get_param;
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
    Vpm vpm = choose_param(get_param(np, vertex_point),
                           get_property_map(CGAL::vertex_point, pmesh));

    unsigned int step = choose_param(get_param(np, number_of_iterations), 1);
    for (int i = 0; i < step; i++) Private::PQQ_1step(pmesh, vpm, mask);
  }

  // 
  template <class PolygonMesh,  class Mask>
  void PTQ(PolygonMesh& pmesh, Mask mask, int step = 1) {
    for (int i = 0; i < step; i++) Private::PTQ_1step(pmesh, get(vertex_point,pmesh), mask);
  }

  template <class PolygonMesh, class Mask, class NamedParameters>
  void PTQ(PolygonMesh& pmesh, Mask mask, const NamedParameters& np) {
    // todo:  static assert that PolygonMesh == Mask::PolygonMesh  
    using boost::choose_param;
    using boost::get_param;
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
    Vpm vpm = choose_param(get_param(np, vertex_point),
                           get_property_map(CGAL::vertex_point, pmesh));

    unsigned int step = choose_param(get_param(np, number_of_iterations), 1);
    for (unsigned int i = 0; i < step; i++) Private::PTQ_1step(pmesh, vpm, mask);
  }


  // 
  template <class PolygonMesh, class Mask>
  void DQQ(PolygonMesh& pmesh, Mask mask, int step = 1) {
    for (int i = 0; i < step; ++i) Private::DQQ_1step(pmesh, get(vertex_point,pmesh), mask);
  }

  template <class PolygonMesh, class Mask, class NamedParameters>
  void DQQ(PolygonMesh& pmesh, Mask mask, const NamedParameters& np) {
    // todo:  static assert that PolygonMesh == Mask::PolygonMesh  
    using boost::choose_param;
    using boost::get_param;
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
    Vpm vpm = choose_param(get_param(np, vertex_point),
                           get_property_map(CGAL::vertex_point, pmesh));

    unsigned int step = choose_param(get_param(np, number_of_iterations), 1);
    for (unsigned int i = 0; i < step; i++) Private::DQQ_1step(p, vpm, mask);
  }

  // 
  template <class PolygonMesh, class Mask>
  void Sqrt3(PolygonMesh& pmesh, Mask mask, int step = 1) {
    for (int i = 0; i < step; i++) Private::Sqrt3_1step(pmesh, get(vertex_point,pmesh), mask);
  }  

  template <class PolygonMesh, class Mask, class NamedParameters>
  void Sqrt3(PolygonMesh& pmesh, Mask mask, const NamedParameters& np) {
    // todo:  static assert that PolygonMesh == Mask::PolygonMesh  
    using boost::choose_param;
    using boost::get_param;
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
    Vpm vpm = choose_param(get_param(np, vertex_point),
                           get_property_map(CGAL::vertex_point, pmesh));

    unsigned int step = choose_param(get_param(np, number_of_iterations), 1);
    for (unsigned int i = 0; i < step; i++) Private::Sqrt3_1step(p, vpm, mask);
  }

  //
  template <class PolygonMesh>
  void CatmullClark_subdivision(PolygonMesh& pmesh, int step = 1) {
    PQQ(pmesh, CatmullClark_mask_3<PolygonMesh>(pmesh,get(vertex_point,pmesh)), step);
  }
  
  template <class PolygonMesh, class NamedParameters>
  void CatmullClark_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
    using boost::choose_param;
    using boost::get_param;
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
    Vpm vpm = choose_param(get_param(np, vertex_point),
                           get_property_map(CGAL::vertex_point, pmesh));

    unsigned int step = choose_param(get_param(np, number_of_iterations), 1);
    CatmullClark_mask_3<PolygonMesh,Vpm> mask(pmesh,vpm);

    for(unsigned int i = 0; i < step; i++) Private::PQQ_1step(pmesh, vpm, mask);
  }
  

  //
  template <class PolygonMesh>
  void Loop_subdivision(PolygonMesh& pmesh, int step = 1) {
    PTQ(pmesh, Loop_mask_3<PolygonMesh>(pmesh,get(vertex_point,pmesh)) , step);
  }

  template <class PolygonMesh, class NamedParameters>
  void Loop_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
    using boost::choose_param;
    using boost::get_param;
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
    Vpm vpm = choose_param(get_param(np, vertex_point),
                           get_property_map(CGAL::vertex_point, pmesh));

    unsigned int step = choose_param(get_param(np, number_of_iterations), 1);
    Loop_mask_3<PolygonMesh,Vpm> mask(pmesh,vpm);

    for(unsigned int i = 0; i < step; i++) Private::PTQ_1step(pmesh, vpm, mask);
  }
  
  //
  template <class PolygonMesh>
  void DooSabin_subdivision(PolygonMesh& pmesh, int step = 1) {
    DQQ(pmesh, DooSabin_mask_3<PolygonMesh>(pmesh,get(vertex_point,pmesh)), step);
  }

  template <class PolygonMesh, class NamedParameters>
  void DooSabin_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
    using boost::choose_param;
    using boost::get_param;
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
    Vpm vpm = choose_param(get_param(np, vertex_point),
                           get_property_map(CGAL::vertex_point, pmesh));

    unsigned int step = choose_param(get_param(np, number_of_iterations), 1);
    DooSabin_mask_3<PolygonMesh,Vpm> mask(pmesh,vpm);

    for(unsigned int i = 0; i < step; i++) Private::DQQ_1step(pmesh, vpm, mask);
  }
  
  //
  template <class PolygonMesh>
  void Sqrt3_subdivision(PolygonMesh& pmesh, int step = 1) {
    Sqrt3(pmesh, Sqrt3_mask_3<PolygonMesh>(pmesh,get(vertex_point,pmesh)), step);
  }

  template <class PolygonMesh, class NamedParameters>
  void Sqrt_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
    using boost::choose_param;
    using boost::get_param;
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
    Vpm vpm = choose_param(get_param(np, vertex_point),
                           get_property_map(CGAL::vertex_point, pmesh));

    unsigned int step = choose_param(get_param(np, number_of_iterations), 1);
    Sqrt3_mask_3<PolygonMesh,Vpm> mask(pmesh,vpm);

    for(unsigned int i = 0; i < step; i++) Private::Sqrt_1step(pmesh, vpm, mask);
  }
  
}

} //namespace CGAL

#endif //CGAL_POLYHEDRON_SUBDIVISION_H_01292002
