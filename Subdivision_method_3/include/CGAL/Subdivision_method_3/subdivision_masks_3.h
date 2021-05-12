// Copyright (c) 2005-2017 GeometryFactory (France).  All Rights Reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Le-Jeng Shiue <Andy.Shiue@gmail.com>
//

#ifndef CGAL_SUBDIVISION_MASKS_3_H
#define CGAL_SUBDIVISION_MASKS_3_H

#include <CGAL/basic.h>
#include <CGAL/Origin.h>

#include <CGAL/circulator.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/property_maps.h>

namespace CGAL {

// ======================================================================
/// The stencil of the Primal-Quadrilateral-Quadrisection
template <class PolygonMesh,
          class VertexPointMap = typename boost::property_map<PolygonMesh,
                                                              vertex_point_t>::type >
class PQQ_stencil_3 {

public:
  typedef PolygonMesh                                               Mesh;
  typedef typename boost::property_map<Mesh, vertex_point_t>::type  Vertex_pmap;

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor       face_descriptor;

  typedef typename boost::property_traits<Vertex_pmap>::value_type  Point;

  typedef typename Kernel_traits<Point>::Kernel                     Kernel;
  typedef typename Kernel::FT                                       FT;
  typedef typename Kernel::Vector_3                                 Vector;

  Mesh* pmesh;
  VertexPointMap vpmap;

public:
  PQQ_stencil_3(Mesh* pmesh)
    : pmesh(pmesh), vpmap(get(vertex_point, pmesh))
  { }

  PQQ_stencil_3(Mesh* pmesh, VertexPointMap vpmap)
    : pmesh(pmesh), vpmap(vpmap)
  { }

  void face_node(face_descriptor, Point&) { }
  void edge_node(halfedge_descriptor, Point&) { }
  void vertex_node(vertex_descriptor, Point&) { }

  void border_node(halfedge_descriptor, Point&, Point&) { }
};

// ======================================================================
/// Bi-linear geometry mask for PQQ, PTQ, and Sqrt(3) schemes
template <class PolygonMesh,
          class VertexPointMap = typename boost::property_map<PolygonMesh,
                                                              vertex_point_t>::type >
class Linear_mask_3 : public PQQ_stencil_3<PolygonMesh, VertexPointMap> {
  typedef PQQ_stencil_3<PolygonMesh, VertexPointMap> Base;
public:
  typedef PolygonMesh                                Mesh;

#ifndef DOXYGEN_RUNNING
  typedef typename Base::vertex_descriptor           vertex_descriptor;
  typedef typename Base::halfedge_descriptor         halfedge_descriptor;
  typedef typename Base::face_descriptor             face_descriptor;

  typedef typename Base::Kernel                      Kernel;

  typedef typename Base::FT                          FT;
  typedef typename Base::Point                       Point;
  typedef typename Base::Vector                      Vector;

  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref;
#endif

public:
  Linear_mask_3(Mesh* pmesh)
    : Base(pmesh, get(vertex_point, pmesh))
  { }

  Linear_mask_3(Mesh* pmesh, VertexPointMap vpmap)
    : Base(pmesh, vpmap)
  { }

  void face_node(face_descriptor facet, Point& pt) {
    int n = 0;
    Point p(0,0,0);
    for(vertex_descriptor vd :
                  vertices_around_face(halfedge(facet, *(this->pmesh)), *(this->pmesh)))
    {
      p = p + ( get(this->vpmap,vd) - ORIGIN);
      ++n;
    }
    pt = ORIGIN + (p - ORIGIN)/FT(n);
  }

  void edge_node(halfedge_descriptor edge, Point& pt) {
    Point_ref p1 = get(this->vpmap, target(edge, *(this->pmesh)));
    Point_ref p2 = get(this->vpmap, source(edge, *(this->pmesh)));
    pt = Point((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2);
  }

  void vertex_node(vertex_descriptor vertex, Point& pt) {
    pt = get(this->vpmap, vertex);
  }

  void border_node(halfedge_descriptor edge, Point& ept, Point& /*vpt*/){
    edge_node(edge, ept);
  }
};

// ======================================================================
/*!
\ingroup PkgSurfaceSubdivisionMethod3Ref

The geometry mask of Catmull-Clark subdivision.

A stencil determines a source neighborhood
whose points contribute to the position of a refined point.
The geometry mask of a stencil specifies
the computation on the nodes of the stencil.
`CatmullClark_mask_3` implements the geometry masks of
Catmull-Clark subdivision on models of `MutableFaceGraph`,
such as `Polyhedron_3` and `Surface_mesh`.

\tparam PolygonMesh must be a model of the concept `MutableFaceGraph`.
\tparam VertexPointMap must be a model of `WritablePropertyMap` with value type `Point_3`

\image html CCBorderMask.svg


\cgalModels `PQQMask_3`

\sa `CGAL::Subdivision_method_3`
*/
template <class PolygonMesh,
          class VertexPointMap = typename boost::property_map<PolygonMesh,
                                                              vertex_point_t>::type >
class CatmullClark_mask_3 : public Linear_mask_3<PolygonMesh, VertexPointMap> {
  typedef Linear_mask_3<PolygonMesh, VertexPointMap> Base;
public:
  typedef PolygonMesh                                Mesh;

#ifndef DOXYGEN_RUNNING
  typedef typename Base::vertex_descriptor           vertex_descriptor;
  typedef typename Base::halfedge_descriptor         halfedge_descriptor;
  typedef typename Base::face_descriptor             face_descriptor;

  typedef typename Base::Kernel                      Kernel;

  typedef typename Base::FT                          FT;
  typedef typename Base::Point                       Point;
  typedef typename Base::Vector                      Vector;

  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref;
#endif

public:
/// \name Creation
/// @{

  /// Constructor.
  /// The default vertex point property map, `get(vertex_point, pmesh)`, is used.
  CatmullClark_mask_3(Mesh* pmesh)
    : Base(pmesh, get(vertex_point, pmesh))
  { }

  /// Constructor with a custom vertex point property map.
  CatmullClark_mask_3(Mesh* pmesh, VertexPointMap vpmap)
    : Base(pmesh, vpmap)
  { }

/// @}

/// \name Stencil functions
/// @{

  /// computes the Catmull-Clark edge-point `pt` of the edge `edge`.
  void edge_node(halfedge_descriptor edge, Point& pt) {
    Point_ref p1 = get(this->vpmap,target(edge, *(this->pmesh)));
    Point_ref p2 = get(this->vpmap,source(edge, *(this->pmesh)));
    Point f1, f2;
    this->face_node(face(edge, *(this->pmesh)), f1);
    this->face_node(face(opposite(edge, *(this->pmesh)), *(this->pmesh)), f2);
    pt = Point((p1[0]+p2[0]+f1[0]+f2[0])/4,
               (p1[1]+p2[1]+f1[1]+f2[1])/4,
               (p1[2]+p2[2]+f1[2]+f2[2])/4 );
  }

  /// computes the Catmull-Clark vertex-point `pt` of the vertex `vertex`.
  void vertex_node(vertex_descriptor vertex, Point& pt) {
    Halfedge_around_target_circulator<Mesh> vcir(vertex, *(this->pmesh));
    int n = static_cast<int>(degree(vertex, *(this->pmesh)));

    FT Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
    Point_ref S = get(this->vpmap,vertex);

    Point q;
    for (int i = 0; i < n; i++, ++vcir) {
      Point_ref p2 = get(this->vpmap, target(opposite(*vcir, *(this->pmesh)), *(this->pmesh)));
      R[0] += (S[0] + p2[0]) / 2;
      R[1] += (S[1] + p2[1]) / 2;
      R[2] += (S[2] + p2[2]) / 2;
      this->face_node(face(*vcir, *(this->pmesh)), q);
      Q[0] += q[0];
      Q[1] += q[1];
      Q[2] += q[2];
    }
    R[0] /= n;    R[1] /= n;    R[2] /= n;
    Q[0] /= n;    Q[1] /= n;    Q[2] /= n;

    pt = Point((Q[0] + 2*R[0] + S[0]*(n-3))/n,
               (Q[1] + 2*R[1] + S[1]*(n-3))/n,
               (Q[2] + 2*R[2] + S[2]*(n-3))/n );
  }

  /// computes the Catmull-Clark edge-point `ept` and the Catmull-Clark
  /// vertex-point `vpt` of the border edge `edge`.
  void border_node(halfedge_descriptor edge, Point& ept, Point& vpt) {
    Point_ref ep1 = get(this->vpmap,target(edge, *(this->pmesh)));
    Point_ref ep2 = get(this->vpmap,target(opposite(edge, *(this->pmesh)), *(this->pmesh)));
    ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

    Halfedge_around_target_circulator<Mesh> vcir(edge, *(this->pmesh));
    Point_ref vp1  = get(this->vpmap,target(opposite(*vcir, *(this->pmesh)), *(this->pmesh)));
    Point_ref vp0  = get(this->vpmap, target(*vcir, *(this->pmesh)));
    --vcir;
    Point_ref vp_1 = get(this->vpmap, target(opposite(*vcir, *(this->pmesh)), *(this->pmesh)));
    vpt = Point((vp_1[0] + 6*vp0[0] + vp1[0])/8,
                (vp_1[1] + 6*vp0[1] + vp1[1])/8,
                (vp_1[2] + 6*vp0[2] + vp1[2])/8 );
  }
/// @}
};

// ======================================================================
/*!
\ingroup PkgSurfaceSubdivisionMethod3Ref

The geometry mask of Loop subdivision.

A stencil determines a source neighborhood
whose points contribute to the position of a refined point.
The geometry mask of a stencil specifies
the computation on the nodes of the stencil.
`Loop_mask_3` implements the geometry masks of
Loop subdivision on a triangulated model of `MutableFaceGraph`,
such as `Polyhedron_3` and `Surface_mesh`.

\tparam PolygonMesh must be a model of the concept `MutableFaceGraph`. Additionally all faces must be triangles.
\tparam VertexPointMap must be a model of `WritablePropertyMap` with value type `Point_3`

\image html LoopBorderMask.png
\image latex LoopBorderMask.png

\cgalModels `PTQMask_3`

\sa `CGAL::Subdivision_method_3`

*/
template <class PolygonMesh,
          class VertexPointMap = typename boost::property_map<PolygonMesh,
                                                              vertex_point_t>::type >
class Loop_mask_3 : public PQQ_stencil_3<PolygonMesh, VertexPointMap> {
  typedef PQQ_stencil_3<PolygonMesh, VertexPointMap>   Base;
public:
  typedef PolygonMesh                                Mesh;

#ifndef DOXYGEN_RUNNING
  typedef typename Base::vertex_descriptor           vertex_descriptor;
  typedef typename Base::halfedge_descriptor         halfedge_descriptor;
  typedef typename Base::face_descriptor             face_descriptor;

  typedef typename Base::Kernel                      Kernel;

  typedef typename Base::FT                          FT;
  typedef typename Base::Point                       Point;
  typedef typename Base::Vector                      Vector;
  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref;
#endif

  typedef Halfedge_around_face_circulator<Mesh> Halfedge_around_facet_circulator;
  typedef Halfedge_around_target_circulator<Mesh> Halfedge_around_vertex_circulator;

public:
/// \name Creation
/// @{

  /// Constructor.
  /// The default vertex point property map, `get(vertex_point, pmesh)`, is used.
  Loop_mask_3(Mesh* pmesh)
    : Base(pmesh, get(vertex_point, pmesh))
  { }

  /// Constructor with a custom vertex point property map.
  Loop_mask_3(Mesh* pmesh, VertexPointMap vpmap)
    : Base(pmesh, vpmap)
  { }

/// @}

/// \name Stencil functions
/// @{

  /// computes the Loop edge-point `pt` of the edge `edge`.
  void edge_node(halfedge_descriptor edge, Point& pt) {
    Point_ref p1 = get(this->vpmap,target(edge, *(this->pmesh)));
    Point_ref p2 = get(this->vpmap,target(opposite(edge, *(this->pmesh)), *(this->pmesh)));
    Point_ref f1 = get(this->vpmap,target(next(edge, *(this->pmesh)), *(this->pmesh)));
    Point_ref f2 = get(this->vpmap,target(next(opposite(edge, *(this->pmesh)), *(this->pmesh)), *(this->pmesh)));

    pt = Point((3*(p1[0]+p2[0])+f1[0]+f2[0])/8,
               (3*(p1[1]+p2[1])+f1[1]+f2[1])/8,
               (3*(p1[2]+p2[2])+f1[2]+f2[2])/8 );
  }

  /// computes the Loop vertex-point `pt` of the vertex `vertex`.
  void vertex_node(vertex_descriptor vertex, Point& pt) {
    Halfedge_around_vertex_circulator vcir(vertex, *(this->pmesh));
    size_t n = circulator_size(vcir);

    FT R[] = {0.0, 0.0, 0.0};
    Point_ref S = get(this->vpmap,vertex);

    for (size_t i = 0; i < n; i++, ++vcir) {
      Point_ref p = get(this->vpmap,target(opposite(*vcir, *(this->pmesh)), *(this->pmesh)));
      R[0] += p[0];         R[1] += p[1];         R[2] += p[2];
    }
    if (n == 6) {
      pt = Point((10*S[0]+R[0])/16, (10*S[1]+R[1])/16, (10*S[2]+R[2])/16);
    } else {
      const FT Cn = (FT) (5.0/8.0 - CGAL::square(3+2*std::cos(2 * CGAL_PI/(double) n))/64.0);
      const FT Sw = (double)n*(1-Cn)/Cn;
      const FT W = (double)n/Cn;
      pt = Point((Sw*S[0]+R[0])/W, (Sw*S[1]+R[1])/W, (Sw*S[2]+R[2])/W);
    }
  }

  //
  //void face_node(face_descriptor facet, Point& pt) {};
  //

  /// computes the Loop edge-point `ept` and the Loop vertex-point `vpt` of the border edge `edge`.
  void border_node(halfedge_descriptor edge, Point& ept, Point& vpt) {
    Point_ref ep1 = get(this->vpmap,target(edge, *(this->pmesh)));
    Point_ref ep2 = get(this->vpmap,target(opposite(edge, *(this->pmesh)), *(this->pmesh)));
    ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

    Halfedge_around_vertex_circulator vcir(edge, *(this->pmesh));
    Point_ref vp1  = get(this->vpmap,target(opposite(*vcir, *(this->pmesh) ), *(this->pmesh)));
    Point_ref vp0  = get(this->vpmap,target(*vcir, *(this->pmesh)));
    --vcir;
    Point_ref vp_1 = get(this->vpmap,target(opposite(*vcir, *(this->pmesh)), *(this->pmesh)));
    vpt = Point((vp_1[0] + 6*vp0[0] + vp1[0])/8,
                (vp_1[1] + 6*vp0[1] + vp1[1])/8,
                (vp_1[2] + 6*vp0[2] + vp1[2])/8 );
  }

/// @}
};


//==========================================================================
/// The stencil of the Dual-Quadrilateral-Quadrisection
template <class PolygonMesh,
          class VertexPointMap = typename boost::property_map<PolygonMesh,
                                                              vertex_point_t>::type >
class DQQ_stencil_3 {
public:
  typedef PolygonMesh                                              Mesh;

  typedef typename boost::property_map<Mesh, vertex_point_t>::type Vertex_pmap;

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor      face_descriptor;

  typedef typename boost::property_traits<Vertex_pmap>::value_type Point;

  typedef typename Kernel_traits<Point>::Kernel                    Kernel;
  typedef typename Kernel::FT                                      FT;
  typedef typename Kernel::Vector_3                                Vector;

  Mesh* pmesh;
  Vertex_pmap vpm;

public:
  DQQ_stencil_3(Mesh* pmesh)
    : pmesh(pmesh), vpm(get(vertex_point, pmesh))
  { }

  DQQ_stencil_3(Mesh* pmesh, VertexPointMap vpmap)
    : pmesh(pmesh), vpm(vpmap)
  { }
};


// ======================================================================
/*!
\ingroup PkgSurfaceSubdivisionMethod3Ref

The geometry mask of Doo-Sabin subdivision.

A stencil determines a source neighborhood
whose points contribute to the position of a refined point.
The geometry mask of a stencil specifies
the computation on the nodes of the stencil.
`DooSabin_mask_3` implements the geometry masks of
Doo-Sabin subdivision on models of `MutableFaceGraph`,
such as `Polyhedron_3` and `Surface_mesh`.

\tparam PolygonMesh must be a model of the concept `MutableFaceGraph`.
\tparam VertexPointMap must be a model of `WritablePropertyMap` with value type `Point_3`

\image html DSCornerMask.png
\image latex DSCornerMask.png

\cgalModels `DQQMask_3`

\sa `CGAL::Subdivision_method_3`

*/
template <class PolygonMesh,
          class VertexPointMap = typename boost::property_map<PolygonMesh,
                                                              vertex_point_t>::type >
class DooSabin_mask_3 : public DQQ_stencil_3<PolygonMesh, VertexPointMap> {
  typedef DQQ_stencil_3<PolygonMesh, VertexPointMap> Base;
public:
  typedef PolygonMesh                                Mesh;

#ifndef DOXYGEN_RUNNING
  typedef typename Base::vertex_descriptor           vertex_descriptor;
  typedef typename Base::halfedge_descriptor         halfedge_descriptor;
  typedef typename Base::face_descriptor             face_descriptor;

  typedef typename Base::Kernel                      Kernel;

  typedef typename Base::FT                          FT;
  typedef typename Base::Point                       Point;
  typedef typename Base::Vector                      Vector;
  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref;
#endif

public:
/// \name Creation
/// @{

  /// Constructor.
  /// The default vertex point property map, `get(vertex_point, pmesh)`, is used.
  DooSabin_mask_3(Mesh* pmesh)
    : Base(pmesh, get(vertex_point, pmesh))
  { }

  /// Constructor with a custom vertex point property map.
  DooSabin_mask_3(Mesh* pmesh, VertexPointMap vpmap)
    : Base(pmesh, vpmap)
  { }

/// @}

/// \name Stencil functions
/// @{

  /// computes the Doo-Sabin point `pt` of the vertex pointed by the halfedge `he`.
  void corner_node(halfedge_descriptor he, Point& pt) {
    size_t n = 0;
    halfedge_descriptor hd = he;
    do{
      hd = next(hd, *(this->pmesh));
      ++n;
    } while(hd != he);

    Vector cv(0,0,0);
    if (n == 4) {
      cv = cv + (get(this->vpm, target(he, *(this->pmesh))) - CGAL::ORIGIN)*9;
      cv = cv + (get(this->vpm, target(next(he, *(this->pmesh)), *(this->pmesh)))-CGAL::ORIGIN)*3;
      cv = cv + (get(this->vpm, target(next(next(he, *(this->pmesh)), *(this->pmesh)), *(this->pmesh)))-CGAL::ORIGIN);
      cv = cv + (get(this->vpm, target(prev(he, *(this->pmesh)), *(this->pmesh)))-CGAL::ORIGIN)*3;
      cv = cv/16;
    }  else {
      FT a;
      for (size_t k = 0; k < n; ++k, he = next(he, *(this->pmesh))) {
        if (k == 0) a = (FT) ((5.0/n) + 1);
        else a = (FT) (3+2*std::cos(2*k*CGAL_PI/n))/n;
        cv = cv + (get(this->vpm, target(he, *(this->pmesh)))-CGAL::ORIGIN)*a;
      }
      cv = cv/4;
    }
    pt = CGAL::ORIGIN + cv;
  }
/// @}
};

// ======================================================================
/*!
\ingroup PkgSurfaceSubdivisionMethod3Ref

The geometry mask of Sqrt(3) subdivision.

A stencil determines a source neighborhood
whose points contribute to the position of a refined point.
The geometry mask of a stencil specifies
the computation on the nodes of the stencil.
`Sqrt3_mask_3` implements the geometry masks of
\f$ \sqrt{3}\f$ subdivision on a triangulated
model of `MutableFaceGraph`,
such as `Polyhedron_3` and `Surface_mesh`.

\tparam PolygonMesh must be a model of the concept `MutableFaceGraph`. Additionally all faces must be triangles.
\tparam VertexPointMap must be a model of `WritablePropertyMap` with value type `Point_3`

\cgalModels `Sqrt3Mask_3`

\sa `CGAL::Subdivision_method_3`

*/
template <class PolygonMesh,
          class VertexPointMap = typename boost::property_map<PolygonMesh,
                                                              vertex_point_t>::type >
class Sqrt3_mask_3 : public Linear_mask_3<PolygonMesh, VertexPointMap> {
  typedef Linear_mask_3<PolygonMesh, VertexPointMap>       Base;
public:
  typedef PolygonMesh                                Mesh;

#ifndef DOXYGEN_RUNNING
  typedef typename Base::vertex_descriptor           vertex_descriptor;
  typedef typename Base::halfedge_descriptor         halfedge_descriptor;
  typedef typename Base::face_descriptor             face_descriptor;

  typedef typename Base::Kernel                      Kernel;

  typedef typename Base::FT                          FT;
  typedef typename Base::Point                       Point;
  typedef typename Base::Vector                      Vector;
#endif

public:
/// \name Creation
/// @{

  /// Constructor.
  /// The default vertex point property map, `get(vertex_point, pmesh)`, is used.
  Sqrt3_mask_3(Mesh* pmesh)
    : Base(pmesh, get(vertex_point, pmesh))
  { }

  /// Constructor with a custom vertex point property map.
  Sqrt3_mask_3(Mesh* pmesh, VertexPointMap vpmap)
    : Base(pmesh, vpmap)
  { }

/// @}

/// \name Stencil functions
/// @{

  /// computes the \f$ \sqrt{3}\f$ vertex-point `pt` of the vertex `vd`.
  void vertex_node(vertex_descriptor vertex, Point& pt) {
    Halfedge_around_target_circulator<Mesh> vcir(vertex, *(this->pmesh));
    const typename boost::graph_traits<Mesh>::degree_size_type n = degree(vertex, *(this->pmesh));

    const FT a = (FT) ((4.0-2.0*std::cos(2.0*CGAL_PI/(double)n))/9.0);

    Vector cv = ((FT)(1.0-a)) * (get(this->vpmap, vertex) - CGAL::ORIGIN);
    for (typename boost::graph_traits<Mesh>::degree_size_type i = 1; i <= n; ++i, --vcir) {
      cv = cv + (a/FT(n))*(get(this->vpmap, target(opposite(*vcir, *(this->pmesh)), *(this->pmesh)))-CGAL::ORIGIN);
    }

    pt = CGAL::ORIGIN + cv;
  }

  /// computes the \f$ \sqrt{3}\f$ edge-points `ept1` and `ept2` of the halfedge `hd`.
  /// The updated point coordinates for the target vertex of `hd` is also computed and put in `vpt`.
  /// \attention Border subdivision only happens every second step of a <em>single</em>
  ///            successive \f$ \sqrt{3}\f$ subdivision (thus requiring a depth larger than 1).
  void border_node(halfedge_descriptor hd, Point& ept1, Point& ept2, Point& vpt) {
    // this function takes the opposite of a BORDER halfedge
    halfedge_descriptor bhd = opposite(hd, *(this->pmesh));
    CGAL_precondition(is_border(bhd, *(this->pmesh)));
    vertex_descriptor prev_s = source(prev(bhd, *(this->pmesh)), *(this->pmesh));
    Vector prev_sv = get(this->vpmap, prev_s) - CGAL::ORIGIN;

    vertex_descriptor s = source(bhd, *(this->pmesh));
    Vector sv = get(this->vpmap, s) - CGAL::ORIGIN;

    vertex_descriptor t = target(bhd, *(this->pmesh));
    Vector tv = get(this->vpmap, t) - CGAL::ORIGIN;

    vertex_descriptor next_t = target(next(bhd, *(this->pmesh)), *(this->pmesh));
    Vector next_tv = get(this->vpmap, next_t) - CGAL::ORIGIN;

    const FT denom = 1./27.;
    ept1 = CGAL::ORIGIN + denom * ( prev_sv + 16.*sv + 10.*tv );
    ept2 = CGAL::ORIGIN + denom * ( 10.*sv + 16.*tv + next_tv );
    vpt = CGAL::ORIGIN + 1./27. * ( 4*prev_sv + 19*sv + 4*tv );
  }
/// @}
};

} // namespace CGAL

#endif // CGAL_SUBDIVISION_MASKS_3_H
