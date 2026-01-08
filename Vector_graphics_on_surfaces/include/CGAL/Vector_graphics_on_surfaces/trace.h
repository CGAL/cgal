// Copyright (c) 2023-2026 GeometryFactory and Claudio Mancinelli.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli and Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_TRACE_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_TRACE_H

#include <CGAL/license/Vector_graphics_on_surfaces.h>

#include <CGAL/Vector_graphics_on_surfaces/locally_shortest_path.h>
#include <CGAL/Vector_graphics_on_surfaces/straightest_geodesic.h>
#include <CGAL/Vector_graphics_on_surfaces/recursive_de_Casteljau.h>

namespace CGAL {
namespace Vector_graphics_on_surfaces {

/*!
 * \ingroup VGSFunctions
 * computes the face location of each vertex of a 2D polygon on `tmesh`.
 * The vertices of the polygon are given as a pair of direction and distance with respect to a point corresponding to `center` on `tmesh`.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam K a model of `Kernel` with `K::FT` being a floating point number type (float or double)
 * \param center the location on `tmesh` used as reference. The y-axis used for coordinates is `halfedge(center.first, tmesh)`.
 * \param directions contains the direction one need to move from `center` to reach each vertex of the polygon.
 * \param lengths the distance one need to move from `center` along the direction at the same position in `directions` to reach each vertex of the polygon.
 * \param tmesh input triangle mesh supporting the vertices of the output polygon
 * \param solver container for the precomputed information. If not initialized, it will be initialized internally.
 * \return a face location for each vertex of the polygon
 * \todo add named parameters
 * \todo polygon orientation is not handled in the function and should be done outside of the function for now
 * \todo offer something better than a 2D vector for the direction
 * \todo directly handle polyline?
 * \todo why the first polygon vertex is duplicated by the function? (most probably for the example but it shouldn't be done here)
 */
template <class K, class TriangleMesh>
std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh,typename K::FT>>
trace_geodesic_polygon(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT> &center,
                       const std::vector<typename K::Vector_2>& directions,
                       const std::vector<typename K::FT>& lengths,
                       const TriangleMesh &tmesh,
                       const Dual_geodesic_solver<typename K::FT>& solver = {})
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  std::size_t n=directions.size();
  std::vector<PMP::Face_location<TriangleMesh,typename K::FT>> result;
  std::vector<PMP::Face_location<TriangleMesh, typename K::FT>> vertices(n);
#ifdef CGAL_DEBUG_BSURF
  std::cout << "trace_geodesic_polygon\n";
#endif
  for(std::size_t i=0;i<n;++i)
  {
    vertices[i]= straightest_geodesic<K>(center,directions[i],lengths[i],tmesh).back();
#ifdef CGAL_DEBUG_BSURF
    std::cout << PMP::construct_point(vertices[i], tmesh) << "\n";
#endif
  }
  std::vector<PMP::Edge_location<TriangleMesh,typename K::FT>> edge_locations;

  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }

  for(std::size_t i=0;i<n;++i)
  {
    edge_locations.clear();
    locally_shortest_path<typename K::FT>(vertices[i],vertices[(i+1)%n],tmesh, edge_locations, *solver_ptr);
    result.push_back(vertices[i]);
    for(auto& el : edge_locations)
    {
      result.push_back(PMP::to_face_location(el, tmesh));
    }
  }
  result.push_back(vertices[0]);

  return result;
}


/*!
 *
 *
 *  Doc Me
 *
 */
template <class K>
struct Default_placement
{
  using Vector_2 = typename K::Vector_2;
  using Point_2 = typename K::Point_2;
  using FT = typename K::FT;

  template <class TriangleMesh, class VPM>
  void register_polygons(const TriangleMesh& /* tmesh */,
                         VPM /* vpm */,
                         const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& /* center */,
                         const FT& /* scaling */,
                         const std::vector<std::vector<Point_2>>& polygons)
  {
    m_polygons = polygons;
    m_polygon_bboxes.reserve(polygons.size());

    for (const std::vector<Point_2>& polygon : polygons)
    {
      m_polygon_bboxes.push_back( bbox_2(polygon.begin(), polygon.end()) );
      m_gbox+=m_polygon_bboxes.back();
    }
  }

  template <class TriangleMesh, class VPM>
  CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>
  operator()(std::size_t pid,
             const TriangleMesh& tmesh,
             VPM vpm,
             const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>& center,
             const FT& scaling,
             const Dual_geodesic_solver<typename K::FT>& solver,
             std::vector<Vector_2>& directions, std::vector<FT>& lens)
  {
    using Impl = internal::Locally_shortest_path_imp<K, TriangleMesh, VPM>;
    using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
    namespace PMP = CGAL::Polygon_mesh_processing;

    const std::vector<Point_2>& polygon = m_polygons.value().get()[pid];

    // Vector_2  initial_dir(0,1);// TODO: fix this case (bad transport?)
    Vector_2  initial_dir(1,0);// TODO: input parameter or 2D PCA of the centers?

    Vector_2 dir( (-m_gbox.xmin()-m_gbox.xmax()+m_polygon_bboxes[pid].xmin()+m_polygon_bboxes[pid].xmax())/2.,
                  (-m_gbox.ymin()-m_gbox.ymax()+m_polygon_bboxes[pid].ymin()+m_polygon_bboxes[pid].ymax())/2. );
    std::vector< PMP::Face_location<TriangleMesh, FT> > spath =
      straightest_geodesic<K>(center, dir, scaling * std::sqrt(dir.squared_length()),tmesh);
    PMP::Face_location<TriangleMesh, FT> polygon_center = spath.back();

    // TODO: avoid using the shortest path and directly use the straightest!
    std::vector<PMP::Edge_location<TriangleMesh, FT>> shortest_path;
      locally_shortest_path<FT>(center, polygon_center, tmesh, shortest_path, solver);

    // update direction
    Vector_2 v = initial_dir;
    for(std::size_t i=0;i<shortest_path.size();++i)
    {
      halfedge_descriptor h_curr = halfedge(shortest_path[i].first, tmesh);
      v = Impl::parallel_transport_f2f(h_curr, v, vpm, tmesh);
    }

    std::vector<std::pair<FT, FT>> polar_coords =
      convert_to_polar_coordinates<K>(polygon,
                                      Point_2((m_polygon_bboxes[pid].xmin()+m_polygon_bboxes[pid].xmax())/2.,
                                              (m_polygon_bboxes[pid].ymin()+m_polygon_bboxes[pid].ymax())/2.));
    if (polygon.front()==polygon.back())
      polar_coords.pop_back();

    double theta = atan2(v.y(),v.x());

    lens.reserve(polar_coords.size());
    directions.reserve(polar_coords.size());

    for (const std::pair<double, double>& coord : polar_coords)
    {
      lens.push_back(scaling * coord.first);
      directions.emplace_back(std::cos(coord.second+theta), std::sin(coord.second+theta));
    }

    return polygon_center;
  }

  std::optional<
    std::reference_wrapper<
      const std::vector<std::vector<Point_2>>
  > > m_polygons;

  std::vector<Bbox_2> m_polygon_bboxes;
  Bbox_2 m_gbox;
};

/*!
 * \ingroup VGSFunctions
 * computes for each vertex of each polygon in `polygons` a face location on `tmesh`, where `center` represents the center of the 2D bounding box of the polygons.
 * This method computes the location of the center of the bounding box of each polygon on the mesh with respect to `center` and calls `trace_geodesic_polygon()` with that center with
 * appropriate directions and distances to have a consistent orientation for the polygons.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam K a model of `Kernel` with `K::FT` being a floating point number type (float or double)
 * \param center the location on `tmesh` corresponding to the center of the 2D bounding box of the polygons.
 * \param polygons 2D polygons
 * \param scaling a scaling factor to scale the polygons on `tmesh` (considering geodesic distances on `tmesh`)
 * \param tmesh input triangle mesh supporting the vertices of the output polygon
 * \param solver container for the precomputed information. If not initialized, it will be initialized internally.
 * \return a face location for each vertex of each polygon
 * \todo add named parameters
 * \todo polygon orientation is not handled in the function and should be done outside of the function for now
 * \todo for better rendering we can group polygons to have one center for the same group of polygon (useful for letters that are not simply connected)
 * \todo what if boundary is reached
 */
template <class K, class TriangleMesh, class Placement = Default_placement<K> >
std::vector<std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh,typename K::FT>>>
trace_geodesic_polygons(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT> &center,
                        const std::vector<std::vector<typename K::Point_2>>& polygons,
                        const typename K::FT scaling,
                        const TriangleMesh &tmesh,
                        const Dual_geodesic_solver<typename K::FT>& solver = {},
                              Placement placement = Placement())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  VPM vpm = get(CGAL::vertex_point, tmesh);

  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }

  std::vector<std::vector<PMP::Face_location<TriangleMesh,typename K::FT>>> result(polygons.size());
  placement.register_polygons(tmesh, vpm, center, scaling, polygons);

  for(std::size_t i=0;i<polygons.size();++i)
  {
    std::vector<typename K::Vector_2> directions;
    std::vector<typename K::FT> lens;
    auto polygon_center = placement(i, tmesh, vpm, center, scaling, *solver_ptr, directions, lens);
    result[i] = trace_geodesic_polygon<K>(polygon_center, directions, lens, tmesh, *solver_ptr);
  }

  return result;
}


template <class TriangleMesh, class K>
struct Placement_along_a_line
  : public Default_placement<K>
{
  using Base = Default_placement<K>;
  using Vector_2 = typename K::Vector_2;
  using Point_2 = typename K::Point_2;
  using FT = typename K::FT;

  template <class VPM>
  void register_polygons(const TriangleMesh& tmesh,
                         VPM vpm,
                         const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& center,
                         const FT& scaling,
                         const std::vector<std::vector<Point_2>>& polygons)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;
    Base::register_polygons(tmesh, vpm, center, scaling, polygons);

    m_c2 = Point_2( (this->m_gbox.xmin()+this->m_gbox.xmax())/2., (this->m_gbox.ymin()+this->m_gbox.ymax())/2. );
    Point_2 left_most(this->m_gbox.xmin(), m_c2.y());
    Point_2 right_most(this->m_gbox.xmax(), m_c2.y());
    double len = (this->m_gbox.xmax()-this->m_gbox.xmin())/2.;

    m_left_path = straightest_geodesic<K>(center, left_most-m_c2, scaling * len, tmesh);
    m_right_path = straightest_geodesic<K>(center, right_most-m_c2, scaling * len, tmesh);

    CGAL_assertion(m_left_path.size() >=2);
    CGAL_assertion(m_right_path.size() >=2);

  // TODO: precompute distances along supporting curve + stop when exceeding max distance needed
  }

  template <class VPM>
  auto get_polygon_center(const TriangleMesh& tmesh,
                          VPM vpm,
                          const std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>>& path,
                          double targetd)
  {
    namespace PMP=CGAL::Polygon_mesh_processing;

    // use left
    double acc=0.;
    std::size_t k=0;
    while(true)
    {
      double len = std::sqrt(squared_distance(PMP::construct_point(path[k], tmesh),
                                              PMP::construct_point(path[k+1], tmesh)));
      acc+=len;
      if (acc == targetd)
      {
        // TODO: if you land here and path[k] is on an input vertex, it might be that loc_k==loc_k1

        double theta=0;
        PMP::Face_location<TriangleMesh, typename K::FT> loc_k=path[k], loc_k1=path[k+1];
        CGAL_assertion_code(bool OK=)
        PMP::locate_in_common_face(loc_k, loc_k1, tmesh);
        CGAL_assertion(OK);
        std::array<Vector_2,3> flat_triangle =
          internal::init_flat_triangle<K>(halfedge(loc_k.first,tmesh),vpm,tmesh);
        Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
        Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

        Vector_2 dir2 = tgt-src;
        theta = atan2(dir2.y(), dir2.x());

        return std::make_pair(path[k+1], theta);
      }

      if (acc > targetd)
      {
        double excess = acc-targetd;

        PMP::Face_location<TriangleMesh, typename K::FT> loc_k=path[k], loc_k1=path[k+1];
        CGAL_assertion_code(bool OK=)
        PMP::locate_in_common_face(loc_k, loc_k1, tmesh);
        CGAL_assertion(OK);

        PMP::Face_location<TriangleMesh, typename K::FT> polygon_center;
        polygon_center.first=loc_k.first;
        double alpha = excess/len;

        for(int ii=0; ii<3;++ii)
          polygon_center.second[ii] = loc_k.second[ii]*alpha+loc_k1.second[ii]*(1.-alpha);

        std::array<Vector_2,3> flat_triangle =
          internal::init_flat_triangle<K>(halfedge(polygon_center.first,tmesh),vpm,tmesh);
        Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
        Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

        Vector_2 dir2 = tgt-src;
        double theta = atan2(dir2.y(), dir2.x());

        return std::make_pair(polygon_center, theta);
      }

      if (++k==path.size()-1)
      {
        PMP::Face_location<TriangleMesh, typename K::FT> loc_k=path[k-1], loc_k1=path[k];
        CGAL_assertion_code(bool OK=)
        PMP::locate_in_common_face(loc_k, loc_k1, tmesh);
        CGAL_assertion(OK);
        std::array<Vector_2,3> flat_triangle =
          internal::init_flat_triangle<K>(halfedge(loc_k.first,tmesh),vpm,tmesh);
        Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
        Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

        Vector_2 dir2 = tgt-src;
        double theta = atan2(dir2.y(), dir2.x());

        return std::make_pair(path.back(), theta);
      }
    }
  }

  template <class VPM>
  CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>
  operator()(std::size_t pid,
             const TriangleMesh& tmesh,
             VPM vpm,
             const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>& center,
             const FT& scaling,
             const Dual_geodesic_solver<typename K::FT>& solver,
             std::vector<Vector_2>& directions, std::vector<FT>& lens)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    PMP::Face_location<TriangleMesh, typename K::FT> polygon_center;
    double xc = (this->m_polygon_bboxes[pid].xmin()+this->m_polygon_bboxes[pid].xmax())/2.;

    double theta;
    if (xc<m_c2.x())
    {
      std::tie(polygon_center, theta)=get_polygon_center(tmesh, vpm, m_left_path, scaling * (m_c2.x()-xc));
      theta=CGAL_PI+theta;
    }
    else
      std::tie(polygon_center, theta)=get_polygon_center(tmesh, vpm, m_right_path, scaling * (xc-m_c2.x()));

    std::vector<PMP::Edge_location<TriangleMesh, typename K::FT>> shortest_path;
      locally_shortest_path<typename K::FT>(center, polygon_center, tmesh, shortest_path, solver);


    std::vector<std::pair<typename K::FT, typename K::FT>> polar_coords =
      convert_to_polar_coordinates<K>(this->m_polygons.value().get()[pid],
                                      typename K::Point_2((this->m_polygon_bboxes[pid].xmin()+this->m_polygon_bboxes[pid].xmax())/2.,
                                                          (this->m_gbox.ymin()+this->m_gbox.ymax())/2.));

    //already duplicated in trace_geodesic_polygon
    if (polar_coords.front()==polar_coords.back()) polar_coords.pop_back();

    lens.reserve(polar_coords.size());
    directions.reserve(polar_coords.size());

    for (const std::pair<double, double>& coord : polar_coords)
    {
      lens.push_back(scaling * coord.first);
      directions.emplace_back(std::cos(coord.second+theta), std::sin(coord.second+theta));
    }

    return polygon_center;
  }

  Point_2 m_c2;
  std::vector< CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT> > m_left_path, m_right_path;
};

/*!
 * \ingroup VGSFunctions
 * computes for each vertex of each polygon in `polygons` a face location on `tmesh`, where `center` represents the center of the 2D bounding box of the polygons.
 * This method starts by considering the segment splitting in two halves along the y-axis the bounding box of the polygons. 2D centers for each polygon are
 * computed on this segment as the intersection with the line splitting the bounding box of the polygon in two halves along the x-axis.
 * The splitting segment is then drawn on `tmesh` and the face location of the 2D centers is found.
 * `trace_geodesic_polygon()` is then called for each polygon and center, with appropriate directions and distances to have a consistent orientation for the polygons.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam K a model of `Kernel` with `K::FT` being a floating point number type (float or double)
 * \param center the location on `tmesh` corresponding to the center of the 2D bounding box of the polygons.
 * \param polygons 2D polygons
 * \param scaling a scaling factor to scale the polygons on `tmesh` (considering geodesic distances on `tmesh`)
 * \param tmesh input triangle mesh supporting the vertices of the output polygon
 * \param solver container for the precomputed information. If not initialized, it will be initialized internally.
 * \return a face location for each vertex of each polygon
 * \todo add named parameters
 * \todo polygon orientation is not handled in the function and should be done outside of the function for now
 * \todo for better rendering we can group polygons to have one center for the same group of polygon (useful for letters that are not simply connected)
 * \todo what if boundary is reached
 */
template <class K, class TriangleMesh>
std::vector<std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh,typename K::FT>>>
trace_geodesic_label(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT> &center,
                     const std::vector<std::vector<typename K::Point_2>>& polygons,
                     const typename K::FT scaling,
                     const TriangleMesh &tmesh,
                     const Dual_geodesic_solver<typename K::FT>& solver = {})
{
  //TODO: missing the direction...
  Placement_along_a_line<TriangleMesh, K> placement;
  return trace_geodesic_polygons<K>(center, polygons, scaling, tmesh, solver, placement);
}

template <class TriangleMesh, class K>
struct Placement_along_a_curve
  : public Default_placement<K>
{
  using Base = Default_placement<K>;
  using Vector_2 = typename K::Vector_2;
  using Point_2 = typename K::Point_2;
  using FT = typename K::FT;

  Placement_along_a_curve(const std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>>& supporting_curve,
                          const FT& padding,
                          bool is_centered)
    : m_supporting_curve(supporting_curve)
    , m_padding(padding)
    , m_is_centered(is_centered)
  {}

  template <class VPM>
  void register_polygons(const TriangleMesh& tmesh,
                         VPM vpm,
                         const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& center,
                         const FT& scaling,
                         const std::vector<std::vector<Point_2>>& polygons)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;
    Base::register_polygons(tmesh, vpm, center, scaling, polygons);

    m_c2 = Point_2( (this->m_gbox.xmin()+this->m_gbox.xmax())/2., (this->m_gbox.ymin()+this->m_gbox.ymax())/2. );
    Point_2 left_most(this->m_gbox.xmin(), m_c2.y());
    Point_2 right_most(this->m_gbox.xmax(), m_c2.y());


    // TODO: filter the path?

    m_support_lengths.resize(m_supporting_curve.size(),0);
    for (std::size_t i=0; i<m_supporting_curve.size()-1; ++i)
      m_support_lengths[i+1]=m_support_lengths[i]+std::sqrt(squared_distance(PMP::construct_point(m_supporting_curve[i], tmesh),
                                                                             PMP::construct_point(m_supporting_curve[i+1], tmesh)));

    CGAL_assertion( (m_support_lengths.back() > m_padding + scaling * (this->m_gbox.xmax()-this->m_gbox.xmin())) ||
                    (m_is_centered && (m_support_lengths.back() > scaling * (this->m_gbox.xmax()-this->m_gbox.xmin()))) );

    m_padding = m_is_centered ? (m_support_lengths.back()-(scaling*(this->m_gbox.xmax()-this->m_gbox.xmin())))/2
                              : m_padding;






    //~ double len = (this->m_gbox.xmax()-this->m_gbox.xmin())/2.;

    //~ m_left_path = straightest_geodesic<K>(center, left_most-m_c2, scaling * len, tmesh);
    //~ m_right_path = straightest_geodesic<K>(center, right_most-m_c2, scaling * len, tmesh);

    //~ CGAL_assertion(m_left_path.size() >=2);
    //~ CGAL_assertion(m_right_path.size() >=2);
  }

  template <class VPM>
  auto get_polygon_center(const TriangleMesh& tmesh,
                          VPM vpm,
                          double targetd)
  {
    namespace PMP=CGAL::Polygon_mesh_processing;
    // use left
    std::size_t k=0;
    while(true) // TODO get rid of the while and std::lower_bound
    {
      double acc = m_support_lengths[k];
      if (acc == targetd)
      {
        // note: if the supporting curve goes though m_supporting_curve[k], that vertex might be duplicated (if it's a path).
        //       But it is not an issue for the code as the 0 length segments do not contribute to the distance and
        //       points of loc_k and loc_k1 are expected to be always different.

        double theta=0;
        //TODO: should pick k-1 if k==m_supporting_curve.size()-1
        PMP::Face_location<TriangleMesh, typename K::FT> loc_k=m_supporting_curve[k], loc_k1=m_supporting_curve[k+1];
        CGAL_assertion_code(bool OK=)
        PMP::locate_in_common_face(loc_k, loc_k1, tmesh);
        CGAL_assertion(OK);
        std::array<Vector_2,3> flat_triangle =
          internal::init_flat_triangle<K>(halfedge(loc_k.first,tmesh),vpm,tmesh);
        Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
        Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

        Vector_2 dir2 = tgt-src;
        theta = atan2(dir2.y(), dir2.x());

        return std::make_pair(m_supporting_curve[k], theta);
      }

      if (acc > targetd)
      {
        double excess = acc-targetd;

        PMP::Face_location<TriangleMesh, typename K::FT> loc_k=m_supporting_curve[k-1], loc_k1=m_supporting_curve[k];
        CGAL_assertion_code(bool OK=)
        PMP::locate_in_common_face(loc_k, loc_k1, tmesh);
        CGAL_assertion(OK);

        PMP::Face_location<TriangleMesh, typename K::FT> polygon_center;
        polygon_center.first=loc_k.first;
        double alpha = excess/(m_support_lengths[k]-m_support_lengths[k-1]);

        for(int ii=0; ii<3;++ii)
          polygon_center.second[ii] = loc_k.second[ii]*alpha+loc_k1.second[ii]*(1.-alpha);

        std::array<Vector_2,3> flat_triangle =
          internal::init_flat_triangle<K>(halfedge(polygon_center.first,tmesh),vpm,tmesh);
        Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
        Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

        Vector_2 dir2 = tgt-src;
        double theta = atan2(dir2.y(), dir2.x());

        return std::make_pair(polygon_center, theta);
      }

      if (++k==m_supporting_curve.size())
      {
        //TODO: shall we throw an exception instead or return false?

        PMP::Face_location<TriangleMesh, typename K::FT> loc_k=m_supporting_curve[k-2], loc_k1=m_supporting_curve[k-1];
        CGAL_assertion_code(bool OK=)
        PMP::locate_in_common_face(loc_k, loc_k1, tmesh);
        CGAL_assertion(OK);
        std::array<Vector_2,3> flat_triangle =
          internal::init_flat_triangle<K>(halfedge(loc_k.first,tmesh),vpm,tmesh);
        Vector_2 src = loc_k.second[0]*flat_triangle[0]+loc_k.second[1]*flat_triangle[1]+loc_k.second[2]*flat_triangle[2];
        Vector_2 tgt = loc_k1.second[0]*flat_triangle[0]+loc_k1.second[1]*flat_triangle[1]+loc_k1.second[2]*flat_triangle[2];

        Vector_2 dir2 = tgt-src;
        double theta = atan2(dir2.y(), dir2.x());

        return std::make_pair(m_supporting_curve.back(), theta);
      }
    }
  }

  template <class VPM>
  CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>
  operator()(std::size_t pid,
             const TriangleMesh& tmesh,
             VPM vpm,
             const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>& /* center */,
             const FT& scaling,
             const Dual_geodesic_solver<typename K::FT>& /* solver */,
             std::vector<Vector_2>& directions, std::vector<FT>& lens)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    double xc = (this->m_polygon_bboxes[pid].xmin()+this->m_polygon_bboxes[pid].xmax())/2.;

    PMP::Face_location<TriangleMesh, FT> polygon_center;
    double theta;
    std::tie(polygon_center, theta)=get_polygon_center(tmesh, vpm, scaling * (xc-this->m_gbox.xmin())+m_padding);

    std::vector<std::pair<FT, FT>> polar_coords =
      convert_to_polar_coordinates<K>(this->m_polygons.value().get()[pid],
                                      Point_2((this->m_polygon_bboxes[pid].xmin()+this->m_polygon_bboxes[pid].xmax())/2.,
                                              (this->m_gbox.ymin()+this->m_gbox.ymax())/2.));
    //already duplicated in trace_geodesic_polygon
    if (polar_coords.front()==polar_coords.back()) polar_coords.pop_back();

    lens.reserve(polar_coords.size());
    directions.reserve(polar_coords.size());

    for (const std::pair<double, double>& coord : polar_coords)
    {
      lens.push_back(scaling * coord.first);
      directions.emplace_back(std::cos(coord.second+theta), std::sin(coord.second+theta));
    }

    return polygon_center;
  }

  Point_2 m_c2;
  std::vector< CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT> > m_left_path, m_right_path;

  const std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>>& m_supporting_curve;
  std::vector<FT> m_support_lengths;
  FT m_padding;
  bool m_is_centered;
};

/*!
 * \ingroup VGSFunctions
 * computes for each vertex of each polygon in `polygons` a face location on `tmesh` along the curve `supporting_curve`.
 * This method starts by considering the segment splitting in two halves along the y-axis the bounding box of the polygons. 2D centers for each polygon are
 * computed on this segment as the intersection with the line splitting the bounding box of the polygon in two halves along the x-axis.
 * The splitting segment is then mapped onto `supporting_curve` by first scaling it using `scaling`, and using `padding` and `is_centered`.
 * Face locations of the 2D center are found on `supporting_curve`.
 * `trace_geodesic_polygon()` is then called for each polygon and center, with appropriate directions and distances to have a consistent orientation for the polygons.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam K a model of `Kernel` with `K::FT` being a floating point number type (float or double)
 * \param supporting_curve a path on `tmesh` that will support the center of the bounding box of each polygon.
 *                         For two consecutive face locations, there must exist a face in `tmesh` containing the two corresponding points.
 * \param polygons 2D polygons
 * \param scaling a scaling factor to scale the polygons on `tmesh` (considering geodesic distances on `tmesh`)
 * \param padding padding applied at the beginning of supporting curve to start the drawing
 * \param is_centered is `true`, `padding` is ignored and the bounding box of the polygon is centered on the supporting curve (in 1D)
 * \param tmesh input triangle mesh supporting the vertices of the output polygon
 * \param solver container for the precomputed information. If not initialized, it will be initialized internally.
 * \return a face location for each vertex of each polygon
 * \todo add named parameters
 * \todo polygon orientation is not handled in the function and should be done outside of the function for now
 * \todo for better rendering we can group polygons to have one center for the same group of polygon (useful for letters that are not simply connected)
 * \todo check padding is ignored if is_centered is used + update the doc if not
 * \todo doc what happens if supporting curve is not long enough + boundary reached
 */
template <class K, class TriangleMesh>
std::vector<std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>>>
trace_geodesic_label_along_curve(const std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>>& supporting_curve,
                                 const std::vector<std::vector<typename K::Point_2>>& polygons,
                                 const typename K::FT scaling,
                                 const typename K::FT padding,
                                 const bool is_centered,
                                 const TriangleMesh &tmesh,
                                 const Dual_geodesic_solver<typename K::FT>& solver = {})
{
  Placement_along_a_curve<TriangleMesh, K> placement(supporting_curve, padding, is_centered);
  auto center = supporting_curve.front(); // arbitrary (not used)
  return trace_geodesic_polygons<K>(center, polygons, scaling, tmesh, solver, placement);
}

/*!
 * \ingroup VGSFunctions
 * computes a path representing a Bézier curve defined by four control points.
 * Control points are defined by the endpoints of straightest geodesic curves starting from `center` along given directions and distances.
 * The iterative de Casteljau subdivision algorithm is applied to create more control points that are then connected with locally shortest paths.
 * The output path is such that for two consecutive face locations, there must exist a face in `tmesh` containing the two corresponding points.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam K a model of `Kernel` with `K::FT` being a floating point number type (float or double)
 * \param center the location on `tmesh` where straightest geodesic for the placement of control points starts. The y-axis used is `halfedge(center.first, tmesh)`.
 * \param directions contains the direction of the straightest geodesic for each control point
 * \param lengths contains the length of the straightest geodesic for each control point
 * \param num_subdiv the number of iterations of the de Casteljau subdivision algorithm
 * \param tmesh input triangle mesh supporting the vertices of the output polygon
 * \param solver container for the precomputed information. If not initialized, it will be initialized internally.
 * \return a face location for each vertex of each polygon
 * \todo add named parameters
 * \todo polygon orientation is not handled in the function and should be done outside of the function for now
 * \todo for better rendering we can group polygons to have one center for the same group of polygon (useful for letters that are not simply connected)
 * \todo what if boundary is reached
 */
template <class K, class TriangleMesh>
std::vector< std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>> >
trace_Bezier_curves(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT> &center,
                    const std::vector<std::array<typename K::Vector_2, 4>>& directions,
                    const std::vector<std::array<typename K::FT, 4>>& lengths,
                    const int num_subdiv,
                    const TriangleMesh &tmesh
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                    , const Dual_geodesic_solver<typename K::FT>& solver = {}
#endif
)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using FT = typename K::FT;

  std::size_t n=directions.size();
  std::vector< std::vector<PMP::Face_location<TriangleMesh, typename K::FT>> > result(n);

#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }
#endif

#ifdef CGAL_DEBUG_BSURF
  std::ofstream debug_cp("/tmp/control_points.xyz");
  std::ofstream debug_ep("/tmp/end_points.xyz");
  debug_cp << std::setprecision(17);
  debug_ep << std::setprecision(17);
#endif
  for (std::size_t i=0; i<n; ++i)
  {
    Bezier_segment<TriangleMesh, FT> control_loc;
    for (int k=0;k<4; ++k)
    {
      control_loc[k] = straightest_geodesic<K>(center,directions[i][k],lengths[i][k],tmesh).back();
    }

    #ifdef CGAL_DEBUG_BSURF
      debug_ep << PMP::construct_point(control_loc[0], tmesh) << "\n";
      debug_ep << PMP::construct_point(control_loc[3], tmesh) << "\n";
      debug_cp << PMP::construct_point(control_loc[1], tmesh) << "\n";
      debug_cp << PMP::construct_point(control_loc[2], tmesh) << "\n";
    #endif

    std::vector<PMP::Face_location<TriangleMesh, FT>> bezier =
      recursive_de_Casteljau(tmesh, control_loc, num_subdiv
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                            , *solver_ptr
#endif
                            );

    result[i].reserve(bezier.size());
    result[i].push_back(bezier[0]);
    for(std::size_t b=0; b<bezier.size()-1; ++b)
    {
      const PMP::Face_location<TriangleMesh, FT>& loc = bezier[b];
      const PMP::Face_location<TriangleMesh, FT>& loc1 = bezier[b+1];

      // connect the two face locations with shortest path is they are in different faces
      if (loc.first!=loc1.first)
      {
        std::vector<PMP::Edge_location<TriangleMesh, FT>> edge_locations;
        locally_shortest_path<FT>(loc, loc1, tmesh, edge_locations, solver);
        result[i].reserve(result[i].size()+edge_locations.size());
        for (const PMP::Edge_location<TriangleMesh, FT>& e : edge_locations)
          result[i].push_back(PMP::to_face_location(e, tmesh));
      }
      result[i].push_back(loc1);
    }
  }

  return result;
}


/*!
 * \ingroup VGSFunctions
 * computes a path representing a Bézier polyline (a sequence of Bézier curves having a common control points,
 * that is the fourth control point of the nth curve is the first control point of the (n+1)th curve).
 * Control points are defined by the endpoints of straightest geodesic curves starting from `center` along given directions and distances.
 * The iterative de Casteljau subdivision algorithm is applied to create more control points that are then connected with locally shortest paths.
 * The output path is such that for two consecutive face locations, there must exist a face in `tmesh` containing the two corresponding points.
 * The first portion of the curve is defined by the four first values in `directions` and `lengths`. The second portion is defined by the fourth value
 * and the next three, and so on until the end is reached. If `is_closed` is true, then the first value will be used with the last three to define the last portion.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam K a model of `Kernel` with `K::FT` being a floating point number type (float or double)
 * \param center the location on `tmesh` where straightest geodesic for the placement of control points starts. The y-axis used is `halfedge(center.first, tmesh)`.
 * \param directions contains the direction of the straightest geodesic for each control point
 * \param lengths contains the length of the straightest geodesic for each control point
 * \param num_subdiv the number of iterations of the de Casteljau subdivision algorithm
 * \param is_closed if true [directions/lengths].front() will be used as additional last point, generating a closed path
 * \param tmesh input triangle mesh supporting the vertices of the output polygon
 * \param solver container for the precomputed information. If not initialized, it will be initialized internally.
 * \return a face location for each vertex of each polygon
 * \todo add named parameters
 * \todo polygon orientation is not handled in the function and should be done outside of the function for now
 * \todo for better rendering we can group polygons to have one center for the same group of polygon (useful for letters that are not simply connected)
 * \todo what if boundary is reached
 */
//TODO: make sure it is consistent with the rest to not duplicate the last point if closed
template <class K, class TriangleMesh>
std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>>
trace_polyline_of_Bezier_curves(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT> &center,
                                const std::vector<typename K::Vector_2>& directions,
                                const std::vector<typename K::FT>& lengths,
                                bool is_closed, // use [directions/lengths].front as last control point?
                                const int num_subdiv,
                                const TriangleMesh &tmesh
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                                , const Dual_geodesic_solver<typename K::FT>& solver = {}
#endif
)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using FT = typename K::FT;

  std::size_t n = (directions.size() - (is_closed?0:1))/3;
  CGAL_assertion( n * 3 + (is_closed?0:1) == directions.size() );

  std::vector<PMP::Face_location<TriangleMesh, typename K::FT>> result;

  // n is the number of quadruple of control points
  // After num_subdiv steps, we have 2^num_subdiv * n quadruples of control points


  // even if closed we will duplicate the last point
  // (this is a lower bound without taking into account shortest path between points)
  result.reserve( (1<<num_subdiv) * 3 + 2 );

#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }
#endif

#ifdef CGAL_DEBUG_BSURF
  std::ofstream debug_cp("/tmp/control_points.xyz");
  std::ofstream debug_ep("/tmp/end_points.xyz");
  debug_cp << std::setprecision(17);
  debug_ep << std::setprecision(17);
#endif

  PMP::Face_location<TriangleMesh, FT> prev_loc = straightest_geodesic<K>(center,directions[0],lengths[0],tmesh).back(),
                                       first_loc = prev_loc;

  for (std::size_t i=0; i<n; ++i)
  {
    Bezier_segment<TriangleMesh, FT> control_loc;
    control_loc[0]=prev_loc;
    for (int k=1;k<4; ++k)
    {
      if (k!=3 || !is_closed || 3*i+k!=directions.size())
        control_loc[k] = straightest_geodesic<K>(center,directions[3*i+k],lengths[3*i+k],tmesh).back();
      else
        control_loc[k] = first_loc;
    }
    prev_loc=control_loc[3];

    #ifdef CGAL_DEBUG_BSURF
      debug_ep << PMP::construct_point(control_loc[0], tmesh) << "\n";
      debug_ep << PMP::construct_point(control_loc[3], tmesh) << "\n";
      debug_cp << PMP::construct_point(control_loc[1], tmesh) << "\n";
      debug_cp << PMP::construct_point(control_loc[2], tmesh) << "\n";
    #endif

    std::vector<PMP::Face_location<TriangleMesh, FT>> bezier =
      recursive_de_Casteljau(tmesh, control_loc, num_subdiv
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                            , *solver_ptr
#endif
                            );

    if (i==0)
      result.push_back(bezier[0]);
    for(std::size_t b=0; b<bezier.size()-1; ++b)
    {
      const PMP::Face_location<TriangleMesh, FT>& loc = bezier[b];
      const PMP::Face_location<TriangleMesh, FT>& loc1 = bezier[b+1];

      // connect the two face locations with shortest path is they are in different faces
      if (loc.first!=loc1.first)
      {
        std::vector<PMP::Edge_location<TriangleMesh, FT>> edge_locations;
        locally_shortest_path<FT>(loc, loc1, tmesh, edge_locations, solver);
        for (const PMP::Edge_location<TriangleMesh, FT>& e : edge_locations)
          result.push_back(PMP::to_face_location(e, tmesh));
      }
      result.push_back(loc1);
    }
  }

  return result;
}

// template <class K, class TriangleMesh>
// std::vector<typename K::Vector_2>
// tangent_path_direction(const std::vector<CGAL::Polygon_mesh_processing::Edge_location<TriangleMesh,typename K::FT>>& path,
//                        const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>& src,
//                        const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>& tgt,
//                        const TriangleMesh &tmesh,const bool initial=true)
// {
//   auto find = [](const std::array<typename K::FT,3> &vec, int x) {
//     for (int i = 0; i < vec.size(); i++)
//       if (vec[i] == x)
//         return i;
//     return -1;
//   };
//  typename K::Vector_2 direction;
//   using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
//   VPM vpm = get(CGAL::vertex_point, tmesh);
//
//  if(initial)
//  {
//
//     halfedge_descriptor h_ref=halfedge(src.first,mesh);
//     std::array<Vector_2, 3> flat_tid=internal::init_flat_triangle<K>(h_ref,vpm,tmesh);
//     if(path.size()==0)
//     {
//       assert(src.first==tgt.first);//TODO:src and tgt may have different faces because we do not update them when cleaning the strip
//       typename K::Vector_2 flat_src=src.second[0]*flat_tid[0]+src.second[1]*flat_tid[1]+src.second[2]*flat_tid[2];
//       typename K::Vector_2 flat_tgt=tgt.second[0]*flat_tid[0]+tgt.second[1]*flat_tid[1]+tgt.second[2]*flat_tid[2];
//       direction=normalize(flat_tgt-flat_src);
//     }else{
//     halfedge_descriptor h_edge=halfedge(path[0].first,tmesh);
//     int k=0;

//     if(h_edge==prev(h_ref,mesh))
//       k=2;
//     else if(h_edge==next(h_ref,mesh))
//       k=1;
//     else
//       assert(h_edge==h_ref);
//     }
//  }
// }

template <class K, class TriangleMesh>
std::vector<typename K::Point_3>
trace_agap_polygon(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT> &center,
                     const std::vector<typename K::Vector_2>& directions,
                     const std::vector<typename K::FT>& lengths,
                     const TriangleMesh &tmesh,
                     const Dual_geodesic_solver<typename K::FT>& solver = {})
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  size_t n=directions.size();
  std::vector<typename K::Point_3> result;
  std::vector<PMP::Face_location<TriangleMesh, typename K::FT>> vertices(n);
  for(std::size_t i=0;i<n;++i)
    vertices[i]= straightest_geodesic<K>(center,directions[i],lengths[i],tmesh).back();

  std::vector<PMP::Edge_location<TriangleMesh,typename K::FT>> edge_locations;

  const Dual_geodesic_solver<typename K::FT>* solver_ptr=&solver;
  Dual_geodesic_solver<typename K::FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, tmesh);
  }

  for(std::size_t i=0;i<n;++i)
  {
    edge_locations.clear();
    locally_shortest_path<typename K::FT>(vertices[i],vertices[(i+1)%n],tmesh, edge_locations, *solver_ptr);
    result.push_back(PMP::construct_point(vertices[i],tmesh));
    for(auto& el : edge_locations)
    {
      result.push_back(PMP::construct_point(el, tmesh));
    }
  }
  result.push_back(PMP::construct_point(vertices[0],tmesh));

  return result;

}

} // namespace Vector_graphics_on_surfaces
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_BSURF_TRACE_H
