// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Hossam Saeed
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H

#include <CGAL/license/Polygon_mesh_processing/interpolated_corrected_curvature_measures.h>

#include <CGAL/assertions.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <Eigen/Eigenvalues>

#include <numeric>
#include <queue>
#include <unordered_set>
#include <functional>

#define EXPANDING_RADIUS_EPSILON 1e-6

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template<typename PolygonMesh, typename GT>
typename GT::FT average_edge_length(const PolygonMesh& pmesh)
{
    const std::size_t n = edges(pmesh).size();
    if (n == 0)
        return 0;

    typename GT::FT avg_edge_length = 0;
    for (auto e : edges(pmesh))
        avg_edge_length += edge_length(e, pmesh);

    avg_edge_length /= n;
    return avg_edge_length;
}

enum Curvature_measure_index {
    MU0_AREA_MEASURE, ///< corrected area density
    MU1_MEAN_CURVATURE_MEASURE, ///< corrected mean curvature density
    MU2_GAUSSIAN_CURVATURE_MEASURE ///< corrected gaussian curvature density
};

template<typename GT>
typename GT::FT interpolated_corrected_area_measure_face(const std::vector<typename GT::Vector_3>& u,
                                                         const std::vector<typename GT::Vector_3>& x = {})
{
    const std::size_t n = x.size();
    CGAL_precondition(u.size() == n);
    CGAL_precondition(n >= 3);

    typename GT::Construct_cross_product_vector_3 cross_product;

    // Triangle: use triangle formula
    if (n == 3)
    {
        const typename GT::Vector_3 um = (u[0] + u[1] + u[2]) / 3.0;
        return 0.5 * um * cross_product(x[1] - x[0], x[2] - x[0]);
    }
    // Quad: use bilinear interpolation formula
    else if (n == 4)
    {
        // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
        // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3
        return (1.0 / 36.0) * (
              (4 * u[0] + 2 * u[1] + 2 * u[3] + u[2]) * cross_product(x[1] - x[0], x[3] - x[0])
            + (2 * u[0] + 4 * u[1] + u[3] + 2 * u[2]) * cross_product(x[1] - x[0], x[2] - x[1])
            + (2 * u[0] + u[1] + 4 * u[3] + 2 * u[2]) * cross_product(x[2] - x[3], x[3] - x[0])
            + (u[0] + 2 * u[1] + 2 * u[3] + 4 * u[2]) * cross_product(x[2] - x[3], x[2] - x[1])
            );
    }
    // N-gon: split into n triangles by polygon center and use triangle formula for each
    else
    {
        typename GT::FT mu0 = 0;

        // getting center of points
        typename GT::Vector_3 xc =
            std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
        xc /= n;

        // getting unit average normal of points
        typename GT::Vector_3 uc =
            std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
        uc /= sqrt(uc * uc);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            mu0 += interpolated_corrected_area_measure_face<GT>(
                std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc},
                std::vector<typename GT::Vector_3> {x[i], x[(i + 1) % n], xc}
            );
        }
        return mu0;
    }
}

template<typename GT>
typename GT::FT interpolated_corrected_mean_curvature_measure_face(const std::vector<typename GT::Vector_3>& u,
                                                                   const std::vector<typename GT::Vector_3>& x = {})
{
    const std::size_t n = x.size();
    CGAL_precondition(u.size() == n);
    CGAL_precondition(n >= 3);

    typename GT::Construct_cross_product_vector_3 cross_product;

    // Triangle: use triangle formula
    if (n == 3)
    {
        const typename GT::Vector_3 um = (u[0] + u[1] + u[2]) / 3.0;

        return 0.5 * um * (cross_product(u[2] - u[1], x[0])
            + cross_product(u[0] - u[2], x[1])
            + cross_product(u[1] - u[0], x[2]));
    }
    // Quad: use bilinear interpolation formula
    else if (n == 4)
    {
        // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
        // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3

        const typename GT::Vector_3 u02 = u[2] - u[0];
        const typename GT::Vector_3 u13 = u[3] - u[1];
        const typename GT::Vector_3 x0_cross = cross_product(u13, x[0]);
        const typename GT::Vector_3 x1_cross = -cross_product(u02, x[1]);
        const typename GT::Vector_3 x3_cross = cross_product(u02, x[3]);
        const typename GT::Vector_3 x2_cross = -cross_product(u13, x[2]);

        return (1.0 / 12.0) * (
              u[0] * (2 * x0_cross - cross_product((u[3] + u[2]), x[1]) + cross_product((u[1] + u[2]), x[3]) + x2_cross)
            + u[1] * (cross_product((u[3] + u[2]), x[0]) + 2 * x1_cross + x3_cross - cross_product((u[0] + u[3]), x[2]))
            + u[3] * (-cross_product((u[1] + u[2]), x[0]) + x1_cross + 2 * x3_cross + cross_product((u[0] + u[1]), x[2]))
            + u[2] * (x0_cross + cross_product((u[0] + u[3]), x[1]) - cross_product((u[0] + u[1]), x[3]) + 2 * x2_cross)
            );
    }
    // N-gon: split into n triangles by polygon center and use triangle formula for each
    else
    {
        typename GT::FT mu1 = 0;

        // getting center of points
        typename GT::Vector_3 xc =
            std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
        xc /= n;

        // getting unit average normal of points
        typename GT::Vector_3 uc =
            std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
        uc /= sqrt(uc * uc);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            mu1 += interpolated_corrected_mean_curvature_measure_face<GT>(
                std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc},
                std::vector<typename GT::Vector_3> {x[i], x[(i + 1) % n], xc}
            );
        }
        return mu1;
    }
}

template<typename GT>
typename GT::FT interpolated_corrected_gaussian_curvature_measure_face(const std::vector<typename GT::Vector_3>& u,
                                                                       const std::vector<typename GT::Vector_3>& x = {})
{
    const std::size_t n = u.size();
    CGAL_precondition(n >= 3);

    typename GT::Construct_cross_product_vector_3 cross_product;

    // Triangle: use triangle formula
    if (n == 3)
    {
        return 0.5 * u[0] * cross_product(u[1], u[2]);
    }
    // Quad: use bilinear interpolation formula
    else if (n == 4)
    {
        // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
        // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3
        return (1.0 / 36.0) * (
              (4 * u[0] + 2 * u[1] + 2 * u[3] + u[2]) * cross_product(u[1] - u[0], u[3] - u[0])
            + (2 * u[0] + 4 * u[1] + u[3] + 2 * u[2]) * cross_product(u[1] - u[0], u[2] - u[1])
            + (2 * u[0] + u[1] + 4 * u[3] + 2 * u[2]) * cross_product(u[2] - u[3], u[3] - u[0])
            + (u[0] + 2 * u[1] + 2 * u[3] + 4 * u[2]) * cross_product(u[2] - u[3], u[2] - u[1])
            );
    }
    // N-gon: split into n triangles by polygon center and use triangle formula for each
    else
    {
        typename GT::FT mu2 = 0;

        // getting unit average normal of points
        typename GT::Vector_3 uc =
            std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
        uc /= sqrt(uc * uc);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            mu2 += interpolated_corrected_gaussian_curvature_measure_face<GT>(
                std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc}
            );
        }
        return mu2;
    }
}

template<typename GT>
std::array<typename GT::FT, 3 * 3> interpolated_corrected_anisotropic_measure_face(const std::vector<typename GT::Vector_3>& u,
                                                                                   const std::vector<typename GT::Vector_3>& x)
{
    const std::size_t n = x.size();
    CGAL_precondition(u.size() == n);
    CGAL_precondition(n >= 3);

    typename GT::Construct_cross_product_vector_3 cross_product;
    std::array<typename GT::FT, 3 * 3> muXY {0};

    // Triangle: use triangle formula
    if (n == 3)
    {
        const typename GT::Vector_3 u01 = u[1] - u[0];
        const typename GT::Vector_3 u02 = u[2] - u[0];
        const typename GT::Vector_3 x01 = x[1] - x[0];
        const typename GT::Vector_3 x02 = x[2] - x[0];
        const typename GT::Vector_3 um = (u[0] + u[1] + u[2]) / 3.0;

        for (std::size_t ix = 0; ix < 3; ix++)
        {
            typename GT::Vector_3 X;
            if (ix == 0)
                X = typename GT::Vector_3(1, 0, 0);
            if (ix == 1)
                X = typename GT::Vector_3(0, 1, 0);
            if (ix == 2)
                X = typename GT::Vector_3(0, 0, 1);

            for (std::size_t iy = 0; iy < 3; iy++)
                muXY[ix * 3 + iy] = 0.5 * um * (cross_product(u02[iy] * X, x01) - cross_product(u01[iy] * X, x02));
        }
    }
    // Quad: use bilinear interpolation formula
    else if (n == 4)
    {
        // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
        // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3
        for (std::size_t ix = 0; ix < 3; ix++)
        {
            typename GT::Vector_3 X;
            if (ix == 0)
                X = typename GT::Vector_3(1, 0, 0);
            if (ix == 1)
                X = typename GT::Vector_3(0, 1, 0);
            if (ix == 2)
                X = typename GT::Vector_3(0, 0, 1);

            const typename GT::Vector_3 u0xX = cross_product(u[0], X);
            const typename GT::Vector_3 u1xX = cross_product(u[1], X);
            const typename GT::Vector_3 u2xX = cross_product(u[2], X);
            const typename GT::Vector_3 u3xX = cross_product(u[3], X);

            for (std::size_t iy = 0; iy < 3; iy++)
                muXY[ix * 3 + iy] = (1.0 / 72.0) * (

                    u[0][iy] * (    u0xX * (    - x[0] - 11 * x[1] + 13 * x[3] -      x[2])
                                  + u1xX * ( -5 * x[0] -  7 * x[1] + 11 * x[3] +      x[2])
                                  + u3xX * (      x[0] -  7 * x[1] + 11 * x[3] -  5 * x[2])
                                  + u2xX * (    - x[0] -  5 * x[1] +  7 * x[3] -      x[2])
                      )
                  + u[1][iy] * (    u0xX * ( 13 * x[0] -      x[1] -  7 * x[3] -  5 * x[2])
                                  + u1xX * ( 17 * x[0] -  5 * x[1] -  5 * x[3] -  7 * x[2])
                                  + u3xX * (  5 * x[0] +      x[1] +      x[3] -  7 * x[2])
                                  + u2xX * (  7 * x[0] -      x[1] +  5 * x[3] - 11 * x[2])
                      )
                  + u[2][iy] * (    u0xX * (-11 * x[0] +  5 * x[1] -      x[3] +  7 * x[2])
                                  + u1xX * (- 7 * x[0] +      x[1] +      x[3] +  5 * x[2])
                                  + u3xX * (- 7 * x[0] -  5 * x[1] -  5 * x[3] + 17 * x[2])
                                  + u2xX * (- 5 * x[0] -  7 * x[1] -      x[3] + 13 * x[2])
                      )
                  + u[3][iy] * (    u0xX * (-     x[0] +  7 * x[1] -  5 * x[3] -      x[2])
                                  + u1xX * (- 5 * x[0] + 11 * x[1] -  7 * x[3] +      x[2])
                                  + u3xX * (      x[0] + 11 * x[1] -  7 * x[3] -  5 * x[2])
                                  + u2xX * (-     x[0] + 13 * x[1] - 11 * x[3] -      x[2])
                      )

                  );
        }
    }
    // N-gon: split into n triangles by polygon center and use triangle formula for each
    else
    {
        // getting center of points
        typename GT::Vector_3 xc =
            std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
        xc /= n;

        // getting unit average normal of points
        typename GT::Vector_3 uc =
            std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
        uc /= sqrt(uc * uc);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            std::array<typename GT::FT, 3 * 3> muXY_curr_triangle =
                interpolated_corrected_anisotropic_measure_face<GT>(
                    std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc},
                    std::vector<typename GT::Vector_3> {x[i], x[(i + 1) % n], xc}
                );

            for (std::size_t ix = 0; ix < 3; ix++)
                for (std::size_t iy = 0; iy < 3; iy++)
                    muXY[ix * 3 + iy] += muXY_curr_triangle[ix * 3 + iy];
        }
    }
    return muXY;
}

template<typename PolygonMesh, typename FaceMeasureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void
    interpolated_corrected_measure_mesh(const PolygonMesh& pmesh,
                                        FaceMeasureMap fmm,
                                        const Curvature_measure_index mu_i,
                                        const NamedParameters& np = parameters::default_values())
{

    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

    typedef dynamic_vertex_property_t<typename GT::Vector_3> Vector_map_tag;
    typedef typename boost::property_map<PolygonMesh, Vector_map_tag>::const_type Default_vector_map;
    typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
        NamedParameters,
        Default_vector_map>::type       VNM;

    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
            get_const_property_map(CGAL::vertex_point, pmesh));

    VNM vnm = choose_parameter(get_parameter(np, internal_np::vertex_normal_map),
                               get(Vector_map_tag(), pmesh));

    if (is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value)
        compute_vertex_normals(pmesh, vnm, np);

    std::function
        <typename GT::FT(const std::vector<typename GT::Vector_3>&, const std::vector<typename GT::Vector_3>&)>
        iccm_function;
    switch (mu_i)
    {
    case MU0_AREA_MEASURE:
        iccm_function = &interpolated_corrected_area_measure_face<GT>;
        break;
    case MU1_MEAN_CURVATURE_MEASURE:
        iccm_function = &interpolated_corrected_mean_curvature_measure_face<GT>;
        break;
    case MU2_GAUSSIAN_CURVATURE_MEASURE:
        iccm_function = &interpolated_corrected_gaussian_curvature_measure_face<GT>;
        break;
    }

    std::vector<typename GT::Vector_3> x;
    std::vector<typename GT::Vector_3> u;

    for (face_descriptor f : faces(pmesh))
    {

        for (vertex_descriptor v : vertices_around_face(halfedge(f, pmesh), pmesh))
        {
            typename GT::Point_3 p = get(vpm, v);
            x.push_back(typename GT::Vector_3(p.x(), p.y(), p.z()));
            u.push_back(get(vnm, v));
        }

        put(fmm, f, iccm_function(u, x));
        x.clear();
        u.clear();
    }
}

template<typename PolygonMesh, typename FaceMeasureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void
    interpolated_corrected_anisotropic_measure_mesh(const PolygonMesh& pmesh,
        FaceMeasureMap fmm,
        const NamedParameters& np = parameters::default_values())
{

    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

    typedef dynamic_vertex_property_t<typename GT::Vector_3> Vector_map_tag;
    typedef typename boost::property_map<PolygonMesh, Vector_map_tag>::const_type Default_vector_map;
    typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
        NamedParameters,
        Default_vector_map>::type       VNM;

    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
            get_const_property_map(CGAL::vertex_point, pmesh));

    VNM vnm = choose_parameter(get_parameter(np, internal_np::vertex_normal_map),
        get(Vector_map_tag(), pmesh));

    if (is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value)
        compute_vertex_normals(pmesh, vnm, np);

    std::vector<typename GT::Vector_3> x;
    std::vector<typename GT::Vector_3> u;

    for (face_descriptor f : faces(pmesh))
    {

        for (vertex_descriptor v : vertices_around_face(halfedge(f, pmesh), pmesh))
        {
            typename GT::Point_3 p = get(vpm, v);
            x.push_back(typename GT::Vector_3(p.x(), p.y(), p.z()));
            u.push_back(get(vnm, v));
        }

        put(fmm, f, interpolated_corrected_anisotropic_measure_face<GT>(u, x));
        x.clear();
        u.clear();

    }
}


//
//template<typename GT>
//typename GT::FT triangle_in_ball_ratio_1(const typename GT::Vector_3 x1,
//                                         const typename GT::Vector_3 x2,
//                                         const typename GT::Vector_3 x3,
//                                         const typename GT::FT r,
//                                         const typename GT::Vector_3 c,
//                                         const std::size_t res = 3)
//{
//    const typename GT::FT R = r * r;
//    const typename GT::FT acc = 1.0 / res;
//    std::size_t samples_in = 0;
//    for (GT::FT alpha = acc / 3; alpha < 1; alpha += acc)
//        for (GT::FT beta = acc / 3; beta < 1 - alpha; beta += acc)
//        {
//            if ((alpha * x1 + beta * x2 + (1 - alpha - beta) * x3 - c).squared_length() < R)
//                samples_in++;
//        }
//    return samples_in / (typename GT::FT)(res * (res + 1) / 2);
//}

template<typename GT>
typename GT::FT face_in_ball_ratio(const std::vector<typename GT::Vector_3>& x,
                                   const typename GT::FT r,
                                   const typename GT::Vector_3 c)
{
    const std::size_t n = x.size();

    // getting center of points
    typename GT::Vector_3 xm =
        std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
    xm /= n;

    typename GT::FT d_min = (xm - c).squared_length();
    typename GT::FT d_max = d_min;

    for (const typename GT::Vector_3 xi : x)
    {
        const typename GT::FT d_sq = (xi - c).squared_length();
        d_max = (std::max)(d_sq, d_max);
        d_min = (std::min)(d_sq, d_min);
    }

    if (d_max <= r * r) return 1.0;
    else if (r * r <= d_min) return 0.0;

    d_max = sqrt(d_max);
    d_min = sqrt(d_min);

    return (r - d_min) / (d_max - d_min);
}

template<typename PolygonMesh, typename FaceMeasureMap, typename VertexMeasureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void expand_interpolated_corrected_measure_vertex(const PolygonMesh& pmesh,
        FaceMeasureMap area_fmm,
        FaceMeasureMap curvature_fmm,
        VertexMeasureMap area_vmm,
        VertexMeasureMap curvature_vmm,
        const typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
        const NamedParameters& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

    const typename GT::FT
        r = choose_parameter(get_parameter(np, internal_np::ball_radius), 0);

    typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
            get_const_property_map(CGAL::vertex_point, pmesh));


    std::queue<face_descriptor> bfs_queue;
    std::unordered_set<face_descriptor> bfs_visited;

    typename GT::Point_3 vp = get(vpm, v);
    typename GT::Vector_3 c = typename GT::Vector_3(vp.x(), vp.y(), vp.z());

    typename GT::FT corrected_mu0 = 0;
    typename GT::FT corrected_mui = 0;

    for (face_descriptor f : faces_around_target(halfedge(v, pmesh), pmesh)) {
        if (f != boost::graph_traits<PolygonMesh>::null_face())
        {
            bfs_queue.push(f);
            bfs_visited.insert(f);
        }
    }
    while (!bfs_queue.empty()) {
        face_descriptor fi = bfs_queue.front();
        bfs_queue.pop();

        // looping over vertices in face to get point coordinates
        std::vector<typename GT::Vector_3> x;
        for (vertex_descriptor vi : vertices_around_face(halfedge(fi, pmesh), pmesh))
        {
            typename GT::Point_3 pi = get(vpm, vi);
            x.push_back(typename GT::Vector_3(pi.x(), pi.y(), pi.z()));
        }

        const typename GT::FT f_ratio = face_in_ball_ratio<GT>(x, r, c);

        if (f_ratio != 0.0)
        {
            corrected_mu0 += f_ratio * get(area_fmm, fi);
            corrected_mui += f_ratio * get(curvature_fmm, fi);
            for (face_descriptor fj : faces_around_face(halfedge(fi, pmesh), pmesh))
            {
                if (bfs_visited.find(fj) == bfs_visited.end() && fj != boost::graph_traits<PolygonMesh>::null_face())
                {
                    bfs_queue.push(fj);
                    bfs_visited.insert(fj);
                }
            }
        }
    }
    put(area_vmm, v, corrected_mu0);
    put(curvature_vmm, v, corrected_mui);
}

template<typename PolygonMesh, typename AreaFaceMeasureMap, typename AnisotropicFaceMeasureMap,
    typename AreaVertexMeasureMap, typename AnisotropicVertexMeasureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void expand_interpolated_corrected_anisotropic_measure_vertex(const PolygonMesh& pmesh,
        AreaFaceMeasureMap area_fmm,
        AnisotropicFaceMeasureMap aniso_fmm,
        AreaVertexMeasureMap area_vmm,
        AnisotropicVertexMeasureMap aniso_vmm,
        const typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
        const NamedParameters& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

    const typename GT::FT
        r = choose_parameter(get_parameter(np, internal_np::ball_radius), 0);

    typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
            get_const_property_map(CGAL::vertex_point, pmesh));


    std::queue<face_descriptor> bfs_queue;
    std::unordered_set<face_descriptor> bfs_visited;

    typename GT::Point_3 vp = get(vpm, v);
    typename GT::Vector_3 c = typename GT::Vector_3(vp.x(), vp.y(), vp.z());

    typename GT::FT corrected_mu0 = 0;
    Eigen::Matrix<typename GT::FT, 3, 3> corrected_muXY = Eigen::Matrix<typename GT::FT, 3, 3>::Zero();

    for (face_descriptor f : faces_around_target(halfedge(v, pmesh), pmesh)) {
        if (f != boost::graph_traits<PolygonMesh>::null_face())
        {
            bfs_queue.push(f);
            bfs_visited.insert(f);
        }
    }
    while (!bfs_queue.empty()) {
        face_descriptor fi = bfs_queue.front();
        bfs_queue.pop();

        // looping over vertices in face to get point coordinates
        std::vector<typename GT::Vector_3> x;
        for (vertex_descriptor vi : vertices_around_face(halfedge(fi, pmesh), pmesh))
        {
            typename GT::Point_3 pi = get(vpm, vi);
            x.push_back(typename GT::Vector_3(pi.x(), pi.y(), pi.z()));
        }

        const typename GT::FT f_ratio = face_in_ball_ratio<GT>(x, r, c);

        if (f_ratio != 0.0)
        {
            corrected_mu0 += f_ratio * get(area_fmm, fi);

            std::array<typename GT::FT, 3 * 3> muXY_face = get(aniso_fmm, fi);

            for (std::size_t ix = 0; ix < 3; ix++)
                for (std::size_t iy = 0; iy < 3; iy++)
                    corrected_muXY(ix, iy) += f_ratio * muXY_face[ix * 3 + iy];

            for (face_descriptor fj : faces_around_face(halfedge(fi, pmesh), pmesh))
            {
                if (bfs_visited.find(fj) == bfs_visited.end() && fj != boost::graph_traits<PolygonMesh>::null_face())
                {
                    bfs_queue.push(fj);
                    bfs_visited.insert(fj);
                }
            }
        }
    }
    put(area_vmm, v, corrected_mu0);
    put(aniso_vmm, v, corrected_muXY);
}

} // namespace internal

/**
* \ingroup PMP_corrected_curvatures_grp
*
* Computes the interpolated corrected mean curvature across the mesh
* and stores it in a veretex property map `vcm`.
*
* @tparam PolygonMesh a model of `FaceListGraph`.
* @tparam VertexCurvatureMap model of `WritablePropertyMap` with
* `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and `GT::FT` as value type.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param pmesh the polygon mesh.
* @param vcm the vertex property map in which the computed mean curvatures are stored.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters".
*
* @see `interpolated_corrected_gaussian_curvature()`
* @see `interpolated_corrected_principal_curvatures()`
*/

template<typename PolygonMesh, typename VertexCurvatureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void interpolated_corrected_mean_curvature(const PolygonMesh& pmesh,
        VertexCurvatureMap vcm,
        const NamedParameters& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::property_map<PolygonMesh,
        CGAL::dynamic_face_property_t<typename GT::FT>>::const_type FaceMeasureMap;

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::property_map<PolygonMesh,
        CGAL::dynamic_vertex_property_t<typename GT::FT>>::const_type VertexMeasureMap;

    typename GT::FT
        r = choose_parameter(get_parameter(np, internal_np::ball_radius), 0);

    if (r == 0)
        r = internal::average_edge_length<PolygonMesh, GT>(pmesh) * EXPANDING_RADIUS_EPSILON;

    FaceMeasureMap mu0_map = get(CGAL::dynamic_face_property_t<typename GT::FT>(), pmesh);
    FaceMeasureMap mu1_map = get(CGAL::dynamic_face_property_t<typename GT::FT>(), pmesh);

    VertexMeasureMap mu0_expand_map = get(CGAL::dynamic_vertex_property_t<typename GT::FT>(), pmesh);
    VertexMeasureMap mu1_expand_map = get(CGAL::dynamic_vertex_property_t<typename GT::FT>(), pmesh);

    internal::interpolated_corrected_measure_mesh(pmesh, mu0_map, internal::MU0_AREA_MEASURE, np);
    internal::interpolated_corrected_measure_mesh(pmesh, mu1_map, internal::MU1_MEAN_CURVATURE_MEASURE, np);

    for (vertex_descriptor v : vertices(pmesh))
    {
        internal::expand_interpolated_corrected_measure_vertex(pmesh, mu0_map, mu1_map, mu0_expand_map, mu1_expand_map, v, np.ball_radius(r));

        typename GT::FT v_mu0 = get(mu0_expand_map, v);
        if (v_mu0 != 0.0)
            put(vcm, v, 0.5 * get(mu1_expand_map, v) / v_mu0);
        else
            put(vcm, v, 0);
    }
}

/**
* \ingroup PMP_corrected_curvatures_grp
*
* Computes the interpolated corrected gaussian curvature across the mesh
* and stores it in a veretex property map `vcm`.
*
* @tparam PolygonMesh a model of `FaceListGraph`.
* @tparam VertexCurvatureMap model of `WritablePropertyMap` with
* `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and `GT::FT` as value type.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param pmesh the polygon mesh.
* @param vcm the vertex property map in which the computed gaussian curvatures are stored.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters".
*
* @see `interpolated_corrected_mean_curvature()`
* @see `interpolated_corrected_principal_curvatures()`
*/
template<typename PolygonMesh, typename VertexCurvatureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void interpolated_corrected_gaussian_curvature(const PolygonMesh& pmesh,
        VertexCurvatureMap vcm,
        const NamedParameters& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::property_map<PolygonMesh,
        CGAL::dynamic_face_property_t<typename GT::FT>>::const_type FaceMeasureMap;

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::property_map<PolygonMesh,
        CGAL::dynamic_vertex_property_t<typename GT::FT>>::const_type VertexMeasureMap;

    typename GT::FT
        r = choose_parameter(get_parameter(np, internal_np::ball_radius), 0);

    if (r == 0)
        r = internal::average_edge_length<PolygonMesh, GT>(pmesh) * EXPANDING_RADIUS_EPSILON;

    FaceMeasureMap mu0_map = get(CGAL::dynamic_face_property_t<typename GT::FT>(), pmesh);
    FaceMeasureMap mu2_map = get(CGAL::dynamic_face_property_t<typename GT::FT>(), pmesh);

    VertexMeasureMap mu0_expand_map = get(CGAL::dynamic_vertex_property_t<typename GT::FT>(), pmesh);
    VertexMeasureMap mu2_expand_map = get(CGAL::dynamic_vertex_property_t<typename GT::FT>(), pmesh);

    internal::interpolated_corrected_measure_mesh(pmesh, mu0_map, internal::MU0_AREA_MEASURE, np);
    internal::interpolated_corrected_measure_mesh(pmesh, mu2_map, internal::MU2_GAUSSIAN_CURVATURE_MEASURE, np);

    for (vertex_descriptor v : vertices(pmesh))
    {
        internal::expand_interpolated_corrected_measure_vertex(pmesh, mu0_map, mu2_map, mu0_expand_map, mu2_expand_map, v, np.ball_radius(r));

        typename GT::FT v_mu0 = get(mu0_expand_map, v);
        if(v_mu0 != 0.0)
            put(vcm, v, get(mu2_expand_map, v) / v_mu0);
        else
            put(vcm, v, 0);
    }
}

/**
* \ingroup PMP_corrected_curvatures_grp
*
* Computes the interpolated corrected principal curvatures across the mesh
* and stores it in a veretex property map `vcm`.
*
* @tparam PolygonMesh a model of `FaceListGraph`.
* @tparam VertexCurvatureMap model of `WritablePropertyMap` with
* `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
* `std::tuple<GT::FT, GT::FT, Eigen::Vector<GT::FT, 3>, Eigen::Vector<GT::FT, 3>>` as value type.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param pmesh the polygon mesh.
* @param vcm the vertex property map in which the computed principal curvatures are stored.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for
*                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_normal_map}
*     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Vector_3` as value type}
*     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, vertex normals will be
*                     computed using compute_vertex_normals()}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @see `interpolated_corrected_mean_curvature()`
* @see `interpolated_corrected_gaussian_curvature()`
*/
template<typename PolygonMesh, typename VertexCurvatureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void interpolated_corrected_principal_curvatures(const PolygonMesh& pmesh,
        VertexCurvatureMap vcm,
        const NamedParameters& np = parameters::default_values())
{
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

    typedef dynamic_vertex_property_t<typename GT::Vector_3> Vector_map_tag;
    typedef typename boost::property_map<PolygonMesh, Vector_map_tag>::const_type Default_vector_map;
    typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
        NamedParameters,
        Default_vector_map>::type       VNM;

    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
            get_const_property_map(CGAL::vertex_point, pmesh));

    VNM vnm = choose_parameter(get_parameter(np, internal_np::vertex_normal_map),
        get(Vector_map_tag(), pmesh));

    typename GT::FT
        r = choose_parameter(get_parameter(np, internal_np::ball_radius), 0);

    if (r == 0)
        r = internal::average_edge_length<PolygonMesh, GT>(pmesh) * EXPANDING_RADIUS_EPSILON;

    if (is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value)
        compute_vertex_normals(pmesh, vnm, np);



    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef std::unordered_map<face_descriptor, typename GT::FT> FaceScalarMeasureMap_tag;
    // using std:: array to store FT values on the 9 combinations of the standard 3D basis
    typedef std::unordered_map<face_descriptor, std::array<typename GT::FT, 3 * 3>> FaceArrayMeasureMap_tag;

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef std::unordered_map<vertex_descriptor, typename GT::FT> VertexScalarMeasureMap_tag;
    // using Eigen matrix to store & Process FT values on the 9 combinations of the standard 3D basis
    typedef std::unordered_map<vertex_descriptor, Eigen::Matrix<typename GT::FT, 3, 3>> VertexMatrixMeasureMap_tag;


    FaceScalarMeasureMap_tag mu0_init;
    boost::associative_property_map<FaceScalarMeasureMap_tag>
        mu0_map(mu0_init);

    FaceArrayMeasureMap_tag muXY_init;
    boost::associative_property_map<FaceArrayMeasureMap_tag>
        muXY_map(muXY_init);

    VertexScalarMeasureMap_tag mu0_expand_init;
    boost::associative_property_map<VertexScalarMeasureMap_tag>
        mu0_expand_map(mu0_expand_init);

    VertexMatrixMeasureMap_tag muXY_expand_init;
    boost::associative_property_map<VertexMatrixMeasureMap_tag>
        muXY_expand_map(muXY_expand_init);

    internal::interpolated_corrected_measure_mesh(pmesh, mu0_map, internal::MU0_AREA_MEASURE, np);
    internal::interpolated_corrected_anisotropic_measure_mesh(pmesh, muXY_map, np);

    for (vertex_descriptor v : vertices(pmesh))
    {
        internal::expand_interpolated_corrected_anisotropic_measure_vertex(pmesh, mu0_map, muXY_map, mu0_expand_map, muXY_expand_map, v, np.ball_radius(r));

        typename GT::FT v_mu0 = get(mu0_expand_map, v);
        Eigen::Matrix<typename GT::FT, 3, 3> v_muXY = get(muXY_expand_map, v);

        typename GT::Vector_3 u_GT = get(vnm, v);

        Eigen::Matrix<typename GT::FT, 3, 1> u(u_GT.x(), u_GT.y(), u_GT.z());

        const typename GT::FT K = 1000 * v_mu0;

        v_muXY = 0.5 * (v_muXY + v_muXY.transpose()) + K * u * u.transpose();

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix <typename GT::FT, 3, 3>> eigensolver;

        eigensolver.computeDirect(v_muXY);

        if (eigensolver.info() != Eigen::Success)
        {
            put(vcm, v, std::make_tuple(
                0,
                0,
                typename GT::Vector_3(0, 0, 0),
                typename GT::Vector_3(0, 0, 0)));
            continue;
        }

        const Eigen::Matrix<typename GT::FT, 3, 1> eig_vals = eigensolver.eigenvalues();
        const Eigen::Matrix<typename GT::FT, 3, 3> eig_vecs = eigensolver.eigenvectors();

        const typename GT::Vector_3 min_eig_vec(eig_vecs(0, 1), eig_vecs(1, 1), eig_vecs(2, 1));
        const typename GT::Vector_3 max_eig_vec(eig_vecs(0, 0), eig_vecs(1, 0), eig_vecs(2, 0));

        put(vcm, v, std::make_tuple(
            (v_mu0 != 0.0) ? -eig_vals[1] / v_mu0 : 0.0,
            (v_mu0 != 0.0) ? -eig_vals[0] / v_mu0 : 0.0,
            min_eig_vec,
            max_eig_vec));
    }
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H