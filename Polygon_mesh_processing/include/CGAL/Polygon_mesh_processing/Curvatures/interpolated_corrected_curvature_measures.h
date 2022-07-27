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
#endif

#include <CGAL/license/Polygon_mesh_processing/interpolated_corrected_curvature_measures.h>

#include <CGAL/assertions.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>

#include <numeric>
#include <queue>
#include <unordered_set>

namespace CGAL {

namespace Polygon_mesh_processing {

/*!
 * \ingroup PMP_corrected_curvatures_grp
 * Enumeration type used to specify which measure of a given face
 * is computed for the interpolated corrected curvature functions
 */
// enum
enum Curvature_measure_index {
    MU0_AREA_MEASURE, ///< corrected area density
    MU1_MEAN_CURVATURE_MEASURE, ///< corrected mean curvature density
    MU2_GAUSSIAN_CURVATURE_MEASURE ///< corrected gaussian curvature density
};

/**
* \ingroup PMP_corrected_curvatures_grp
*
* computes the interpolated corrected measure of specific triangle.
*
* @tparam GT is the geometric traits class.
*
* @param x0 is the position of vertex #0.
* @param x1 is the position of vertex #1.
* @param x2 is the position of vertex #2.
* @param u0 is the normal vector of vertex #0.
* @param u1 is the normal vector of vertex #1.
* @param u2 is the normal vector of vertex #2.
* @param mu_i an enum for choosing between computing the area measure,
* the mean curvature measure, or the gaussian curvature measure.
*
* @return a scalar of type `GT::FT`.
* This is the value of the interpolated corrected measure of the given triangle.
*
* @see `interpolated_corrected_measure_face()`
* @see `interpolated_corrected_measure_quad()`
*/
template<typename GT>
typename GT::FT interpolated_corrected_measure_triangle(const typename GT::Vector_3 x0,
                                                        const typename GT::Vector_3 x1,
                                                        const typename GT::Vector_3 x2,
                                                        const typename GT::Vector_3 u0,
                                                        const typename GT::Vector_3 u1,
                                                        const typename GT::Vector_3 u2,
                                                        const Curvature_measure_index mu_i)
{
    typename GT::Construct_cross_product_3 cross_product;
    switch (mu_i)
    {
    case MU0_AREA_MEASURE:
    {
        const typename GT::Vector_3 um = (u0 + u1 + u2) / 3.0;

        return 0.5 * um * cross_product(x1 - x0, x2 - x0);
    }
    case MU1_MEAN_CURVATURE_MEASURE:
    {
        const typename GT::Vector_3 um = (u0 + u1 + u2) / 3.0;

        return 0.5 * um * (cross_product(u2 - u1, x0)
            + cross_product(u0 - u2, x1)
            + cross_product(u1 - u0, x2));
    }
    case MU2_GAUSSIAN_CURVATURE_MEASURE:

        return 0.5 * u0 * cross_product(u1, u2);

    default: return 0;
    }
}

/**
* \ingroup PMP_corrected_curvatures_grp
*
* computes the interpolated corrected measure of specific quad
* Note that the vertices 0 to 3 are ordered like this \n
* v0  _  v1 \n
* v2 |_| v3
*
* @tparam GT is the geometric traits class.
*
* @param x0 is the position of vertex #0.
* @param x1 is the position of vertex #1.
* @param x2 is the position of vertex #2.
* @param x3 is the position of vertex #3.
* @param u0 is the normal vector of vertex #0.
* @param u1 is the normal vector of vertex #1.
* @param u2 is the normal vector of vertex #2.
* @param u3 is the normal vector of vertex #3.
* @param mu_i an enum for choosing between computing the area measure,
* the mean curvature measure, or the gaussian curvature measure.
*
* @return a scalar of type `GT::FT`.
* This is the value of the interpolated corrected measure of the given triangle.
*
* @see `interpolated_corrected_measure_face()`
* @see `interpolated_corrected_measure_triangle()`
*/
template<typename GT>
typename GT::FT interpolated_corrected_measure_quad(const typename GT::Vector_3 x0,
                                                    const typename GT::Vector_3 x1,
                                                    const typename GT::Vector_3 x2,
                                                    const typename GT::Vector_3 x3,
                                                    const typename GT::Vector_3 u0,
                                                    const typename GT::Vector_3 u1,
                                                    const typename GT::Vector_3 u2,
                                                    const typename GT::Vector_3 u3,
                                                    const Curvature_measure_index mu_i)
{
    // x0  _  x1
    // x2 |_| x3
    typename GT::Construct_cross_product_3 cross_product;
    switch (mu_i)
    {
    case MU0_AREA_MEASURE:

        return (1 / 36.0) * (
              (4 * u0 + 2 * u1 + 2 * u2 + u3) * cross_product(x1 - x0, x2 - x0)
            + (2 * u0 + 4 * u1 + u2 + 2 * u3) * cross_product(x1 - x0, x3 - x1)
            + (2 * u0 + u1 + 4 * u2 + 2 * u3) * cross_product(x3 - x2, x2 - x0)
            + (u0 + 2 * u1 + 2 * u2 + 4 * u3) * cross_product(x3 - x2, x3 - x1)
            );

    case MU1_MEAN_CURVATURE_MEASURE:
    {
        const typename GT::Vector_3 u03 = u3 - u0;
        const typename GT::Vector_3 u12 = u2 - u1;
        const typename GT::Vector_3 x0_cross = cross_product(u12, x0);
        const typename GT::Vector_3 x1_cross = -cross_product(u03, x1);
        const typename GT::Vector_3 x2_cross = cross_product(u03, x2);
        const typename GT::Vector_3 x3_cross = -cross_product(u12, x3);


        return (1 / 12.0) * (
              u0 * (2 * x0_cross - cross_product((u2 + u3), x1) + cross_product((u1 + u3), x2) + x3_cross)
            + u1 * (cross_product((u2 + u3), x0) + 2 * x1_cross + x2_cross - cross_product((u0 + u2), x3))
            + u2 * (-cross_product((u1 + u3), x0) + x1_cross + 2 * x2_cross + cross_product((u0 + u1), x3))
            + u3 * (x0_cross + cross_product((u0 + u2), x1) - cross_product((u0 + u1), x2) + 2 * x3_cross)
            );
    }
    case MU2_GAUSSIAN_CURVATURE_MEASURE:

        return (1 / 36.0) * (
              (4 * u0 + 2 * u1 + 2 * u2 + u3) * cross_product(u1 - u0, u2 - u0)
            + (2 * u0 + 4 * u1 + u2 + 2 * u3) * cross_product(u1 - u0, u3 - u1)
            + (2 * u0 + u1 + 4 * u2 + 2 * u3) * cross_product(u3 - u2, u2 - u0)
            + (u0 + 2 * u1 + 2 * u2 + 4 * u3) * cross_product(u3 - u2, u3 - u1)
            );

    default: return 0;
    }
}


/**
* \ingroup PMP_corrected_curvatures_grp
*
* computes the interpolated corrected measure of specific face.
*
* @tparam GT is the geometric traits class.
*
* @param x is a vector of the vertex positions of the face.
* @param u is a vector of the vertex nomrals of the face.
* @param mu_i an enum for choosing between computing the area measure,
* the mean curvature measure, or the gaussian curvature measure.
*
* @return a scalar of type `GT::FT`.
* This is the value of the interpolated corrected measure of the given face.
*
* @see `interpolated_corrected_measure_triangle()`
* @see `interpolated_corrected_measure_quad()`
* @see `interpolated_corrected_measure_mesh()`
*/
template<typename GT>
typename GT::FT interpolated_corrected_measure_face(const std::vector<typename GT::Vector_3>& x,
                                                    const std::vector<typename GT::Vector_3>& u,
                                                    const Curvature_measure_index mu_i)
{
    std::size_t n = x.size();
    CGAL_precondition(u.size() == n);
    CGAL_precondition(n >= 3);

    // Triangle: use triangle formulas
    if (n == 3)
        return interpolated_corrected_measure_triangle<GT>(x[0], x[1], x[2],
            u[0], u[1], u[2], mu_i);

    // Quad: use bilinear interpolation formulas
    else if (n == 4)
        // x[0]  _  x[1]  --->  x0  _  x1   (reason for changing order)
        // x[3] |_| x[2]  --->  x2 |_| x3
        return interpolated_corrected_measure_quad<GT>(x[0], x[1], x[3], x[2],
            u[0], u[1], u[3], u[2], mu_i);

    // N-gon: split into n triangles by barycenter and use triangle formulas for each
    else {
        typename GT::FT mu0 = 0;

        // getting center of points
        typename GT::Vector_3 xm =
            std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
        xm /= n;

        // getting unit average normal of points
        typename GT::Vector_3 um =
            std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
        um /= sqrt(um * um);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            mu0 += interpolated_corrected_measure_triangle<GT>(x[i], x[(i + 1) % n], xm,
                u[i], u[(i + 1) % n], um, mu_i);
        }
        return mu0;
    }
}


/**
* \ingroup PMP_corrected_curvatures_grp
*
* computes the interpolated corrected curvature measure on each face of the mesh
*
* @tparam PolygonMesh a model of `FaceGraph`
* @tparam FaceMeasureMap a a model of `WritablePropertyMap` with
* `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and `GT::FT` as value type.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param pmesh the polygon mesh
* @param fmm (face measure map) the property map used for storing the computed interpolated corrected measure
* @param mu_i an enum for choosing between computing
*             the area measure, the mean curvature measure or the gaussian curvature measure
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
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
*     \cgalParamExtra{If this parameter is omitted, vertex normals should be
*                     computed inside the function body.}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @see `interpolated_corrected_measure_face()`
*/
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

    for (face_descriptor f : faces(pmesh))
    {
        std::vector<typename GT::Vector_3> x;
        std::vector<typename GT::Vector_3> u;

        for (vertex_descriptor v : vertices_around_face(halfedge(f, pmesh), pmesh))
        {
            typename GT::Point_3 p = get(vpm, v);
            x.push_back(typename GT::Vector_3(p.x(), p.y(), p.z()));
            u.push_back(get(vnm, v));
        }

        put(fmm, f, interpolated_corrected_measure_face<GT>(x, u, mu_i));
    }
}

//
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
typename GT::FT face_in_ball_ratio_2(const std::vector<typename GT::Vector_3>& x,
                                     const typename GT::FT r,
                                     const typename GT::Vector_3 c)
{
    std::size_t n = x.size();

    // getting center of points
    typename GT::Vector_3 xm =
        std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
    xm /= n;

    typename GT::FT d_min = (xm - c).squared_length();
    typename GT::FT d_max = d_min;

    for (const typename GT::Vector_3 xi : x)
    {
        const typename GT::FT d_sq = (xi - c).squared_length();
        d_max = std::max(d_sq, d_max);
        d_min = std::min(d_sq, d_min);
    }

    if (d_max <= r * r) return 1.0;
    else if (r * r <= d_min) return 0.0;

    d_max = sqrt(d_max);
    d_min = sqrt(d_min);

    return (r - d_min) / (d_max - d_min);
}

template<typename PolygonMesh, typename FaceMeasureMap, typename VertexCurvatureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void expand_interpolated_corrected_measure_vertex(const PolygonMesh& pmesh,
        FaceMeasureMap fmm,
        VertexCurvatureMap vcm,
        const typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
        const NamedParameters& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

    const typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
        r = choose_parameter(get_parameter(np, internal_np::ball_radius), 0.01);

    if (r < 0.000001)
        return;

    typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
            get_const_property_map(CGAL::vertex_point, pmesh));


    std::queue<face_descriptor> bfs_q;
    std::unordered_set<face_descriptor> bfs_v;

    typename GT::Point_3 vp = get(vpm, v);
    typename GT::Vector_3 c = typename GT::Vector_3(vp.x(), vp.y(), vp.z());

    typename GT::FT corrected_mui = 0;

    for (face_descriptor f : faces_around_target(halfedge(v, pmesh), pmesh)) {
        if (f != boost::graph_traits<PolygonMesh>::null_face())
        {
            bfs_q.push(f);
            bfs_v.insert(f);
        }
    }
    while (!bfs_q.empty()) {
        face_descriptor fi = bfs_q.front();
        bfs_q.pop();

        // looping over vertices in face to get point coordinates
        std::vector<typename GT::Vector_3> x;
        for (vertex_descriptor vi : vertices_around_face(halfedge(fi, pmesh), pmesh))
        {
            typename GT::Point_3 pi = get(vpm, vi);
            x.push_back(typename GT::Vector_3(pi.x(), pi.y(), pi.z()));
        }

        const typename GT::FT f_ratio = face_in_ball_ratio_2<GT>(x, r, c);

        if (f_ratio > 0.000001)
        {
            corrected_mui += f_ratio * get(fmm, fi);
            for (face_descriptor fj : faces_around_face(halfedge(fi, pmesh), pmesh))
            {
                if (bfs_v.find(fj) == bfs_v.end() && fj != boost::graph_traits<PolygonMesh>::null_face())
                {
                    bfs_q.push(fj);
                    bfs_v.insert(fj);
                }
            }
        }
    }

    put(vcm, v, corrected_mui);
}

//template<typename PolygonMesh, typename FaceMeasureMap,
//    typename NamedParameters = parameters::Default_named_parameters>
//    void expand_interpolated_corrected_measure_mesh(const PolygonMesh& pmesh,
//        FaceMeasureMap fmm,
//        const NamedParameters& np = parameters::default_values())
//{
//    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
//    for (face_descriptor f : faces(pmesh))
//        expand_interpolated_corrected_measure_face(pmesh, fmm, f, np);
//}

template<typename PolygonMesh, typename VertexCurvatureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void interpolated_corrected_mean_curvature(const PolygonMesh& pmesh,
        VertexCurvatureMap vcm,
        const NamedParameters& np = parameters::default_values())
{
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef std::unordered_map<face_descriptor, typename GT::FT> FaceMeasureMap_tag;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef std::unordered_map<vertex_descriptor, typename GT::FT> VertexMeasureMap_tag;

    FaceMeasureMap_tag mu0_init, mu1_init;
    boost::associative_property_map<FaceMeasureMap_tag>
        mu0_map(mu0_init), mu1_map(mu1_init);

    VertexMeasureMap_tag mu0_expand_init, mu1_expand_init;
    boost::associative_property_map<VertexMeasureMap_tag>
        mu0_expand_map(mu0_expand_init), mu1_expand_map(mu1_expand_init);

    interpolated_corrected_measure_mesh(pmesh, mu0_map, MU0_AREA_MEASURE);
    interpolated_corrected_measure_mesh(pmesh, mu1_map, MU1_MEAN_CURVATURE_MEASURE);

    for (vertex_descriptor v : vertices(pmesh))
    {
        expand_interpolated_corrected_measure_vertex(pmesh, mu0_map, mu0_expand_map, v, np);
        expand_interpolated_corrected_measure_vertex(pmesh, mu1_map, mu1_expand_map, v, np);

        typename GT::FT v_mu0 = get(mu0_expand_map, v);
        if (v_mu0 > 0.000001)
            put(vcm, v, 0.5 * get(mu1_expand_map, v) / v_mu0);
        else
            put(vcm, v, 0);
    }
}

template<typename PolygonMesh, typename VertexCurvatureMap,
    typename NamedParameters = parameters::Default_named_parameters>
    void interpolated_corrected_gaussian_curvature(const PolygonMesh& pmesh,
        VertexCurvatureMap vcm,
        const NamedParameters& np = parameters::default_values())
{
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef std::unordered_map<face_descriptor, typename GT::FT> FaceMeasureMap_tag;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef std::unordered_map<vertex_descriptor, typename GT::FT> VertexMeasureMap_tag;

    FaceMeasureMap_tag mu0_init, mu2_init;
    boost::associative_property_map<FaceMeasureMap_tag>
        mu0_map(mu0_init), mu2_map(mu2_init);

    VertexMeasureMap_tag mu0_expand_init, mu2_expand_init;
    boost::associative_property_map<VertexMeasureMap_tag>
        mu0_expand_map(mu0_expand_init), mu2_expand_map(mu2_expand_init);

    interpolated_corrected_measure_mesh(pmesh, mu0_map, MU0_AREA_MEASURE);
    interpolated_corrected_measure_mesh(pmesh, mu2_map, MU2_GAUSSIAN_CURVATURE_MEASURE);

    for (vertex_descriptor v : vertices(pmesh))
    {
        expand_interpolated_corrected_measure_vertex(pmesh, mu0_map, mu0_expand_map, v, np);
        expand_interpolated_corrected_measure_vertex(pmesh, mu2_map, mu2_expand_map, v, np);

        typename GT::FT v_mu0 = get(mu0_expand_map, v);
        if(v_mu0 > 0.000001)
            put(vcm, v, get(mu2_expand_map, v) / v_mu0);
        else
            put(vcm, v, 0);
    }
}



}
}
