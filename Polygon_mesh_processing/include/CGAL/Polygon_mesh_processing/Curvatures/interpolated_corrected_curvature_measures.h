#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H
#endif

#include <numeric>

#include <CGAL/assertions.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/property_map.h>


namespace CGAL {

namespace Polygon_mesh_processing {

/*!
 * \ingroup PMP_corrected_curvatures_grp
 * Enumeration type used to specify which measure is computed for the
 * interpolated corrected curvature functions
 */
// enum
enum Measure_index {
    MU0_AREA_MEASURE, ///< corrected area density of the given face
    MU1_MEAN_CURVATURE_MEASURE, ///< corrected mean curvature density of the given face
    MU2_GAUSSIAN_CURVATURE_MEASURE ///< corrected gaussian curvature density of the given face
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
* @return a scalar of type GT::FT. This is the value of the interpolated corrected measure of the given triangle.
*
* @see `interpolated_corrected_measure_face()`
* @see `interpolated_corrected_measure_quad()`
*/
template<typename GT>
typename GT::FT interpolated_corrected_measure_triangle(const typename GT::Vector_3 x0, const typename GT::Vector_3 x1, const typename GT::Vector_3 x2,
    const typename GT::Vector_3 u0, const typename GT::Vector_3 u1, const typename GT::Vector_3 u2, const Measure_index mu_i)
{
    switch (mu_i)
    {
    case MU0_AREA_MEASURE:
    {
        const typename GT::Vector_3 um = (u0 + u1 + u2) / 3.0;

        return 0.5 * um * CGAL::cross_product(x1 - x0, x2 - x0);
    }
    case MU1_MEAN_CURVATURE_MEASURE:
    {
        const typename GT::Vector_3 um = (u0 + u1 + u2) / 3.0;

        return 0.5 * um * (CGAL::cross_product(u2 - u1, x0)
            + CGAL::cross_product(u0 - u2, x1)
            + CGAL::cross_product(u1 - u0, x2));
    }
    case MU2_GAUSSIAN_CURVATURE_MEASURE:

        return 0.5 * u0 * CGAL::cross_product(u1, u2);

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
* @return a scalar of type GT::FT. This is the value of the interpolated corrected measure of the given triangle.
*
* @see `interpolated_corrected_measure_face()`
* @see `interpolated_corrected_measure_triangle()`
*/
template<typename GT>
typename GT::FT interpolated_corrected_measure_quad(const typename GT::Vector_3 x0, const typename GT::Vector_3 x1, const typename GT::Vector_3 x2, const typename GT::Vector_3 x3,
    const typename GT::Vector_3 u0, const typename GT::Vector_3 u1, const typename GT::Vector_3 u2, const typename GT::Vector_3 u3, const Measure_index mu_i)
{
    // x0  _  x1
    // x2 |_| x3

    switch (mu_i)
    {
    case MU0_AREA_MEASURE:

        return (1 / 36.0) * ((4 * u0 + 2 * u1 + 2 * u2 + u3) * CGAL::cross_product(x1 - x0, x2 - x0)
            + (2 * u0 + 4 * u1 + u2 + 2 * u3) * CGAL::cross_product(x1 - x0, x3 - x1)
            + (2 * u0 + u1 + 4 * u2 + 2 * u3) * CGAL::cross_product(x3 - x2, x2 - x0)
            + (u0 + 2 * u1 + 2 * u2 + 4 * u3) * CGAL::cross_product(x3 - x2, x3 - x1));

    case MU1_MEAN_CURVATURE_MEASURE:
    {
        const typename GT::Vector_3 u03 = u3 - u0;
        const typename GT::Vector_3 u12 = u2 - u1;
        const typename GT::Vector_3 x0_cross = CGAL::cross_product(u12, x0);
        const typename GT::Vector_3 x1_cross = -CGAL::cross_product(u03, x1);
        const typename GT::Vector_3 x2_cross = CGAL::cross_product(u03, x2);
        const typename GT::Vector_3 x3_cross = -CGAL::cross_product(u12, x3);


        return (1 / 12.0) * (u0 * (2 * x0_cross - CGAL::cross_product((u2 + u3), x1) + CGAL::cross_product((u1 + u3), x2) + x3_cross)
            + u1 * (CGAL::cross_product((u2 + u3), x0) + 2 * x1_cross + x2_cross - CGAL::cross_product((u0 + u2), x3))
            + u2 * (CGAL::cross_product(-(u1 + u3), x0) + x1_cross + 2 * x2_cross + CGAL::cross_product((u0 + u1), x3))
            + u3 * (x0_cross + CGAL::cross_product((u0 + u2), x1) - CGAL::cross_product((u0 + u1), x2) + 2 * x3_cross));
    }
    case MU2_GAUSSIAN_CURVATURE_MEASURE:

        return (1 / 36.0) * ((4 * u0 + 2 * u1 + 2 * u2 + u3) * CGAL::cross_product(u1 - u0, u2 - u0)
            + (2 * u0 + 4 * u1 + u2 + 2 * u3) * CGAL::cross_product(u1 - u0, u3 - u1)
            + (2 * u0 + u1 + 4 * u2 + 2 * u3) * CGAL::cross_product(u3 - u2, u2 - u0)
            + (u0 + 2 * u1 + 2 * u2 + 4 * u3) * CGAL::cross_product(u3 - u2, u3 - u1));

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
* @return a scalar of type GT::FT. This is the value of the interpolated corrected measure of the given face.
*
* @see `interpolated_corrected_measure_triangle()`
* @see `interpolated_corrected_measure_quad()`
* @see `interpolated_corrected_measure_mesh()`
*/
template<typename GT>
typename GT::FT interpolated_corrected_measure_face(const std::vector<typename GT::Vector_3>& x, const std::vector<typename GT::Vector_3>& u, const Measure_index mu_i)
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

        // getting barycenter of points
        typename GT::Vector_3 xm = std::accumulate(x.begin(), x.end(), GT::Vector_3(0, 0, 0));
        xm /= n;

        // getting unit average normal of points
        typename GT::Vector_3 um = std::accumulate(u.begin(), u.end(), GT::Vector_3(0, 0, 0));
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
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param pmesh the polygon mesh
* @param mu_i an enum for choosing between computing the area measure, the mean curvature measure or the gaussian curvature measure
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_normal_map}
*     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Vector_3` as value type}
*     \cgalParamDefault{`TODO`}
*     \cgalParamExtra{If this parameter is omitted, vertex normals should be computed inside the function body.}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @return a vector of the computed measure on all faces. The return type is a std::vector<GT::FT>.
* GT is the type of the Geometric Triats deduced from the PolygonMesh and the NamedParameters arguments
* This is to be changed later to a property_map<face_descriptor, GT::FT>.
*
* @see `interpolated_corrected_measure_face()`
*/
template<typename PolygonMesh,
    typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
    std::vector<FT>
#else
    std::vector<typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT>
#endif
    interpolated_corrected_measure_mesh(
        const PolygonMesh& pmesh,
        const Measure_index mu_i,
        NamedParameters& np = parameters::default_values())
{
    typedef GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
    typedef dynamic_vertex_property_t<GT::Vector_3> Vector_map_tag;
    typedef typename boost::property_map<PolygonMesh, Vector_map_tag>::type Default_vector_map;
    typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
        NamedParameters,
        Default_vector_map>::type       VNM;

    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;


    typedef boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
            get_const_property_map(CGAL::vertex_point, pmesh));

    VNM vnm = choose_parameter<Default_vector_map>(
        get_parameter(np, internal_np::vertex_normal_map), Vector_map_tag(), pmesh);

    if (is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value)
        compute_vertex_normals(pmesh, vnm, np);

    std::vector<GT::FT> mu_i_map;


    for (face_descriptor f : faces(pmesh))
    {
        halfedge_descriptor h_start = pmesh.halfedge(f);
        halfedge_descriptor h_iter = h_start;

        std::vector<GT::Vector_3> x;
        std::vector<GT::Vector_3> u;

        // looping over vertices in face
        do {
            vertex_descriptor v = source(h_iter, pmesh);
            GT::Point_3 p = get(vpm, v);
            x.push_back(GT::Vector_3(p.x(),p.y(),p.z()));
            u.push_back(get(vnm, v));
            h_iter = next(h_iter, pmesh);
        } while (h_iter != h_start);


        mu_i_map.push_back(interpolated_corrected_measure_face<GT>(x, u, mu_i));
    }
    return mu_i_map;
}
}
}