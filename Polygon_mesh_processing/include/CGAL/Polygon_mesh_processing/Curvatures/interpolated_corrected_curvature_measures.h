#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <numeric>

#include <CGAL/assertions.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/property_map.h>


namespace CGAL {

    namespace Polygon_mesh_processing {

        typedef Exact_predicates_inexact_constructions_kernel Epic;

        // enum to specify which measure is computed
        enum Measure_index {
            MU0_AREA_MEASURE,
            MU1_MEAN_CURVATURE_MEASURE,
            MU2_GAUSSIAN_CURVATURE_MEASURE
        };

        Epic::FT interpolated_mu_i_triangle(const Epic::Vector_3 x0, const Epic::Vector_3 x1, const Epic::Vector_3 x2,
            const Epic::Vector_3 u0, const Epic::Vector_3 u1, const Epic::Vector_3 u2, const Measure_index mu_i)
        {
            Epic::Vector_3 um;
            switch (mu_i)
            {
            case MU0_AREA_MEASURE:

                um = (u0 + u1 + u2) / 3.0;

                return 0.5 * um * CGAL::cross_product(x1 - x0, x2 - x0);

            case MU1_MEAN_CURVATURE_MEASURE:

                um = (u0 + u1 + u2) / 3.0;

                return 0.5 * um * (CGAL::cross_product(u2 - u1, x0)
                    + CGAL::cross_product(u0 - u2, x1)
                    + CGAL::cross_product(u1 - u0, x2));

            case MU2_GAUSSIAN_CURVATURE_MEASURE:

                return 0.5 * u0 * CGAL::cross_product(u1, u2);

            default: return 0;
            }
        }

        Epic::FT interpolated_mu_i_quad(const Epic::Vector_3 x0, const Epic::Vector_3 x1, const Epic::Vector_3 x2, const Epic::Vector_3 x3,
            const Epic::Vector_3 u0, const Epic::Vector_3 u1, const Epic::Vector_3 u2, const Epic::Vector_3 u3, const Measure_index mu_i)
        {
            /// x0  _  x1
            /// x2 |_| x3

            switch (mu_i)
            {
            case MU0_AREA_MEASURE:

                return (1 / 36.0) * ((4 * u0 + 2 * u1 + 2 * u2 + u3) * CGAL::cross_product(x1 - x0, x2 - x0)
                    + (2 * u0 + 4 * u1 + u2 + 2 * u3) * CGAL::cross_product(x1 - x0, x3 - x1)
                    + (2 * u0 + u1 + 4 * u2 + 2 * u3) * CGAL::cross_product(x3 - x2, x2 - x0)
                    + (u0 + 2 * u1 + 2 * u2 + 4 * u3) * CGAL::cross_product(x3 - x2, x3 - x1));

            case MU1_MEAN_CURVATURE_MEASURE:
            {
                const Epic::Vector_3 u03 = u3 - u0;
                const Epic::Vector_3 u12 = u2 - u1;
                const Epic::Vector_3 x0_cross = CGAL::cross_product(u12, x0);
                const Epic::Vector_3 x1_cross = -CGAL::cross_product(u03, x1);
                const Epic::Vector_3 x2_cross = CGAL::cross_product(u03, x2);
                const Epic::Vector_3 x3_cross = -CGAL::cross_product(u12, x3);


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

        Epic::FT interpolated_mu_i_face(const std::vector<Epic::Vector_3>& x, const std::vector<Epic::Vector_3>& u, const Measure_index mu_i)
        {
            std::size_t n = x.size();
            CGAL_precondition(u.size() == n);
            CGAL_precondition(n >= 3);

            // Triangle: use triangle formulas
            if (n == 3)
                return interpolated_mu_i_triangle(x[0], x[1], x[2],
                    u[0], u[1], u[2], mu_i);

            // Quad: use bilinear interpolation formulas
            else if (n == 4)
                /// x[0]  _  x[1]  --->  x0  _  x1   (reason for changing order)
                /// x[3] |_| x[2]  --->  x2 |_| x3
                return interpolated_mu_i_quad(x[0], x[1], x[3], x[2],
                    u[0], u[1], u[3], u[2], mu_i);

            // N-gon: split into n triangles by barycenter and use triangle formulas for each
            else {
                Epic::FT mu0 = 0;

                // getting barycenter of points
                Epic::Vector_3 xm = std::accumulate(x.begin(), x.end(), Epic::Vector_3(0, 0, 0));
                xm /= n;

                // getting unit average normal of points
                Epic::Vector_3 um = std::accumulate(u.begin(), u.end(), Epic::Vector_3(0, 0, 0));
                um /= sqrt(um * um);

                // summing each triangle's measure after triangulation by barycenter split.
                for (std::size_t i = 0; i < n; i++)
                {
                    mu0 += interpolated_mu_i_triangle(x[i], x[(i + 1) % n], xm,
                        u[i], u[(i + 1) % n], um, mu_i);
                }
                return mu0;
            }
        }

        
        /// TODO:
        /// 1- Handle if VNM is not given
        /// 2- use GT instead of Epic

        template<typename PolygonMesh,
            typename VertexNormalMap,
            typename NamedParameters = parameters::Default_named_parameters>
            std::vector<Epic::FT>
            interpolated_corrected_measure_i(
                const PolygonMesh& pmesh,
                const Measure_index mu_i,
                VertexNormalMap vnm,
                const NamedParameters& np = parameters::default_values())
        {
            using parameters::choose_parameter;
            using parameters::get_parameter;

            typedef boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
            typedef boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
            typedef boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

            typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
                vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
                    get_const_property_map(CGAL::vertex_point, pmesh));

            //std::unordered_map<vertex_descriptor, Epic::Vector_3> vnm_init;

            //boost::associative_property_map
            //    <std::unordered_map<boost::graph_traits<PolygonMesh>::vertex_descriptor, Epic::Vector_3>> vnm(vnm_init);

            ////if (!vnm)
            // {
            //    compute_vertex_normals(pmesh, vnm, np);
            //}

            std::vector<Epic::FT> mu_i_map;


            for (face_descriptor f : faces(pmesh))
            {
                halfedge_descriptor h_start = pmesh.halfedge(f);
                halfedge_descriptor h_iter = h_start;

                std::vector<Epic::Vector_3> x;
                std::vector<Epic::Vector_3> u;

                // looping over vertices in face
                do {
                    vertex_descriptor v = source(h_iter, pmesh);
                    Epic::Point_3 p = get(vpm, v);
                    x.push_back(Epic::Vector_3(p.x(),p.y(),p.z()));
                    u.push_back(get(vnm, v));
                    h_iter = next(h_iter, pmesh);
                } while (h_iter != h_start);


                mu_i_map.push_back(interpolated_mu_i_face(x, u, mu_i));
            }
            return mu_i_map;
        }

    }
}

