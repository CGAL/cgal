#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Variational_shape_approximation.h>
#include <boost/functional/hash_fwd.hpp>
#include <boost/property_map/vector_property_map.hpp>

#include <CGAL/Surface_mesh_approximation/L21_metric_plane_proxy.h>
#include <CGAL/Variational_shape_approximation.h>

#include <limits>

#include <CGAL/Surface_mesh.h>

namespace PMP = CGAL::Polygon_mesh_processing;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using SurfaceMesh = CGAL::Surface_mesh<K::Point_3>;


typedef boost::property_map<SurfaceMesh, boost::vertex_point_t>::type VertexPointMap;
typedef typename CGAL::Kernel_traits<
    typename boost::property_traits<VertexPointMap>::value_type>::Kernel GeomTraits;
typedef CGAL::Surface_mesh_approximation::
    L21_metric_plane_proxy<SurfaceMesh, VertexPointMap, GeomTraits>
        L21_Metric;
typedef CGAL::Variational_shape_approximation<SurfaceMesh, VertexPointMap, L21_Metric>
    L21_MeshApproximation;

struct PolygonSoup
{
    using PointType = K::Point_3;
    std::vector<PointType> points;
    std::vector<std::array<std::size_t, 3>> faces;
};

auto meshArea = [](const SurfaceMesh& mesh) {
    double area = 0;
    for (const auto& f : mesh.faces())
    {
        std::vector<K::Point_3> pts;
        for (const auto& v : mesh.vertices_around_face(mesh.halfedge(f)))
        {
            pts.push_back(mesh.point(v));
        }
        K::Triangle_3 tri(pts[0], pts[1], pts[2]);
        area += std::sqrt(tri.squared_area());
    }
    return area;
};

auto soupArea = [](const PolygonSoup& soup) {
    double area = 0;
    for (const auto& f : soup.faces)
    {
        std::vector<K::Point_3> pts = {soup.points[f[0]], soup.points[f[1]], soup.points[f[2]]};
        K::Triangle_3 tri(pts[0], pts[1], pts[2]);
        area += std::sqrt(tri.squared_area());
    }
    return area;
};
auto soupToMesh = [](const PolygonSoup& soup) {
    SurfaceMesh mesh;
    std::vector<SurfaceMesh::Vertex_index> vertices;
    for (const auto& p : soup.points)
    {
        vertices.push_back(mesh.add_vertex(p));
    }
    for (const auto& f : soup.faces)
    {
        mesh.add_face(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
    }
    return mesh;
};

int main(int, char*[])
{
    std::array<double,6> subdivisionRatios          = { 0.10 , 2.0 , 2.0 , 10.0 , 10.0 , 10.00 };
    std::array<double, 6> boundarySubdivisionRatios = { 0.01 , 2.0 , 0.1 , 10.0 ,  1.0 ,  0.01 };

    for (std::size_t i = 0; i < subdivisionRatios.size(); ++i)
    {
        const auto subdivisionRatio = subdivisionRatios[i];
        const auto boundarySubdivisionRatio = boundarySubdivisionRatios[i];
        std::string filename = "./VSA-" + std::to_string(subdivisionRatio) + "-" + std::to_string(boundarySubdivisionRatio);
        SurfaceMesh mesh;

        std::ifstream in("data/patch.ply");
        CGAL::IO::read_PLY(in, mesh);
        std::cout << vertices(mesh).size() << std::endl;
        const int maxProxies = 10;
        const int defaultIterations = 20;

        VertexPointMap vpmap = get(boost::vertex_point, mesh);
        L21_Metric metric(mesh, vpmap);
        L21_MeshApproximation approx(mesh, vpmap, metric);

        const auto seedingParams = //
            CGAL::parameters::seeding_method(
                CGAL::Surface_mesh_approximation::Seeding_method::INCREMENTAL)
                .max_number_of_proxies(maxProxies)
                .min_error_drop(0.001)
                .number_of_relaxations(5);
        approx.initialize_seeds(seedingParams);

        approx.run(defaultIterations);

        const auto meshParams = CGAL::parameters::subdivision_ratio(subdivisionRatio)
                                    .boundary_subdivision_ratio(boundarySubdivisionRatio)
                                    .with_dihedral_angle(false)
                                    .optimize_boundary_anchor_location(false)
                                    .optimize_anchor_location(true)
                                    .pca_plane(false);
        const bool isOutputManifold = approx.extract_mesh(meshParams);
        if (!isOutputManifold)
        {
            throw std::runtime_error("not manifold");
        }

        PolygonSoup soup;
        const auto outputParams = CGAL::parameters::anchors(std::back_inserter(soup.points))
                                      .triangles(std::back_inserter(soup.faces));
        approx.output(outputParams);

        std::cerr << " Mesh area = " << meshArea(mesh) << " Soup area = " << soupArea(soup)
                  << std::endl;
        {
            std::ofstream f(filename + "soup.ply");
            CGAL::IO::write_PLY(f, soupToMesh(soup));
        }
    }

    return 0;
}
