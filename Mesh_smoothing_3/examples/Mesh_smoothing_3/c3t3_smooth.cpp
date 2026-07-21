// Image
#include <CGAL/ImageIO.h>
#include <CGAL/Image_3.h>
// Domain
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Detect_features_in_image.h>
// C3t3
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_3/C3T3_helpers.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_segment_primitive_3.h>

#include <CGAL/Kernel/global_functions.h>

#include <map>
#include <set>

#include <Mesh_smoothing_3/Mesh_smoothing_3.h>

// Parallel tag
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// Domain
typedef CGAL::Labeled_mesh_domain_3<K> Image_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Image_domain> Mesh_domain;
// Triangulations
typedef CGAL::Mesh_vertex_base_3<K, Mesh_domain> Vertex_base;
typedef CGAL::Compact_mesh_cell_base_3<K, Mesh_domain> Cell_base;
typedef std::array<std::size_t, 2> Vertex_info;
typedef std::array<std::size_t, 2> Cell_info;
typedef CGAL::Triangulation_vertex_base_with_info_3<Vertex_info, K, Vertex_base> Vertex_base_with_info;
typedef CGAL::Triangulation_cell_base_with_info_3<Cell_info, K, Cell_base> Cell_base_with_info;
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag, Vertex_base_with_info, Cell_base_with_info>::type Tr;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// C3t3
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;


#if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(6, 2, 0)
    using Image_internal = CGAL::_image;
    auto image_WK_FIXED = CGAL::WK_FIXED;
    auto image_SGN_UNSIGNED = CGAL::SGN_UNSIGNED;
// Code for CGAL 6.0 or newer
#else
    using Image_internal = _image;
    auto image_WK_FIXED = WK_FIXED;
    auto image_SGN_UNSIGNED = SGN_UNSIGNED;
#endif


// AABB traits for C3t3 facets and edges
struct Facet_to_triangle_property_map
{
    typedef C3t3::Facet key_type;
    typedef K::Triangle_3 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;
};

inline Facet_to_triangle_property_map::reference
get(const Facet_to_triangle_property_map&, const C3t3::Facet& f)
{
    const Tr::Cell_handle cell = f.first;
    const int i = f.second;
    const K::Point_3& p0 = cell->vertex((i + 1) % 4)->point().point();
    const K::Point_3& p1 = cell->vertex((i + 2) % 4)->point().point();
    const K::Point_3& p2 = cell->vertex((i + 3) % 4)->point().point();
    return K::Triangle_3(p0, p1, p2);
}

struct Facet_to_point_property_map
{
    typedef C3t3::Facet key_type;
    typedef K::Point_3 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;
};

inline Facet_to_point_property_map::reference
get(const Facet_to_point_property_map&, const C3t3::Facet& f)
{
    const Tr::Cell_handle cell = f.first;
    const int i = f.second;
    return cell->vertex((i + 1) % 4)->point().point();
}

struct Edge_to_segment_property_map
{
    typedef C3t3::Edge key_type;
    typedef K::Segment_3 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;
};

inline Edge_to_segment_property_map::reference
get(const Edge_to_segment_property_map&, const C3t3::Edge& e)
{
    const Tr::Cell_handle cell = e.first;
    return K::Segment_3(cell->vertex(e.second)->point().point(), cell->vertex(e.third)->point().point());
}

struct Edge_to_point_property_map
{
    typedef C3t3::Edge key_type;
    typedef K::Point_3 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;
};

inline Edge_to_point_property_map::reference
get(const Edge_to_point_property_map&, const C3t3::Edge& e)
{
    const Tr::Cell_handle cell = e.first;
    return cell->vertex(e.second)->point().point();
}

typedef CGAL::AABB_primitive<C3t3::Facet,
                             Facet_to_triangle_property_map,
                             Facet_to_point_property_map,
                             CGAL::Tag_false,
                             CGAL::Tag_false> Facet_primitive;
typedef CGAL::AABB_traits_3<K, Facet_primitive> Facet_traits;
typedef CGAL::AABB_tree<Facet_traits> Facet_tree;

typedef CGAL::AABB_primitive<C3t3::Edge,
                             Edge_to_segment_property_map,
                             Edge_to_point_property_map,
                             CGAL::Tag_false,
                             CGAL::Tag_false> Edge_primitive;
typedef CGAL::AABB_traits_3<K, Edge_primitive> Edge_traits;
typedef CGAL::AABB_tree<Edge_traits> Edge_tree;

// To avoid verbose function and named parameters call
namespace params = CGAL::parameters;


// image of type 'unsigned char' initialised with '0'
CGAL::Image_3 create_cgal_image(const std::size_t& xdim, const std::size_t& ydim, const std::size_t& zdim,
                                double vx = 1.0, double vy = 1.0, double vz = 1.0)
{
    Image_internal* im = _createImage(xdim, ydim, zdim, 1,
                                    vx, vy, vz, 1,
                                    image_WK_FIXED, image_SGN_UNSIGNED);
    std::fill_n(static_cast<unsigned char*>(im->data), xdim*ydim*zdim, 0);
    return CGAL::Image_3(im);
}

void add_rectangle_in_image(const K::Point_3& rectangle_min, // [0..1]^3
                            const K::Point_3& rectangle_max, // [0..1]^3
                            const unsigned char& label,
                            CGAL::Image_3 &image)
{
    using CGAL::IMAGEIO::static_evaluate;
    Image_internal* im = image.image();
    const std::size_t min_x = rectangle_min.x() * im->xdim;
    const std::size_t min_y = rectangle_min.y() * im->ydim;
    const std::size_t min_z = rectangle_min.z() * im->zdim;
    const std::size_t max_x = rectangle_max.x() * im->xdim;
    const std::size_t max_y = rectangle_max.y() * im->ydim;
    const std::size_t max_z = rectangle_max.z() * im->zdim;
    for (std::size_t k = min_z; k < max_z; k++)
        for (std::size_t j = min_y; j < max_y; j++)
            for (std::size_t i = min_x; i < max_x; i++)
                static_evaluate<unsigned char>(im, i, j, k) = label;
}

void deform_c3t3_smooth_fold(C3t3& c3t3, const K::FT& angle, const K::FT& nb_planes)
{
    typedef C3t3::Point WPoint_3;
    K::Point_3 center = CGAL::midpoint(K::Point_3(c3t3.bbox().min(0), c3t3.bbox().min(1), c3t3.bbox().min(2)),
                                       K::Point_3(c3t3.bbox().max(0), c3t3.bbox().max(1), c3t3.bbox().max(2)));
    K::Plane_3 start_plane(center, K::Vector_3(1,0,0));
    Tr& triangulation = c3t3.triangulation();
    std::list<Tr::Vertex*> vertices;
    for (Tr::Vertex& vertex : triangulation.tds().vertices())
    {
        vertices.emplace_back(&vertex);
    }
    const K::FT angle_smooth_length = (c3t3.bbox().max(0) - c3t3.bbox().min(0)) * 0.1875; // * 0.25;
    const K::FT angle_smooth_length_increment = angle_smooth_length / nb_planes;
    const K::FT angle_increment = angle/nb_planes;
    const K::FT ca_i = cos(angle_increment);
    const K::FT sa_i = sin(angle_increment);
    K::Point_3 plane_center = center;
    plane_center -= K::Vector_3(angle_smooth_length * 0.5, 0, 0);
    for (std::size_t p = 0; p < nb_planes; p++)
    {
        const K::Vector_3 plane_normal(cos(p*angle_increment), sin(p*angle_increment), 0);
        const K::Plane_3 plane(plane_center, plane_normal);
        const K::FT plane_x = plane_center[0];
        const K::FT plane_y = plane_center[1];
        const K::FT plane_z = plane_center[2];
        for (std::list<Tr::Vertex*>::iterator it = vertices.begin(); it != vertices.end();)
        {
            Tr::Vertex& vertex = **it;
            const WPoint_3& point = vertex.point();
            if (plane.has_on_positive_side(point.point()))
            {
                const K::FT px = point.point()[0] - plane_x;
                const K::FT py = point.point()[1] - plane_y;
                const K::FT pz = point.point()[2] - plane_z;
                const K::Point_3 rotated_point(px*ca_i-py*sa_i + plane_x, px*sa_i+py*ca_i + plane_y, pz + plane_z);
                vertex.set_point(WPoint_3(rotated_point, point.weight()));
                it++;
            }
            else
            {
                it = vertices.erase(it);
            }
        }
        plane_center += angle_smooth_length_increment * K::Vector_3(cos((p+1) * angle_increment), sin((p+1) * angle_increment), 0);
    }
}


int main(int argc, char* argv[])
{
    CGAL::Image_3 image = create_cgal_image(64, 64, 64);
    add_rectangle_in_image(K::Point_3(0.2, 0.45, 0.2), K::Point_3(0.5, 0.55, 0.5), 1, image);
    add_rectangle_in_image(K::Point_3(0.5, 0.45, 0.2), K::Point_3(0.8, 0.55, 0.5), 2, image);
    add_rectangle_in_image(K::Point_3(0.2, 0.45, 0.5), K::Point_3(0.5, 0.55, 0.8), 3, image);
    add_rectangle_in_image(K::Point_3(0.5, 0.45, 0.5), K::Point_3(0.8, 0.55, 0.8), 4, image);
    _writeImage(image.image(), "input_image.inr");

    // Make c3t3
    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image,
        params::features_detector = CGAL::Mesh_3::Detect_features_in_image());
    CGAL::Bbox_3 bbox = domain.bbox();
    double diag = CGAL::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                             CGAL::square(bbox.ymax() - bbox.ymin()) +
                             CGAL::square(bbox.zmax() - bbox.zmin()));
    double sizing_default = diag * 0.01;
    Mesh_criteria criteria(params::edge_size = sizing_default,
        params::facet_angle = 30,
        params::facet_size = sizing_default,
        params::facet_distance = sizing_default / 10,
        params::facet_topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH,
        params::cell_radius_edge_ratio = 1.5,
        params::cell_size = 0
    );
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        params::no_exude(),
                                        params::no_perturb());
    // Output
    CGAL::dump_c3t3(c3t3, "c3t3_initial");

    // Make deformed c3t3
    C3t3 c3t3_deformed(static_cast<const C3t3>(c3t3));
    deform_c3t3_smooth_fold(c3t3_deformed, 3.141592635 * 0.5, 100);

    C3t3 c3t3_deformed_ref(static_cast<const C3t3>(c3t3));
    deform_c3t3_smooth_fold(c3t3_deformed_ref, 3.141592635 * 0.5, 100);
    CGAL::dump_c3t3(c3t3_deformed, "c3t3_deformed");

    
    // AABB trees of the surface and curve features of the c3t3
    std::map<C3t3::Surface_patch_index, std::vector<C3t3::Facet>> facets_by_patch;
    for (auto f : c3t3_deformed_ref.facets_in_complex())
    {
        facets_by_patch[c3t3_deformed_ref.surface_patch_index(f)].push_back(f);
    }
    
    std::map<C3t3::Curve_index, std::vector<C3t3::Edge>> edges_by_curve;
    for (auto e : c3t3_deformed_ref.edges_in_complex())
    {
        edges_by_curve[c3t3_deformed_ref.curve_index(e)].push_back(e);
    }

    std::map<C3t3::Surface_patch_index, Facet_tree> facet_trees;
    for (auto& kv : facets_by_patch)
    {
        facet_trees.try_emplace(kv.first, kv.second.begin(), kv.second.end());
    }

    std::map<C3t3::Curve_index, Edge_tree> edge_trees;
    for (auto& kv : edges_by_curve)
    {
        edge_trees.try_emplace(kv.first, kv.second.begin(), kv.second.end());
    }

    std::cout << "Built " << facet_trees.size() << " facet AABB trees and " << edge_trees.size() << " edge AABB trees." << std::endl;


    // Smoothing

    Mesh_smoothing_3::C3t3_smoother smoother(c3t3_deformed);
    for (auto c : c3t3.vertices_in_complex())
    {
        smoother.set_vertex_lock(c, true);
    }

    smoother.set_boundary_query([&](K::Point_3 const &pt, C3t3::Surface_patch_index id, double radius) -> std::tuple<K::Point_3, K::Vector_3, double> {
        auto res = facet_trees.at(id).closest_point_and_primitive(pt);
        K::Point_3 closest_point = res.first;
        const auto triangle = c3t3_deformed_ref.triangulation().triangle(res.second);
        K::Vector_3 normal = CGAL::unit_normal(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));
        return std::make_tuple(closest_point, normal, 1.);
    });

    smoother.set_curves_query([&](K::Point_3 const &pt, C3t3::Curve_index id, double radius) -> std::tuple<K::Point_3, K::Vector_3, double> {
        auto res = edge_trees.at(id).closest_point_and_primitive(pt);
        K::Point_3 closest_point = res.first;
        const auto segment = c3t3_deformed_ref.triangulation().segment(res.second);
        K::Vector_3 direction = (segment.target() - segment.source());
        if (direction.squared_length() > 1e-8)
        {
            direction /= CGAL::sqrt(direction.squared_length());
        }
        return std::make_tuple(closest_point, direction, 1.);
    });

    smoother.set_verbose();
    smoother.set_max_number_of_iteration(100);
    smoother.run();

    CGAL::dump_c3t3(c3t3_deformed, "c3t3_smoothed");

    return 0;
}