#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Octree_domain.h>
#include <CGAL/Octree_wrapper.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <math.h>

#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Vector_3 Vector_3;
typedef typename Kernel::Point_3 Point_3;

typedef std::vector<Point_3> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

typedef CGAL::Isosurfacing::Octree_wrapper<Kernel> Octree_wrapper_;
typedef CGAL::Isosurfacing::Octree_domain<Kernel> Octree_domain_;

Kernel::FT sphere_function(const Point_3& point) {
    return std::sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z());
}

struct Refine_one_eighth {
    std::size_t min_depth_;
    std::size_t max_depth_;

    std::size_t octree_dim_;

    Octree_wrapper_::Uniform_coords uniform_coordinates(const Octree_wrapper_::Octree::Node& node) const {
        auto coords = node.global_coordinates();
        const std::size_t depth_factor = std::size_t(1) << (max_depth_ - node.depth());
        for (int i = 0; i < Octree_wrapper_::Octree::Node::Dimension::value; ++i) {
            coords[i] *= (uint32_t)depth_factor;
        }

        return coords;
    }

    Refine_one_eighth(std::size_t min_depth, std::size_t max_depth) : min_depth_(min_depth), max_depth_(max_depth) {
        octree_dim_ = std::size_t(1) << max_depth_;
    }

    bool operator()(const Octree_wrapper_::Octree::Node& n) const {
        // n.depth()
        if (n.depth() < min_depth_) {
            return true;
        }
        if (n.depth() == max_depth_) {
            return false;
        }

        auto leaf_coords = uniform_coordinates(n);

        if (leaf_coords[0] >= octree_dim_ / 2) {
            return false;
        }
        if (leaf_coords[1] >= octree_dim_ / 2) {
            return false;
        }
        if (leaf_coords[2] >= octree_dim_ / 2) {
            // return false;
        }
        return true;
    }
};

int main() {
    Octree_wrapper_ octree_wrap({-1, -1, -1, 1, 1, 1});

    Refine_one_eighth split_predicate(4, 6);
    octree_wrap.refine(split_predicate);

    Octree_domain_ octree_domain(octree_wrap);

    auto lam = [&](const Octree_domain_::Vertex_handle& v) {
        Point_3 p = octree_domain.position(v);
        const auto val = sphere_function(p);
        Vector_3 gradient = p - CGAL::ORIGIN;
        gradient = gradient / std::sqrt(gradient.squared_length());
        octree_wrap.value(v) = val;
        octree_wrap.gradient(v) = gradient;
    };
    octree_domain.iterate_vertices(lam, CGAL::Sequential_tag());

    Point_range points;
    Polygon_range polygons;

    CGAL::Isosurfacing::make_quad_mesh_using_dual_contouring(octree_domain, 0.8, points, polygons);

    CGAL::IO::write_OFF("result.off", points, polygons);
}
