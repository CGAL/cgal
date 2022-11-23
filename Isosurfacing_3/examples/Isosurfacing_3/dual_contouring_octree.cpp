#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Implicit_octree_domain.h>
#include <CGAL/Octree_wrapper.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <math.h>

#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Point_3 Point;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

typedef CGAL::Isosurfacing::Octree_wrapper<Kernel> Octree_wrapper_;

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
    const CGAL::Bbox_3 bbox(-1, -1, -1, 1, 1, 1);
    std::shared_ptr<Octree_wrapper_> octree_wrap = std::make_shared<Octree_wrapper_>(bbox);

    Refine_one_eighth split_predicate(3, 4);
    octree_wrap->refine(split_predicate);

    auto sphere_function = [&](const Point& p) { return std::sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z()); };

    auto sphere_gradient = [&](const Point& p) {
        const Vector g = p - CGAL::ORIGIN;
        return g / std::sqrt(g.squared_length());
    };

    auto domain = CGAL::Isosurfacing::create_implicit_octree_domain(octree_wrap, sphere_function, sphere_gradient);

    Point_range points;
    Polygon_range polygons;

    CGAL::Isosurfacing::dual_contouring(domain, 0.8, points, polygons);

    CGAL::IO::write_OFF("result.off", points, polygons);

    return 0;
}
