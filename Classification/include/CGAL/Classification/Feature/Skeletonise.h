#ifndef CGAL_CLASSIFICATION_FEATURE_SKELETONISE_H
#define CGAL_CLASSIFICATION_FEATURE_SKELETONISE_H

#include <cmath>
#include <utility>

#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Classification/Point_set_neighborhood.h>

#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <boost/graph/undirected_dfs.hpp>

namespace CGAL {
namespace   Classification {
namespace     Feature {

template <typename Point>
std::pair<std::vector<Point>, std::vector<float> > convert_polylines_to_point_cloud(const std::vector<std::vector<Point>>& polylines, int num_samples, const std::vector<float>& cycles) {
    std::vector<Point> point_cloud;
    std::vector<float> v1;
    int k = 0;
    for (const auto& polyline : polylines) {
        int num_points = polyline.size();

        if (num_points < 2) {
            continue;
        }

        // Go through each segment of the polyline
        for (int i = 0; i < num_points - 1; i++) {
            const Point& prev_point = polyline[i];
            const Point& next_point = polyline[i + 1];

            // Sample points along the segment and calculate distances in the cloud point
            Point current_point = prev_point;
            for (int j = 0; j < num_samples; j++) {
                float t = static_cast<float>(j + 1) / (num_samples + 1);
                Point interpolated_point;
                interpolated_point = Point((1 - t) * prev_point.x() + t * next_point.x(),
                    (1 - t) * prev_point.y() + t * next_point.y(),
                    (1 - t) * prev_point.z() + t * next_point.z());

                current_point = interpolated_point;
                            
                point_cloud.push_back(interpolated_point);
                v1.push_back(cycles.at(k));
            }
        }
        k++;
    }
    return std::pair<std::vector<Point>, std::vector<float> >(point_cloud, v1);
}
            

template <typename Skeleton_vertex>
struct detect_loops : public boost::dfs_visitor<> {
    detect_loops(int& n_cycles) : m_cycles(n_cycles) {}

    template < class Vertex, class Graph >
    void start_vertex(Vertex u, const Graph&) {
        m_source = u;
    }
    template < class Edge, class Graph >
    void back_edge(Edge e, const Graph& g) {
        Skeleton_vertex td = target(e, g);
        if (td == m_source) {
            m_cycles++;
        }
    }

protected:
    Skeleton_vertex m_source;
    int& m_cycles;
};

template <typename Skeleton_vertex>
struct detect_branches : public boost::dfs_visitor<> {
    detect_branches(int& branches, const std::map<Skeleton_vertex, int>& cycles) : m_branches(branches), m_cycles(cycles) {}

    template < class Vertex, class Graph >
    void start_vertex(Vertex u, const Graph&) {
        m_source = u;
    }
    template < class Edge, class Graph >
    void tree_edge(Edge e, const Graph& g) {
        Skeleton_vertex sd = source(e, g);
        Skeleton_vertex td = target(e, g);
        if (sd == m_source) {
            m_branches++;
            m_branch = true;
        }

        if (m_cycles.at(td) > 0) {
            if (m_branch) m_branches--;
            m_branch = false;
        }
    }

protected:
    Skeleton_vertex m_source;
    int& m_branches;
    bool m_branch = true;
    const std::map<Skeleton_vertex, int>& m_cycles;
};

template <typename GeomTraits, typename PointRange, typename PointMap, typename Mesh>
class Skeleton : public CGAL::Classification::Feature_base_multi_dim
{
    using MeshKernel = typename Mesh::Point::R::Kernel;
    using MeshPoint = typename Mesh::Point;
    using Mesh_to_point = CGAL::Cartesian_converter<MeshKernel, GeomTraits>;

    typedef typename GeomTraits::Point_3 Point;
    typedef std::vector<Point> Point_container;
    typedef CGAL::Identity_property_map<Point> Pmap;

    typedef CGAL::Classification::Point_set_neighborhood<GeomTraits, Point_container, Pmap> Neighborhood;
    typedef typename Neighborhood::K_neighbor_query KNeighborQuery;

    typedef CGAL::Mean_curvature_flow_skeletonization<Mesh> Skeletonization;
    typedef typename Skeletonization::Skeleton              Skele;
    typedef typename Skele::vertex_descriptor               Skeleton_vertex;
    typedef typename Skele::edge_descriptor                 Skeleton_edge;

    typedef typename boost::property_traits<boost::associative_property_map<std::map<Skeleton_vertex, boost::default_color_type> > >::value_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;

    std::vector<float> m_cycle_values;
    std::vector<float> m_branch_values;
public:
    Skeleton(const PointRange& input, PointMap point_map, const Mesh& wrap, const int num_samples = 10)
    {
        this->set_name("skele");
        m_cycle_values.resize(input.number_of_points(), 0.f);
        m_branch_values.resize(input.number_of_points(), 0.f);

        Mesh_to_point to_point;

        // Extract mean curvature flow skeleton
        Skele skeleton;
        CGAL::extract_mean_curvature_flow_skeleton(wrap, skeleton);

        // Convert polylines to point cloud  
        std::vector<std::vector<Point>> polylines;
        for (const Skeleton_edge& e : CGAL::make_range(edges(skeleton))) {
            Point s1 = to_point(skeleton[source(e, skeleton)].point);
            Point t1 = to_point(skeleton[target(e, skeleton)].point);
            polylines.push_back(std::vector<Point>({ s1,t1 }));
        }

        // Find cycles in skeleton graph
        std::map<Skeleton_edge, boost::default_color_type> edge_color;
        auto ecmap = boost::make_assoc_property_map(edge_color);
        std::map<Skeleton_vertex, boost::default_color_type> vertex_color;
        auto vcmap = boost::make_assoc_property_map(vertex_color);

        std::map<Skeleton_vertex, int> cycle_map;
        for (Skeleton_vertex v : CGAL::make_range(vertices(skeleton))) {
            int n_cycles = 0;
            detect_loops<Skeleton_vertex> vis(n_cycles);

            // mark all vertices/edges as unvisited
            typename Skele::vertex_iterator ui, ui_end;
            for (boost::tie(ui, ui_end) = vertices(skeleton); ui != ui_end; ++ui) {
                boost::put(vcmap, *ui, Color::white());
                vis.initialize_vertex(*ui, skeleton);
            }
            typename Skele::edge_iterator ei, ei_end;
            for (boost::tie(ei, ei_end) = boost::edges(skeleton); ei != ei_end; ++ei)
                boost::put(ecmap, *ei, Color::white());

            vis.start_vertex(v, skeleton);
            boost::detail::undir_dfv_impl(skeleton, v, vis, vcmap, ecmap);

            cycle_map[v] = n_cycles;
        }

        // Find branches in skeleton graph
        std::map<Skeleton_vertex, int> branch_map;
        for (Skeleton_vertex v : CGAL::make_range(vertices(skeleton))) {
            int n_branches = 0;
            detect_branches<Skeleton_vertex> vis(n_branches, cycle_map);

            // mark all vertices/edges as unvisited
            typename Skele::vertex_iterator ui, ui_end;
            for (boost::tie(ui, ui_end) = vertices(skeleton); ui != ui_end; ++ui) {
                boost::put(vcmap, *ui, Color::white());
                vis.initialize_vertex(*ui, skeleton);
            }
            typename Skele::edge_iterator ei, ei_end;
            for (boost::tie(ei, ei_end) = boost::edges(skeleton); ei != ei_end; ++ei)
                boost::put(ecmap, *ei, Color::white());

            vis.start_vertex(v, skeleton);
            boost::detail::undir_dfv_impl(skeleton, v, vis, vcmap, ecmap);

            branch_map[v] = n_branches;
        }

        std::vector<float> cycles_edges;
        for (const Skeleton_edge& e : CGAL::make_range(edges(skeleton))) {
            int s = cycle_map.at(source(e, skeleton));
            int t = cycle_map.at(target(e, skeleton));
            float a = float(s + t) / 2;
            cycles_edges.push_back(a);
        }

        std::vector<float> branches_edges;
        for (const Skeleton_edge& e : CGAL::make_range(edges(skeleton))) {
            int s = branch_map.at(source(e, skeleton));
            int t = branch_map.at(target(e, skeleton));
            float a = float(s + t) / 2;
            branches_edges.push_back(a);
        }

        std::vector<Point> pc;
        std::vector<float> cycle_vec;
        std::tie(pc, cycle_vec) = convert_polylines_to_point_cloud<Point>(polylines, num_samples, cycles_edges);
        std::vector<float> branch_vec = convert_polylines_to_point_cloud<Point>(polylines, num_samples, branches_edges).second;
        /*
        {
            typedef PointRange::Property_map<float> Fmap;
            PointRange pts;
            for (int i = 0; i < pc.size(); ++i) {
                pts.insert(pc[i]);
            }
            Fmap cyclemap = pts.add_property_map<float>("scalar_Cycle", 0).first;
            for (int i = 0; i < pts.number_of_points(); ++i) {
                cyclemap[i] = cycle_vec.at(i);
            }
            Fmap branchmap = pts.add_property_map<float>("scalar_Branch", 0).first;
            for (int i = 0; i < pts.number_of_points(); ++i) {
                branchmap[i] = branch_vec.at(i);
            }

            CGAL::IO::write_PLY("skel_points.ply", pts, CGAL::parameters::use_binary_mode(true));
        }
        */

        // Construct a tree for efficient nearest neighbor queries
        Pmap point_map2 = Pmap();
        Neighborhood neighborhood(pc, point_map2);
        KNeighborQuery neighbor_query = neighborhood.k_neighbor_query(1);
        for (std::size_t i = 0; i < input.number_of_points(); ++i) {
            std::vector<std::size_t> neighbors;

            const Point& point = get(input.point_map(), *(input.begin() + i));
            neighbor_query(point, std::back_inserter(neighbors));

            std::size_t idx = neighbors[0];
            m_cycle_values[i] = cycle_vec.at(idx);
            m_branch_values[i] = branch_vec.at(idx);
        }
    }

    float value(std::size_t pt_index, std::size_t dim_index)
    {
        if (dim_index == 0)
            return m_cycle_values[pt_index];
        else if (dim_index == 1)
            return m_branch_values[pt_index];
        return 0;
    }
};


} // namespace Feature
} // namespace Classification
} // namespace CGAL
#endif CGAL_CLASSIFICATION_FEATURE_SKELETONISE_H