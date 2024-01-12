

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#define LOG() std::cout<<__FILE__<<":"<<__LINE__<< std::endl;

#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Classification/ETHZ/Random_forest_classifier.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Classification/Point_set_neighborhood.h>
#include <CGAL/IO/write_ply_points.h>
#include <utility>

#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

// Boost directed graph cycle detection example https://www.boost.org/doc/libs/1_82_0/libs/graph/doc/file_dependency_example.html#sec:cycles
        // Stack Overflow BGL directed graph cycle detection example https://stackoverflow.com/questions/41027377/does-undirected-dfs-detect-all-cycle-of-graph
#include <boost/graph/undirected_dfs.hpp>

namespace CGAL {
namespace   Classification {
namespace     Feature {


class Skeleton_branches : public CGAL::Classification::Feature_base {
    std::vector<int> s0;
public:
    Skeleton_branches(std::vector<int>& s0) : s0(s0) {
        this->set_name("skele_branches");
    }

    float value(std::size_t pt_index) {
        return (float)s0[pt_index];
    }
};

class Skeleton_cycles : public CGAL::Classification::Feature_base {
    std::vector<int> s1;
public:
    Skeleton_cycles(std::vector<int>& s1) : s1(s1) {
        this->set_name("skele_cycles");
    }

    float value(std::size_t pt_index) {
        return (float)s1[pt_index];
    }
};


template <typename Point>
std::pair<std::vector<Point>, std::vector<float> > convert_polylines_to_point_cloud(const std::vector<std::vector<Point>>& polylines, int num_samples,const std::vector<float>& cycles) {
    std::vector<Point> point_cloud;
    std::vector<float> v1;
    int k = 0;
    for (const auto& polyline : polylines) {
        int num_points = polyline.size();

        if (num_points < 2) {
            continue;
        }

        // Parcourir les segments de la polyline
        for (int i = 0; i < num_points - 1; i++) {
            const Point& prev_point = polyline[i];
            const Point& next_point = polyline[i + 1];

            // Échantillonner des points le long du segment et calculer les distances dans le point cloud
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
            

struct detect_loops : public boost::dfs_visitor<> {
    detect_loops(int& n_cycles) : _n_cycles(n_cycles) {}

    template < class Vertex, class Graph >
    void start_vertex(Vertex u, const Graph&) {
        _source = u;
    }
    //template < class Edge, class Graph > void tree_edge(Edge e, const Graph& g) { std::cout << "    (te) " << source(e, g) << " -- " << target(e, g) << std::endl; }
    template < class Edge, class Graph >
    void back_edge(Edge e, const Graph& g) {
        //std::cout << "    (be) " << source(e, g) << " -- " << target(e, g) << std::endl;
        int _target = target(e, g);
        if (_target == _source) {
            _n_cycles++;
        }
    }

protected:
    int _source;
    int& _n_cycles;
};

template <typename Skeleton_vertex>
struct detect_branches : public boost::dfs_visitor<> {
    detect_branches(int& branches, const std::map<Skeleton_vertex, int>& cycles) : _branches(branches), _cycles(cycles) {}

    template < class Vertex, class Graph >
    void start_vertex(Vertex u, const Graph&) {
        _source = u;
    }
    template < class Edge, class Graph >
    void tree_edge(Edge e, const Graph& g) {
        Skeleton_vertex sd = source(e, g);
        Skeleton_vertex td = target(e, g);
        //std::cout << "    (te) " << source(e, g) << " -- " << target(e, g) << " cycles(s): " << _cycles.at(sd) << " cycles(t): " << _cycles.at(td) << ", branch: " << _branch << ", branches: " << _branches << std::endl;
        if (_cycles.at(td) > 0 && _cycles.at(sd) > 0) {
            _branches = false;
        }
        if (_branch) {
            if (sd == _source) {
                _branches++;
            }
        }
    }

protected:
    int _source;
    int& _branches;
    bool _branch = true;
    const std::map<Skeleton_vertex, int>& _cycles;
};


template <typename GeomTraits, typename PointRange, typename PointMap>
std::tuple< std::vector<int>, std::vector<int> > compute_skeleton_features(const PointRange& input, PointMap point_map, const Mesh& wrap) {
    using Image_float = CGAL::Classification::Image<float>;

    using FloatMap = typename PointRange::template Property_map<float>;
    std::vector<typename FloatMap::value_type> values;

    typedef GeomTraits::Point_3 Point;
    typedef std::vector<Point> Point_container;
    typedef CGAL::Identity_property_map<Point> Pmap;

    namespace Classification = CGAL::Classification;
    typedef CGAL::Classification::Point_set_neighborhood<GeomTraits, Point_container, Pmap> Neighborhood;
    typedef Neighborhood::K_neighbor_query KNeighborQuery;

    typedef CGAL::Mean_curvature_flow_skeletonization<Mesh> Skeletonization;
    typedef Skeletonization::Skeleton                       Skele;
    typedef Skele::vertex_descriptor                     Skeleton_vertex;
    typedef Skele::edge_descriptor                       Skeleton_edge;

                
    // Step 2: Extract mean curvature flow skeleton
    // Skeleton_graph skeleton_graph = CGAL::extract_mean_curvature_flow_skeleton(wrap);

    Skele skeleton;
    CGAL::extract_mean_curvature_flow_skeleton(wrap, skeleton);
    std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
    std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";


    // Step 3: Split the skeleton into polylines
    //Skeleton_polylines skeleton_polylines = CGAL::split_graph_into_polylines(skeleton_graph);

    // Step 4: Convert polylines to point cloud  
        //adjacent_vertices()
    std::vector<std::vector<Point>> polylines;
    for (Skeleton_edge e : CGAL::make_range(edges(skeleton))) {
        Point_Alpha s = skeleton[source(e, skeleton)].point;
        Point_Alpha t = skeleton[target(e, skeleton)].point;
        //std::cout << "2 " << s << " " << t << "\n";
        Point s1(s.x(), s.y(), s.z());
        Point t1(t.x(), t.y(), t.z());
        polylines.push_back(std::vector<Point>({ s1,t1 }));
    }

    // find cycles in skeleton graph
    std::map<Skeleton_edge, boost::default_color_type> edge_color;
    auto ecmap = boost::make_assoc_property_map(edge_color);
    std::map<Skeleton_vertex, boost::default_color_type> vertex_color;
    auto vcmap = boost::make_assoc_property_map(vertex_color);

    typedef typename boost::property_traits<boost::associative_property_map<std::map<Skeleton_vertex, boost::default_color_type> > >::value_type ColorValue;
    typedef boost::color_traits< ColorValue > Color;

    std::map<Skeleton_vertex, int> cycle_map;
    for (Skeleton_vertex v : CGAL::make_range(vertices(skeleton))) {
        int n_cycles = 0;
        detect_loops vis(n_cycles);

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
        //boost::undirected_dfs(skeleton, boost::root_vertex(v).visitor(vis).edge_color_map(ecmap).vertex_color_map(vcmap));

        //std::cout << "v: " << v << " no. of cycles: " << n_cycles << std::endl;
        cycle_map[v] = n_cycles;
    }

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

        //std::cout << "(branch alg) Exporing vertex: " << v << std::endl;
        boost::detail::undir_dfv_impl(skeleton, v, vis, vcmap, ecmap);

        //std::cout << "v: " << v << " no. of branches: " << n_branches << std::endl;
        branch_map[v] = n_branches;
    }

    {
        typedef PointRange::Property_map<float> Fmap;
        PointRange pts1;
        for (Skeleton_vertex v : CGAL::make_range(vertices(skeleton)))
            pts1.insert(Point(skeleton[v].point.x(), skeleton[v].point.y(), skeleton[v].point.z()));
        Fmap cyclemap = pts1.add_property_map<float>("scalar_Cycle", 0).first;
        for (int i = 0; i < pts1.number_of_points(); ++i) {
            cyclemap[i] = (float)cycle_map.at(i);
        }
        Fmap branchmap = pts1.add_property_map<float>("scalar_Branch", 0).first;
        for (int i = 0; i < pts1.number_of_points(); ++i) {
            branchmap[i] = (float)branch_map.at(i);
        }

        std::ofstream f("skel_points123.ply");
        f.precision(18);
        f << pts1;
    }

    std::vector<float> cycles_edges;
    for (Skeleton_edge e : CGAL::make_range(edges(skeleton))) {
        int s = cycle_map.at(source(e, skeleton));
        int t = cycle_map.at(target(e, skeleton));
        float a = float(s + t) / 2;
        cycles_edges.push_back(a);
    }

    std::vector<float> branches_edges;
    for (Skeleton_edge e : CGAL::make_range(edges(skeleton))) {
        int s = branch_map.at(source(e, skeleton));
        int t = branch_map.at(target(e, skeleton));
        float a = float(s + t) / 2;
        branches_edges.push_back(a);
    }

    int num_samples = 10;
    //std::vector<Point> pc = convert_polylines_to_point_cloud<Point>(polylines, num_samples);
    std::pair<std::vector<Point>, std::vector<float> > pc_cycle_pair = convert_polylines_to_point_cloud<Point>(polylines, num_samples, branches_edges);
    std::pair<std::vector<Point>, std::vector<float> > pc_cycle_pair2 = convert_polylines_to_point_cloud<Point>(polylines, num_samples, cycles_edges);
    std::vector<Point> pc = pc_cycle_pair.first;
    std::vector<float> cycle_vec = pc_cycle_pair2.second;
    std::vector<float> branch_vec = pc_cycle_pair.second;
    std::cout << " conversion done" << std::endl;
    // Step 5: Construct a tree for efficient nearest neighbor queries
        // Tree tree = construct_tree(skeleton_points);


    std::vector<int> cycle_values;
    std::vector<int> branch_values;
    cycle_values.resize(input.size(), 0.f);
    branch_values.resize(input.size(), 0.f);

    Pmap point_map2 = Pmap();
    Neighborhood neighborhood(pc, point_map2);
    int n = 1; // number of neighbours
    KNeighborQuery neighbor_query = neighborhood.k_neighbor_query(n);
    for (std::size_t i = 0; i < input.number_of_points(); ++i) {
        std::vector<std::size_t> neighbors;
                            
        const Point& point = get(input.point_map(), *(input.begin() + i));
        neighbor_query(point, std::back_inserter(neighbors));
        //Point& closest_point = get(point_map2, *(pc.begin() + neighbors[0]));
        //float d = CGAL::sqrt((closest_point - point).squared_length());

        size_t idx = neighbors[0];
        cycle_values[i] = cycle_vec.at(idx);
        branch_values[i] = branch_vec.at(idx);
    }

    return std::make_tuple(cycle_values,branch_values);
}

        
}}}

