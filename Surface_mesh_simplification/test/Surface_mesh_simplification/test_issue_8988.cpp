#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE 0

#include <iostream>
#include <chrono>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Face_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <iomanip>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

// Stats structure to track simplification progress
struct Stats
{
    std::size_t collected = 0;
    std::size_t processed = 0;
    std::size_t collapsed = 0;
    std::size_t non_collapsable = 0;
    std::size_t cost_uncomputable = 0;
    std::size_t placement_uncomputable = 0;
};

// Visitor class to track and display simplification progress
struct My_visitor : SMS::Edge_collapse_visitor_base<Mesh>
{
    My_visitor(Stats* s, size_t initial_face_count, size_t target_face_count, unsigned int update_frequency = 1000) :
        stats(s),
        initial_face_count(initial_face_count),
        target_face_count(target_face_count),
        total_face_count_to_process(initial_face_count - target_face_count),
        update_frequency(update_frequency),
        collect_counter(0),
        process_counter(0),
        last_update_time(std::chrono::steady_clock::now()) {}

    void OnCollected(const SMS::Edge_profile<Mesh>&, const std::optional<double>&)
    {
        ++(stats->collected);
        ++collect_counter;

        auto now = std::chrono::steady_clock::now();
        if (collect_counter >= update_frequency ||
            std::chrono::duration_cast<std::chrono::milliseconds>(now - last_update_time).count() > 250) {
            std::cout << "\rEdges collected: " << stats->collected << std::flush;
            collect_counter = 0;
            last_update_time = now;
        }
    }

    void OnSelected(const SMS::Edge_profile<Mesh>& mesh,
                    std::optional<double> cost,
                    std::size_t,
                    std::size_t)
    {
        if (stats->processed == 0) {
            std::cout << "\n"
                      << std::flush;
        }
        ++(stats->processed);
        ++process_counter;
        if (!cost)
            ++(stats->cost_uncomputable);
        auto now = std::chrono::steady_clock::now();
        if (process_counter >= update_frequency ||
            std::chrono::duration_cast<std::chrono::milliseconds>(now - last_update_time).count() > 250) {
            size_t current_face_count = mesh.surface_mesh().number_of_faces();
            double processed_ratio = static_cast<double>(initial_face_count - current_face_count) / total_face_count_to_process;
            std::cout << "\rFaces left: " << current_face_count << " ("
                      << std::fixed << std::setprecision(2) << processed_ratio * 100.0 << "%)"
                      << std::flush;
            process_counter = 0;
            last_update_time = now;
        }
    }

    void OnCollapsing(const SMS::Edge_profile<Mesh>&,
                      std::optional<Point_3> placement)
    {
        if (!placement)
            ++(stats->placement_uncomputable);
    }

    void OnNonCollapsable(const SMS::Edge_profile<Mesh>&)
    {
        ++(stats->non_collapsable);
    }

    void OnCollapsed(const SMS::Edge_profile<Mesh>&, Mesh::Vertex_index)
    {
        ++(stats->collapsed);
    }

    Stats* stats;
    size_t initial_face_count;
    size_t target_face_count;
    size_t total_face_count_to_process;
    unsigned int update_frequency;
    unsigned int collect_counter;
    unsigned int process_counter;
    std::chrono::steady_clock::time_point last_update_time;
};

// Function to create a large planar mesh
Mesh create_planar_mesh(int rows, int cols)
{
    Mesh mesh;
    std::vector<Mesh::Vertex_index> vertices;
    for (int i = 0; i <= rows; ++i) {
        for (int j = 0; j <= cols; ++j) {
            vertices.push_back(mesh.add_vertex(Point_3(j, i, 0)));
        }
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            Mesh::Vertex_index v1 = vertices[i * (cols + 1) + j];
            Mesh::Vertex_index v2 = vertices[i * (cols + 1) + j + 1];
            Mesh::Vertex_index v3 = vertices[(i + 1) * (cols + 1) + j];
            Mesh::Vertex_index v4 = vertices[(i + 1) * (cols + 1) + j + 1];
            mesh.add_face(v1, v2, v4);
            mesh.add_face(v1, v4, v3);
        }
    }
    return mesh;
}

int main()
{
    int rows = 1500;
    int cols = 1500;
    Mesh mesh = create_planar_mesh(rows, cols);
    CGAL::IO::write_polygon_mesh("input_"+std::to_string(cols)+".off", mesh, CGAL::parameters::stream_precision(17));
    std::cout << "Initial mesh: " << mesh.number_of_faces() << " faces." << std::endl;

    size_t target_face_count = 15000;
    SMS::Face_count_stop_predicate<Mesh> stop(target_face_count);

    // Apply simplification using Garland-Heckbert policies
    std::cout << "Starting simplification with Garland-Heckbert policy..." << std::endl;
    auto t_start = std::chrono::steady_clock::now();

    Stats stats;
    unsigned int update_frequency = std::max(1000u, static_cast<unsigned int>(mesh.number_of_edges() / 20));
    size_t initial_face_count = mesh.number_of_faces();
    My_visitor visitor(&stats, initial_face_count, target_face_count, update_frequency);

    typedef SMS::GarlandHeckbert_policies<Mesh, Kernel> GHPolicies;
    GHPolicies gh_policies(mesh);
    auto gh_cost = gh_policies.get_cost();
    auto gh_placement = gh_policies.get_placement();

    int r = SMS::edge_collapse(mesh,
                               stop,
                               CGAL::parameters::get_cost(gh_cost)
                                   .get_placement(gh_placement)
                                   .visitor(visitor));

    auto t_end = std::chrono::steady_clock::now();
    double elapsed_sec = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();

    std::cout << "\nSimplification finished." << std::endl;
    std::cout << "Edges removed: " << r << std::endl;
    std::cout << "Final mesh: " << mesh.number_of_faces() << " faces." << std::endl;
    std::cout << "Execution time: " << elapsed_sec << " seconds." << std::endl;

    // Display statistics
    std::cout << "\nEdges collected: " << stats.collected
              << "\nEdges processed: " << stats.processed
              << "\nEdges collapsed: " << stats.collapsed
              << std::endl
              << "\nEdges not collapsed due to topological constraints: " << stats.non_collapsable
              << "\nEdge not collapsed due to cost computation constraints: " << stats.cost_uncomputable
              << "\nEdge not collapsed due to placement computation constraints: " << stats.placement_uncomputable
              << std::endl;

    CGAL::IO::write_polygon_mesh("out_"+std::to_string(cols)+".off", mesh, CGAL::parameters::stream_precision(17));

    return 0;
}