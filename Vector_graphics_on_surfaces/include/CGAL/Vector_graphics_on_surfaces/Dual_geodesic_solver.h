// Copyright (c) 2023-2026 GeometryFactory and Claudio Mancinelli.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli and Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_DUAL_GEODESIC_SOLVER_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_DUAL_GEODESIC_SOLVER_H

#include <CGAL/license/Vector_graphics_on_surfaces.h>

#include <CGAL/Vector_graphics_on_surfaces/internal/utils.h>

#ifdef CGAL_BSURF_USE_DIJKSTRA_SP
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <CGAL/boost/graph/Dual.h>
#endif

namespace CGAL {
namespace Vector_graphics_on_surfaces {

#ifndef DOXYGEN_RUNNING
template <class FT>
struct Dual_geodesic_solver;
#endif

namespace internal {

template <class K, class TriangleMesh, class VertexPointMap, class VertexIndexMap, class FaceIndexMap>
struct Geodesic_circle_impl
{
  using face_descriptor =
      typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using vertex_descriptor =
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor =
      typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  using Point_2 = typename K::Point_2;
  using Point_3 = typename K::Point_3;
  using Vector_2 = typename K::Vector_2;
  using Vector_3 = typename K::Vector_3;
  using FT = typename K::FT;

  using Face_location = Polygon_mesh_processing::Face_location<TriangleMesh, FT>;
  using Edge_location = Polygon_mesh_processing::Edge_location<TriangleMesh, FT>;

  struct geodesic_solver {
    struct graph_edge {
      int node = -1;
      double len=DBL_MAX;
    };
    std::vector<std::vector<graph_edge>> graph;
  };

  static
  void connect_nodes(geodesic_solver &solver,
                    const vertex_descriptor& a,
                    const vertex_descriptor& b,
                    const VertexIndexMap& vidmap,
                    const FT& len)
  {
    // TODO: avoid cast
    unsigned int vida=get(vidmap,a);
    unsigned int vidb=get(vidmap,b);
    solver.graph[vida].push_back({static_cast<int>(vidb), len});
    solver.graph[vidb].push_back({static_cast<int>(vida), len});

  }

  static
  double opposite_nodes_arc_length(const VertexPointMap &vpm,
                                   const vertex_descriptor& a,
                                   const vertex_descriptor& c,
                                   const vertex_descriptor& b,
                                   const vertex_descriptor& d)
  {
    // Triangles (a, b, d) and (b, d, c) are connected by (b, d) edge
    // Nodes a and c must be connected.

    Vector_3 ba = get(vpm,a) - get(vpm,b);
    Vector_3 bc = get(vpm,c) - get(vpm,b);
    Vector_3 bd = get(vpm,d) - get(vpm,b);

    Vector_3 ba_norm=ba/sqrt(ba.squared_length());
    Vector_3 bd_norm=bd/sqrt(bd.squared_length());
    Vector_3 bc_norm=bc/sqrt(bc.squared_length());

    double cos_alpha = ba_norm * bd_norm;
    double cos_beta  = bc_norm * bd_norm;
    double sin_alpha = sqrt((std::max)(0.0, 1 - cos_alpha * cos_alpha));
    double sin_beta  = sqrt((std::max)(0.0, 1 - cos_beta * cos_beta));


    // cos(alpha + beta)
    double cos_alpha_beta = cos_alpha * cos_beta - sin_alpha * sin_beta;
    CGAL_assertion(cos_alpha_beta>-1);
    if (cos_alpha_beta <= -1) return DBL_MAX;

    // law of cosines (generalized Pythagorean theorem)
    double len = ba.squared_length() + bc.squared_length() -
                 sqrt(ba.squared_length()) * sqrt(bc.squared_length()) * 2 * cos_alpha_beta;

     CGAL_assertion(len>0);
    if (len <= 0)
      return DBL_MAX;
    else
      return sqrt(len);
  }

  static
  void connect_opposite_nodes(geodesic_solver& solver,
                              const VertexPointMap &vpm,
                              const TriangleMesh &mesh,
                              const VertexIndexMap& vidmap,
                              const vertex_descriptor& a,
                              const vertex_descriptor& b,
                              const halfedge_descriptor& h)
  {
    vertex_descriptor v0 =target(next(h,mesh),mesh);
    vertex_descriptor v1= target(next(opposite(h,mesh),mesh),mesh);
    auto len = opposite_nodes_arc_length(vpm, v0,v1,a,b);
     //std::cout<<len<<std::endl;
    connect_nodes(solver, v0, v1, vidmap, len);
  }

  static
  geodesic_solver
  make_geodesic_solver(const VertexPointMap &vpm,
                       const VertexIndexMap& vidmap,
                       const TriangleMesh &mesh)
  {
    geodesic_solver solver;
    solver.graph.resize(vertices(mesh).size());
    for (face_descriptor f : faces(mesh)) {
      halfedge_descriptor h=halfedge(f,mesh);
      for(std::size_t i=0;i<3;++i)
    {
      vertex_descriptor a=source(h,mesh);
      vertex_descriptor b=target(h,mesh);

      if(a<b)
        {
          FT len=sqrt(squared_distance(get(vpm,a),get(vpm,b)));
          //std::cout<<len<<std::endl;
          CGAL_assertion(len<DBL_MAX && len>0);
          connect_nodes(solver,a,b,vidmap,len);
        }
        face_descriptor nei=face(opposite(h,mesh),mesh);
        if(f<nei)
            connect_opposite_nodes(solver,vpm,mesh,vidmap,a,b,h);

        h=next(h,mesh);


      }


    }
    return solver;
  }
  // geodesic_solver
  // extended_geodesic_solver(const VertexPointMap& vpm,
  //                          const VertexIndexMap& vidmap,
  //                          const TriangleMesh& mesh,
  //                          const dual_geodesic_solver& dual_solver,
  //                          const int k)
  // {
  //   geodesic_solver solver;
  //   solver.graph.resize(vertices(mesh).size());
  //   for(auto& v:vertices(mesh))
  //   {

  //   }
  // }
  static
  Dual_geodesic_solver<FT>
  make_dual_geodesic_solver(const VertexPointMap &vpm,
                            const FaceIndexMap& tidmap,
                            const TriangleMesh &mesh)
  {
    auto compute_dual_weights=[&mesh,&vpm](const halfedge_descriptor& h)
    {
      std::array<Vector_2,3> flat_tid = internal::init_flat_triangle<K>(h,vpm,mesh);
      std::array<Vector_2,3> flat_nei = internal::unfold_face<K>(h,vpm,mesh,flat_tid);

      Vector_2 c0=0.33*(flat_tid[0]+flat_tid[1]+flat_tid[2]);
      Vector_2 c1=0.33*(flat_nei[0]+flat_nei[1]+flat_nei[2]);

      return sqrt((c1 - c0).squared_length());
    };

    Dual_geodesic_solver<FT> solver;
    solver.graph.resize(faces(mesh).size());
    for (auto f : faces(mesh)) {
      halfedge_descriptor h=halfedge(f,mesh);
      int entry=get(tidmap,f);
      for (auto i = 0; i < 3; ++i) {
        solver.graph[entry][i].node = get(tidmap,face(opposite(h,mesh),mesh));
        solver.graph[entry][i].len = compute_dual_weights(h);
        h=next(h,mesh);
      }
    }
    return solver;
  }

  // `update` is a function that is executed during expansion, every time a node
  // is put into queue. `exit` is a function that tells whether to expand the
  // current node or perform early exit.
  template <typename Update, typename Stop, typename Exit>
  static
  void visit_geodesic_graph(std::vector<double> &field, const geodesic_solver &solver,
                            const std::vector<int> &sources, Update &&update,
                            Stop &&stop, Exit &&exit)
  {
    /*
      This algorithm uses the heuristic Small Label First and Large Label Last
      https://en.wikipedia.org/wiki/Shortest_Path_Faster_Algorithm

      Large Label Last (LLL): When extracting nodes from the queue, pick the
      front one. If it weights more than the average weight of the queue, put
      on the back and check the next node. Continue this way.
      Sometimes average_weight is less than every value due to floating point
      errors (doesn't happen with double precision).

      Small Label First (SLF): When adding a new node to queue, instead of
      always pushing it to the end of the queue, if it weights less than the
      front node of the queue, it is put on front. Otherwise the node is put at
      the end of the queue.
    */

    auto in_queue = std::vector<bool>(solver.graph.size(), false);

    // Cumulative weights of elements in queue. Used to keep track of the
    // average weight of the queue.
    auto cumulative_weight = 0.0;

    // setup queue
    auto queue = std::deque<int>();
    for (auto source : sources) {
      in_queue[source] = true;
      cumulative_weight += field[source];
      queue.push_back(source);
    }

    while (!queue.empty()) {
      auto node = queue.front();
      auto average_weight = (cumulative_weight / queue.size());

      // Large Label Last (see comment at the beginning)
      for (std::size_t tries = 0; tries < queue.size() + 1; tries++) {
        if (field[node] <= average_weight)
          break;
        queue.pop_front();
        queue.push_back(node);
        node = queue.front();
      }

      // Remove node from queue.
      queue.pop_front();
      in_queue[node] = false;
      cumulative_weight -= field[node];

      // Check early exit condition.
      if (exit(node))
        break;
      if (stop(node))
        continue;

      for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
        // Distance of neighbor through this node
        double new_distance = field[node] + solver.graph[node][i].len;

        auto neighbor = solver.graph[node][i].node;

        auto old_distance = field[neighbor];
        if (new_distance >= old_distance)
          continue;

        if (in_queue[neighbor]) {
          // If neighbor already in queue, don't add it.
          // Just update cumulative weight.
          cumulative_weight += new_distance - old_distance;
        } else {
          // If neighbor not in queue, add node to queue using Small Label
          // First (see comment at the beginning).
          if (queue.empty() || (new_distance < field[queue.front()]))
            queue.push_front(neighbor);
          else
            queue.push_back(neighbor);

          // Update queue information.
          in_queue[neighbor] = true;
          cumulative_weight += new_distance;
        }

        // Update distance of neighbor.
        field[neighbor] = new_distance;
        update(node, neighbor, new_distance);
      }
    }
  }

  template <typename Update, typename Stop, typename Exit, typename FT>
  static
  void visit_dual_geodesic_graph(std::vector<double> &field,
                                const Dual_geodesic_solver<FT> &solver,
                                const std::vector<int> &sources,
                                Update &&update,
                                Stop &&stop,
                                Exit &&exit)
  {
    /*
      This algorithm uses the heuristic Small Label First and Large Label Last
      https://en.wikipedia.org/wiki/Shortest_Path_Faster_Algorithm

      Large Label Last (LLL): When extracting nodes from the queue, pick the
      front one. If it weights more than the average weight of the queue, put
      on the back and check the next node. Continue this way.
      Sometimes average_weight is less than every value due to floating point
      errors (doesn't happen with double precision).

      Small Label First (SLF): When adding a new node to queue, instead of
      always pushing it to the end of the queue, if it weights less than the
      front node of the queue, it is put on front. Otherwise the node is put at
      the end of the queue.
    */

    auto in_queue = std::vector<bool>(solver.graph.size(), false);

    // Cumulative weights of elements in queue. Used to keep track of the
    // average weight of the queue.
    auto cumulative_weight = 0.0;

    // setup queue
    auto queue = std::deque<int>();
    for (auto source : sources) {
      in_queue[source] = true;
      cumulative_weight += field[source];
      queue.push_back(source);
    }

    while (!queue.empty()) {
      auto node = queue.front();
      auto average_weight = (cumulative_weight / queue.size());

      // Large Label Last (see comment at the beginning)
      for (std::size_t tries = 0; tries < queue.size() + 1; tries++) {
        if (field[node] <= average_weight)
          break;
        queue.pop_front();
        queue.push_back(node);
        node = queue.front();
      }

      // Remove node from queue.
      queue.pop_front();
      in_queue[node] = false;
      cumulative_weight -= field[node];

      // Check early exit condition.
      if (exit(node))
        break;
      if (stop(node))
        continue;

      for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
        // Distance of neighbor through this node
        auto new_distance = field[node] + solver.graph[node][i].len;
        auto neighbor = solver.graph[node][i].node;

        auto old_distance = field[neighbor];
        if (new_distance >= old_distance)
          continue;

        if (in_queue[neighbor]) {
          // If neighbor already in queue, don't add it.
          // Just update cumulative weight.
          cumulative_weight += new_distance - old_distance;
        } else {
          // If neighbor not in queue, add node to queue using Small Label
          // First (see comment at the beginning).
          if (queue.empty() || (new_distance < field[queue.front()]))
            queue.push_front(neighbor);
          else
            queue.push_back(neighbor);

          // Update queue information.
          in_queue[neighbor] = true;
          cumulative_weight += new_distance;
        }

        // Update distance of neighbor.
        field[neighbor] = new_distance;
        update(node, neighbor, new_distance);
      }
    }
  }

  static
  std::vector<double>
  compute_geodesic_distances(const geodesic_solver &solver,
                             const VertexIndexMap& vidmap,
                             const std::vector<std::pair<vertex_descriptor, double>> &sources_and_dist)
  {
    auto update = [](int /* node */, int /* neighbor */, double /* new_distance */) {};
    auto stop = [](int /* node */) { return false; };
    auto exit = [](int /* node */) { return false; };

    auto distances = std::vector<double>{};
    distances.assign(solver.graph.size(), DBL_MAX);
    std::vector<int>sources_id((sources_and_dist.size()));
    for (std::size_t i = 0; i < sources_and_dist.size(); ++i) {
      sources_id[i] = get(vidmap,sources_and_dist[i].first);
      distances[sources_id[i]] = sources_and_dist[i].second;
    }
    visit_geodesic_graph(distances, solver, sources_id, update, stop, exit);

    return distances;
  }

  static
  std::vector<double> solve_with_targets(const geodesic_solver& solver,
                                        const VertexIndexMap& vidmap,
                                        const std::vector<std::pair<vertex_descriptor, double>> &sources_and_dist,
                                        const std::vector<std::pair<vertex_descriptor, double>> &targets_and_dist)
  {
    auto update       = [](int node, int neighbor, float new_distance) {};
    auto stop         = [](int node) { return false; };
    double max_distance = DBL_MAX;
    std::vector<int> exit_verts(targets_and_dist.size());
    for (auto i = 0; i < targets_and_dist.size(); ++i) {
        exit_verts[i]=get(vidmap,targets_and_dist[i].first);
    }
    auto exit = [&exit_verts](int node) {
      auto it = find(exit_verts.begin(), exit_verts.end(), node);
      if (it != exit_verts.end())
        exit_verts.erase(it);

    if (exit_verts.empty())
        return true;

      return false;
    };

    auto distances  = std::vector<double>(solver.graph.size(), DBL_MAX);
    std::vector<int>sources_id((sources_and_dist.size()));
    for (auto i = 0; i < sources_and_dist.size(); ++i) {
      sources_id[i] = get(vidmap,sources_and_dist[i].first);
      distances[sources_id[i]] = sources_and_dist[i].second;
    }

    visit_geodesic_graph(distances, solver, sources_id, update, stop, exit);
    return distances;
  }

  //compute the length between two opposite vertices by flattening
  //TODO: handle concave configurations
  static
  double length_by_flattening(const VertexPointMap &vpm,
                            const TriangleMesh &mesh,
                            const halfedge_descriptor& h)
  {
    std::array<Vector_2,3> flat_tid = internal::init_flat_triangle<K>(h,vpm,mesh);
    std::array<Vector_2,3> flat_nei = internal::unfold_face<K>(h,vpm,mesh,flat_tid);
    return sqrt((flat_tid[2]-flat_nei[2]).squared_length());
  }

  static
  Point_3 eval_position(const VertexPointMap &vpm,
                        const TriangleMesh &mesh,
                        const Face_location& p)
  {
    halfedge_descriptor h=halfedge(p.first,mesh);
    return CGAL::barycenter(get(vpm, source(h,mesh)), p.second[0], get(vpm, target(h,mesh)), p.second[1], get(vpm, target(next(h,mesh),mesh)), p.second[2]);
  }

  // compute the distance between a point p and some vertices around him
  // TODO: consider to take more vertices (increase accuracy)
  // TODO: handle concave configurations
  static
  std::vector<std::pair<vertex_descriptor, double>>
  nodes_around_point(const VertexPointMap &vpm,
                    const TriangleMesh &mesh,
                    const Face_location& p)
  {
    auto get_vid=[&mesh](const int k,const face_descriptor& tid)
    {
      halfedge_descriptor h=halfedge(tid,mesh);
      switch(k)
      {
        case 0:
          return source(h,mesh);
        case 1:
          return target(h,mesh);
        default:
          return target(next(h,mesh),mesh);
      }
    };

    std::vector<std::pair<vertex_descriptor, double>> nodes;
    nodes.reserve(6);
    auto [is_vert,offset]=internal::point_is_vert<K>(p);
    if (is_vert) {
      vertex_descriptor vid = get_vid(offset,p.first);
      nodes.push_back({vid, 0});
    } else {
      face_descriptor tid = p.first;
      Point_3 pos = eval_position(vpm, mesh, p);
      halfedge_descriptor h=halfedge(tid,mesh);
      for (auto i = 0; i < 3; ++i) {
        vertex_descriptor p0 = source(h,mesh);
        //connect to current vertex
        double d = sqrt(squared_distance(get(vpm,p0),pos));
        nodes.push_back(std::make_pair(p0, d));

        //connecting to opposite vertex w.r.t to current halfedge
        vertex_descriptor opp = target(next(opposite(h,mesh),mesh),mesh);
        double l =
            length_by_flattening(vpm,mesh,h);
        nodes.push_back(std::make_pair(opp, l));
        h=next(h,mesh);
      }
    }

    return nodes;
  }

  //compute geodesic distance field from p
  //TODO: can be easily extended to more than one source
  static
  std::vector<double>
  compute_geodesic_distances(const geodesic_solver& solver,
                             const VertexPointMap& vpm,
                             const VertexIndexMap& vim,
                             const TriangleMesh &mesh,
                             const Face_location& p)
  {
    std::vector<std::pair<vertex_descriptor,double>> source_nodes=nodes_around_point(vpm,mesh,p);

    return compute_geodesic_distances(solver, vim, source_nodes);
  }

  //compute the geodesic distance field from src, and stop the propagation
  //once tgt is reached
  static
  std::vector<double>
  compute_pruned_geodesic_distances(const geodesic_solver &solver,
                                    const VertexPointMap &vpm,
                                    const VertexIndexMap& vidmap,
                                    const TriangleMesh &mesh,
                                    const Face_location& src,
                                    const Face_location& tgt)
  {
    std::vector<std::pair<vertex_descriptor, double>>
    source_nodes = nodes_around_point(vpm,mesh,src);

    std::vector<std::pair<vertex_descriptor, double>>
    target_nodes = nodes_around_point(vpm,mesh,tgt);

    return solve_with_targets(solver, source_nodes, target_nodes);
  }

  template <class FT>
  static
  std::vector<halfedge_descriptor>
  strip_on_dual_graph(const Dual_geodesic_solver<FT> &solver,
                      const TriangleMesh &mesh,
                      const int src,
                      const int tgt)
  {
    if (src == tgt)
      return {};

    auto common_halfedge = [&mesh](face_descriptor f1, face_descriptor f2) {
      halfedge_descriptor h = halfedge(f1, mesh);
      for (int i = 0; i < 3; ++i) {
        if (face(opposite(h, mesh), mesh) == f2)
          return h;
        h = next(h, mesh);
      }
      CGAL_assertion(!"faces do no share a common edge");
      return halfedge_descriptor();
    };
    // initialize once for all and sparsely cleanup at the end of every solve
    std::vector<int> parents(solver.graph.size(), -1);
    std::vector<double> field(solver.graph.size(), DBL_MAX);
    std::vector<face_descriptor> id_to_face_map(faces(mesh).begin(), faces(mesh).end());

    field[src]=0.0;
    std::vector<int> sources = {src};
    auto update = [&parents](int node, int neighbor, double) {
      parents[neighbor] = node;
    };
    auto stop = [](int) { return false; };
    auto exit = [tgt](int node) { return node==tgt; };

    visit_dual_geodesic_graph(field,solver, sources, update, stop, exit);

    // extract_strip
    std::vector<halfedge_descriptor> strip;
    int node = tgt;
    CGAL_assertion(parents[tgt] != -1);
    //update the result using id_to_face_map
    strip.reserve((int)std::sqrt(parents.size()));
    while (parents[node] != -1) {
      strip.push_back(common_halfedge(id_to_face_map[node],id_to_face_map[parents[node]]));
      node = parents[node];
    }
    std::reverse(strip.begin(),strip.end());
    return strip;
  }
};

#ifdef CGAL_BSURF_USE_DIJKSTRA_SP
  class Dijkstra_end_exception : public std::exception
  {
    const char* what() const throw ()
    {
      return "Dijkstra shortest path: reached the target vertex.";
    }
  };

  template <class Graph>
  class Stop_at_target_Dijkstra_visitor : boost::default_dijkstra_visitor
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

    vertex_descriptor destination_vd;

  public:
    Stop_at_target_Dijkstra_visitor(vertex_descriptor destination_vd)
      : destination_vd(destination_vd)
    { }

    void initialize_vertex(const vertex_descriptor& /*s*/, const Graph& /*mesh*/) const { }
    void examine_vertex(const vertex_descriptor& /*s*/, const Graph& /*mesh*/) const { }
    void examine_edge(const edge_descriptor& /*e*/, const Graph& /*mesh*/) const { }
    void edge_relaxed(const edge_descriptor& /*e*/, const Graph& /*mesh*/) const { }
    void discover_vertex(const vertex_descriptor& /*s*/, const Graph& /*mesh*/) const { }
    void edge_not_relaxed(const edge_descriptor& /*e*/, const Graph& /*mesh*/) const { }
    void finish_vertex(const vertex_descriptor &vd, const Graph& /* mesh*/) const
    {
      if(vd == destination_vd)
        throw Dijkstra_end_exception();
    }
  };
#endif

} // end of internal namespace

/*!
 * \ingroup VGSMiscellaneous
 * Geodesic solver class used to store precomputed information of a given mesh for
 * approximate geodesic computation. Those information are computed by the function
 * `init_geodesic_dual_solver()`.
 * \tparam FT floating point number type (float or double)
 */
template <class FT>
struct Dual_geodesic_solver
{
  struct Edge {
    int node = -1;
    FT len = DBL_MAX;
  };
  std::vector<std::array<Edge, 3>> graph = {};
};

/*!
 * \ingroup VGSMiscellaneous
 * fills `solver` for a given mesh `tmesh`. It is the user responsibility to
 * call again this function if `tmesh` or the points of its vertices are modified.
 * If `solver` was used in a previous call to this function, information will be overwritten.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam FT floating point number type (float or double)
 * \param solver the container for the precomputed information
 * \param tmesh triangle mesh to be considered for the precomputations
 * \todo add named parameters
 * \todo make sure solver.graph is cleared before filling it
 */
template <class FT, class TriangleMesh>
void init_geodesic_dual_solver(Dual_geodesic_solver<FT>& solver, const TriangleMesh& tmesh)
{
  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using K =  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;
  VPM vpm = get(CGAL::vertex_point, tmesh);
  typedef typename GetInitializedFaceIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type FIM;
  typedef typename GetInitializedVertexIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type VIM;
  const FIM fim = get_initialized_face_index_map(tmesh, parameters::default_values());

  using Impl2 = typename internal::Geodesic_circle_impl<K, TriangleMesh, VPM, VIM, FIM>;
  solver=Impl2::make_dual_geodesic_solver(vpm, fim, tmesh);
}

/*!
 * \ingroup VGSMiscellaneous
 *  computes the approximate geodesic distances of `center` to all vertices of the mesh.
 \pre `tmesh` must consist of a single connected component
 * and put the distance in `distance_map`
 */
template <class FT, class TriangleMesh, class VertexDistanceMap>
void approximate_geodesic_distance_field(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& center,
                                         VertexDistanceMap distance_map,
                                         const TriangleMesh& tmesh)
{
  // TODO: the solver could be init once and used several times for different centers
  //       in particular, it can be tweaked to compute the Voronoi diagram of the initial centers
  //       or geodesic furthest point sampling.
  // TODO: add a parameter for the link size to increase to precision of the approximation of the distance
  //       that is you increase the size of the neighborhood of each vertex and you connect in the graph each vertex to its neighbors
  //       with weight being the geodesic shortest path.
  //       the more neighbors you have, the better is the approximation
  // TODO: graph construction can be done in parallel
  // TODO: concave flattening should be handled to improve the approximation of the distance
  //       (shortest path is not a line in that case)

  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  VPM vpm = get(CGAL::vertex_point, tmesh);
  using K =  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;

  typedef typename GetInitializedVertexIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type VIM;
  const VIM vim = get_initialized_vertex_index_map(tmesh, parameters::default_values());
  typedef typename GetInitializedFaceIndexMap<TriangleMesh, parameters::Default_named_parameters>::const_type FIM;

  using Impl = typename internal::Geodesic_circle_impl<K, TriangleMesh, VPM, VIM, FIM>;

  typename Impl::geodesic_solver solver = Impl::make_geodesic_solver(vpm, vim,tmesh);
  std::vector<double> distances = Impl::compute_geodesic_distances(solver, vpm, vim, tmesh, center);

  for (typename Impl::vertex_descriptor v : vertices(tmesh))
  {
    put(distance_map, v, distances[get(vim,v)]);
  }
}

} } // end of CGAL::Vector_graphics_on_surfaces

#endif // CGAL_POLYGON_MESH_PROCESSING_BSURF_DUAL_GEODESIC_SOLVER_H