#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Vector_graphics_on_surfaces/locally_shortest_path.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Union_find.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>

#include <boost/graph/adjacency_list.hpp>

#include <nanosvg.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
typedef PMP::Face_location<Mesh, double>                      Face_location;
typedef PMP::Edge_location<Mesh, double>                      Edge_location;


int main(int argc, char** argv)
{
  std::string mesh_filename = (argc > 1) ? std::string(argv[1])
    : CGAL::data_file_path("meshes/elephant.off");

  std::string svg_filename = (argc > 2) ? std::string(argv[2])
    : CGAL::data_file_path("polylines_2/nano.svg");

  const double expected_diag = (argc>3) ? atof(argv[3]) : 0.6; // user parameter for scaling
  int face_index = (argc>4) ? atoi(argv[4]) : 2154;
  double b0 = (argc>5) ? atof(argv[5]) : 0.3;
  double b1 = (argc>6) ? atof(argv[6]) : 0.3;
  double b2 = (argc>7) ? atof(argv[7]) : 0.4;

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(mesh_filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input mesh." << std::endl;
    return 1;
  }

  NSVGimage* g_image = nsvgParseFromFile(svg_filename.c_str(), "px", 96.0f);
  if (g_image == NULL) {
    printf("Could not open SVG image.\n");
    return 1;
  }

  // extract control points
  std::vector< std::array<K::Point_2, 4> > bezier_control_points;
  CGAL::Bbox_2 bb2;

  // in SVG's the y axis points downward, so we must take the opposite y coordinates
  for (NSVGshape* shape = g_image->shapes; shape != NULL; shape = shape->next)
  {
    for (NSVGpath* path = shape->paths; path != NULL; path = path->next)
    {
      CGAL::Bbox_2 path_bbox(path->bounds[0], -path->bounds[1],
                             path->bounds[2], -path->bounds[3]);
      bb2+=path_bbox;

      float* pts=path->pts;
      int npts=path->npts;

      for (int i=0; i<npts-1; i += 3)
      {
        bezier_control_points.emplace_back();
        float* p = &pts[i*2];
        bezier_control_points.back()[0]=K::Point_2(p[0],-p[1]);
        bezier_control_points.back()[1]=K::Point_2(p[2],-p[3]);
        bezier_control_points.back()[2]=K::Point_2(p[4],-p[5]);
        bezier_control_points.back()[3]=K::Point_2(p[6],-p[7]);
      }
    }
  }

  nsvgDelete(g_image);

  std::cout << "#Bezier curves read: " << bezier_control_points.size() << "\n";

  // convert control points to polar coordinates
  typename K::Point_2 center_2((bb2.xmax()+bb2.xmin())/2., (bb2.ymax()+bb2.ymin())/2.);
  double diag = std::sqrt( CGAL::square(bb2.xmin()-bb2.xmax()) + CGAL::square(bb2.xmin()-bb2.xmax()) );
  const double scaling = expected_diag/diag;

  //TODO: do the scaling at read time!

  PMP::Dual_geodesic_solver<double> solver;
  PMP::init_geodesic_dual_solver(solver, mesh);

  std::size_t nb_faces = faces(mesh).size();
  Mesh::Face_index f = *std::next(faces(mesh).begin(), (face_index)%nb_faces);
  Face_location center(f, CGAL::make_array(b0,b1,b2));

  std::size_t num_subdiv = 3;

  // option 1: trace bezier curves with a single center
  {
    CGAL::Real_timer time;
    time.start();
    std::vector<std::array<K::Vector_2, 4>> directions;
    std::vector<std::array<K::FT, 4>> lengths;
    directions.reserve(bezier_control_points.size());
    lengths.reserve(bezier_control_points.size());

    for (const std::array<K::Point_2, 4>& bezier  : bezier_control_points)
    {
      std::vector<std::pair<double, double>> polar_coords =
        PMP::convert_to_polar_coordinates<K>(bezier, center_2);

      directions.emplace_back();
      lengths.emplace_back();

      assert(polar_coords.size()==4);

      for (int i=0;i<4; ++i)
      {
        lengths.back()[i] = scaling * polar_coords[i].first;
        directions.back()[i]=K::Vector_2(std::cos(polar_coords[i].second), std::sin(polar_coords[i].second));
      }
    }

    std::vector< std::vector<Face_location> > res =
      PMP::trace_bezier_curves<K>(center, directions, lengths, num_subdiv, mesh, solver);

    // write result
    std::ofstream out("svg_option1.polylines.txt");
    out << std::setprecision(17);
    for (const auto& b : res)
    {
      std::vector<K::Point_3> poly;
      poly.reserve(b.size());
      PMP::convert_path_to_polyline(b, mesh, std::back_inserter(poly));


      out << poly.size();
      for (const K::Point_3& p : poly)
        out << " " << p;
      out << "\n";
    }
    time.stop();
    std::cout << "option 1 took: " << time.time() << "\n";
  }

  // option 2: regroup control polygon and trace bezier curves with a center per group
  {
    CGAL::Real_timer time;
    time.start();

    typedef boost::adjacency_list< boost::vecS, boost::vecS,
                                   boost::undirectedS,
                                   K::Point_2> Graph;
    using Graph_vertex = Graph::vertex_descriptor;

    Graph graph;

    std::vector<Graph_vertex> vrts;
    std::map<K::Point_2, Graph_vertex> pt_map;
    Graph_vertex null_vertex = boost::graph_traits<Graph>::null_vertex();

    // build bezier control segment graph
    for (std::size_t eid=0; eid<bezier_control_points.size(); ++eid)
    {
      const std::array<K::Point_2, 4>& bezier = bezier_control_points[eid];

      auto insert_res_0 = pt_map.emplace(bezier[0], null_vertex);
      if (insert_res_0.second)
      {
        vrts.push_back(add_vertex(graph));
        graph[vrts.back()]=bezier[0];
        insert_res_0.first->second = vrts.back();
      }

      auto insert_res_3 = pt_map.emplace(bezier[3], null_vertex);
      if (insert_res_3.second)
      {
        vrts.push_back(add_vertex(graph));
        graph[vrts.back()]=bezier[3];
        insert_res_3.first->second = vrts.back();
      }

      Graph_vertex v1 = add_vertex(graph);
      Graph_vertex v2 = add_vertex(graph);
      graph[v1]=bezier[1];
      graph[v2]=bezier[2];

      add_edge(insert_res_0.first->second, v1, graph);
      add_edge(v1, v2, graph);
      add_edge(v2, insert_res_3.first->second,  graph);
    }

    // visitor for split_graph_into_polylines that will also take care of
    // the drawing and writing in the output file
    struct Collect_visitor
    {
      const Graph& graph;
      std::vector<std::vector<K::Point_2>>& polylines_2;
      std::vector<Graph_vertex> current_polyline;

      Collect_visitor(const Graph& g, std::vector<std::vector<K::Point_2>>& poly2)
        : graph(g)
        , polylines_2(poly2)
      {}

      void start_new_polyline()
      {
        current_polyline.clear();
      }

      void add_node(Graph_vertex v)
      {
        current_polyline.push_back(v);
      }

      void end_polyline()
      {
        // recover polyline of control points
        polylines_2.emplace_back();
        polylines_2.back().reserve(current_polyline.size());
        for (Graph_vertex v  : current_polyline)
          polylines_2.back().push_back(graph[v]);
      }
    };


    std::vector<std::vector<K::Point_2>> polylines_2;
    Collect_visitor svisitor(graph, polylines_2);
    CGAL::split_graph_into_polylines(graph, svisitor);
    std::size_t nb_poly = polylines_2.size();

    // regroup polylines if they have overlapping bboxes
    std::vector<CGAL::Bbox_2> bboxes_2(nb_poly);
    using UF = CGAL::Union_find<std::size_t>;
    std::vector<UF::handle> uf_handles(nb_poly);
    UF uf;
    for(std::size_t i=0;i<nb_poly; ++i)
    {
      bboxes_2[i]=CGAL::bbox_2(polylines_2[i].begin(), polylines_2[i].end());
      uf_handles[i]=uf.make_set(i);
    }

    // TODO: use box intersection d?
    for (std::size_t i=0; i<nb_poly-1; ++i)
    {
      for (std::size_t j=i+1; j<nb_poly; ++j)
      {
        if (CGAL::do_overlap(bboxes_2[i], bboxes_2[j]))
          uf.unify_sets(uf_handles[i], uf_handles[j]);
      }
    }

    std::vector<std::size_t> partition_id(nb_poly,-1);
    std::vector< std::vector<std::size_t> > partition;
    partition.reserve(uf.number_of_sets());

    for (std::size_t i=0; i<nb_poly; ++i)
    {
      std::size_t master_i = *uf.find(uf_handles[i]);
      if (partition_id[master_i]==std::size_t(-1))
      {
        partition_id[master_i]=partition.size();
        partition.emplace_back();
      }
      partition[partition_id[master_i]].push_back(i);
    }

    // now draw each group of polyline
    std::vector<std::vector<Face_location>> polygons_3;
    std::ofstream out("svg_option2.polylines.txt");
    out << std::setprecision(17);
    bool all_closed=true;
    for (const std::vector<std::size_t>& pids : partition)
    {
      // center of the polyline
      CGAL::Bbox_2 poly_bb;
      for (std::size_t i : pids)
        poly_bb+=bboxes_2[i];
      K::Point_2 poly_center((poly_bb.xmin()+poly_bb.xmax())/2., (poly_bb.ymin()+poly_bb.ymax())/2.);

      // path from the general center to the polyline center
      K::Vector_2 dir(center_2, poly_center);
      std::vector<Face_location> path_2_center =
        PMP::straightest_geodesic<K>(center, dir, scaling * std::sqrt(dir.squared_length()), mesh);
      Face_location poly_center_loc = path_2_center.back();

      // update initial angle from the global center to the polyline center
      using Impl=PMP::internal::Locally_shortest_path_imp<K, Mesh, decltype(mesh.points())>;
      K::Vector_2 initial_dir(1,0);// TODO: this is arbitray and could be a user input or from PCA...


      // TODO: avoid using the shortest path and directly use the straightest!
      std::vector<Edge_location> shortest_path;
        PMP::locally_shortest_path<K::FT>(center, poly_center_loc, mesh, shortest_path, solver);

      typename K::Vector_2 v = initial_dir;
      // TODO: the code uses the straightest path but since the code expects Edge_location,
      // it is not working directly. We need to extract the edge encoded by the face location
      // for(std::size_t i=1;i<path_2_center.size();++i)
      // {
      //   Mesh::Halfedge_index h_curr = halfedge(path_2_center[i].first, mesh);
      //   v = Impl::parallel_transport_f2f(h_curr, v, mesh.points(), mesh);
      // }
      for(std::size_t i=0;i<shortest_path.size();++i)
      {
        Mesh::Halfedge_index h_curr = halfedge(shortest_path[i].first, mesh);
        v = Impl::parallel_transport_f2f(h_curr, v, mesh.points(), mesh);
      }
      double theta = atan2(v.y(),v.x());

      for (std::size_t pid : pids)
      {
        std::vector<K::Vector_2> directions;
        std::vector<K::FT> lengths;
        directions.reserve(polylines_2[pid].size());
        lengths.reserve(polylines_2[pid].size());

        // polar coordinates from the polyline center and convertion to lengths and directions
        std::vector<std::pair<double, double>> polar_coords =
          PMP::convert_to_polar_coordinates<K>(polylines_2[pid], poly_center);

        for (const std::pair<double, double>& polar_coord : polar_coords)
        {
          lengths.push_back(scaling * polar_coord.first);
          directions.push_back(K::Vector_2(std::cos(polar_coord.second+theta), std::sin(polar_coord.second+theta)));
        }

        bool is_closed=polylines_2[pid].front()==polylines_2[pid].back();
        if (is_closed)
        {
          lengths.pop_back();
          directions.pop_back();
        }
        else
          all_closed=false;


        std::vector<Face_location> res =
          PMP::trace_polyline_of_bezier_curves<K>(poly_center_loc, directions, lengths,
                                                  is_closed, // use [directions/lengths].front as last control point?
                                                  num_subdiv, mesh, solver);
        // write result
        std::vector<K::Point_3> poly3;
        poly3.reserve(res.size());
        PMP::convert_path_to_polyline(res, mesh, std::back_inserter(poly3));

        out << poly3.size();
        for (const K::Point_3& p : poly3)
          out << " " << p;
        out << "\n";

        polygons_3.push_back(std::move(res));
      }
    }
    out.close();

    time.stop();
    std::cout << "option 2 took: " << time.time() << "\n";


    if (all_closed)
    {
      std::cout << "Now carve the input mesh\n";
      // copy/pasted from trace_polygon_example.cpp
      // now refine the input mesh
      std::vector<Mesh::Halfedge_index> cst_hedges;
      auto vnm = mesh.add_property_map<Mesh::Vertex_index, K::Vector_3>("vnm", K::Vector_3(0,0,0)).first;
      auto fnm = mesh.add_property_map<Mesh::Face_index, K::Vector_3>("fnm", K::Vector_3(0,0,0)).first;
      using VNM = decltype(vnm);

      PMP::compute_normals(mesh, vnm, fnm);
      PMP::refine_mesh_along_paths<K>(polygons_3, mesh, vnm, fnm, std::back_inserter(cst_hedges));

      std::ofstream("mesh_refined.off") << std::setprecision(17) << mesh;

      std::ofstream cst_edges("refinement_edges.polylines.txt");
      cst_edges.precision(17);
      for (Mesh::Halfedge_index h : cst_hedges)
        cst_edges << "2 " << mesh.point(source(h,mesh)) << " " << mesh.point(target(h,mesh)) << "\n";


      auto ecm = mesh.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
      for (Mesh::Halfedge_index h : cst_hedges)
        ecm[edge(h, mesh)]=true;


      // face index for doing flood fill and mark inside-out
      Mesh::Face_index out_face(2612);
      std::vector<int> in_out(num_faces(mesh), -1);

      bool inorout=false;
      std::vector<Mesh::Face_index> queue, next_queue;
      queue.push_back(out_face);

      while(!queue.empty())
      {
        Mesh::Face_index f = queue.back();
        queue.pop_back();
        if (in_out[f]==-1)
        {
          in_out[f]=inorout?1:0;
          Mesh::Halfedge_index h=halfedge(f, mesh);
          for (int i=0; i<3; ++i)
          {
            Mesh::Face_index nf = face(opposite(h, mesh), mesh);
            if (nf!=boost::graph_traits<Mesh>::null_face() && in_out[nf]==-1)
            {
              if (ecm[edge(h,mesh)])
                next_queue.push_back(nf);
              else
                queue.push_back(nf);
            }
            h=next(h, mesh);
          }
        }
        if (queue.empty())
        {
          queue.swap(next_queue);
          inorout=!inorout;
        }
      }

      struct Visitor
        : public PMP::Corefinement::Default_visitor<Mesh>
      {
        VNM vnm;
        Visitor(VNM vnm) : vnm(vnm) {}

        std::vector<std::pair<Mesh::Halfedge_index, Mesh::Halfedge_index> > hedge_map;
        void after_edge_duplicated(Mesh::Halfedge_index h, Mesh::Halfedge_index new_hedge, const Mesh&)
        {
          hedge_map.emplace_back(h, new_hedge);
        }

        void after_vertex_copy(Mesh::Vertex_index v, const Mesh&, Mesh::Vertex_index nv, const Mesh&)
        {
          put(vnm, nv, get(vnm, v));
        }
      };
      Visitor visitor(vnm);

      PMP::internal::split_along_edges(mesh, ecm, mesh.points(), visitor);

      double delta = -0.0005;
      for (const auto& ph : visitor.hedge_map)
      {
        Mesh::Halfedge_index h1 = ph.first;
        Mesh::Halfedge_index h2 = ph.second;
        if (is_border(h1, mesh)) h1=opposite(h1, mesh);
        if (is_border(h2, mesh)) h2=opposite(h2, mesh);
        Mesh::Halfedge_index h = in_out[face(h1, mesh)]==1 ? h1 : h2;

        Mesh::Vertex_index v = target(h, mesh);
        K::Vector_3 n = get(vnm, v);
        mesh.point(v) = mesh.point(v)+delta*n;
      }

      // interior vertices
      for (Mesh::Vertex_index v : vertices(mesh))
      {
        bool skip=false;
        Mesh::Halfedge_index h = halfedge(v, mesh);
        for (Mesh::Halfedge_index h : CGAL::halfedges_around_target(v, mesh))
        {
          if (is_border(h, mesh) || in_out[face(h, mesh)]==0)
          {
            skip=true;
            break;
          }
        }
        if (!skip)
        {
          K::Vector_3 n = get(vnm, v);
          mesh.point(v) = mesh.point(v)+delta*n;
        }
      }

      std::vector<Mesh::Halfedge_index> b1(visitor.hedge_map.size());
      std::vector<Mesh::Halfedge_index> b2(visitor.hedge_map.size());
      for (std::size_t i=0; i<visitor.hedge_map.size(); ++i)
      {
        Mesh::Halfedge_index h1 = visitor.hedge_map[i].first;
        Mesh::Halfedge_index h2 = visitor.hedge_map[i].second;
        if (is_border(h1, mesh)) h1=opposite(h1, mesh);
        if (is_border(h2, mesh)) h2=opposite(h2, mesh);
        if (in_out[face(h1, mesh)]==1) std::swap(h1,h2);

        b1[i]=opposite(h1, mesh);
        b2[i]=opposite(h2, mesh);

      }

      PMP::extrude_impl::create_strip(b1, b2, mesh);

      std::ofstream("mesh_refined_split.off") << std::setprecision(17) << mesh;
    }
  }

  return 0;
}
