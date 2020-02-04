#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>

#include <CGAL/Classification.h>

#include <CGAL/Random.h>

#include <fstream>
#include <queue>

#include "include/Color_ramp.h"

///////////////////////////////////////////////////////////////////
//! [TIN DS]

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Projection_traits = CGAL::Projection_traits_xy_3<Kernel>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Segment_3 = Kernel::Segment_3;

// Triangulated Irregular Network
using TIN = CGAL::Delaunay_triangulation_2<Projection_traits>;

//! [TIN DS]
///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
//! [TIN_with_info DS]

// Triangulated Irregular Network (with info)
using Point_set = CGAL::Point_set_3<Point_3>;
using Vbi = CGAL::Triangulation_vertex_base_with_info_2 <Point_set::Index, Projection_traits>;
using Fbi = CGAL::Triangulation_face_base_with_info_2<int, Projection_traits>;
using TDS = CGAL::Triangulation_data_structure_2<Vbi, Fbi>;
using TIN_with_info = CGAL::Delaunay_triangulation_2<Projection_traits, TDS>;

//! [TIN_with_info DS]
///////////////////////////////////////////////////////////////////

namespace Classification = CGAL::Classification;

#ifdef CGAL_LINKED_WITH_TBB
using Concurrency_tag = CGAL::Parallel_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

///////////////////////////////////////////////////////////////////
//! [Contouring functions]

bool face_has_isovalue (TIN_with_info::Face_handle fh, double isovalue)
{
  bool above = false, below = false;
  for (int i = 0; i < 3; ++ i)
  {
    if (fh->vertex(i)->point().z() > isovalue)
      above = true;
    if (fh->vertex(i)->point().z() < isovalue)
      below = true;
  }

  return (above && below);
}

Segment_3 isocontour_in_face (TIN_with_info::Face_handle fh, double isovalue)
{
  Point_3 source;
  Point_3 target;
  bool source_found = false;

  for (int i = 0; i < 3; ++ i)
  {
    Point_3 p0 = fh->vertex((i+1) % 3)->point();
    Point_3 p1 = fh->vertex((i+2) % 3)->point();

    if ((p0.z() - isovalue) * (p1.z() - isovalue) > 0)
      continue;

    double zbottom = p0.z();
    double ztop = p1.z();
    if (zbottom > ztop)
    {
      std::swap (zbottom, ztop);
      std::swap (p0, p1);
    }
    
    double ratio = (isovalue - zbottom) / (ztop - zbottom);
    
    Point_3 p = CGAL::barycenter (p0, (1 - ratio), p1,ratio);
    
    if (source_found)
      target = p;
    else
    {
      source = p;
      source_found = true;
    }
  }

  return Segment_3 (source, target);
}


//! [Contouring functions]
///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
//! [Contouring visitor]

template <typename Graph>
class Polylines_visitor
{
private:
  Graph& graph;
  std::vector<std::vector<Point_3> >& polylines;

public:

  Polylines_visitor (Graph& graph, std::vector<std::vector<Point_3> >& polylines)
    : polylines (polylines), graph(graph) { }
  
  void start_new_polyline()
  {
    polylines.push_back (std::vector<Point_3>());
  }
  
  void add_node (typename Graph::vertex_descriptor vd)
  {
    polylines.back().push_back (graph[vd]);
  }
  
  void end_polyline()
  {
    // filter small polylines
    if (polylines.back().size() < 50)
      polylines.pop_back();
  }
};

//! [Contouring visitor]
///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
//! [CDT]

namespace PS = CGAL::Polyline_simplification_2;
using CDT_vertex_base = PS::Vertex_base_2<Projection_traits>;
using CDT_face_base = CGAL::Constrained_triangulation_face_base_2<Projection_traits>;
using CDT_TDS = CGAL::Triangulation_data_structure_2<CDT_vertex_base, CDT_face_base>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<Projection_traits, CDT_TDS>;
using CTP = CGAL::Constrained_triangulation_plus_2<CDT>;

//! [CDT]
///////////////////////////////////////////////////////////////////

int main (int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " points.ply" << std::endl;
    return EXIT_FAILURE;
  }

  ///////////////////////////////////////////////////////////////////
  //! [Init TIN]

  // Read points
  std::ifstream ifile (argv[1], std::ios_base::binary);
  CGAL::Point_set_3<Point_3> points;
  ifile >> points;
  std::cerr << points.size() << " point(s) read" << std::endl;

  // Create TIN
  TIN tin (points.points().begin(), points.points().end());

  //! [Init TIN]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Save TIN]

  CGAL::Surface_mesh<Point_3> tin_mesh;
  CGAL::copy_face_graph (tin, tin_mesh);
  std::ofstream tin_ofile ("tin.ply", std::ios_base::binary);
  CGAL::set_binary_mode (tin_ofile);
  CGAL::write_ply (tin_ofile, tin_mesh);

  //! [Save TIN]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [TIN_with_info]

  auto idx_to_point_with_info
    = [&](const Point_set::Index& idx) -> std::pair<Point_3, Point_set::Index>
      {
        return std::make_pair (points.point(idx), idx);
      };

  TIN_with_info tin_with_info
    (boost::make_transform_iterator (points.begin(), idx_to_point_with_info),
     boost::make_transform_iterator (points.end(), idx_to_point_with_info));
     
  //! [TIN_with_info]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Components]

  double spacing = CGAL::compute_average_spacing<Concurrency_tag>(points, 6);
  spacing *= 2;

  auto face_height
    = [&](const TIN_with_info::Face_handle fh) -> double
      {
        double out = 0.;
        for (int i = 0; i < 3; ++ i)
          out = (std::max) (out, CGAL::abs(fh->vertex(i)->point().z() - fh->vertex((i+1)%3)->point().z()));
        return out;
      };

  // Initialize faces info
  for (TIN_with_info::Face_handle fh : tin_with_info.all_face_handles())
    if (tin_with_info.is_infinite(fh) || face_height(fh) > spacing) // Filtered faces are given info() = -2
      fh->info() = -2;
    else // Pending faces are given info() = -1;
      fh->info() = -1;

  // Flooding algorithm
  std::vector<int> component_size;
  for (TIN_with_info::Face_handle fh : tin_with_info.finite_face_handles())
  {
    if (fh->info() != -1)
      continue;

    std::queue<TIN_with_info::Face_handle> todo;
    todo.push(fh);

    int size = 0;
    while (!todo.empty())
    {
      TIN_with_info::Face_handle current = todo.front();
      todo.pop();

      if (current->info() != -1)
        continue;
      current->info() = component_size.size();
      ++ size;
      
      for (int i = 0; i < 3; ++ i)
        todo.push (current->neighbor(i));
    }

    component_size.push_back (size);
  }

  std::cerr << component_size.size() << " connected component(s) found" << std::endl;

  //! [Components]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Save TIN with info]

  using Mesh = CGAL::Surface_mesh<Point_3>;
  
  Mesh tin_colored_mesh;
  
  Mesh::Property_map<Mesh::Face_index, CGAL::Color>
    color_map = tin_colored_mesh.add_property_map<Mesh::Face_index, CGAL::Color>("f:color").first;
  
  CGAL::copy_face_graph (tin_with_info, tin_colored_mesh,
                         CGAL::parameters::face_to_face_output_iterator
                         (boost::make_function_output_iterator
                          ([&](const std::pair<TIN_with_info::Face_handle, Mesh::Face_index>& ff)
                           {
                             // Color unassigned faces grey
                             if (ff.first->info() < 0)
                               color_map[ff.second] = CGAL::Color(128, 128, 128);
                             else
                             {
                               // Random color seeded by the component ID
                               CGAL::Random r (ff.first->info());
                               color_map[ff.second] = CGAL::Color (r.get_int(64, 192),
                                                                   r.get_int(64, 192),
                                                                   r.get_int(64, 192));
                             }
                           })));
                               
  std::ofstream tin_colored_ofile ("colored_tin.ply", std::ios_base::binary);
  CGAL::set_binary_mode (tin_colored_ofile);
  CGAL::write_ply (tin_colored_ofile, tin_colored_mesh);

  //! [Save TIN with info]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Filtering]

  int min_size = 100000;

  std::vector<TIN_with_info::Vertex_handle> to_remove;
  std::ofstream dbg ("dbg.xyz");
  dbg.precision(18);

  for (TIN_with_info::Vertex_handle vh : tin_with_info.finite_vertex_handles())
  {
    TIN_with_info::Face_circulator circ = tin_with_info.incident_faces (vh),
      start = circ;

    // Remove a vertex if it's only adjacent to components smaller than threshold
    bool keep = false;
    do
    {
      if (circ->info() >= 0 && component_size[std::size_t(circ->info())] > min_size)
      {
        keep = true;
        break;
      }
    }
    while (++ circ != start);

    if (!keep)
    {
      to_remove.push_back (vh);
      dbg << vh->point() << std::endl;
    }
  }

  std::cerr << to_remove.size() << " vertices(s) will be removed after filtering" << std::endl;
  for (TIN_with_info::Vertex_handle vh : to_remove)
    tin_with_info.remove (vh);
  
  //! [Filtering]
  ///////////////////////////////////////////////////////////////////

  // Save as Mesh
  CGAL::Surface_mesh<Point_3> dsm_mesh;
  CGAL::copy_face_graph (tin_with_info, dsm_mesh);
  std::ofstream dsm_ofile ("dsm.ply", std::ios_base::binary);
  CGAL::set_binary_mode (dsm_ofile);
  CGAL::write_ply (dsm_ofile, dsm_mesh);

  ///////////////////////////////////////////////////////////////////
  //! [Rastering]

  CGAL::Bbox_3 bbox = CGAL::bbox_3 (points.points().begin(), points.points().end());
  
  // Generate raster image 1920-pixels large
  std::size_t width = 1920;
  std::size_t height = std::size_t((bbox.ymax() - bbox.ymin()) * 1920 / (bbox.xmax() - bbox.xmin()));

  std::cerr << "Rastering with resolution " << width << "x" << height << std::endl;

  // Use PPM format (Portable PixMap) for simplicity
  std::ofstream raster_ofile ("raster.ppm", std::ios_base::binary);
  raster_ofile << "P6" << std::endl << width << " " << height << std::endl << 255 << std::endl;

  TIN_with_info::Face_handle location;

  // Use rainbow color ramp output
  Color_ramp color_ramp;
  
  for (std::size_t y = 0; y < height; ++ y)
    for (std::size_t x = 0; x < width; ++ x)
    {
      Point_3 query (bbox.xmin() + x * (bbox.xmax() - bbox.xmin()) / double(width),
                     bbox.ymin() + y * (bbox.ymax() - bbox.ymin()) / double(height),
                     0); // not relevant for location in 2D
      
      location = tin_with_info.locate (query, location);

      std::array<unsigned char, 3> colors { 0, 0, 0 };
      
      if (!tin_with_info.is_infinite(location))
      {
        std::array<double, 3> barycentric_coordinates
          = CGAL::Polygon_mesh_processing::barycentric_coordinates
          (Point_2 (location->vertex(0)->point().x(), location->vertex(0)->point().y()),
           Point_2 (location->vertex(1)->point().x(), location->vertex(1)->point().y()),
           Point_2 (location->vertex(2)->point().x(), location->vertex(2)->point().y()),
           Point_2 (query.x(), query.y()),
           Kernel());

        double height_at_query
          = (barycentric_coordinates[0] * location->vertex(0)->point().z()
             + barycentric_coordinates[1] * location->vertex(1)->point().z()
             + barycentric_coordinates[2] * location->vertex(2)->point().z());

        double height_ratio = (height_at_query - bbox.zmin()) / (bbox.zmax() - bbox.zmin());
        colors = color_ramp.get(height_ratio);
      }
      raster_ofile.write ((char*)(&colors), 3);

    }
  
  //! [Rastering]
  ///////////////////////////////////////////////////////////////////

  // Smooth heights with 5 successive Gaussian filters
  double gaussian_variance = 4 * spacing * spacing;
  for (std::size_t i = 0; i < 5; ++ i)
    for (TIN_with_info::Vertex_handle vh : tin_with_info.finite_vertex_handles())
    {
      double z = vh->point().z();
      double total_weight = 1;

      TIN_with_info::Vertex_circulator circ = tin_with_info.incident_vertices (vh),
        start = circ;

      do
      {
        if (!tin_with_info.is_infinite(circ))
        {
          double sq_dist = CGAL::squared_distance (vh->point(), circ->point());

          double weight = std::exp(- sq_dist / gaussian_variance);
          z += weight * circ->point().z();
          total_weight += weight;
        }
      }
      while (++ circ != start);

      z /= total_weight;

      vh->point() = Point_3 (vh->point().x(), vh->point().y(), z);
    }
  
  ///////////////////////////////////////////////////////////////////
  //! [Contouring extraction]

  std::array<double, 50> isovalues; // Contour 50 isovalues
  for (std::size_t i = 0; i < isovalues.size(); ++ i)
    isovalues[i] = bbox.zmin() + ((i+1) * (bbox.zmax() - bbox.zmin()) / (isovalues.size() - 2));

  // First find on each face if they are crossed by some isovalues and
  // extract segments in a graph
  using Segment_graph = boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, Point_3>;
  Segment_graph graph;
  using Map_p2v = std::map<Point_3, Segment_graph::vertex_descriptor>;
  Map_p2v map_p2v;
  for (TIN_with_info::Face_handle vh : tin_with_info.finite_face_handles())
    for (double iv : isovalues)
      if (face_has_isovalue (vh, iv))
      {
        Segment_3 segment = isocontour_in_face (vh, iv);
        for (const Point_3& p : { segment.source(), segment.target() })
        {
          Map_p2v::iterator iter;
          bool inserted;
          std::tie (iter, inserted) = map_p2v.insert (std::make_pair (p, Segment_graph::vertex_descriptor()));
          if (inserted)
          {
            iter->second = boost::add_vertex (graph);
            graph[iter->second] = p;
          }
        }
        boost::add_edge (map_p2v[segment.source()], map_p2v[segment.target()], graph);
      }

  //! [Contouring extraction]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Contouring split]

  // Split segments into polylines
  std::vector<std::vector<Point_3> > polylines;
  Polylines_visitor<Segment_graph> visitor (graph, polylines);
  CGAL::split_graph_into_polylines (graph, visitor);

  std::cerr << polylines.size() << " polylines computed, with "
            << map_p2v.size() << " vertices in total" << std::endl;

  // Output to polyline file
  std::ofstream contour_ofile ("contour.polylines.txt");
  contour_ofile.precision(18);
  for (const std::vector<Point_3>& poly : polylines)
  {
    contour_ofile << poly.size();
    for (const Point_3& p : poly)
      contour_ofile << " " << p;
    contour_ofile << std::endl;
  }

  //! [Contouring split]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Contouring simplify]

  // Construct constrained Delaunay triangulation with polylines as constraints
  CTP ctp;
  for (const std::vector<Point_3>& poly : polylines)
    ctp.insert_constraint (poly.begin(), poly.end());

  // Simplification algorithm with limit on distance
  PS::simplify (ctp, PS::Squared_distance_cost(), PS::Stop_above_cost_threshold (16 * spacing * spacing));

  // Output to file
  std::ofstream simplified_ofile ("simplified.polylines.txt");
  simplified_ofile.precision(18);
  std::size_t nb_vertices = 0;
  for (CTP::Constraint_id cid : ctp.constraints())
  {
    simplified_ofile << ctp.vertices_in_constraint (cid).size();
    nb_vertices += ctp.vertices_in_constraint (cid).size();
    for (CTP::Vertex_handle vh : ctp.vertices_in_constraint(cid))
      simplified_ofile << " " << vh->point();
    simplified_ofile << std::endl;
  }

  std::cerr << nb_vertices << " vertices remaining after simplification ("
            << 100. * (nb_vertices / double(map_p2v.size())) << "%)" << std::endl;
  

  //! [Contouring simplify]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Classification]

  // Get training from input
  Point_set::Property_map<int> training_map;
  bool training_found;
  std::tie (training_map, training_found) = points.property_map<int>("training");
  
  if (training_found)
  {
    std::cerr << "Classifying ground/vegetation/building" << std::endl;
    
    // Create labels
    Classification::Label_set labels ({ "ground", "vegetation", "building" });

    // Generate features
    Classification::Feature_set features;
    Classification::Point_set_feature_generator<Kernel, Point_set, Point_set::Point_map>
      generator (points, points.point_map(), 5); // 5 scales
    
#ifdef CGAL_LINKED_WITH_TBB
    // If TBB is used, features can be computed in parallel
    features.begin_parallel_additions();
    generator.generate_point_based_features (features);
    features.end_parallel_additions();
#else
    generator.generate_point_based_features (features);
#endif

    // Train a random forest classifier
    Classification::ETHZ::Random_forest_classifier classifier (labels, features);
    classifier.train (points.range(training_map));

    // Classify with graphcut regularization
    Point_set::Property_map<int> label_map = points.add_property_map<int>("labels").first;
    Classification::classify_with_graphcut<Concurrency_tag>
      (points, points.point_map(), labels, classifier,
       generator.neighborhood().k_neighbor_query(12), // regularize on 12-neighbors graph
       0.5f, // graphcut weight
       12, // Subdivide to speed-up process
       label_map);

    // Evaluate
    std::cerr << "Mean IoU on training data = "
              << Classification::Evaluation(labels,
                                            points.range(training_map),
                                            points.range(label_map)).mean_intersection_over_union() << std::endl;

    // Save the classified point set
    std::ofstream classified_ofile ("classified.ply");
    CGAL::set_binary_mode (classified_ofile);
    classified_ofile << points;
  }
  
  //! [Classification]
  ///////////////////////////////////////////////////////////////////
  

  return EXIT_SUCCESS;
}
  
