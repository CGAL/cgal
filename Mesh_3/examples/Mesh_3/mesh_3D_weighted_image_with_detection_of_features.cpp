#define CGAL_MESH_3_VERBOSE 1
#define CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG

#include <vector>
#include <iostream>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

/// [Domain definition]
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/generate_label_weights.h>

#include <CGAL/Mesh_3/triple_lines_extraction/combinations.h>
#include <CGAL/Mesh_3/triple_lines_extraction/cases_table.h>
#include <CGAL/Mesh_3/triple_lines_extraction/triple_lines.h>
#include <CGAL/Mesh_3/triple_lines_extraction/cube_isometries.h>
#include <CGAL/Mesh_3/triple_lines_extraction/coordinates.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Image_domain = CGAL::Labeled_mesh_domain_3<K>;
using Mesh_domain = CGAL::Mesh_domain_with_polyline_features_3<Image_domain>;
/// [Domain definition]

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

/// [Add 1D features]
#include "read_polylines.h"
#include <CGAL/Mesh_3/polylines_to_protect.h> // undocumented header

typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef std::vector<Point_3>       Polyline_type;
typedef std::vector<Polyline_type> Polylines;


#include <boost/range/join.hpp>


// Protect the intersection of the object with the box of the image,
// by declaring 1D-features. Note that `CGAL::polylines_to_protect` is
// not documented.
bool add_1D_features(const CGAL::Image_3& image,
                     Mesh_domain& domain)
{
  typedef unsigned char Word_type;
  CGAL::Mesh_3::Triple_line_extractor<Point_3> lines;

  Polylines features_inside;

  const double vx = image.vx();
  const double vy = image.vy();
  const double vz = image.vz();
  const double dist_bound = (std::min)(vx, (std::min)(vy, vz)) / 256;
  const double sq_dist_bound = dist_bound * dist_bound;

  const std::size_t xdim = image.xdim();
  const std::size_t ydim = image.ydim();
  const std::size_t zdim = image.zdim();

  const double tx = image.tx();
  const double ty = image.ty();
  const double tz = image.tz();

  using CGAL::IMAGEIO::static_evaluate;

  typedef CGAL::Delaunay_triangulation_3<K> Del;
  Del triangulation;
  Del::Cell_handle start_cell;

  for (std::size_t k = 0, end_k = zdim - 1; k < end_k; ++k)
    for (std::size_t j = 0, end_j = ydim - 1; j < end_j; ++j)
      for (std::size_t i = 0, end_i = xdim - 1; i < end_i; ++i)
      {
        const K::Vector_3 translation{ i * vx + tx,
                                       j * vy + ty,
                                       k * vz + tz };
        const Cube cube = {
          static_evaluate<unsigned char>(image.image(), i  , j  , k),
          static_evaluate<unsigned char>(image.image(), i + 1, j  , k),
          static_evaluate<unsigned char>(image.image(), i  , j + 1, k),
          static_evaluate<unsigned char>(image.image(), i + 1, j + 1, k),
          static_evaluate<unsigned char>(image.image(), i  , j  , k + 1),
          static_evaluate<unsigned char>(image.image(), i + 1, j  , k + 1),
          static_evaluate<unsigned char>(image.image(), i  , j + 1, k + 1),
          static_evaluate<unsigned char>(image.image(), i + 1, j + 1, k + 1),
        }; /// TODO: optimize the access to the image data
        bool monocolor = cube[0] == cube[1];
        for (int i = 2; i < 8; ++i) monocolor = monocolor && (cube[0] == cube[i]);
        if (monocolor) continue;

        std::array<int, 8> inv_color_transformation{ INT_MIN, INT_MIN, INT_MIN, INT_MIN,
                                                     INT_MIN, INT_MIN, INT_MIN, INT_MIN };
        std::array<int, 256> color_transformation;
        std::fill(color_transformation.begin(), color_transformation.end(), INT_MIN);
        int nb_color = 0;
        for (int i = 0; i < 8; ++i) {
          if (color_transformation[cube[i]] == INT_MIN) {
            color_transformation[cube[i]] = nb_color;
            inv_color_transformation[nb_color] = cube[i];
            ++nb_color;
          }
        }
        if (nb_color > 3) {
          CGAL_warning_msg(nb_color > 3, "voxel with more than 3 colors");
          continue;
        }
        Cube reference_cube = {
          (unsigned char)(color_transformation[cube[0]]),
          (unsigned char)(color_transformation[cube[1]]),
          (unsigned char)(color_transformation[cube[2]]),
          (unsigned char)(color_transformation[cube[3]]),
          (unsigned char)(color_transformation[cube[4]]),
          (unsigned char)(color_transformation[cube[5]]),
          (unsigned char)(color_transformation[cube[6]]),
          (unsigned char)(color_transformation[cube[7]]),
        };
        auto case_it = find_case(cases, reference_cube);
        using std::end;
        const bool case_found = (case_it != end(cases));
        if (case_found) reference_cube = combinations[(*case_it)[8]];
        else {
          //std::cerr << "Warning: case not found: " << reference_cube << '\n';
          CGAL_error();
        };
#ifdef CGAL_DEBUG_TRIPLE_LINES
        std::cerr << "Cube " << cube << std::endl;
        std::cerr << "reference cube " << reference_cube << std::endl;
        std::cerr << "  with transformation " << cube_isometries[(*case_it)[9]] << "\n";
#endif // CGAL_DEBUG_TRIPLE_LINES
        auto fct_it = lines.create_polylines_fcts.find(reference_cube);
        if (fct_it != lines.create_polylines_fcts.end()) {
#ifdef CGAL_DEBUG_TRIPLE_LINES
          std::cerr << "Using the function of " << Cube(fct_it->first) << "\n";
#endif // CGAL_DEBUG_TRIPLE_LINES
          Polylines cube_features = (fct_it->second)(10);
          if (case_found) {
            const Permutation& transformation = cube_isometries[(*case_it)[9]];
            Coordinates a1 = coordinates[transformation[0]];
            Coordinates u = coordinates[transformation[1]] - a1;
            Coordinates v = coordinates[transformation[2]] - a1;
            Coordinates w = coordinates[transformation[4]] - a1;
            const Point_3  pa{ a1[0], a1[1], a1[2] };
            const Vector_3 vu{ u[0], u[1], u[2] };
            const Vector_3 vv{ v[0], v[1], v[2] };
            const Vector_3 vw{ w[0], w[1], w[2] };
#ifdef CGAL_DEBUG_TRIPLE_LINES
            std::cerr << "pa: " << pa << "\n";
            std::cerr << "vu: " << vu << "\n";
            std::cerr << "vv: " << vv << "\n";
            std::cerr << "vw: " << vw << "\n";
#endif // CGAL_DEBUG_TRIPLE_LINES
            for (auto& polyline : cube_features) {
              for (auto& point : polyline) {
                point = pa
                      + point.x() * vu
                      + point.y() * vv
                      + point.z() * vw;
                point = { vx * point.x(),
                          vy * point.y(),
                          vz * point.z(), };
                point = point + translation;
              }
              for (int i = 0; i < 2; ++i) {
                K::Point_3& extremity = (i == 0) ? polyline.front() : polyline.back();
                Del::Vertex_handle vh = triangulation.nearest_vertex(extremity, start_cell);
                if (Del::Vertex_handle() != vh) {
                  if (squared_distance(vh->point(), extremity) < sq_dist_bound) {
                    extremity = vh->point();
                  }
                }
                vh = triangulation.insert(extremity, start_cell);
                start_cell = vh->cell();
              }
              features_inside.push_back(std::move(polyline));
            } // end loop on polylines
          } // end case where the transformation is not the identity
        } // end if the reference_cube has polylines
      }

  // call the split_graph_into_polylines, to create long polylines from the
  // short polylines that were generated per voxel.
  Polylines new_polylines_inside;
  CGAL::polylines_to_protect<Point_3>(new_polylines_inside,
                                      features_inside.begin(),
                                      features_inside.end());

  std::vector<std::vector<Point_3> > polylines_on_bbox;
  CGAL::polylines_to_protect<Point_3, Word_type>(image, polylines_on_bbox,
                                                 new_polylines_inside.begin(),
                                                 new_polylines_inside.end());

  domain.add_features(polylines_on_bbox.begin(), polylines_on_bbox.end());

  // It is very important that the polylines from the file `lines_fname`
  // contain only polylines in the inside of the box of the image.
  domain.add_features(new_polylines_inside.begin(), new_polylines_inside.end());

  std::ofstream output_polylines("out-generated.polylines.txt");
  output_polylines.precision(17);
  for (auto poly : boost::range::join(polylines_on_bbox, new_polylines_inside)) {
    output_polylines << poly.size();
    for (auto p : poly) output_polylines << " " << p;
    output_polylines << std::endl;
  }

  return true;
}
/// [Add 1D features]


int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("images/liver.inr.gz");// 420.inr");
  // Loads image
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }

  /// [Domain creation]
  const float sigma = 10.f;
  CGAL::Image_3 img_weights =
    CGAL::Mesh_3::generate_label_weights(image, sigma);

  CGAL::Mesh_3::postprocess_weights_for_feature_protection(image, img_weights);

  Mesh_domain domain
    = Mesh_domain::create_labeled_image_mesh_domain(image,
                                                   weights = img_weights);
  /// [Domain creation]

  /// Declare 1D-features, see above [Call add_1D_features]
  if(!add_1D_features(image, domain)) {
    return EXIT_FAILURE;
  }
  /// [Call add_1D_features]

  CGAL::Bbox_3 bbox = domain.bbox();
  double diag = CGAL::sqrt((bbox.xmax() - bbox.xmin()) * (bbox.xmax() - bbox.xmin()) + (bbox.ymax() - bbox.ymin()) * (bbox.ymax() - bbox.ymin()) + (bbox.zmax() - bbox.zmin()) * (bbox.zmax() - bbox.zmin()));
  double sizing_default = diag * 0.05;

  /// Note that `edge_size` is needed with 1D-features [Mesh criteria]
  Mesh_criteria criteria(edge_size = sizing_default,
      facet_angle = 30, facet_size = sizing_default, facet_distance = sizing_default / 10,
      facet_topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH,
      cell_radius_edge_ratio = 0, cell_size = 0
  );
  /// [Mesh criteria]

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      CGAL::parameters::no_exude(),
                                      CGAL::parameters::no_perturb(),
                                      CGAL::parameters::mesh_3_dump());

  // Output
  CGAL::dump_c3t3(c3t3, "out");

  return 0;
}
