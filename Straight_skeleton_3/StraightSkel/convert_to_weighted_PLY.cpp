#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <iostream>
#include <fstream>

namespace PMP = ::CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;
using Vector = K::Vector_3;

using Mesh = CGAL::Surface_mesh<Point>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

bool assign_weights(Mesh& sm,
                    const char* weights_filename)
{
  // do not change the pmap's name without changing it in 'PLYFile.cpp'
  auto res = sm.add_property_map<face_descriptor, double>("f:weight");
  if(!res.second) {
    std::cerr << "Error: failed to add property map?" << std::endl;
    return false;
  }

  if(!weights_filename) {
    std::cout << "no input weights, all weights are set to '1'." << std::endl;
    for(face_descriptor f : faces(sm))
      put(res.first, f, 1.);

    return true;
  }

  auto fwm = res.first;

  std::ifstream weights_in(weights_filename) ;
  std::string x1, x2, y1, y2, z1, z2;
  double vx1, vx2, vy1, vy2, vz1, vz2;
  vx1 = vx2 = vy1 = vy2 = vz1 = vz2 = 0;

  if(!(weights_in >> x1 >> vx1
                  >> x2 >> vx2
                  >> y1 >> vy1
                  >> y2 >> vy2)) {
    std::cerr << "Error: failed to read weights" << std::endl;
    return false;
  }

  if(weights_in >> z1 >> vz1
                >> z2 >> vz2) {
    std::cout << "bottom & top info is present" << std::endl;
  } else {
    if(vx2 != vy2) {
      std::cerr << "Error: unknown z speeds, and x-speed and y-speed differ..." << std::endl;
      return false;
    }

    // assign the uniform speed to the "up" direction
    vz2 = vx2;
  }

  if(vx1 < 0 || vx2 < 0 || vy1 < 0 || vy2 < 0 || vz1 < 0 || vz2 < 0) {
    std::cerr << "Error: negative weights?" << std::endl;
    return false;
  }

  double eps_weight = std::numeric_limits<double>::max();
  if(vx1 > 0.) eps_weight = (std::min)(eps_weight, vx1);
  if(vx2 > 0.) eps_weight = (std::min)(eps_weight, vx2);
  if(vy1 > 0.) eps_weight = (std::min)(eps_weight, vy1);
  if(vy2 > 0.) eps_weight = (std::min)(eps_weight, vy2);
  if(vz1 > 0.) eps_weight = (std::min)(eps_weight, vz1);
  if(vz2 > 0.) eps_weight = (std::min)(eps_weight, vz2);

  if(eps_weight == 0.) {
    std::cerr << "Error: all weights to zero" << std::endl;
    return false;
  }

  eps_weight = 1e-2 * eps_weight;

  if(vx1 == 0.) { std::cout << "vx1 to eps weight" << std::endl; vx1 = eps_weight; }
  if(vx2 == 0.) { std::cout << "vx2 to eps weight" << std::endl; vx2 = eps_weight; }
  if(vy1 == 0.) { std::cout << "vy1 to eps weight" << std::endl; vy1 = eps_weight; }
  if(vy2 == 0.) { std::cout << "vy2 to eps weight" << std::endl; vy2 = eps_weight; }
  if(vz1 == 0.) { std::cout << "vz1 to eps weight" << std::endl; vz1 = eps_weight; }
  if(vz2 == 0.) { std::cout << "vz2 to eps weight" << std::endl; vz2 = eps_weight; }

  const Vector east  {  1,  0,  0 }; // x2
  const Vector south {  0, -1,  0 }; // y1
  const Vector west  { -1,  0,  0 }; // x1
  const Vector north {  0,  1,  0 }; // y2
  const Vector up    {  0,  0,  1 }; // z2
  const Vector down  {  0,  0, -1 }; // z1

  for(face_descriptor f : faces(sm))
  {
    double weight = 1.;

    // internal stuff, we don't need to normalize and introduce inexactness
    // @todo do that in exact?
    Vector v = CGAL::NULL_VECTOR;
    PMP::internal::sum_normals<Point>(sm, f, get(CGAL::vertex_point, sm), v, K());
    FT sq_n = v.squared_length();

    if(v.x() == 0 && v.y() == 0) {
      if(v.z() > 0)
        weight = vz2;
      else
        weight = vz1;
    } else {
      if(v.x() >= 0) {
        const FT sq_cos = CGAL::square(CGAL::scalar_product(v, west)) / sq_n;
        if(v.y() >= 0) {
          // north east quadrant
          weight = vy2 * (1 - sq_cos) + vx2 * sq_cos;
        } else {
          // south east quadrant
          weight = vy1 * (1 - sq_cos) + vx2 * sq_cos;
        }
      } else { // x < 0
        const FT sq_cos = CGAL::square(CGAL::scalar_product(v, east)) / sq_n;
        if(v.y() >= 0) {
          // north west quadrant
          weight = vy2 * (1 - sq_cos) + vx1 * sq_cos;
        } else {
          // south west quadrant
          weight = vy1 * (1 - sq_cos) + vx1 * sq_cos;
        }
      }
    }

    std::cout << "face: " << f << " weight: " << weight << std::endl;
    put(fwm, f, weight);
  }

  std::cout << "EWSN weights: " << vx1 << " " << vx2 << " " << vy1 << " " << vy2 << std::endl;

  return true;
}

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: "
              << argv[0] << "\n"
              << "\tinput_filename\n"
              << "\t[output_filename] (PLY)\n"
              << "\t[weights_filename.txt]\n"
              << std::endl;
    std::cerr << "Input format: any format readable with CGAL::IO::read_polygon_mesh()" << std::endl;
    std::cerr << "Output format: PLY" << std::endl;
    std::cerr << "Weight format:" << std::endl;
    std::cerr << "  x1: val" << std::endl;
    std::cerr << "  x2: val" << std::endl;
    std::cerr << "  y1: val" << std::endl;
    std::cerr << "  y2: val" << std::endl;
    std::cerr << "  bottom: val (optional line)" << std::endl;
    std::cerr << "  top: val (optional line)" << std::endl;
    return EXIT_FAILURE;
  }

  const char* input_filename = argv[1];
  const char* output_filename = (argc > 2) ? argv[2] : "out.ply";
  const char* weights_filename = (argc > 3) ? argv[3] : nullptr;

  std::cout << "in: " << input_filename << std::endl;
  std::cout << "out: " << output_filename << std::endl;
  std::cout << "weight: " << weights_filename << std::endl;

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(input_filename, sm)) {
    std::cerr << "Error: failed to read input" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh: " << num_vertices(sm) << " NV " << num_faces(sm) << " NF" << std::endl;

  if(!assign_weights(sm, weights_filename))
    return EXIT_FAILURE;

  std::ofstream out(output_filename);
  if(!out) {
    std::cerr << "Error: failed to create output file" << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::IO::write_PLY(out, sm,
                      CGAL::parameters::use_binary_mode(false)
                                       .stream_precision(17));


  return EXIT_SUCCESS;
}
