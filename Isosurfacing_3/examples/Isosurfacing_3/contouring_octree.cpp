#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>
#include <CGAL/Isosurfacing_3/Octree_partition.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <cmath>
#include <iostream>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Vector = typename Kernel::Vector_3;
using Point = typename Kernel::Point_3;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

using Octree = CGAL::Octree<Kernel, std::vector<Point> >;
using Values = CGAL::Isosurfacing::Value_function_3<Octree>;
using Gradients = CGAL::Isosurfacing::Gradient_function_3<Octree>;
using MC_Domain = CGAL::Isosurfacing::Marching_cubes_domain_3<Octree, Values>;
using Domain = CGAL::Isosurfacing::Dual_contouring_domain_3<Octree, Values, Gradients>;

namespace IS = CGAL::Isosurfacing;

auto sphere_function = [](const Point& p) -> FT
{
  return std::sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
};

auto sphere_gradient = [](const Point& p) -> Vector
{
  const Vector g = p - CGAL::ORIGIN;
  return g / std::sqrt(g.squared_length());
};

auto blobby_function = [](const Point& p) -> FT
{
  return std::exp(-1.5 * ((p.x() - 0.2) * (p.x() - 0.2) + (p.y() - 0.2) * (p.y() - 0.2) + (p.z() - 0.2) * (p.z() - 0.2))) +
         std::exp(-1.5 * ((p.x() + 0.2) * (p.x() + 0.2) + (p.y() + 0.2) * (p.y() + 0.2) + (p.z() + 0.2) * (p.z() + 0.2))) +
         std::exp(-1.5 * ((p.x() - 0.4) * (p.x() - 0.4) + (p.y() + 0.4) * (p.y() + 0.4) + (p.z() - 0.4) * (p.z() - 0.4))) +
         std::exp(-6 * ((p.x() - 0.1) * (p.x() - 0.1) + (p.y() - 0.1) * (p.y() - 0.1))) +
         std::exp(-6 * ((p.y() + 0.1) * (p.y() + 0.1) + (p.z() + 0.1) * (p.z() + 0.1))) +
         std::exp(-6 * ((p.x() + 0.1) * (p.x() + 0.1) + (p.z() - 0.1) * (p.z() - 0.1))) -
         0.3;
};

auto blobby_gradient = [](const Point& p) -> Vector
{
  const FT g1 = -3 * std::exp(-1.5 * ((p.x() - 0.2) * (p.x() - 0.2) + (p.y() - 0.2) * (p.y() - 0.2) + (p.z() - 0.2) * (p.z() - 0.2)));
  const FT g2 = -3 * std::exp(-1.5 * ((p.x() + 0.2) * (p.x() + 0.2) + (p.y() + 0.2) * (p.y() + 0.2) + (p.z() + 0.2) * (p.z() + 0.2)));
  const FT g3 = -3 * std::exp(-1.5 * ((p.x() - 0.4) * (p.x() - 0.4) + (p.y() + 0.4) * (p.y() + 0.4) + (p.z() - 0.4) * (p.z() - 0.4)));
  const FT g4 = -12 * std::exp(-6 * ((p.x() - 0.1) * (p.x() - 0.1) + (p.y() - 0.1) * (p.y() - 0.1)));
  const FT g5 = -12 * std::exp(-6 * ((p.y() + 0.1) * (p.y() + 0.1) + (p.z() + 0.1) * (p.z() + 0.1)));
  const FT g6 = -12 * std::exp(-6 * ((p.x() + 0.1) * (p.x() + 0.1) + (p.z() - 0.1) * (p.z() - 0.1)));

  return Vector(g1 * (p.x() - 0.2) + g2 * (p.x() + 0.2) + g3 * (p.x() - 0.4) + g4 * (p.x() - 0.1) + g6 * (p.x() + 0.1),
                g1 * (p.y() - 0.2) + g2 * (p.y() + 0.2) + g3 * (p.y() + 0.4) + g4 * (p.y() - 0.1) + g5 * (p.y() + 0.1),
                g1 * (p.z() - 0.2) + g2 * (p.z() + 0.2) + g3 * (p.z() - 0.4) + g5 * (p.z() + 0.1) + g6 * (p.z() - 0.1));
};

// This is a naive refinement that is adapted to the isosurface:
// This refines:
// - at the minimum till minimum depth
// - at the maximum till maximum depth
// - we split if the the isovalue goes through the voxel, i.e. if not all vertices of the cell
//   are on the same side of the isosurface defined by a function
// It's not a great refinement technique because the surface can enter and leave a cell
// without involving the cell's vertex. In practice, that means a hole if at nearby adjacent
// cells the voxels did get refined and registered the surface.
struct Refine_around_isovalue
{
  std::size_t min_depth_;
  std::size_t max_depth_;
  std::function<FT(const Point&)> function_;
  FT isovalue_;

  Refine_around_isovalue(std::size_t min_depth,
                         std::size_t max_depth,
                         std::function<FT(const Point&)> function,
                         FT isovalue)
    : min_depth_(min_depth),
      max_depth_(max_depth),
      function_(function),
      isovalue_(isovalue)
  {}

  bool operator()(const Octree::Node_index& ni, const Octree& octree) const
  {
    // Ensure minimum depth refinement
    if (octree.depth(ni) < min_depth_)
      return true;

    // Stop refinement at maximum depth
    if (octree.depth(ni) >= max_depth_)
      return false;

    // Get the bounding box of the node
    auto bbox = octree.bbox(ni);

    // Evaluate the function at the corners of the bounding box
    std::array<FT, 8> corner_values;
    int index = 0;
    for (FT x : {bbox.xmin(), bbox.xmax()})
      for (FT y : {bbox.ymin(), bbox.ymax()})
        for (FT z : {bbox.zmin(), bbox.zmax()})
          corner_values[index++] = function_(Point(x, y, z));

    // Check if the function values cross the isovalue
    bool has_positive = false, has_negative = false;
    for (const auto& value : corner_values)
    {
      if (value > isovalue_)
        has_positive = true;
      if (value < isovalue_)
        has_negative = true;
      if (has_positive && has_negative)
        return true; // Refine if the isosurface intersects the voxel
    }

    return false; // No refinement needed
  }
};

// This is a refinement that is NOT adapted to the isosurface
struct Refine_one_eighth
{
  std::size_t min_depth_;
  std::size_t max_depth_;

  std::size_t octree_dim_;

  Refine_one_eighth(std::size_t min_depth,
                    std::size_t max_depth)
    : min_depth_(min_depth),
      max_depth_(max_depth)
  {
    octree_dim_ = std::size_t(1) << max_depth_;
  }

  Octree::Global_coordinates uniform_coordinates(const Octree::Node_index& node_index, const Octree& octree) const
  {
    auto coords = octree.global_coordinates(node_index);
    const std::size_t depth_factor = std::size_t(1) << (max_depth_ - octree.depth(node_index));
    for(int i=0; i < 3; ++i)
      coords[i] *= uint32_t(depth_factor);

    return coords;
  }

  bool operator()(const Octree::Node_index& ni, const Octree& octree) const
  {
    if(octree.depth(ni) < min_depth_)
      return true;

    if(octree.depth(ni) == max_depth_)
      return false;

    auto leaf_coords = uniform_coordinates(ni, octree);

    if(leaf_coords[0] >= octree_dim_ / 2)
      return false;

    if(leaf_coords[1] >= octree_dim_ / 2)
      return false;

    if(leaf_coords[2] >= octree_dim_ / 2)
      return false;

    return true;
  }
};

template <typename Splitter>
void run_DC_octree(const CGAL::Bbox_3 bbox,
                   const Splitter& split_predicate,
                   const std::function<FT(const Point&)> function,
                   const std::function<Vector(const Point&)> gradient,
                   const FT isovalue,
                   const std::string& name)
{
  std::vector<Point> bbox_points { { bbox.xmin(), bbox.ymin(), bbox.zmin() },
                                   { bbox.xmax(), bbox.ymax(), bbox.zmax() } };

  Octree octree(bbox_points);
  octree.refine(split_predicate);


  auto leaf_range = octree.traverse(CGAL::Orthtrees::Leaves_traversal<Octree>(octree));
  std::size_t leaf_counter = std::distance(leaf_range.begin(), leaf_range.end());

  std::cout << "octree has " << leaf_counter << " leaves" << std::endl;

  // fill up values and gradients
  Values values { function, octree };
  Gradients gradients { gradient, octree };
  Domain domain { octree, values, gradients };

  // output containers
  Point_range points;
  Polygon_range triangles;

  std::cout << "Running Dual Contouring with isovalue = " << isovalue << std::endl;

  // run Dual Contouring
  IS::dual_contouring<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles,
                                                       CGAL::parameters::do_not_triangulate_faces(true)
                                                                        .constrain_to_cell(false));

  std::cout << "Output #vertices (DC): " << points.size() << std::endl;
  std::cout << "Output #triangles (DC): " << triangles.size() << std::endl;

  std::ofstream oo("octree_DC_" + name + ".polylines.txt");
  oo.precision(17);
  octree.dump_to_polylines(oo);

  CGAL::IO::write_polygon_soup("DC_" + name + ".off", points, triangles);
}

template <typename Splitter>
void run_MC_octree(const CGAL::Bbox_3 bbox,
                   const Splitter& split_predicate,
                   const std::function<FT(const Point&)> function,
                   const FT isovalue,
                   const std::string& name)
{
  std::vector<Point> bbox_points { { bbox.xmin(), bbox.ymin(), bbox.zmin() },
  { bbox.xmax(), bbox.ymax(), bbox.zmax() } };

  Octree octree(bbox_points);
  octree.refine(split_predicate);

  Values values { function, octree };

  Point_range points;
  Polygon_range triangles;
  MC_Domain mcdomain { octree, values };

  std::cout << "Running MC" << std::endl;

  CGAL::Isosurfacing::marching_cubes<CGAL::Parallel_if_available_tag>(mcdomain, isovalue, points, triangles);

  std::cout << "Output #vertices (MC): " << points.size() << std::endl;
  std::cout << "Output #triangles (MC): " << triangles.size() << std::endl;

  std::ofstream oo("octree_MC_" + name + ".polylines.txt");
  oo.precision(17);
  octree.dump_to_polylines(oo);

  CGAL::IO::write_polygon_soup("MC_" + name + ".off", points, triangles);
}

// Whether you are using MC, TMC, or DC, there is no guarantee for an octree:
// it should behave well if your nodes are split with a uniform size around the surface,
// but it is sure to produce cracks if you have varying depths around the isosurface.
int main(int argc, char** argv)
{
  const FT isovalue = (argc > 1) ? std::stod(argv[1]) : 0.3;

  const CGAL::Bbox_3 bbox { -1., -1., -1.,  1., 1., 1. };

  Refine_one_eighth one_eight_splitter(3, 5);
  run_DC_octree(bbox, one_eight_splitter, sphere_function, sphere_gradient, isovalue, "one_eight");

  // This is
  Refine_around_isovalue isovalue_splitter(3, 5, sphere_function, isovalue);
  run_DC_octree(bbox, isovalue_splitter, sphere_function, sphere_gradient, isovalue, "sphere_adapted");

  Refine_around_isovalue isvalue_splitter_2(3, 5, blobby_function, isovalue);
  run_DC_octree(bbox, isvalue_splitter_2, blobby_function, blobby_gradient, isovalue, "blobby_adapted");

  // to illustrate cracks
  run_MC_octree(bbox, isovalue_splitter, sphere_function, isovalue, "adapted");

  run_MC_octree(bbox, one_eight_splitter, sphere_function, isovalue, "one_eight");

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
