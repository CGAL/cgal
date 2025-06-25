#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>

#include <CGAL/Isosurfacing_3/internal/tables.h>
#include <CGAL/Isosurfacing_3/internal/Cell_type.h>

#include <openvdb/openvdb.h>
#include <openvdb/io/File.h>
#include <openvdb/io/Stream.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/RayIntersector.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <vector>
#include <string>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

namespace IS = CGAL::Isosurfacing;

class OpenVDB_Edge_intersection
{
public:
  OpenVDB_Edge_intersection(openvdb::FloatGrid::Ptr grid,
                            const float isovalue)
  : _grid(grid),
    _intersector(*grid, isovalue)
  {}

public:
  template <typename Domain> // == Isosurfacing_domain_3 or similar
  bool operator()(const typename Domain::Geom_traits::Point_3& p_0,
                  const typename Domain::Geom_traits::Point_3& p_1,
                  const typename Domain::Geom_traits::FT /*val_0*/,
                  const typename Domain::Geom_traits::FT /*val_1*/,
                  const Domain& /*domain*/,
                  const typename Domain::Geom_traits::FT/* isovalue*/,
                  typename Domain::Geom_traits::Point_3& p) const
  {
    const openvdb::Vec3R eye(p_0.x(), p_0.y(), p_0.z());
    const openvdb::Vec3R direction = openvdb::Vec3R(p_1.x() - p_0.x(),
                                                    p_1.y() - p_0.y(),
                                                    p_1.z() - p_0.z());
    const openvdb::math::Ray<openvdb::Real> ray(eye, direction.unit());
    openvdb::Vec3R I;
    openvdb::Vec3R normal;
    openvdb::Real t;
    if (_intersector.intersectsWS(ray, I, normal, t) && t > 0.0)
    {
      // Only consider intersections within the segment bounds.
      const openvdb::Vec3R vI = I - eye;
      if (vI.lengthSqr() <= direction.lengthSqr())
      {
        p = typename Domain::Geom_traits::Point_3(I.x(), I.y(), I.z());
        return true;
      }
    }
    return false;
  }

private:
  openvdb::FloatGrid::Ptr _grid;
  openvdb::tools::LevelSetRayIntersector<openvdb::FloatGrid> _intersector;
};

class OpenVDB_partition
{
public:
  using Geom_traits = Kernel;

  OpenVDB_partition(openvdb::FloatGrid::Ptr grid)
    : _grid(grid)
    , _gt()
  {}

  const Geom_traits& geom_traits() const { return _gt; }

  openvdb::FloatGrid::Ptr grid() const { return _grid; }

private:
  openvdb::FloatGrid::Ptr _grid;
  Geom_traits _gt;
};


namespace CGAL {
namespace Isosurfacing {

template<>
struct partition_traits<OpenVDB_partition>
{
public:
  using vertex_descriptor = openvdb::Coord; //indices (i,j,k)
  // identifies a cell by its corner vertex with the smallest (i, j, k) index
  using cell_descriptor = openvdb::Coord;
  // identifies an edge by its starting vertex (i, j, k) and the direction x -> 0, y -> 1, z -> 2
  struct edge_descriptor
  {
    vertex_descriptor source;
    int direction;

    bool operator==(const edge_descriptor& other) const {
      return source == other.source && direction == other.direction;
    }
  };

  static constexpr Cell_type CELL_TYPE = CUBICAL_CELL;
  static constexpr std::size_t VERTICES_PER_CELL = 8;
  static constexpr std::size_t EDGES_PER_CELL = 12;

  using Edge_vertices = std::array<vertex_descriptor, 2>;
  using Cells_incident_to_edge = std::array<cell_descriptor, 4>;
  using Cell_vertices = std::array<vertex_descriptor, VERTICES_PER_CELL>;
  using Cell_edges = std::array<edge_descriptor, EDGES_PER_CELL>;

  static Point point(const vertex_descriptor& v,
                     const OpenVDB_partition& partition)
  {
    const auto& grid = partition.grid();
    const openvdb::Vec3R vr = grid->indexToWorld(v);
    return Point(vr.x(), vr.y(), vr.z());
  }

  static Edge_vertices incident_vertices(const edge_descriptor& e, const OpenVDB_partition&)
  {
    vertex_descriptor src = e.source;// start vertex
    vertex_descriptor tgt = src;// end vertex
    tgt[e.direction] += 1;// one position further in the direction of the edge
    return {src, tgt};
  }

  static Cells_incident_to_edge incident_cells(const edge_descriptor& e, const OpenVDB_partition&)
  {
    // lookup the neighbor cells relative to the edge
    const int local = internal::Cube_table::edge_store_index[e.direction];
    auto neighbors = internal::Cube_table::edge_to_voxel_neighbor[local];

    Cells_incident_to_edge cite;
    for (std::size_t i = 0; i < 4; ++i) { //for each cell
      cite[i] = cell_descriptor(e.source.x() + neighbors[i][0],
                                e.source.y() + neighbors[i][1],
                                e.source.z() + neighbors[i][2]);
    }
    return cite;
  }

  static Cell_vertices cell_vertices(const cell_descriptor& c, const OpenVDB_partition&)
  {
    Cell_vertices cv;
    for (std::size_t i = 0; i < cv.size(); ++i)
    { //for each vertex in cell_vertices (8)
      cv[i] = vertex_descriptor(
        c.x() + internal::Cube_table::local_vertex_position[i][0],
        c.y() + internal::Cube_table::local_vertex_position[i][1],
        c.z() + internal::Cube_table::local_vertex_position[i][2]);
    }
    return cv;
  }

  static Cell_edges cell_edges(const cell_descriptor& c, const OpenVDB_partition&)
  {
    Cell_edges ce;
    for (std::size_t i = 0; i < ce.size(); ++i) //for each edge
    {
      // lookup the relative edge indices and offset them by the cell position
      ce[i] = edge_descriptor{
        vertex_descriptor(c.x() + internal::Cube_table::global_edge_id[i][0],
                          c.y() + internal::Cube_table::global_edge_id[i][1],
                          c.z() + internal::Cube_table::global_edge_id[i][2]),
        internal::Cube_table::global_edge_id[i][3]};
    }
    return ce;
  }

  template <typename ConcurrencyTag = CGAL::Sequential_tag, typename Functor>
  static void for_each_vertex(Functor& f,
                              const OpenVDB_partition& partition,
                              ConcurrencyTag)
  {
    const auto& grid = partition.grid();
    for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter.test(); ++iter) {
      f(iter.getCoord());
    }
  }

  template <typename ConcurrencyTag = CGAL::Sequential_tag, typename Functor>
  static void for_each_edge(Functor& f,
                            const OpenVDB_partition& partition,
                            ConcurrencyTag)
  {
    const auto grid = partition.grid();
    for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter.test(); ++iter)
    {
      // all three edges starting at vertex (i, j, k)
      f(edge_descriptor{iter.getCoord(), 0});
      f(edge_descriptor{iter.getCoord(), 1});
      f(edge_descriptor{iter.getCoord(), 2});
    }
  }

  template <typename ConcurrencyTag = CGAL::Sequential_tag, typename Functor>
  static void for_each_cell(Functor& f,
                            const OpenVDB_partition& partition,
                            ConcurrencyTag)
  {
    const auto grid = partition.grid();
    for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter.test(); ++iter) {
      f(iter.getCoord());
    }
  }
};

}
}

namespace std{
  template<>
  struct hash<CGAL::Isosurfacing::partition_traits<OpenVDB_partition>::edge_descriptor>
  {
    using edge_descriptor
      = CGAL::Isosurfacing::partition_traits<OpenVDB_partition>::edge_descriptor;
    using vertex_descriptor
      = CGAL::Isosurfacing::partition_traits<OpenVDB_partition>::vertex_descriptor;

    auto operator()(const edge_descriptor& e) const -> size_t
    {
      return std::hash<vertex_descriptor>{}(e.source) ^ std::hash<int>{}(e.direction);
    }
  };
}  // namespace std

class OpenVDB_value_field
{
public:
  using FT = FT;
  using Point_3 = Point;
  using vertex_descriptor = openvdb::Coord;

  OpenVDB_value_field(openvdb::FloatGrid::Ptr grid)
    : _grid(grid)
    , _accessor(grid->getAccessor())
  {}

  FT operator()(const Point_3& p) const
  {
    // calculate coordinates of min index
    openvdb::Vec3d ijk_p = _grid->worldToIndex(openvdb::Vec3d(p.x(), p.y(), p.z()));
    openvdb::Coord ijk(static_cast<int>(ijk_p.x()),
                       static_cast<int>(ijk_p.y()),
                       static_cast<int>(ijk_p.z()));

    openvdb::Coord ijk_z  (ijk[0] + 0, ijk[1] + 0, ijk[2] + 1);
    openvdb::Coord ijk_y  (ijk[0] + 0, ijk[1] + 1, ijk[2] + 0);
    openvdb::Coord ijk_yz (ijk[0] + 0, ijk[1] + 1, ijk[2] + 1);
    openvdb::Coord ijk_x  (ijk[0] + 1, ijk[1] + 0, ijk[2] + 0);
    openvdb::Coord ijk_xz (ijk[0] + 1, ijk[1] + 0, ijk[2] + 1);
    openvdb::Coord ijk_xy (ijk[0] + 1, ijk[1] + 1, ijk[2] + 0);
    openvdb::Coord ijk_xyz(ijk[0] + 1, ijk[1] + 1, ijk[2] + 1);

    auto min_p = _grid->indexToWorld(ijk);
    auto ijk_z_p = _grid->indexToWorld(ijk_z);
    auto ijk_y_p = _grid->indexToWorld(ijk_y);
    auto ijk_yz_p = _grid->indexToWorld(ijk_yz);
    auto ijk_x_p = _grid->indexToWorld(ijk_x);
    auto ijk_xz_p = _grid->indexToWorld(ijk_xz);
    auto ijk_xy_p = _grid->indexToWorld(ijk_xy);
    auto ijk_xyz_p = _grid->indexToWorld(ijk_xyz);

    std::cout << "p: " << p << std::endl;
    std::cout << "ijk_p: " << ijk_p << std::endl;
    std::cout << "min_p: " << min_p << std::endl;
    std::cout <<  _grid->indexToWorld(ijk) << std::endl;
    std::cout <<  _grid->indexToWorld(ijk_z) << std::endl;
    std::cout <<  _grid->indexToWorld(ijk_y) << std::endl;
    std::cout <<  _grid->indexToWorld(ijk_yz) << std::endl;
    std::cout <<  _grid->indexToWorld(ijk_x) << std::endl;
    std::cout <<  _grid->indexToWorld(ijk_xz) << std::endl;
    std::cout <<  _grid->indexToWorld(ijk_xy) << std::endl;
    std::cout <<  _grid->indexToWorld(ijk_xyz) << std::endl;

    // interpolation factors between 0 and 1
    const FT f_i = CGAL::abs((p.x() - min_p.x()) / (ijk_x_p.x() - min_p.x())); // @fixme? negative voxels
    const FT f_j = CGAL::abs((p.y() - min_p.y()) / (ijk_y_p.y() - min_p.y()));
    const FT f_k = CGAL::abs((p.z() - min_p.z()) / (ijk_z_p.z() - min_p.z()));
    std::cout << "f_i = " << f_i << std::endl;
    std::cout << "f_j = " << f_j << std::endl;
    std::cout << "f_k = " << f_k << std::endl;
    CGAL_assertion(f_i >= 0 && f_j >= 0 && f_k >= 0 && f_i <= 1. && f_j <= 1. && f_k <= 1.);

    // read the value at all 8 corner points
    const FT g000 = _accessor.getValue(ijk);
    const FT g001 = _accessor.getValue(ijk_z);
    const FT g010 = _accessor.getValue(ijk_y);
    const FT g011 = _accessor.getValue(ijk_yz);
    const FT g100 = _accessor.getValue(ijk_x);
    const FT g101 = _accessor.getValue(ijk_xz);
    const FT g110 = _accessor.getValue(ijk_xy);
    const FT g111 = _accessor.getValue(ijk_xyz);

    // interpolate along all axes by weighting the corner points
    const FT lambda000 = (FT(1) - f_i) * (FT(1) - f_j) * (FT(1) - f_k);
    const FT lambda001 = (FT(1) - f_i) * (FT(1) - f_j) * f_k;
    const FT lambda010 = (FT(1) - f_i) * f_j * (FT(1) - f_k);
    const FT lambda011 = (FT(1) - f_i) * f_j * f_k;
    const FT lambda100 = f_i * (FT(1) - f_j) * (FT(1) - f_k);
    const FT lambda101 = f_i * (FT(1) - f_j) * f_k;
    const FT lambda110 = f_i * f_j * (FT(1) - f_k);
    const FT lambda111 = f_i * f_j * f_k;

    // add weighted corners
    return g000 * lambda000 + g001 * lambda001 +
           g010 * lambda010 + g011 * lambda011 +
           g100 * lambda100 + g101 * lambda101 +
           g110 * lambda110 + g111 * lambda111;
  }

  FT operator()(const vertex_descriptor& v) const
  {
    return _accessor.getValue(v);
  }

private:
  openvdb::FloatGrid::Ptr _grid;
  openvdb::FloatGrid::Accessor _accessor;
};

using Partition = OpenVDB_partition;
using Values = OpenVDB_value_field;
using Edge_intersection = OpenVDB_Edge_intersection;

void run_marching_cubes(openvdb::FloatGrid::Ptr grid,
                        const float isovalue,
                        const bool use_tcm)
{
  using Domain = IS::Marching_cubes_domain_3<Partition, Values, Edge_intersection>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Marching Cubes with isovalue = " << isovalue
    << ", tcm = " << std::boolalpha << use_tcm << std::endl;

  Partition partition(grid);
  Values values(grid);
  Edge_intersection intersection_oracle(grid, isovalue);
  Domain domain = IS::create_marching_cubes_domain_3(partition, values, intersection_oracle);

  // run marching cubes isosurfacing
  Point_range points;
  Polygon_range triangles;

  IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles,
    CGAL::parameters::use_topologically_correct_marching_cubes(use_tcm));

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;
  const std::string outfile = use_tcm
                        ? "marching_cubes_openvdb_tcm.off"
                        : "marching_cubes_openvdb.off";
  CGAL::IO::write_polygon_soup(outfile, points, triangles);
}

void run_dual_contouring(openvdb::FloatGrid::Ptr grid,
                         const FT isovalue)
{
  using Gradients = CGAL::Isosurfacing::Finite_difference_gradient_3<Kernel>;
  using Domain = IS::Dual_contouring_domain_3<Partition, Values, Gradients, Edge_intersection>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring with isovalue = " << isovalue << std::endl;

  // fill up values and gradients
  Partition partition(grid);
  Values values(grid);
  const FT step = grid->voxelSize().length() * 0.01; // finite difference step
  Gradients gradients{ values, step };
  Edge_intersection intersection_oracle(grid, isovalue);

  Domain domain = IS::create_dual_contouring_domain_3(partition, values, gradients, intersection_oracle);

  Point_range points;
  Polygon_range triangles;

  // run dual contouring isosurfacing
  IS::dual_contouring<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles,
                                                       CGAL::parameters::do_not_triangulate_faces(true));

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;
  CGAL::IO::write_polygon_soup("dual_contouring_openvdb.off", points, triangles);
}


int main(int argc, char** argv)
{
  const std::string input_filename = argv[1];

  std::ifstream ifs(input_filename, std::ios::binary);
  if (!ifs)
    return EXIT_FAILURE;

  openvdb::initialize();
  openvdb::io::Stream stream(ifs);
  openvdb::GridPtrVecPtr grids = stream.getGrids();
  openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(grids->front());

  run_marching_cubes(grid, 0./*isovalue*/, true/*topologically correct MC*/);
  run_marching_cubes(grid, 0./*isovalue*/, false/*topologically correct MC*/);
  run_dual_contouring(grid, 0./*isovalue*/);

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}