// Domain-decomposition benchmark for CGAL tetrahedral remeshing.
//
// Strategy:
//   1. Partition the input mesh's DOMAIN cells into spatial buckets (uniform
//      DIV_X x DIV_Y x DIV_Z grid, half-open bboxes). A domain cell is assigned
//      to bucket b iff ALL 4 of its vertices lie in b's bbox; domain cells whose
//      vertices span >=2 buckets are "straddlers" and are frozen (selected in no
//      bucket). Cells OUTSIDE the domain (convex-hull filler, subdomain 0) are
//      never bucketed and never remeshed.
//   2. Give bucket b its own copy of the triangulation with only b selected.
//      NOTE: this is a deep copy of the WHOLE mesh, not a standalone bucket T3 --
//      see extract_bucket_copy for why a true extraction is not possible here.
//   3. Remesh each bucket T3 independently and in parallel (8 buckets, 1
//      thread each emergently via TBB), with remesh_boundaries(false) so the
//      bucket boundary (interfaces + skin) is frozen. Each bucket gets its
//      own spatial lock grid so the 8 isolated T3s do not false-contend on a
//      shared static grid.
//   4. Merge: reassemble a whole-mesh T3 from the 8 remeshed buckets plus the
//      untouched frozen straddlers. Frozen boundary vertices are bit-identical
//      across buckets and the original, so coordinate-dedup re-glues the
//      interfaces exactly (in-memory Medit roundtrip).
//   5. Re-partition and repeat. There is NO final global pass: the seams are
//      covered by interface migration, because a cell frozen on a bucket
//      interface in one round sits in a bucket INTERIOR in the next.
//
// The partition is produced by ParMETIS_V3_AdaptiveRepart, warm-started from
// the previous round's partition (see assign_buckets_parmetis), with serial
// METIS seeding the very first round.
//
// Based on 9ac7560 (gsoc2025 baseline with the refactored Elementary_remesher).

#define USE_REFACTORED_TETRAHEDRAL_REMESHING
// Freeze the partition interface independently of remesh_boundaries, so the
// buckets can run with protect_boundaries = false and actually remesh the
// MODEL SURFACE while their shared interface stays bit-identical.
//#define CGAL_TETRAHEDRAL_REMESHING_PROTECT_SUBDOMAIN_INTERFACES
#define CGAL_CONCURRENT_TETRAHEDRAL_REMESHING

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "benchmark_refactored_tetrahedral_remeshing_macros_config.h"
#include "benchmark_tetrahedral_remeshing_common.h"
#include "mesh_quality.h"

#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Real_timer.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>

#ifdef CGAL_DD_USE_METIS
#include <metis.h>
#endif

#ifdef CGAL_DD_USE_PARMETIS
// parmetis.h pulls in metis.h (idx_t / real_t / METIS_OK), so this header is
// self-sufficient even when CGAL_DD_USE_METIS is undefined.
//
// Skip the deprecated MPI C++ bindings: mpi.h drags them in by default and
// they are not in the library MPI::MPI_C links, so every one of their inline
// members becomes an undefined reference at link time. Only the C API is used
// here.
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#define MPI_NO_CPPBIND 1
#include <mpi.h>
#include <parmetis.h>
#endif

#if defined(CGAL_DD_USE_METIS) || defined(CGAL_DD_USE_PARMETIS)
#include <unordered_map>
#endif

#include <tbb/global_control.h>
#include <tbb/task_group.h>

#include <nlohmann/json.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <tuple>
#include <memory>
#include <atomic>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <set>

#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
  #define Concurrency_tag CGAL::Parallel_tag
#else
  #define Concurrency_tag CGAL::Sequential_tag
#endif

using K    = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb   = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
using Cb   = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
using T3   = CGAL::Triangulation_3<K, CGAL::Triangulation_data_structure_3<Vb, Cb, Concurrency_tag>>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<T3, int, int>;
using LockDS = T3::Lock_data_structure;
using Point_3 = K::Point_3;

// A Sequential_tag triangulation, used for the lock-free bucket phase.
// Adaptive_remesher_type_generator dispatches ExecutionPolicy on
// Tr::Concurrency_tag (elementary_remesh_impl.h:48-50): a Sequential_tag Tr
// selects ElementaryOperationExecutionSequential, so the elementary operations
// run single-threaded on one worker and lock_zone() is never called -- no
// spatial lock grid at all. Parallelism then comes purely from the outer
// tbb::task_group over buckets, whose triangulations are disjoint objects.
//
// (Contrast the "par" bucket path, which gives each bucket its own LockDS and
// runs the full Parallel_tag remesher inside the task_group -- so buckets are
// disjoint but each bucket still pays for the lock grid internally.)
using T3seq = CGAL::Triangulation_3<K, CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag>>;

// Number of parts. ParMmg partitions the DUAL GRAPH, so every element belongs
// to exactly one part and the frozen entity is the interface FACE set, not a
// layer of whole cells. We keep that property with a geometric partition by
// cell CENTROID: every domain cell is owned by exactly one bucket, so no cell
// is left unremeshed and the "straddler" concept disappears entirely.
static constexpr int NUM_BUCKETS = 4;

// subdomain_index encoding. The input C3t3 triangulates the whole convex hull:
// cells OUTSIDE the domain already carry subdomain_index == 0, and must keep it
// so hull filler is never promoted into a bucket domain.
//
//   SD_EXTERIOR        convex-hull filler, never bucketed, never remeshed
//   1 .. NUM_BUCKETS   owning bucket of a domain cell
static constexpr int SD_EXTERIOR = 0;

// ---------------------------------------------------------------------------
// A tetrahedron as 4 points + a subdomain ref (written to the .mesh file).
// ---------------------------------------------------------------------------
struct TetSoup
{
  std::array<Point_3, 4> pts;
  int ref;
};

// A surface triangle plus its patch index, written to the .mesh Triangles
// section so the merge preserves the surface patches (and therefore the
// feature edges derived from them).
struct TriSoup
{
  std::array<Point_3, 3> pts;
  int ref;
};

// ---------------------------------------------------------------------------
// write_tet_soup_medit: write a list of tetrahedra to a stream in Medit .mesh
// format, deduplicating vertices by EXACT coordinate (bit-identical doubles).
// Coordinates are written with 17 significant digits so that any frozen
// boundary vertex round-trips bit-identically and re-dedups across buckets.
// ---------------------------------------------------------------------------
static void write_tet_soup_medit(std::ostream& os,
                                 const std::vector<TetSoup>& tets,
                                 const std::vector<TriSoup>& tris = {})
{
  std::map<std::tuple<double, double, double>, int> vmap; // coord -> 1-based idx
  std::vector<std::tuple<double, double, double>> verts;
  std::vector<std::array<int, 4>> tet_idx;
  std::vector<int> tet_ref;
  tet_idx.reserve(tets.size());
  tet_ref.reserve(tets.size());

  auto get_idx = [&](const Point_3& p) -> int {
    const auto key = std::make_tuple(CGAL::to_double(p.x()),
                                     CGAL::to_double(p.y()),
                                     CGAL::to_double(p.z()));
    auto [it, ins] = vmap.emplace(key, static_cast<int>(vmap.size()) + 1);
    if (ins) verts.push_back(key);
    return it->second;
  };

  for (const auto& t : tets)
  {
    std::array<int, 4> idx;
    for (int i = 0; i < 4; ++i) idx[i] = get_idx(t.pts[i]);
    tet_idx.push_back(idx);
    tet_ref.push_back(t.ref);
  }

  // Resolve triangle vertex indices BEFORE emitting anything: get_idx can insert
  // into vmap, which would invalidate a vertex count already written out.
  std::vector<std::array<int, 3>> tri_idx;
  tri_idx.reserve(tris.size());
  for (const TriSoup& t : tris)
    tri_idx.push_back({get_idx(t.pts[0]), get_idx(t.pts[1]), get_idx(t.pts[2])});

  // Optionally renumber the vertices in Hilbert order.
  //
  // Vertex indices above are assigned at FIRST ENCOUNTER while walking the tet
  // list, so the emitted numbering inherits whatever order the cells happened to
  // be collected in. A rebuilt bucket mesh therefore gets a different vertex and
  // cell order than `copy_tds` would have preserved, and collapse is an ORDERED
  // operation: it stable_sorts candidates by length, so the many near-equal
  // edges of a uniform-target remesh keep their collection order. That is the
  // remaining candidate for true extraction's late-phase divergence (it lags
  // whole-copy only once the run turns collapse-dominant). Spatial sort replaces
  // the incidental order with a deliberately coherent one.
  if (std::getenv("DD_SPATIAL_SORT") && std::string(std::getenv("DD_SPATIAL_SORT")) == "1")
  {
    std::vector<Point_3> pts;
    pts.reserve(verts.size());
    for (const auto& v : verts)
      pts.emplace_back(std::get<0>(v), std::get<1>(v), std::get<2>(v));

    std::vector<std::ptrdiff_t> order(pts.size());
    std::iota(order.begin(), order.end(), std::ptrdiff_t{0});
    using Sort_traits = CGAL::Spatial_sort_traits_adapter_3<K, Point_3*>;
    CGAL::spatial_sort(order.begin(), order.end(), Sort_traits(pts.data()));

    std::vector<int> remap(verts.size() + 1, 0);        // old 1-based -> new 1-based
    std::vector<std::tuple<double, double, double>> sorted;
    sorted.reserve(verts.size());
    for (std::size_t i = 0; i < order.size(); ++i)
    {
      sorted.push_back(verts[static_cast<std::size_t>(order[i])]);
      remap[static_cast<std::size_t>(order[i]) + 1] = static_cast<int>(i) + 1;
    }
    verts.swap(sorted);
    for (auto& t : tet_idx) for (int k = 0; k < 4; ++k) t[k] = remap[t[k]];
    for (auto& t : tri_idx) for (int k = 0; k < 3; ++k) t[k] = remap[t[k]];
  }

  os << std::setprecision(17);
  os << "MeshVersionFormatted 1\nDimension 3\n";
  os << "Vertices\n" << verts.size() << "\n";
  for (const auto& v : verts)
    os << std::get<0>(v) << ' ' << std::get<1>(v) << ' ' << std::get<2>(v) << " 0\n";
  if (!tris.empty())
  {
    os << "Triangles\n" << tri_idx.size() << "\n";
    for (std::size_t i = 0; i < tri_idx.size(); ++i)
      os << tri_idx[i][0] << ' ' << tri_idx[i][1] << ' '
         << tri_idx[i][2] << ' ' << tris[i].ref << "\n";
  }
  os << "Tetrahedra\n" << tet_idx.size() << "\n";
  for (std::size_t i = 0; i < tet_idx.size(); ++i)
    os << tet_idx[i][0] << ' ' << tet_idx[i][1] << ' '
       << tet_idx[i][2] << ' ' << tet_idx[i][3] << ' ' << tet_ref[i] << "\n";
  os << "End\n";
}

// ---------------------------------------------------------------------------
// read_t3_from_medit: read a .mesh from a stream into a fresh T3 via the
// CGAL Medit reader (build_triangulation_from_file preserves the exact tet
// connectivity and caps the boundary with infinite cells).
// ---------------------------------------------------------------------------
static bool read_t3_from_medit(std::istream& is, T3& tr)
{
  return CGAL::IO::read_MEDIT(is, tr);
}

// ---------------------------------------------------------------------------
// BucketGrid: an axis-aligned partition whose cut planes can ROTATE and SHIFT
// between iterations.
//
// This is ParMmg's interface migration. ParMmg re-partitions between adaptation
// iterations so that elements frozen on one iteration's parallel interface land
// in a partition INTERIOR on the next and finally get remeshed. Freezing the
// same region every iteration is what starves it: the band arrives at the end
// still in the split-dominant regime the bulk has long left.
//
// Rotating which pair of axes is cut moves the interfaces decisively; the shift
// term offsets the cut planes so a repeated axis set does not reproduce the
// same interface.
// ---------------------------------------------------------------------------
struct BucketGrid
{
  int    div[3];
  double lo[3], hi[3];
  double shift[3];   // cut-plane offset, in fractions of one bucket width
};

// All three sets give NUM_BUCKETS parts, so the bucket arrays stay fixed size.
static const int AXIS_SETS[3][3] = { {2,2,1}, {1,2,2}, {2,1,2} };

// Cut-plane offsets, in fractions of ONE BUCKET WIDTH (that is the unit of the
// `t` coordinate in axis_idx, not a fraction of the bbox).
//
// The seam only has to move further than the frozen band is thick -- a few cell
// diameters -- for last round's interface cells to land in a bucket interior.
// It must NOT move far enough to unbalance the partition, because the bucket
// phase is bounded by the LARGEST bucket, so imbalance is paid directly in
// wall-clock.
//
// 0.12 bucket widths = 6% of the extent for div=2, roughly 5 cell diameters at
// this target size, and moves each cut from 50% to 56% (a 1.27x worst-bucket
// penalty across two cut axes). The previous value of 0.5 moved the cut to 75%
// of the extent, giving a 9:3:3:1 split with one bucket EMPTY -- which made the
// bucket phase sequential and left DD slower than the sequential remesher.
// Alternating the sign keeps the seams oscillating about the midpoint instead
// of drifting one way. Period is 9 rounds (3 axis sets x 3 offsets).
static const double SHIFTS[3] = { 0.0, 0.12, -0.12 };

static BucketGrid grid_for_iteration(const CGAL::Bbox_3& bb, int it)
{
  BucketGrid g;
  const int* d = AXIS_SETS[it % 3];
  const double lo[3] = { bb.xmin(), bb.ymin(), bb.zmin() };
  const double hi[3] = { bb.xmax(), bb.ymax(), bb.zmax() };
  const double s = SHIFTS[(it / 3) % 3];
  for (int a = 0; a < 3; ++a)
  {
    g.div[a] = d[a]; g.lo[a] = lo[a]; g.hi[a] = hi[a];
    g.shift[a] = (d[a] > 1) ? s : 0.0;
  }
  return g;
}

// Half-open axis index; div == 1 collapses the axis (no cutting plane on it).
static int axis_idx(double v, const BucketGrid& g, int a)
{
  const int d = g.div[a];
  if (d <= 1) return 0;
  const double ext = g.hi[a] - g.lo[a];
  if (ext <= 0.0) return 0;
  const double t = (v - g.lo[a]) / ext * d - g.shift[a];
  return std::clamp(static_cast<int>(std::floor(t)), 0, d - 1);
}

static int bucket_of_point(const Point_3& p, const BucketGrid& g)
{
  const int ix = axis_idx(CGAL::to_double(p.x()), g, 0);
  const int iy = axis_idx(CGAL::to_double(p.y()), g, 1);
  const int iz = axis_idx(CGAL::to_double(p.z()), g, 2);
  return ix + g.div[0] * (iy + g.div[1] * iz);
}

// Owning bucket of a cell = bucket of its centroid. Every cell is owned, so the
// interface is a set of FACES rather than a frozen layer of cells.
template <typename CellHandle>
static int bucket_of_cell(CellHandle c, const BucketGrid& g)
{
  double x = 0, y = 0, z = 0;
  for (int i = 0; i < 4; ++i)
  {
    const auto& p = c->vertex(i)->point();
    x += CGAL::to_double(p.x());
    y += CGAL::to_double(p.y());
    z += CGAL::to_double(p.z());
  }
  return bucket_of_point(Point_3(0.25 * x, 0.25 * y, 0.25 * z), g);
}

#if defined(CGAL_DD_USE_METIS) || defined(CGAL_DD_USE_PARMETIS)
// ---------------------------------------------------------------------------
// build_dual_graph: the CSR dual graph of the DOMAIN cells, built the way
// ParMmg builds it (PMMG_graph_meshElts2metis). Nodes are domain cells; an
// edge joins two cells sharing a facet. Infinite and exterior neighbours are
// skipped, so the graph covers exactly the cells that get bucketed.
//
// Shared by the serial-METIS and the ParMETIS repartitioners: both need the
// same graph and only the partitioning call differs.
// ---------------------------------------------------------------------------
static void
build_dual_graph(C3t3& c3t3,
                 std::vector<T3::Cell_handle>& cells,
                 std::vector<idx_t>& xadj,
                 std::vector<idx_t>& adjncy)
{
  T3& tr = c3t3.triangulation();

  std::unordered_map<T3::Cell_handle, idx_t> id;
  cells.clear();
  cells.reserve(tr.number_of_finite_cells());
  for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
  {
    if (cit->subdomain_index() == SD_EXTERIOR) continue;
    id.emplace(cit, static_cast<idx_t>(cells.size()));
    cells.push_back(cit);
  }

  const idx_t n = static_cast<idx_t>(cells.size());
  xadj.assign(static_cast<std::size_t>(n) + 1, 0);
  adjncy.clear();
  adjncy.reserve(static_cast<std::size_t>(n) * 4);
  for (idx_t i = 0; i < n; ++i)
  {
    for (int f = 0; f < 4; ++f)
    {
      const auto it = id.find(cells[i]->neighbor(f));
      if (it != id.end()) adjncy.push_back(it->second);   // skips infinite/exterior
    }
    xadj[i + 1] = static_cast<idx_t>(adjncy.size());
  }
}

// ---------------------------------------------------------------------------
// commit_partition: write part[] back onto the mesh as subdomain indices
// (bucket b -> b+1) and return the per-bucket cell counts.
// ---------------------------------------------------------------------------
static std::array<std::size_t, NUM_BUCKETS>
commit_partition(const std::vector<T3::Cell_handle>& cells,
                 const std::vector<idx_t>& part)
{
  std::array<std::size_t, NUM_BUCKETS> counts{};
  for (std::size_t i = 0; i < cells.size(); ++i)
  {
    const int b = static_cast<int>(part[i]);
    cells[i]->set_subdomain_index(b + 1);
    ++counts[b];
  }
  return counts;
}
#endif

#ifdef CGAL_DD_USE_METIS
// ---------------------------------------------------------------------------
// assign_buckets_metis: partition the dual graph FROM SCRATCH with serial
// METIS, the way ParMmg does it (see parmmg/src/metis_pmmg.c).
//
// Partitioning the dual graph both balances the parts AND minimises the edge
// cut, i.e. the number of shared faces -- which here is the frozen interface.
// The geometric midpoint split does neither: it was 5.8:1 imbalanced on the
// very first round.
//
// The options are now tuned for SPEED rather than cut quality, because this is
// only the SEED partition for the ParMETIS adaptive repartitioner (and the
// fallback when ParMETIS is unavailable):
//   - METIS_OPTION_CONTIG is NOT set. Enforcing connected parts is the most
//     expensive option METIS has, and it is the one ParMmg pays for. A
//     disconnected part costs us only a slightly larger frozen interface.
//   - Kway at every nparts (ParMmg switches to recursive bisection below 8).
//     Kway refines the k-way partition directly; recursive bisection coarsens
//     and uncoarsens log2(nparts) times over.
//   - NITER=1, NCUTS=1: one refinement sweep, one trial partition.
//   - UFACTOR=50, i.e. stop refining at 5% imbalance instead of 3%.
//
// Set DD_METIS_QUALITY=1 to restore the slow ParMmg-faithful settings.
//
// The seed is varied per round so successive partitions differ and the
// interface migrates even when the mesh changes little.
// ---------------------------------------------------------------------------
static std::array<std::size_t, NUM_BUCKETS>
assign_buckets_metis(C3t3& c3t3, int round)
{
  std::vector<T3::Cell_handle> cells;
  std::vector<idx_t> xadj, adjncy;
  build_dual_graph(c3t3, cells, xadj, adjncy);

  idx_t n = static_cast<idx_t>(cells.size());
  idx_t ncon = 1, nparts = NUM_BUCKETS, objval = 0;
  std::vector<idx_t> part(static_cast<std::size_t>(n), 0);

  const char* qenv = std::getenv("DD_METIS_QUALITY");
  const bool quality = (qenv && std::string(qenv) == "1");

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_SEED] = round;
  if (quality)
    options[METIS_OPTION_CONTIG] = 1;
  else
  {
    options[METIS_OPTION_NITER]   = 1;
    options[METIS_OPTION_NCUTS]   = 1;
    options[METIS_OPTION_UFACTOR] = 50;   // imbalance = 1 + x/1000 => 5%
  }

  const int ier = (quality && nparts < 8)
    ? METIS_PartGraphRecursive(&n, &ncon, xadj.data(), adjncy.data(),
                               NULL, NULL, NULL, &nparts, NULL, NULL,
                               options, &objval, part.data())
    : METIS_PartGraphKway(&n, &ncon, xadj.data(), adjncy.data(),
                          NULL, NULL, NULL, &nparts, NULL, NULL,
                          options, &objval, part.data());
  if (ier != METIS_OK) benchmarking::fatal_error("METIS partitioning failed.");

  std::cout << "[metis] " << (quality ? "quality" : "fast")
            << "  edgecut=" << objval << "  graph nodes=" << n
            << "  edges=" << adjncy.size() / 2 << "\n";
  return commit_partition(cells, part);
}
#endif

#ifdef CGAL_DD_USE_PARMETIS
// ---------------------------------------------------------------------------
// assign_buckets_parmetis: REPARTITION the dual graph with ParMETIS's adaptive
// repartitioner, warm-started from the partition already on the mesh.
//
// This is the difference between partitioning and REpartitioning. METIS_Part*
// always starts from nothing: it coarsens the whole graph, partitions the
// coarsest level, and uncoarsens back -- every single round, on a mesh that
// keeps growing. ParMETIS_V3_AdaptiveRepart instead takes the CURRENT
// partition in part[] and only repairs it, diffusing load across the existing
// part boundaries. After one remeshing iteration the mesh has changed only
// locally (splits and collapses near the target length), so the previous
// partition is nearly right and the repair is much cheaper than a from-scratch
// partition. That is what makes it the fastest option here.
//
// The warm start needs the previous partition to SURVIVE the merge, which is
// why the merge writes each bucket's cells with ref = b+1 instead of a uniform
// 1 (see the merge loop): the rebuilt mesh carries its bucket tags.
//
// Everything runs in ONE MPI rank -- the whole graph is local, vtxdist = {0,n}
// -- because the bucket parallelism in this benchmark is TBB, not MPI.
// ParMETIS is used purely for its adaptive algorithm, not for distribution.
// nparts (NUM_BUCKETS) != npes (1), so the coupling must be
// PARMETIS_PSR_UNCOUPLED: part[] is an arbitrary initial partition, not one
// implied by the rank layout.
//
// ipc2redist trades edge cut against data redistribution. Small values pin the
// partition in place -- cheapest, but then the frozen interface never migrates
// and those cells are never remeshed, which is the whole point of
// re-partitioning. Large values let it move freely at the price of a bigger
// repair. Default 1000 (cut-dominant, interface free to migrate); override
// with DD_PARMETIS_IPC2REDIST.
//
// ---------------------------------------------------------------------------
// DOES NOT WORK AT ONE MPI RANK. Kept for the record, opt-in via
// DD_PARTITIONER=parmetis, NOT the default.
//
// ParMETIS 4.0.3 aborts inside an internal MPI_Recv with MPI_ERR_RANK --
// an abort, not a return code -- whenever the warm start is imbalanced enough
// that it actually has to REBALANCE. Its rebalancing path addresses parts by
// MPI rank, so with nparts=4 and npes=1 it receives from ranks 1..3, which do
// not exist. Established with a standalone probe against this exact library:
//
//   balanced warm start, any N        -> OK (and it changes nothing)
//   imbalanced warm start >= ~1.2x    -> abort
//   ubvec raised above the imbalance  -> OK, partition returned UNCHANGED
//   ParMETIS_V3_RefineKway instead    -> OK, partition returned UNCHANGED
//   empty part in the warm start      -> abort
//   isolated vertices / disconnected  -> OK
//
// Every configuration that survives is one where ParMETIS does no rebalancing,
// and rebalancing after a remeshing iteration is the entire reason to call it.
// Making this work needs nparts == npes real MPI ranks, i.e. an MPI-parallel
// benchmark rather than this TBB one with its shared-memory merge.
//
// It would also be optimising almost nothing: measured on bunny00, 3 rounds,
// serial METIS from scratch every round costs 0.34s of a 37.4s run -- 0.9%.
// The bucket remeshing is 86% and the merge 12.6%.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
static std::array<std::size_t, NUM_BUCKETS>
assign_buckets_parmetis(C3t3& c3t3, int round)
{
  std::vector<T3::Cell_handle> cells;
  std::vector<idx_t> xadj, adjncy;
  build_dual_graph(c3t3, cells, xadj, adjncy);

  const idx_t n = static_cast<idx_t>(cells.size());

  // Warm start: the bucket tag each cell already carries. An index outside
  // 1..NUM_BUCKETS means the cell was never bucketed, so fold it into bucket 0
  // and let the repartitioner move it.
  std::vector<idx_t> part(static_cast<std::size_t>(n), 0);
  std::array<std::size_t, NUM_BUCKETS> warm{};
  std::size_t untagged = 0;
  for (idx_t i = 0; i < n; ++i)
  {
    const int sd = static_cast<int>(cells[i]->subdomain_index());
    if (sd >= 1 && sd <= NUM_BUCKETS) part[i] = sd - 1;
    else                              ++untagged;
    ++warm[static_cast<std::size_t>(part[i])];
  }

  // ParMETIS_V3_AdaptiveRepart ABORTS -- MPI_ERR_RANK out of an internal
  // MPI_Recv, not a return code -- if any part of the warm start is EMPTY.
  // Verified against 4.0.3 with a standalone probe: an empty part kills it,
  // while a one-vertex part, a scattered partition, isolated vertices and a
  // disconnected graph are all fine. So an empty part must never reach it.
  //
  // That also covers the first round, where nothing is tagged yet and every
  // cell would land in bucket 0: there is nothing to repair from anyway, and
  // diffusing a partition into existence from a single part is both slow and
  // badly balanced. Fall back to a from-scratch METIS partition in either case.
  const bool empty_part =
    std::find(warm.begin(), warm.end(), std::size_t{0}) != warm.end();

  std::cout << "[parmetis] warm start:";
  for (int b = 0; b < NUM_BUCKETS; ++b) std::cout << "  b" << b << ":" << warm[b];
  if (untagged) std::cout << "  (untagged cells folded into b0: " << untagged << ")";
  std::cout << "\n";

  if (empty_part)
  {
#ifdef CGAL_DD_USE_METIS
    std::cout << "[parmetis] empty part in the warm start"
                 " -- partitioning from scratch with serial METIS\n";
    return assign_buckets_metis(c3t3, round);
#else
    benchmarking::fatal_error("ParMETIS repartition needs a non-empty seed"
                              " partition, but this binary was built without"
                              " serial METIS to produce one.");
    return std::array<std::size_t, NUM_BUCKETS>{};
#endif
  }

  idx_t vtxdist[2] = { 0, n };
  idx_t wgtflag = 0, numflag = 0, ncon = 1, nparts = NUM_BUCKETS, edgecut = 0;
  // vsize is the per-vertex redistribution cost. Every cell costs the same to
  // move here, so a uniform 1 is exact.
  std::vector<idx_t>  vsize(static_cast<std::size_t>(n), 1);
  std::vector<real_t> tpwgts(static_cast<std::size_t>(nparts) * ncon,
                             real_t(1.0) / real_t(nparts));
  std::vector<real_t> ubvec(static_cast<std::size_t>(ncon), real_t(1.05));

  real_t ipc2redist = real_t(1000.0);
  if (const char* e = std::getenv("DD_PARMETIS_IPC2REDIST"))
    ipc2redist = static_cast<real_t>(std::atof(e));

  idx_t options[4];
  options[0] = 1;                      // 1 = honour options[1..3]
  options[1] = 0;                      // no debug output
  options[2] = round;                  // seed, varied per round
  options[3] = PARMETIS_PSR_UNCOUPLED; // nparts != npes; part[] is arbitrary

  MPI_Comm comm = MPI_COMM_WORLD;
  const int ier = ParMETIS_V3_AdaptiveRepart(
      vtxdist, xadj.data(), adjncy.data(),
      NULL /*vwgt*/, vsize.data(), NULL /*adjwgt*/,
      &wgtflag, &numflag, &ncon, &nparts,
      tpwgts.data(), ubvec.data(), &ipc2redist,
      options, &edgecut, part.data(), &comm);
  if (ier != METIS_OK)
    benchmarking::fatal_error("ParMETIS_V3_AdaptiveRepart failed.");

  // How far the interface actually migrated. If this is ~0 the partition is
  // frozen and the seam cells will never be remeshed -- raise ipc2redist.
  std::size_t moved = 0;
  for (idx_t i = 0; i < n; ++i)
    if (part[i] != static_cast<idx_t>(cells[i]->subdomain_index()) - 1) ++moved;

  std::cout << "[parmetis] adaptive repart  edgecut=" << edgecut
            << "  graph nodes=" << n << "  edges=" << adjncy.size() / 2
            << "  ipc2redist=" << ipc2redist
            << "  cells moved=" << moved
            << " (" << (n ? 100.0 * double(moved) / double(n) : 0.0) << "%)\n";
  return commit_partition(cells, part);
}
#endif

// ---------------------------------------------------------------------------
// report_interface_migration: does the frozen interface actually MOVE between
// rounds, or does every partition cut the mesh in the same place?
//
// This is the load-bearing assumption of the whole schedule. Cells on the
// interface are frozen for that round and are only ever remeshed because a
// LATER round puts them in a bucket interior. If successive cuts land in the
// same place -- and they might, since a min-cut is driven by the geometry's
// narrow necks, not by the seed -- then a fixed band of the mesh is never
// remeshed at all and more iterations cannot fix it.
//
// Comparing cell identities across rounds is impossible: the merge rebuilds
// the triangulation from a tet soup, so handles and ordering are new every
// round. Compare the interface's LOCATION instead: voxel-hash the centroids of
// the interface facets on a fixed grid (bbox diagonal / 100) and report the
// Jaccard overlap with the previous round's occupied voxels. 1.0 means the cut
// came back in exactly the same place; 0.0 means it moved entirely.
// ---------------------------------------------------------------------------
static void report_interface_migration(const C3t3& c3t3,
                                       const CGAL::Bbox_3& bbox,
                                       int round)
{
  static std::set<std::array<int, 3>> prev;

  const double diag = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin())
                              + CGAL::square(bbox.ymax() - bbox.ymin())
                              + CGAL::square(bbox.zmax() - bbox.zmin()));
  const double h = (diag > 0) ? diag / 100.0 : 1.0;

  const T3& tr = c3t3.triangulation();
  std::set<std::array<int, 3>> cur;
  std::size_t n_iface = 0;
  for (auto fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)
  {
    const auto f  = *fit;
    const auto mf = tr.mirror_facet(f);
    if (tr.is_infinite(f.first) || tr.is_infinite(mf.first)) continue;
    const int a = static_cast<int>(f.first->subdomain_index());
    const int b = static_cast<int>(mf.first->subdomain_index());
    if (a == b || a == SD_EXTERIOR || b == SD_EXTERIOR) continue;   // not an interface
    ++n_iface;

    double x = 0, y = 0, z = 0;
    for (int i = 0; i < 4; ++i)
    {
      if (i == f.second) continue;
      const auto& p = f.first->vertex(i)->point();
      x += CGAL::to_double(p.x());
      y += CGAL::to_double(p.y());
      z += CGAL::to_double(p.z());
    }
    cur.insert({ int(std::floor((x / 3 - bbox.xmin()) / h)),
                 int(std::floor((y / 3 - bbox.ymin()) / h)),
                 int(std::floor((z / 3 - bbox.zmin()) / h)) });
  }

  std::cout << "[iface " << round << "] frozen facets=" << n_iface
            << "  voxels=" << cur.size();
  if (round > 0 && !prev.empty())
  {
    std::vector<std::array<int, 3>> inter;
    std::set_intersection(prev.begin(), prev.end(), cur.begin(), cur.end(),
                          std::back_inserter(inter));
    const std::size_t uni = prev.size() + cur.size() - inter.size();
    std::cout << "  overlap_with_prev=" << (uni ? double(inter.size()) / double(uni) : 0.0)
              << " (1.0 = cut did not move)";
  }
  std::cout << "\n";
  prev.swap(cur);
}

static std::array<std::size_t, NUM_BUCKETS>
assign_buckets(C3t3& c3t3, const BucketGrid& g)
{
  T3& tr = c3t3.triangulation();
  std::array<std::size_t, NUM_BUCKETS> counts{};
  for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
  {
    // Read the ORIGINAL index before overwriting: without this, hull filler
    // becomes bucket domain and the parallel phase remeshes the convex hull.
    if (cit->subdomain_index() == SD_EXTERIOR) continue;
    const int b = bucket_of_cell(cit, g);
    cit->set_subdomain_index(b + 1);
    ++counts[b];
  }
  return counts;
}


// ---------------------------------------------------------------------------
// collect_selected_tets: collect cells with the given subdomain index from a
// T3 as a tet soup with ref=1 (used to read back a remeshed bucket's selected
// cells for the merge).
// ---------------------------------------------------------------------------
template <typename Tri>
static std::vector<TetSoup> collect_selected_tets(const Tri& tr, int sd)
{
  std::vector<TetSoup> tets;
  for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
  {
    if (cit->subdomain_index() != sd) continue;
    TetSoup t;
    for (int i = 0; i < 4; ++i) t.pts[i] = cit->vertex(i)->point();
    t.ref = 1;
    tets.push_back(std::move(t));
  }
  return tets;
}

// ---------------------------------------------------------------------------
// collect_surface_tris: surface facets of the cells carrying subdomain `sd`,
// with their patch index. A facet is on the surface when its mirror cell is
// infinite or carries a different subdomain index.
// ---------------------------------------------------------------------------
template <typename Tri>
static std::vector<TriSoup> collect_surface_tris(const Tri& tr, int sd)
{
  std::vector<TriSoup> tris;
  for (auto fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)
  {
    const auto f = *fit;
    const auto mf = tr.mirror_facet(f);
    const bool in  = (!tr.is_infinite(f.first)  && f.first->subdomain_index()  == sd);
    const bool min_ = (!tr.is_infinite(mf.first) && mf.first->subdomain_index() == sd);
    if (in == min_) continue;               // not a surface of subdomain sd
    const auto& c = in ? f : mf;
    const auto& other = in ? mf : f;
    if (!tr.is_infinite(other.first)) continue;   // partition interface, not surface
    const int patch = c.first->surface_patch_index(c.second);
    if (patch == 0) continue;               // untagged: not a model surface facet
    TriSoup t;
    int k = 0;
    for (int i = 0; i < 4; ++i)
      if (i != c.second) t.pts[k++] = c.first->vertex(i)->point();
    t.ref = patch;
    tris.push_back(std::move(t));
  }
  return tris;
}

// ---------------------------------------------------------------------------
// compute_t3_bbox: bbox of a triangulation's finite vertices (the lock grid
// needs a Bbox_3; Triangulation_3 has no public bbox() in this CGAL version).
// ---------------------------------------------------------------------------
static CGAL::Bbox_3 compute_t3_bbox(const T3& tr)
{
  double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30, zmin = 1e30, zmax = -1e30;
  bool any = false;
  for (auto vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit)
  {
    any = true;
    const auto p = vit->point();
    xmin = std::min(xmin, CGAL::to_double(p.x()));
    xmax = std::max(xmax, CGAL::to_double(p.x()));
    ymin = std::min(ymin, CGAL::to_double(p.y()));
    ymax = std::max(ymax, CGAL::to_double(p.y()));
    zmin = std::min(zmin, CGAL::to_double(p.z()));
    zmax = std::max(zmax, CGAL::to_double(p.z()));
  }
  if (!any) return CGAL::Bbox_3(0, 0, 0, 1, 1, 1);
  return CGAL::Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
}

// ---------------------------------------------------------------------------
// extract_bucket_copy: build an isolated, VALID T3 for remeshing bucket b by
// deep-copying the whole input triangulation and relabelling cells so that
// only bucket b is selected (subdomain 1) and everything else (straddlers +
// other buckets + exterior hull filler) is frozen (subdomain 0).
//
// We do NOT extract a small standalone bucket T3 (tet-soup roundtrip): the
// all-4-vertices bucket boundary is a jagged internal cut, and capping it
// with infinite cells produces non-manifold edges that read_MEDIT rejects.
// Deep-copying the full mesh keeps the bucket/straddler interface INTERNAL
// (both sides present) so the copy is a valid T3; protect_boundaries then
// freezes that interface. The copies all share the input bbox, so the
// remesher's shared static lock grid covers them and the 8 buckets occupy
// disjoint octants (no false lock contention).
// ---------------------------------------------------------------------------
static T3 extract_bucket_copy(const C3t3& c3t3, int b)
{
  const T3& src = c3t3.triangulation();
  T3 copy(src); // deep copy (points/cells/connectivity preserved exactly)
  const int sd = b + 1;
  for (auto cit = copy.finite_cells_begin(); cit != copy.finite_cells_end(); ++cit)
    cit->set_subdomain_index(cit->subdomain_index() == sd ? 1 : 0);
  return copy;
}

// ---------------------------------------------------------------------------
// extract_bucket_copy_seq: same, but into a Sequential_tag triangulation.
//
// T3 and T3seq are distinct C++ types (Concurrency_tag is a TDS template
// parameter), so `T3seq copy(src)` will not compile. TDS_3::copy_tds is
// templated over the SOURCE tds type and takes vertex/cell converters, which
// gives an exact cross-type connectivity copy with no serialization roundtrip.
//
// The converters need only carry points and subdomain_index: the remesher's
// init_c3t3() calls rescan_after_load_of_triangulation() and re-derives surface
// patches, complex edges and vertex dimensions from the subdomain indices
// (tetrahedral_adaptive_remeshing_impl.h:375-400). The "par" path's deep copy
// preserves more than that, but the extra state is recomputed anyway.
// ---------------------------------------------------------------------------
struct Seq_vertex_converter
{
  // Remeshing_vertex_base_3 does not forward a Point constructor, so
  // default-construct and set_point (vertex dimension/index are re-derived by
  // the remesher's init_c3t3()).
  T3seq::Vertex operator()(const T3::Vertex& v) const
  {
    T3seq::Vertex w;
    w.set_point(v.point());
    return w;
  }
  void operator()(const T3::Vertex&, T3seq::Vertex&) const {}
};
struct Seq_cell_converter
{
  // Carrying subdomain_index alone loses every surface patch index, and with
  // them every 1-D feature edge: init_c3t3 derives features from
  // nb_incident_surface_patches(e) > 1, so a uniform patch yields none. On
  // bunny00 that is the difference between 113,112 feature edges and zero.
  T3seq::Cell operator()(const T3::Cell& c) const
  {
    T3seq::Cell d;
    d.set_subdomain_index(c.subdomain_index());
    for (int i = 0; i < 4; ++i)
      d.set_surface_patch_index(i, c.surface_patch_index(i));
    return d;
  }
  void operator()(const T3::Cell&, T3seq::Cell&) const {}
};

// ---------------------------------------------------------------------------
// extract_bucket_true_seq: build a STANDALONE triangulation for bucket b,
// containing only its own cells plus a one-ring halo -- not a copy of the whole
// mesh.
//
// This is the fix for the real performance ceiling. extract_bucket_copy_seq
// hands every bucket the ENTIRE mesh with a quarter selected, so all four
// threads traverse all N cells (finite_cells / finite_edges / finite_facets)
// while operating on only N/4. Amdahl against the chunked-sequential reference
// put that redundant traversal at ~66% of the work, which caps the whole design
// at ~1.34x no matter how good the partition is -- exactly why METIS's perfect
// balance bought only 8.9%.
//
// Layout of the extracted mesh:
//   core cells (owned by b)  -> subdomain 1, selected, remeshed
//   halo cells (one ring)    -> subdomain 0, out of complex, frozen context
//   everything else          -> absent
//
// The halo is taken by VERTEX adjacency, not facet adjacency, so every core
// vertex keeps its complete star; operations such as collapse walk the full
// 1-ring of a vertex and would otherwise see a truncated star at the border.
// The halo is one cell thick, so it costs O(n^(2/3)) against the core's O(n).
//
// The core/halo interface is a finite subdomain jump (1 vs 0), so it is a
// boundary and protect_boundaries freezes it -- the same contract the buckets
// already relied on. The halo's OUTER boundary faces infinite cells and both
// sides read subdomain 0, so it is not tagged as surface and nothing operates
// there. The real model surface (core against infinite) still is tagged, and
// its patch indices are carried through the Triangles section as before.
// ---------------------------------------------------------------------------
static bool extract_bucket_true_seq(const C3t3& c3t3, int b, T3seq& out,
                                    int halo_rings)
{
  const T3& tr = c3t3.triangulation();
  const int sd = b + 1;

  std::vector<T3::Cell_handle> core;
  std::set<T3::Vertex_handle> core_verts;
  for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
  {
    if (cit->subdomain_index() != sd) continue;
    core.push_back(cit);
    for (int i = 0; i < 4; ++i) core_verts.insert(cit->vertex(i));
  }
  if (core.empty()) return false;

  std::set<T3::Cell_handle> core_set(core.begin(), core.end());
  std::set<T3::Cell_handle> region(core.begin(), core.end());

  // Grow the halo ring by ring. One ring gives every CORE vertex a complete
  // star, which is the minimum for collapse to walk a full 1-ring. It is not
  // necessarily enough for quality: flip evaluates the whole edge star and the
  // collapse/smooth predicates weigh cells beyond the immediate ring. With the
  // whole-mesh copy those cells were present (frozen); with a thin halo they
  // are absent, and the <5deg sliver count rose 12x. Hence the ring count is a
  // knob (DD_HALO_RINGS).
  std::set<T3::Vertex_handle> front = core_verts;
  std::vector<T3::Cell_handle> inc;
  for (int r = 0; r < halo_rings; ++r)
  {
    std::set<T3::Vertex_handle> next;
    for (T3::Vertex_handle v : front)
    {
      inc.clear();
      tr.incident_cells(v, std::back_inserter(inc));
      for (T3::Cell_handle c : inc)
      {
        if (tr.is_infinite(c)) continue;
        if (region.insert(c).second && r + 1 < halo_rings)
          for (int i = 0; i < 4; ++i) next.insert(c->vertex(i));
      }
    }
    if (next.empty()) break;
    front.swap(next);
  }

  // Make the region's boundary manifold before handing it to the reader.
  //
  // A METIS part plus its one-ring can still PINCH: on bunny00's first round
  // 12 of 29,901 boundary edges were used by 4 boundary facets instead of 2,
  // and build_triangulation_from_file rejects that outright. This is the exact
  // failure the original design cited when it chose to copy the whole mesh
  // instead of extracting -- but it is 0.04% of the boundary, not a reason to
  // abandon extraction.
  //
  // Repair by absorbing every cell incident to an offending edge into the
  // region. The added cells are halo (frozen), so this costs a few cells and
  // changes nothing about what gets remeshed. Absorbing can expose new pinches,
  // hence the loop.
  auto vpair = [](T3::Vertex_handle a, T3::Vertex_handle b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
  };
  for (int pass = 0; pass < 12; ++pass)
  {
    std::map<std::pair<T3::Vertex_handle, T3::Vertex_handle>, int> ecount;
    for (T3::Cell_handle c : region)
    {
      for (int i = 0; i < 4; ++i)
      {
        const T3::Cell_handle nb = c->neighbor(i);
        if (!tr.is_infinite(nb) && region.find(nb) != region.end()) continue;
        T3::Vertex_handle vs[3];
        int k = 0;
        for (int j = 0; j < 4; ++j) if (j != i) vs[k++] = c->vertex(j);
        ++ecount[vpair(vs[0], vs[1])];
        ++ecount[vpair(vs[0], vs[2])];
        ++ecount[vpair(vs[1], vs[2])];
      }
    }
    std::vector<std::pair<T3::Vertex_handle, T3::Vertex_handle>> bad;
    for (const auto& kv : ecount) if (kv.second != 2) bad.push_back(kv.first);
    if (bad.empty()) break;

    for (const auto& e : bad)
    {
      T3::Cell_handle c; int i, j;
      if (!tr.is_edge(e.first, e.second, c, i, j)) continue;
      auto circ = tr.incident_cells(T3::Edge(c, i, j));
      const auto end = circ;
      do { if (!tr.is_infinite(circ)) region.insert(circ); } while (++circ != end);
    }
  }

  std::set<T3::Cell_handle> halo;
  for (T3::Cell_handle c : region)
    if (core_set.find(c) == core_set.end()) halo.insert(c);

  std::vector<TetSoup> tets;
  tets.reserve(core.size() + halo.size());
  for (T3::Cell_handle c : core)
  {
    TetSoup t;
    for (int i = 0; i < 4; ++i) t.pts[i] = c->vertex(i)->point();
    t.ref = 1;
    tets.push_back(std::move(t));
  }
  for (T3::Cell_handle c : halo)
  {
    TetSoup t;
    for (int i = 0; i < 4; ++i) t.pts[i] = c->vertex(i)->point();
    t.ref = 0;                       // frozen context, out of the complex
    tets.push_back(std::move(t));
  }

  // Carry the model-surface patches of the WHOLE REGION -- core AND halo.
  //
  // Emitting only the core's surface leaves the halo's surface facets untagged
  // (halo reads subdomain 0 and so does the infinite cell beyond it, so the
  // reader never marks them), which makes the model surface DISCONTINUOUS
  // exactly where it crosses the core/halo interface. That is where the extra
  // slivers appeared: <5deg went 0.005% -> 0.062% against the whole-mesh copy,
  // and widening the halo changed nothing (halo_rings 1 and 2 gave byte-identical
  // meshes), which ruled out "not enough context" and pointed here instead.
  std::vector<TriSoup> tris;
  for (T3::Cell_handle c : region)
  {
    for (int i = 0; i < 4; ++i)
    {
      const int patch = c->surface_patch_index(i);
      if (patch == 0) continue;
      const T3::Cell_handle nb = c->neighbor(i);
      if (!tr.is_infinite(nb) && region.find(nb) != region.end() && nb < c)
        continue;                    // emit each interior patch facet once
      TriSoup t;
      int k = 0;
      for (int j = 0; j < 4; ++j)
        if (j != i) t.pts[k++] = c->vertex(j)->point();
      t.ref = patch;
      tris.push_back(std::move(t));
    }
  }

  std::stringstream ss;
  write_tet_soup_medit(ss, tets, tris);
  out.clear();
  const bool ok = CGAL::IO::read_MEDIT(ss, out);
  std::cout << "[extract b" << b << "] core=" << core.size()
            << " halo=" << halo.size() << " tris=" << tris.size()
            << " read_MEDIT=" << (ok ? "ok" : "FAILED")
            << " cells_out=" << out.number_of_finite_cells() << std::endl;
  if (!ok || out.number_of_finite_cells() == 0)
    return false;
  return true;
}

static void extract_bucket_copy_seq(const C3t3& c3t3, int b, T3seq& out)
{
  const T3& src = c3t3.triangulation();
  out.clear();
  Seq_vertex_converter vc;
  Seq_cell_converter   cc;
  out.set_infinite_vertex(
    out.tds().copy_tds(src.tds(), src.infinite_vertex(), vc, cc));

  const int sd = b + 1;
  for (auto cit = out.finite_cells_begin(); cit != out.finite_cells_end(); ++cit)
    cit->set_subdomain_index(cit->subdomain_index() == sd ? 1 : 0);
}


// ---------------------------------------------------------------------------
// rebuild_complex_from_subdomains: re-derive the 2-D/1-D complex from the
// cells' subdomain indices.
//
// mark_seam_band/restore_bulk move cells in and out of the complex by writing
// subdomain_index() directly, which bypasses add_to_complex/remove_from_complex.
// cells_in_complex() filters live on subdomain_index so cell metrics stay
// correct, but facets_in_complex()/edges_in_complex() filter on CACHED state
// (the per-cell surface_patch_index bits and the edges_ map). After the seam
// pass those caches describe the band/bulk interface -- the labelling that was
// in force while the remesher ran -- not the model surface, so every facet and
// feature-edge metric is measured on the wrong set.
//
// rescan_after_load_of_triangulation() alone does NOT fix this: it only
// recounts using is_in_complex(), which reads the same stale bits. The flags
// have to be cleared and re-derived, which is what this does. Mirrors the
// tagging in Adaptive_remesher::init_c3t3 (tetrahedral_adaptive_remeshing_impl.h
// :387-445): a facet is on the surface iff its two incident cells disagree on
// subdomain index (an infinite mirror cell reads as the default index, so the
// domain's outer boundary is included), and an edge is a feature iff it is
// incident to >2 subdomains, >1 surface patch, or >2 complex facets.
// ---------------------------------------------------------------------------
static void rebuild_complex_from_subdomains(C3t3& c3t3)
{
  using namespace CGAL::Tetrahedral_remeshing;
  T3& tr = c3t3.triangulation();

  // The per-facet surface bits and patch indices live on the CELLS, which come
  // straight from read_MEDIT (the Triangles section restores them), so they are
  // already correct and must not be touched. In particular do NOT call
  // remove_from_complex(f) to "clear" them first: it resets the patch index to
  // the default, which erases exactly the information the feature edges are
  // derived from.
  //
  // What does need rebuilding after swapping a triangulation into the C3t3 is
  // the 1-D complex (edges_ is keyed on the OLD triangulation's handles) and
  // the cached counts.
  std::vector<T3::Edge> stale_edges(c3t3.edges_in_complex_begin(),
                                    c3t3.edges_in_complex_end());
  for (const T3::Edge& e : stale_edges) c3t3.remove_from_complex(e);

  for (const T3::Edge& e : tr.finite_edges())
  {
    if (nb_incident_subdomains(e, c3t3) > 2
     || nb_incident_surface_patches(e, c3t3) > 1
     || nb_incident_complex_facets(e, c3t3) > 2)
      c3t3.add_to_complex(e, 1);
  }

  c3t3.rescan_after_load_of_triangulation();
}


// ---------------------------------------------------------------------------
// export_mesh: dump a triangulation to a Medit .mesh file for visual
// inspection. Templated so it accepts both the Parallel_tag T3 used by the
// "par" bucket path and the Sequential_tag T3seq used by the "seq" path.
// ---------------------------------------------------------------------------
template<typename Tr>
static void export_mesh(const Tr& tr, const std::string& dir, const std::string& name)
{
  std::filesystem::create_directories(dir);
  std::ofstream os(dir + "/" + name + ".mesh");
  CGAL::IO::write_MEDIT(os, tr);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
#ifdef CGAL_DD_USE_PARMETIS
// ParMETIS is an MPI library even when it is driven, as here, with a single
// rank: the whole graph is local and the parallelism is TBB. MPI_Init is still
// mandatory before any ParMETIS_* call. The guard's destructor covers every
// normal return path out of main.
struct Mpi_session
{
  Mpi_session(int& argc, char**& argv) { MPI_Init(&argc, &argv); }
  ~Mpi_session() { MPI_Finalize(); }
};
#endif

int main(int argc, char** argv)
{
#ifdef CGAL_DD_USE_PARMETIS
  Mpi_session mpi_session(argc, argv);
#endif
  using namespace benchmarking;
  using nlohmann::json;
  json results_json;
  std::cout << std::setprecision(17);
  std::cerr << std::setprecision(17);

  if (argc < 7 || argc > 10)
  {
    fatal_error(std::string("Usage: ") + argv[0] +
                " <input_mesh> <num_iterations> <remeshing_target_edge_factor>"
                " <smooth_constrained_edges> <num_threads> <results_json_path>"
                " [bucket_exec: par|seq  (default seq)]"
                " [final_pass: none  (the global passes have been removed)]"
                " [final_pass_iterations (ignored)]");
  }
  const std::string input              = argv[1];
  const int         num_iterations     = std::stoi(argv[2]);
  const double      target_edge_factor = std::stod(argv[3]);
  const bool        smooth_constrained = std::stoi(argv[4]) != 0;
  const int         num_threads        = std::stoi(argv[5]);
  const std::string results_json_path  = argv[6];
  const std::string bucket_exec        = (argc >= 8) ? argv[7] : "seq";
  // The final whole-mesh passes have been REMOVED (see the end of main). The
  // argument is still accepted so existing driver scripts keep working, but
  // only 'none' is valid.
  const std::string final_pass         = (argc >= 9) ? argv[8] : "none";
  if (bucket_exec != "par" && bucket_exec != "seq")
    fatal_error("bucket_exec must be 'par' or 'seq'.");
  if (final_pass != "none")
    fatal_error("final_pass: the global passes have been removed from this"
                " benchmark; only 'none' is accepted.");

  std::filesystem::create_directories(
    std::filesystem::path(std::filesystem::absolute(results_json_path)).parent_path());

  if (num_threads <= 0) fatal_error("num_threads must be positive.");

  // DD_EXPORT_DIR overrides the default snapshot directory, and DD_EXPORT_TAG
  // prefixes every dumped mesh, so a run's mesh series is self-describing on
  // disk. DD_EXPORT_SUBMESHES additionally dumps every bucket's standalone
  // submesh (one file per partition, per iteration) right after extraction --
  // off by default since it is one file per bucket per iteration.
  std::string mesh_export_dir = "dd_mesh_exports";
  if (const char* env = std::getenv("DD_EXPORT_DIR")) mesh_export_dir = env;
  std::string export_tag = "dd";
  if (const char* env = std::getenv("DD_EXPORT_TAG")) export_tag = env;
  const bool export_submeshes = std::getenv("DD_EXPORT_SUBMESHES") != nullptr;

  #if TBB_VERSION_MAJOR >= 2018
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);
  #else
    tbb::task_scheduler_init sched(num_threads);
  #endif

  std::cout << "Domain-decomp (ParMmg-style: centroid partition + interface"
            << " migration, " << NUM_BUCKETS << " buckets) (threads: "
            << num_threads << ")\n";
  append_run_info(results_json, "Technique", "DomainDecomp-ExtractRemeshMerge");
  append_run_info(results_json, "Num_threads", num_threads);

  // --- Load mesh -----------------------------------------------------------
  C3t3 c3t3;
  T3& tr = c3t3.triangulation();
  {
    std::ifstream is(input, std::ios_base::in);
    if (!CGAL::IO::read_MEDIT(is, tr))
      fatal_error(std::string("Could not read input mesh '") + input + "'");
  }
  std::cout << "Number of vertices: " << tr.number_of_vertices() << "\n";
  std::cout << "Number of cells: "    << tr.number_of_cells()    << "\n";

  const std::string input_name = std::filesystem::path(input).stem().string();
  write_triangulation_info(results_json, tr, input_name);

  const double avg_edge_length = compute_average_edge_length(tr);
  if (avg_edge_length <= 0.0) fatal_error("Could not compute average edge length.");
  const double target_edge_length = avg_edge_length * target_edge_factor;
  append_run_info(results_json, "Edge Length", target_edge_length);

  // --- Strip the finite hull filler ---------------------------------------
  // is_protected_interface() distinguishes the partition interface from the
  // model surface by "the mirror cell is infinite". That only holds once the
  // convex-hull filler is gone: while it is present the model surface is a
  // finite/finite subdomain jump, indistinguishable from an interface, and
  // would be frozen exactly as before. One tet-soup roundtrip establishes the
  // invariant for iteration 0; every later iteration inherits it from the merge.
  {
    std::vector<TetSoup> dom;
    dom.reserve(static_cast<std::size_t>(tr.number_of_finite_cells()));
    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
    {
      if (cit->subdomain_index() == SD_EXTERIOR) continue;
      TetSoup t;
      for (int i = 0; i < 4; ++i) t.pts[i] = cit->vertex(i)->point();
      t.ref = 1;
      dom.push_back(std::move(t));
    }
    std::vector<TriSoup> dom_tris;
    for (auto fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)
    {
      const auto f = *fit;
      const int patch = f.first->surface_patch_index(f.second);
      if (patch == 0) continue;
      TriSoup t;
      int k = 0;
      for (int i = 0; i < 4; ++i)
        if (i != f.second) t.pts[k++] = f.first->vertex(i)->point();
      t.ref = patch;
      dom_tris.push_back(std::move(t));
    }
    T3 stripped;
    std::stringstream ss;
    write_tet_soup_medit(ss, dom, dom_tris);
    if (!read_t3_from_medit(ss, stripped) || stripped.number_of_finite_cells() == 0)
      fatal_error("Hull-filler strip: read_MEDIT failed or empty.");
    std::cout << "[strip] hull filler removed: " << dom.size()
              << " domain cells kept, " << stripped.number_of_vertices()
              << " vertices, " << dom_tris.size() << " surface triangles carried\n";
    c3t3.triangulation().swap(stripped);
  }

  // Dump the stripped mesh: this, not the raw input, is what the adaptation
  // loop actually starts from, so it is the only fair input for a sequential
  // baseline. Feeding the baseline the raw input instead would have it remesh
  // the convex-hull filler as well, and give it a different complex.
  export_mesh(c3t3.triangulation(), mesh_export_dir, export_tag + "_stripped");

  // --- ParMmg-style adaptation loop ---------------------------------------
  // Each outer iteration: re-partition (rotated / shifted cut planes) -> remesh
  // every bucket concurrently for ONE iteration with its interface frozen ->
  // merge. Because the partition moves, cells frozen on an interface in this
  // iteration sit in a bucket INTERIOR in the next one, so no region is
  // systematically starved. That is the point of the restructuring: the
  // previous design froze the same band every time, which left it stranded in
  // the split-dominant regime the bulk had already left.
  //
  // NOTE (unfixed): remesh_boundaries(false) freezes the whole complex
  // boundary, i.e. the bucket interface AND the model surface. Migration moves
  // the interface but never the model surface, so the surface is still never
  // remeshed.
  const CGAL::Bbox_3 mesh_bbox = compute_t3_bbox(tr);
  double t_extract_total = 0, t_bucket_total = 0, t_merge_total = 0;
  // The partition alone, carved out of t_extract: "segmentation" and "build the
  // per-bucket meshes" are different costs and only the first is intrinsic to
  // the decomposition -- the copy is an artefact of how buckets are extracted.
  double t_partition_total = 0;

  // DD_PARTITIONER selects how the mesh is split each round:
  //   metis     serial METIS, partitioned from scratch every round (default)
  //   grid      the original geometric rotated/shifted grid split
  //   parmetis  ParMETIS adaptive REpartition, warm-started from the previous
  //             round's partition. BROKEN AT ONE MPI RANK -- opt-in only, see
  //             assign_buckets_parmetis. Kept so the finding is reproducible.
  enum class Partitioner { Grid, Metis, ParMetis };
#if defined(CGAL_DD_USE_METIS)
  Partitioner partitioner = Partitioner::Metis;
#else
  Partitioner partitioner = Partitioner::Grid;
#endif
  if (const char* env = std::getenv("DD_PARTITIONER"))
  {
    const std::string v = env;
    if      (v == "grid")     partitioner = Partitioner::Grid;
    else if (v == "metis")    partitioner = Partitioner::Metis;
    else if (v == "parmetis")
    {
      partitioner = Partitioner::ParMetis;
      std::cout << "[partitioner] WARNING: ParMETIS aborts (MPI_ERR_RANK) as"
                   " soon as it has real rebalancing to do at one MPI rank."
                   " See assign_buckets_parmetis.\n";
    }
    else fatal_error("DD_PARTITIONER must be 'grid', 'metis' or 'parmetis'.");
  }
#ifndef CGAL_DD_USE_METIS
  if (partitioner == Partitioner::Metis)
    fatal_error("DD_PARTITIONER=metis but this binary was built without METIS.");
#endif
#ifndef CGAL_DD_USE_PARMETIS
  if (partitioner == Partitioner::ParMetis)
    fatal_error("DD_PARTITIONER=parmetis but this binary was built without ParMETIS.");
#endif
  std::cout << "[partitioner] "
            << (partitioner == Partitioner::ParMetis
                  ? "parmetis (adaptive repartition, warm-started)"
                  : partitioner == Partitioner::Metis
                  ? "metis (dual graph, from scratch)"
                  : "geometric grid") << "\n";

  // DD_EXTRACT=true builds standalone bucket meshes (core + one-ring halo);
  // anything else keeps the whole-mesh copy.
  bool true_extract = false;
  if (const char* env = std::getenv("DD_EXTRACT"))
    true_extract = (std::string(env) == "true");
  int halo_rings = 2;
  if (const char* env = std::getenv("DD_HALO_RINGS")) halo_rings = std::max(1, std::atoi(env));
  std::cout << "[extract] " << (true_extract ? "true (core + halo)" : "whole-mesh copy")
            << (true_extract ? "  halo_rings=" + std::to_string(halo_rings) : std::string())
            << "\n";

  // DD_BUCKET_ITERS is how many remeshing iterations each bucket runs per outer
  // pass. The original schedule was 1 iteration inside a num_iterations-long
  // outer loop that re-partitions and re-merges every pass (ParMmg-style
  // interface migration). Setting it to N with num_iterations=1 gives the
  // "decompose once, remesh each region N times, merge once" schedule instead.
  int bucket_iters = 1;
  if (const char* env = std::getenv("DD_BUCKET_ITERS"))
    bucket_iters = std::max(1, std::atoi(env));
  std::cout << "[schedule] outer_passes=" << num_iterations
            << "  bucket_iters_per_pass=" << bucket_iters << "\n";

  for (int it = 0; it < num_iterations; ++it)
  {
    const BucketGrid grid = grid_for_iteration(mesh_bbox, it);

    CGAL::Real_timer t_extract; t_extract.start();
    CGAL::Real_timer t_partition; t_partition.start();
    std::array<std::size_t, NUM_BUCKETS> counts{};
    switch (partitioner)
    {
#ifdef CGAL_DD_USE_PARMETIS
      case Partitioner::ParMetis: counts = assign_buckets_parmetis(c3t3, it); break;
#endif
#ifdef CGAL_DD_USE_METIS
      case Partitioner::Metis:    counts = assign_buckets_metis(c3t3, it);    break;
#endif
      default:                    counts = assign_buckets(c3t3, grid);        break;
    }
    t_partition.stop();
    t_partition_total += t_partition.time();
    const std::size_t owned =
      std::accumulate(counts.begin(), counts.end(), std::size_t{0});
    const std::size_t cmin = *std::min_element(counts.begin(), counts.end());
    const std::size_t cmax = *std::max_element(counts.begin(), counts.end());
    // cmax/cmin reported 0 whenever a bucket was EMPTY, which is exactly the
    // pathological case. What actually bounds the bucket phase is the largest
    // bucket against a perfectly balanced share, so report that.
    const double balance = (owned > 0)
      ? double(cmax) / (double(owned) / NUM_BUCKETS)
      : 1.0;

    if (partitioner == Partitioner::ParMetis)
      std::cout << "[iter " << it << "] parmetis";
    else if (partitioner == Partitioner::Metis)
      std::cout << "[iter " << it << "] metis";
    else
      std::cout << "[iter " << it << "] grid " << grid.div[0] << "x" << grid.div[1]
                << "x" << grid.div[2] << " shift=" << grid.shift[0];
    std::cout
              << "  owned=" << owned
              << "  max/mean=" << balance
              << "  max/min=" << (cmin ? double(cmax) / double(cmin) : -1.0);
    for (int b = 0; b < NUM_BUCKETS; ++b)
      std::cout << "  b" << b << ":" << counts[b];
    std::cout << "\n";

    report_interface_migration(c3t3, mesh_bbox, it);

    std::vector<T3>    bucket_t3s;
    std::vector<T3seq> bucket_t3s_seq;
    if (bucket_exec == "par") bucket_t3s.resize(NUM_BUCKETS);
    else                      bucket_t3s_seq.resize(NUM_BUCKETS);

    // Extract the buckets concurrently. The four extractions are independent:
    // each only READS the shared source triangulation and writes its own output
    // mesh. `copy_tds` takes the source by const reference, and `incident_cells`
    // tracks visited vertices in a FUNCTION-LOCAL set (the vertex bases have no
    // `visited` member, so the marking specialisation of Visited_vertex is not
    // selected) -- so no shared state is mutated.
    //
    // VTune showed why this matters: the main thread carried more CPU (328s)
    // than any worker, because extraction and merge ran single-threaded while
    // three workers idled.
    {
      std::atomic<bool> extract_ok{true};
      tbb::task_group tg_ex;
      for (int b = 0; b < NUM_BUCKETS; ++b)
      {
        tg_ex.run([&, b]() {
          if (bucket_exec == "par") bucket_t3s[b] = extract_bucket_copy(c3t3, b);
          else if (true_extract)
          {
            if (!extract_bucket_true_seq(c3t3, b, bucket_t3s_seq[b], halo_rings))
              extract_ok = false;
          }
          else                      extract_bucket_copy_seq(c3t3, b, bucket_t3s_seq[b]);
        });
      }
      tg_ex.wait();
      if (!extract_ok)
        fatal_error("True extraction failed for a bucket (empty or non-manifold).");
    }
    t_extract.stop();
    t_extract_total += t_extract.time();
    if (bucket_exec == "seq")
    {
      std::size_t tot = 0;
      for (int b = 0; b < NUM_BUCKETS; ++b) tot += bucket_t3s_seq[b].number_of_finite_cells();
      std::cout << "[iter " << it << "] extracted cells total=" << tot
                << " (mesh=" << c3t3.triangulation().number_of_finite_cells()
                << ", ratio=" << double(tot) / double(c3t3.triangulation().number_of_finite_cells())
                << ")\n";
    }

    // Each bucket is a standalone, independently valid triangulation (see
    // extract_bucket_true_seq / extract_bucket_copy): dump it as its own mesh
    // file rather than only the global merged mesh, so the partition itself
    // can be inspected and each submesh loaded on its own.
    if (export_submeshes)
    {
      for (int b = 0; b < NUM_BUCKETS; ++b)
      {
        const std::string name = export_tag + "_iter" + std::to_string(it)
                                + "_bucket" + std::to_string(b);
        if (bucket_exec == "par") export_mesh(bucket_t3s[b], mesh_export_dir, name);
        else                      export_mesh(bucket_t3s_seq[b], mesh_export_dir, name);
      }
    }

    // Each "par" copy needs its OWN lock grid: the copies share the input bbox,
    // so one shared grid would route every bucket's operations onto the same
    // grid cells. The "seq" path needs no grid at all -- a Sequential_tag
    // triangulation selects the sequential executor, so lock_zone() is never
    // called and parallelism comes purely from the task_group.
    std::vector<std::unique_ptr<LockDS>> bucket_locks;
    if (bucket_exec == "par")
    {
      bucket_locks.resize(NUM_BUCKETS);
      for (int b = 0; b < NUM_BUCKETS; ++b)
      {
        bucket_locks[b] = std::make_unique<LockDS>(mesh_bbox, 8);
        bucket_t3s[b].set_lock_data_structure(bucket_locks[b].get());
      }
    }

    CGAL::Real_timer t_bucket; t_bucket.start();
    {
      tbb::task_group tg;
      for (int b = 0; b < NUM_BUCKETS; ++b)
      {
        if (bucket_exec == "par")
          tg.run([&bucket_t3s, b, target_edge_length, smooth_constrained, bucket_iters]() {
            CGAL::tetrahedral_isotropic_remeshing(
              bucket_t3s[b], target_edge_length,
              CGAL::parameters::number_of_iterations(bucket_iters)
                               .remesh_boundaries(false)
                               .smooth_constrained_edges(smooth_constrained));
          });
        else
          tg.run([&bucket_t3s_seq, b, target_edge_length, smooth_constrained, bucket_iters]() {
            try {
              CGAL::tetrahedral_isotropic_remeshing(
                bucket_t3s_seq[b], target_edge_length,
                CGAL::parameters::number_of_iterations(bucket_iters)
                                 .remesh_boundaries(false)
                                 .smooth_constrained_edges(smooth_constrained));
            } catch (const std::exception& ex) {
              std::cerr << "[bucket " << b << "] EXCEPTION: " << ex.what()
                        << std::endl;
              throw;
            } catch (...) {
              std::cerr << "[bucket " << b << "] UNKNOWN EXCEPTION" << std::endl;
              throw;
            }
          });
      }
      tg.wait();
    }
    t_bucket.stop();
    t_bucket_total += t_bucket.time();

    // Merge: every domain cell is owned by exactly one bucket, so the union of
    // the buckets' selected cells is the whole domain -- there is no separate
    // straddler set to splice back in.
    CGAL::Real_timer t_merge; t_merge.start();
    std::vector<TetSoup> merge_tets;
    std::vector<TriSoup> merge_tris;
    merge_tets.reserve(static_cast<std::size_t>(
      c3t3.triangulation().number_of_finite_cells()) + 64);
    for (int b = 0; b < NUM_BUCKETS; ++b)
    {
      auto bt = (bucket_exec == "par") ? collect_selected_tets(bucket_t3s[b], 1)
                                       : collect_selected_tets(bucket_t3s_seq[b], 1);
      // Carry the bucket tag into the merged mesh (ref = b+1, not a uniform 1)
      // so the NEXT round's ParMETIS repartition can warm-start from THIS
      // round's partition. Every domain cell still gets a non-zero subdomain
      // index, so nothing that only separates domain from SD_EXTERIOR changes.
      for (TetSoup& t : bt) t.ref = b + 1;
      merge_tets.insert(merge_tets.end(),
                        std::make_move_iterator(bt.begin()),
                        std::make_move_iterator(bt.end()));
      auto bs = (bucket_exec == "par") ? collect_surface_tris(bucket_t3s[b], 1)
                                       : collect_surface_tris(bucket_t3s_seq[b], 1);
      merge_tris.insert(merge_tris.end(),
                        std::make_move_iterator(bs.begin()),
                        std::make_move_iterator(bs.end()));
    }
    bucket_t3s.clear();
    bucket_t3s_seq.clear();
    bucket_locks.clear();

    T3 merged_t3;
    {
      std::stringstream ss;
      write_tet_soup_medit(ss, merge_tets, merge_tris);
      if (!read_t3_from_medit(ss, merged_t3) || merged_t3.number_of_finite_cells() == 0)
        fatal_error("Merge: read_MEDIT failed or empty.");
    }
    {
      std::set<int> patches;
      std::size_t tagged = 0;
      for (auto fit = merged_t3.finite_facets_begin();
           fit != merged_t3.finite_facets_end(); ++fit)
      {
        const int pi = fit->first->surface_patch_index(fit->second);
        if (pi != 0) { patches.insert(pi); ++tagged; }
      }
      std::cout << "[iter " << it << "] merged mesh: " << tagged
                << " tagged facets, " << patches.size() << " distinct patches\n";
    }
    c3t3.triangulation().swap(merged_t3);
    t_merge.stop();
    t_merge_total += t_merge.time();

    std::cout << "[iter " << it << "] surface tris carried: " << merge_tris.size() << "\n";
    std::cout << "[iter " << it << "] extract " << t_extract.time()
              << "s  bucket " << t_bucket.time()
              << "s  merge " << t_merge.time()
              << "s  -> " << c3t3.triangulation().number_of_finite_cells()
              << " cells, " << c3t3.triangulation().number_of_vertices()
              << " vertices\n";
  }

  // The merge tags every cell with its bucket (for the ParMETIS warm start), so
  // the last round's bucket interfaces would otherwise read as real subdomain
  // boundaries and be counted as surface. Collapse the buckets back into one
  // domain before re-deriving the complex.
  {
    T3& last_tr = c3t3.triangulation();
    for (auto cit = last_tr.finite_cells_begin();
         cit != last_tr.finite_cells_end(); ++cit)
      if (cit->subdomain_index() != SD_EXTERIOR) cit->set_subdomain_index(1);
  }

  // The swap leaves the C3t3's cached facet/edge complex describing the
  // pre-merge mesh, so re-derive it before measuring quality.
  rebuild_complex_from_subdomains(c3t3);

  // The merged mesh BEFORE any global pass. This is what the decomposition on
  // its own produces, seams and all, and it is the baseline the global passes
  // are judged against. Dumped OUTSIDE t_global so the mesh I/O never enters
  // the reported time.
  export_mesh(c3t3.triangulation(), mesh_export_dir, export_tag + "_merged_global00");

#ifdef _WIN32
  FILETIME _ft_create, _ft_exit, _cpu_kernel_before, _cpu_user_before;
  GetProcessTimes(GetCurrentProcess(), &_ft_create, &_ft_exit,
                  &_cpu_kernel_before, &_cpu_user_before);
#endif
  // --- Final global passes: REMOVED ---------------------------------------
  // The whole-mesh sequential passes that used to run here were the seam
  // recovery. They are gone by design: ParMmg does not run one either, and now
  // that the partition is REpartitioned every round the interface migrates, so
  // a cell frozen on a seam in one round sits in a bucket interior in the next
  // and gets remeshed there. Keeping a global pass would also mask exactly the
  // effect this benchmark exists to measure. t_global therefore stays at ~0 and
  // Final_Pass_Time is reported as 0 so the results schema is unchanged.
  CGAL::Real_timer t_global; t_global.start();
  t_global.stop();
#ifdef _WIN32
  FILETIME _cpu_kernel_after, _cpu_user_after;
  GetProcessTimes(GetCurrentProcess(), &_ft_create, &_ft_exit,
                  &_cpu_kernel_after, &_cpu_user_after);
  auto _ft_to_s = [](const FILETIME& ft) -> double {
    return static_cast<double>((static_cast<ULONGLONG>(ft.dwHighDateTime) << 32)
                               | ft.dwLowDateTime) / 1e7;
  };
  double cpu_time = (_ft_to_s(_cpu_user_after) - _ft_to_s(_cpu_user_before))
                  + (_ft_to_s(_cpu_kernel_after) - _ft_to_s(_cpu_kernel_before));
#else
  double cpu_time = 0.0;
#endif
  std::cout << "[final pass] (" << final_pass << ") " << t_global.time() << "s\n";
  {
    const T3& ftr = c3t3.triangulation();
    std::size_t n_dom = 0;
    for (auto cit = ftr.finite_cells_begin(); cit != ftr.finite_cells_end(); ++cit)
      if (cit->subdomain_index() != SD_EXTERIOR) ++n_dom;
    std::cout << "[final] " << ftr.number_of_finite_cells() << " finite cells, "
              << n_dom << " domain cells, "
              << ftr.number_of_vertices() << " vertices\n";
    std::cout << "[complex] cells=" << c3t3.number_of_cells_in_complex()
              << " facets=" << c3t3.number_of_facets_in_complex()
              << " edges=" << c3t3.number_of_edges_in_complex() << "\n";
    append_metric_result(results_json, "Performance", "Final_Domain_Cells", "Value",
                         static_cast<double>(n_dom));
  }

  const double total_time = t_extract_total + t_bucket_total + t_merge_total + t_global.time();

  // --- Results -------------------------------------------------------------
  append_metric_result(results_json, "Performance", "Total_Time",        "Value", total_time);
  append_metric_result(results_json, "Performance", "Extract_Time",      "Value", t_extract_total);
  // Partition_Time is INSIDE Extract_Time; Extraction_Only_Time is the rest of it
  // (building each bucket's mesh). Kept separate so a runtime can be quoted as
  // partition + bucket remeshing + global passes, with the copy and the merge --
  // both pure decomposition overhead -- excluded.
  append_metric_result(results_json, "Performance", "Partition_Time",    "Value", t_partition_total);
  append_metric_result(results_json, "Performance", "Extraction_Only_Time", "Value",
                       t_extract_total - t_partition_total);
  append_metric_result(results_json, "Performance", "Remesh_Only_Time",  "Value",
                       t_partition_total + t_bucket_total + t_global.time());
  append_metric_result(results_json, "Performance", "Bucket_Phase_Time", "Value", t_bucket_total);
  append_metric_result(results_json, "Performance", "Merge_Time",        "Value", t_merge_total);
  append_metric_result(results_json, "Performance", "Final_Pass_Time",   "Value", t_global.time());
  append_metric_result(results_json, "Performance", "CPU_Time",          "Value", cpu_time);

  generate_quality_metrics(c3t3, results_json);
  append_metric_result(results_json, "Performance", "Memory", "Value",
                       CGAL::Memory_sizer().virtual_size() >> 20);
  append_execution_status(results_json, "success");
  write_results_json(results_json, results_json_path);
  return 0;
}
