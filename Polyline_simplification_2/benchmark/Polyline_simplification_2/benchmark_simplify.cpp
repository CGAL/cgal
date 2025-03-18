#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/IO/WKT.h>

#include <benchmark/benchmark.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>

namespace PS = CGAL::Polyline_simplification_2;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polygon_2 = CGAL::Polygon_2<K>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

using Vb = PS::Vertex_base_2<K>;
using Fb = CGAL::Constrained_triangulation_face_base_2<K>;
using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, CGAL::Exact_predicates_tag>;
using CT = CGAL::Constrained_triangulation_plus_2<CDT>;
using Stop = PS::Stop_below_count_ratio_threshold;
using Cost = PS::Squared_distance_cost;

static void BM_Simplify(benchmark::State& state, std::string filename) {
  using Point_2 = K::Point_2;
  using MultiPoint = std::vector<Point_2>;

  using LineString = std::vector<Point_2>;
  using MultiLineString = std::deque<LineString>;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
  using MultiPolygon = std::deque<Polygon_with_holes_2>;

  std::ifstream ifs(filename);
  MultiPoint points;
  MultiLineString polylines;
  MultiPolygon polygons;
  if(!CGAL::IO::read_WKT(ifs, points, polylines, polygons) && false) {
     state.SkipWithError("Cannot read file " + filename);
    return;
  }

  state.counters["#points"] = points.size();
  state.counters["#polylines"] = polylines.size();
  state.counters["#polygons"] = polygons.size();

  CT ct;
  for(const auto& point : points) {
    ct.insert(point);
  }
  for(const auto& polyline : polylines) {
    ct.insert_constraint(polyline);
  }
  for(const auto& polygon_with_holes : polygons) {
    const Polygon_2& outer_polygon = polygon_with_holes.outer_boundary();
    ct.insert_constraint(outer_polygon);
    for(Polygon_with_holes_2::Hole_const_iterator it = polygon_with_holes.holes_begin();
        it != polygon_with_holes.holes_end(); ++it)
    {
      const Polygon_2& hole = *it;
      ct.insert_constraint(hole);
    }
  }

  state.counters["nb of constraints"] = ct.number_of_constraints();
  state.counters["nb of vertices"] = ct.number_of_vertices();
  state.counters["nb of sub-constraints"] = ct.number_of_subconstraints();

  for([[maybe_unused]] auto _ : state) {
    state.PauseTiming();
    CT ct_copy = ct; // Copy the object `ct` in the loop
    state.ResumeTiming();
    PS::simplify(ct_copy, Cost(), Stop(0.5));
  }
}

int main(int argc, char** argv) {
  std::string filename = CGAL::data_file_path("wkt/norway-MP.wkt");
  if(argc > 1) {
    std::string_view arg1{argv[1]};
    if(arg1.size() < 2 || arg1[0] != '-' || arg1[1] != '-') {
      --argc;
      ++argv;
      filename = arg1;
    }
  }
  benchmark::RegisterBenchmark("simplify file " + filename, BM_Simplify, filename);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
}
