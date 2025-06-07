#ifndef CGAL_MESH_3_BENCHMARK_MESH_3_MESH_QUALITY_H
#define CGAL_MESH_3_BENCHMARK_MESH_3_MESH_QUALITY_H

#include <CGAL/number_utils.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <nlohmann/json.hpp>
#include "benchmark_tetrahedral_remeshing_common.h"

#include <limits>
#include <vector>

using benchmarking::append_metric_result;

template <typename TriangleMesh,
          typename NamedParameters = CGAL::parameters::Default_named_parameters>
bool has_degenerate_faces(const TriangleMesh& mesh,
                          const NamedParameters& np = CGAL::parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  for(auto f : faces(mesh))
    if(CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(f, mesh, np))
      return true;

  return false;
}

struct Surface_quality
{
  double minimum_edge_length;
  double mean_edge_length;
  double maximum_edge_length;

  double minimum_area;
  double mean_area;
  double maximum_area;
  double total_area;

  double minimum_angle;
  double maximum_angle;
};

struct Volume_quality
{
  double minimum_volume;
  double mean_volume;
  double maximum_volume;
  double total_volume;

  double minimum_dihedral_angle;
  double mean_dihedral_angle;
  double maximum_dihedral_angle;

  double smallest_edge_radius_ratio;
  double smallest_radius_radius_ratio;
  double biggest_v_sma;
};

template <typename Remeshing_triangulation>
void generate_surface_quality_metrics(const Remeshing_triangulation& tr,
                                     Surface_quality& surface_quality)
{
  surface_quality.minimum_edge_length = -1.;
  surface_quality.mean_edge_length = -1.;
  surface_quality.maximum_edge_length = -1.;

  surface_quality.minimum_area = -1.;
  surface_quality.mean_area = -1.;
  surface_quality.maximum_area = -1.;

  surface_quality.minimum_angle = -1.;
  surface_quality.maximum_angle = -1.;

  surface_quality.total_area = -1.;

  // edge length
  if(tr.number_of_finite_edges() > 0)
  {
    std::vector<double> edge_lengths;

    for(auto eit = tr.finite_edges_begin(); eit != tr.finite_edges_end(); ++eit)
    {
      auto c = eit->first;
      int i = eit->second;
      int j = eit->third;
      const auto& pi = c->vertex(i)->point();
      const auto& pj = c->vertex(j)->point();

      const double edge_length = CGAL::approximate_sqrt(CGAL::squared_distance(pi, pj));
      edge_lengths.push_back(edge_length);
    }

    surface_quality.minimum_edge_length = *std::min_element(edge_lengths.begin(), edge_lengths.end());
    surface_quality.maximum_edge_length = *std::max_element(edge_lengths.begin(), edge_lengths.end());
    surface_quality.mean_edge_length = std::accumulate(edge_lengths.begin(), edge_lengths.end(), 0.) / edge_lengths.size();
  }

  // area
  if(tr.number_of_finite_facets() > 0)
  {
    std::vector<double> areas;

    for(auto fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)
    {
      auto c = fit->first;
      int s = fit->second;


      const auto& pi = c->vertex((s+1)%4)->point();
      const auto& pj = c->vertex((s+2)%4)->point();
      const auto& pk = c->vertex((s+3)%4)->point();

      const double area = CGAL::approximate_sqrt(CGAL::squared_area(pi, pj, pk));
      areas.push_back(area);
    }

    surface_quality.minimum_area = *std::min_element(areas.begin(), areas.end());
    surface_quality.maximum_area = *std::max_element(areas.begin(), areas.end());
    surface_quality.total_area = std::accumulate(areas.begin(), areas.end(), 0.);
    surface_quality.mean_area = surface_quality.total_area / areas.size();
  }

  // angle
  if(tr.number_of_finite_facets() > 0)
  {
    std::vector<double> angles;

    for(auto fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)
    {
      auto c = fit->first;
      int s = fit->second;
      std::array<int, 3> indices = {(s+1)%4, (s+2)%4, (s+3)%4};

      for(int o=0; o<3; ++o)
      {
        const auto& ei = c->vertex(indices[o])->point();
        const auto& ej = c->vertex(indices[(o+1)%3])->point();
        const auto& ek = c->vertex(indices[(o+2)%3])->point();

        const double angle = CGAL::approximate_angle(ei, ej, ek);
        angles.push_back(angle);
      }
    }

    surface_quality.minimum_angle = *std::min_element(angles.begin(), angles.end());
    surface_quality.maximum_angle = *std::max_element(angles.begin(), angles.end());
  }
}

template <typename Remeshing_triangulation>
void generate_volume_quality_metrics(const Remeshing_triangulation& tr,
                                     Volume_quality& volume_quality)
{
  volume_quality.minimum_dihedral_angle = (std::numeric_limits<double>::max)();
  volume_quality.mean_dihedral_angle = 0;
  volume_quality.maximum_dihedral_angle = 0;

  volume_quality.minimum_volume = (std::numeric_limits<double>::max)();
  volume_quality.mean_volume = 0;
  volume_quality.maximum_volume = 0;

  volume_quality.total_volume = 0;

  volume_quality.smallest_edge_radius_ratio = (std::numeric_limits<double>::max)();
  volume_quality.smallest_radius_radius_ratio = (std::numeric_limits<double>::max)();
  volume_quality.biggest_v_sma = 0;

  // volume
  if(tr.number_of_finite_cells() > 0)
  {
    std::vector<double> volumes;

    for(auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
    {
      const auto& p0 = cit->vertex(0)->point();
      const auto& p1 = cit->vertex(1)->point();
      const auto& p2 = cit->vertex(2)->point();
      const auto& p3 = cit->vertex(3)->point();

      const double volume = CGAL::volume(p0, p1, p2, p3);
      volumes.push_back(volume);
    }

    volume_quality.minimum_volume = *std::min_element(volumes.begin(), volumes.end());
    volume_quality.maximum_volume = *std::max_element(volumes.begin(), volumes.end());
    volume_quality.total_volume = std::accumulate(volumes.begin(), volumes.end(), 0.);
    volume_quality.mean_volume = volume_quality.total_volume / volumes.size();
  }

  // dihedral angle
  if(tr.number_of_finite_cells() > 0)
  {
    std::vector<double> dihedral_angles;

    for(auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
    {
      std::array<typename Remeshing_triangulation::Point, 4> pts;
      for(int i=0; i<4; ++i)
        pts[i] = cit->vertex(i)->point();

      constexpr std::array<int, 24> dhis = { 0, 1, 2, 3,
                                             2, 0, 1, 3,
                                             0, 3, 1, 2,
                                             2, 1, 3, 0,
                                             3, 1, 0, 2,
                                             3, 2, 1, 0 };
      for(int dhi=0; dhi<6; ++dhi)
      {
        double dh_angle = CGAL::approximate_dihedral_angle(pts[dhis[4*dhi]], pts[dhis[4*dhi+1]],
                                                           pts[dhis[4*dhi+2]], pts[dhis[4*dhi+3]]);
        dihedral_angles.push_back(std::abs(dh_angle));
      }
    }

    volume_quality.minimum_dihedral_angle = *std::min_element(dihedral_angles.begin(), dihedral_angles.end());
    volume_quality.maximum_dihedral_angle = *std::max_element(dihedral_angles.begin(), dihedral_angles.end());
    volume_quality.mean_dihedral_angle = std::accumulate(dihedral_angles.begin(), dihedral_angles.end(), 0.) / dihedral_angles.size();
  }

  // smallest edge-radius ratio
  if(tr.number_of_finite_cells() > 0)
  {
    for(auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
    {
      const auto& p0 = cit->vertex(0)->point();
      const auto& p1 = cit->vertex(1)->point();
      const auto& p2 = cit->vertex(2)->point();
      const auto& p3 = cit->vertex(3)->point();
      double v = std::abs(CGAL::volume(p0, p1, p2, p3));
      double circumradius = CGAL::approximate_sqrt(CGAL::squared_radius(p0, p1, p2, p3));

      double edges[6];
      edges[0] = std::sqrt(CGAL::squared_distance(p0, p1));
      edges[1] = std::sqrt(CGAL::squared_distance(p0, p2));
      edges[2] = std::sqrt(CGAL::squared_distance(p0, p3));
      edges[3] = std::sqrt(CGAL::squared_distance(p2, p1));
      edges[4] = std::sqrt(CGAL::squared_distance(p2, p3));
      edges[5] = std::sqrt(CGAL::squared_distance(p1, p3));

      double min_edge = edges[0];
      for(int i=1; i<6; ++i)
      {
       if(edges[i]<min_edge)
         min_edge=edges[i];
      }

      double sumar = CGAL::approximate_sqrt(CGAL::squared_area(p0,p1,p2))
                   + CGAL::approximate_sqrt(CGAL::squared_area(p1,p2,p3))
                   + CGAL::approximate_sqrt(CGAL::squared_area(p2,p3,p0))
                   + CGAL::approximate_sqrt(CGAL::squared_area(p3,p1,p0));

      double inradius = 3 * v / sumar;

      // sqrt(6)/4 so that the perfect tet ratio is 1
      double smallest_edge_radius = std::sqrt(6) * min_edge / ( 4. * circumradius);
      if(smallest_edge_radius < volume_quality.smallest_edge_radius_ratio)
        volume_quality.smallest_edge_radius_ratio = smallest_edge_radius;

      // 3 so that the perfect tet ratio is 1 instead of 1/3
      double smallest_radius_radius = 3 * inradius / circumradius;
      if(smallest_radius_radius < volume_quality.smallest_radius_radius_ratio)
        volume_quality.smallest_radius_radius_ratio = smallest_radius_radius;

      // 6*sqrt(2) so that the perfect tet ratio is 1 instead
      double biggest_v_sma_cube = 6 * std::sqrt(2) * v / std::pow(min_edge, 3);
      if(biggest_v_sma_cube > volume_quality.biggest_v_sma)
        volume_quality.biggest_v_sma = biggest_v_sma_cube;
    }
  }
}

template <typename Remeshing_triangulation>
void generate_quality_metrics(const Remeshing_triangulation& tr, nlohmann::json& results_json)
{
  Surface_quality surface_quality;
  Volume_quality volume_quality;

  generate_surface_quality_metrics(tr, surface_quality);
  generate_volume_quality_metrics(tr, volume_quality);

  // Edge Length
  append_metric_result(results_json, "Quality", "Edge_Length", "Minimum", surface_quality.minimum_edge_length);
  append_metric_result(results_json, "Quality", "Edge_Length", "Mean", surface_quality.mean_edge_length);
  append_metric_result(results_json, "Quality", "Edge_Length", "Maximum", surface_quality.maximum_edge_length);

  // Facet Area
  append_metric_result(results_json, "Quality", "Facet_Area", "Minimum", surface_quality.minimum_area);
  append_metric_result(results_json, "Quality", "Facet_Area", "Mean", surface_quality.mean_area);
  append_metric_result(results_json, "Quality", "Facet_Area", "Maximum", surface_quality.maximum_area);
  append_metric_result(results_json, "Quality", "Facet_Area", "Total", surface_quality.total_area);

  // Facet Angle
  append_metric_result(results_json, "Quality", "Facet_Angle", "Minimum", surface_quality.minimum_angle);
  append_metric_result(results_json, "Quality", "Facet_Angle", "Maximum", surface_quality.maximum_angle);

  // Cell Volume
  append_metric_result(results_json, "Quality", "Cell_Volume", "Minimum", volume_quality.minimum_volume);
  append_metric_result(results_json, "Quality", "Cell_Volume", "Mean", volume_quality.mean_volume);
  append_metric_result(results_json, "Quality", "Cell_Volume", "Maximum", volume_quality.maximum_volume);
  append_metric_result(results_json, "Quality", "Cell_Volume", "Total", volume_quality.total_volume);

  // Dihedral Angle
  append_metric_result(results_json, "Quality", "Dihedral_Angle", "Minimum", volume_quality.minimum_dihedral_angle);
  append_metric_result(results_json, "Quality", "Dihedral_Angle", "Mean", volume_quality.mean_dihedral_angle);
  append_metric_result(results_json, "Quality", "Dihedral_Angle", "Maximum", volume_quality.maximum_dihedral_angle);

  // Other quality metrics
  append_metric_result(results_json, "Quality", "Edge_Radius_Ratio", "Minimum", volume_quality.smallest_edge_radius_ratio);
  append_metric_result(results_json, "Quality", "Radius_Radius_Ratio", "Minimum", volume_quality.smallest_radius_radius_ratio);
  append_metric_result(results_json, "Quality", "V_SMA", "Maximum", volume_quality.biggest_v_sma);
}

#endif // CGAL_MESH_3_BENCHMARK_MESH_3_MESH_QUALITY_H
