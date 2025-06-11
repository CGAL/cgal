#ifndef ZONE_EXP_H
#define ZONE_EXP_H
#include "CGAL/Arr_trapezoid_ric_point_location.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_2/Arr_compute_zone_visitor.h>
#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/draw_arrangement_2.h>

void test_zone() {
  using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Segment_traits = CGAL::Arr_segment_traits_2<Exact_kernel>;
  using Point_2 = Segment_traits::Point_2;
  using Arrangement = CGAL::Arrangement_2<Segment_traits>;
  using X_monotone_curve_2 = Segment_traits::X_monotone_curve_2;
  using Traits = CGAL::Arr_linear_traits_2<Exact_kernel>;
  using Face_const_handle = Arrangement::Face_const_handle;
  using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
  using Vertex_const_handle = Arrangement::Vertex_const_handle;
  using Face_handle = Arrangement::Face_handle;
  using Halfedge_handle = Arrangement::Halfedge_handle;
  using Vertex_handle = Arrangement::Vertex_handle;
  using Halfedge_const_handle = Arrangement::Halfedge_const_handle;

  Arrangement arr;
  auto traits = arr.traits();

  auto cst_x_curve = traits->construct_x_monotone_curve_2_object();
  auto tri = {
      cst_x_curve({-1, 0}, {1, 0}),
      cst_x_curve({1, 0}, {0, 3}),
      cst_x_curve({0, 3}, {-1, 0}),
  };
  CGAL::insert(arr, tri.begin(), tri.end());

  // point location
  auto pl = CGAL::Arr_trapezoid_ric_point_location<Arrangement>(arr);

  // test zone construction
  using Feature = std::variant<Vertex_handle, Halfedge_handle, Face_handle>;
  using Feature_vector = std::vector<Feature>;
  using Feature_vector_out_iter = std::back_insert_iterator<Feature_vector>;
  using Zone_visitor = CGAL::Arr_compute_zone_visitor<Arrangement, Feature_vector_out_iter>;

  Feature_vector zone_features;
  Feature_vector_out_iter out_iter(zone_features);
  Zone_visitor zone_visitor(out_iter);
  auto cv = cst_x_curve({-1, 0}, {1, 0});
  CGAL::Arrangement_zone_2<Arrangement, Zone_visitor> zone(arr, &zone_visitor);

  zone.init(cv, pl);
  zone.compute_zone();

  for(const auto& feature : zone_features) {
    if(auto* hf = std::get_if<Halfedge_handle>(&feature)) {
      std::cout << "Halfedge: " << (*hf)->curve() << std::endl;
      Halfedge_const_handle che(*hf);
      std::cout << "  const curve" << che->curve() << std::endl;
    } else if(auto* vh = std::get_if<Vertex_handle>(&feature)) {
      std::cout << "Vertex: " << (*vh)->point() << std::endl;
    } else if(auto* fh = std::get_if<Face_handle>(&feature)) {
      std::cout << "Face: " << std::endl;
    }
  }

  CGAL::draw(arr);
}
#endif // ZONE_EXP_H