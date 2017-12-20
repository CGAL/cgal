#ifndef CGAL_OPTIMAL_TRANSPORTATION_RECONSTRUCTION_KERNEL_2_H_
#define CGAL_OPTIMAL_TRANSPORTATION_RECONSTRUCTION_KERNEL_2_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>

#include <CGAL/OTR_2/Sample.h>
#include <utility>      // std::pair
#include <vector>

//Qt
#include <QColor>

class GlViewer;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::FT FT;

typedef std::pair<Point, FT> PointMassPair;
typedef std::vector<PointMassPair> PointMassList;
typedef PointMassList::const_iterator InputIterator;

typedef CGAL::value_type_traits<InputIterator>::type MassPoint;

typedef CGAL::First_of_pair_property_map<PointMassPair> Point_property_map;
typedef CGAL::Second_of_pair_property_map<PointMassPair> Mass_property_map;

typedef CGAL::Optimal_transportation_reconstruction_2<K, Point_property_map,
    Mass_property_map> Otr_2;

class Optimal_transportation_reconstruction_kerneled_2:
    public Otr_2 {

public:

  template <typename InputRange>
  Optimal_transportation_reconstruction_kerneled_2(const InputRange& input_range,
      Point_property_map point_pmap, Mass_property_map mass_pmap) :
      Otr_2(input_range, point_pmap,
          mass_pmap),
      viewer(NULL){
  }

  Optimal_transportation_reconstruction_kerneled_2() :
    Otr_2() {
  }

  // RENDER //
  void print_stats() const;

  QColor get_color(float value) const;

  void draw_point(const Point& point);

  void draw_segment(const Point& s, const Point& t);

  void draw_edge(const Edge& edge);

  void draw_face(Face_handle face);

  void draw_edge_with_arrow(const Point& s, const Point& t);

  void draw_vertices(const float point_size, const float red,
      const float green, const float blue);

  void draw_edges(const float line_width, const float red, const float green,
      const float blue);

  void draw_footpoints(const float line_width, const float red,
      const float green, const float blue);

  void draw_mesh_footpoints(const Triangulation& mesh, const float line_width,
      const float red, const float green, const float blue);

  void draw_edge_footpoints(const Triangulation& mesh, const Edge& edge,
      const float red, const float green, const float blue);

  void draw_pedges(const float line_width);

  void draw_one_pedge(const Edge& edge, const FT value, const FT min_value,
      const FT max_value, const float line_width);

  void draw_costs(const float line_width, const bool view_ghost);

  void draw_one_cost(const Edge& edge, const FT min_value, const FT max_value,
      const bool view_ghost);

  void draw_relevance(const float line_width, const int nb);

  void draw_bins(const float thickness);

  void draw_bins_plan0(const Edge& edge);

  void draw_bins_plan1(const Edge& edge);

  void draw_relocation();

  void draw_bezier_curves(const unsigned int nb);

  void draw_one_bezier_curve(const Edge& edge, const unsigned int nb);

  bool locate_edge(const Point& query, Edge& edge);

  void draw_one_ring(const float point_size, const float line_width,
      const Point& query);

  void draw_mesh_one_ring(const float point_size, const float line_width,
      const Triangulation& mesh, const Edge& edge);

  void draw_blocking_edges(const float point_size, const float line_width,
      const Point& query);

  void draw_mesh_blocking_edges(const float point_size,
      const float line_width, const Triangulation& mesh,
      const Edge& edge);

  void draw_collapsible_edge(const float point_size, const float line_width,
      const Point& query);

  void draw_cost_stencil(const float point_size, const float line_width,
      const Point& query);

  void draw_remove_queue_stencil(const float point_size,
      const float line_width, const Point& query);

  void draw_push_queue_stencil(const float point_size, const float line_width,
      const Point& query);

  void draw_bg_faces(const Triangulation& mesh, const float red,
      const float green, const float blue, const float alpha);

  void draw_bg_edges(const Triangulation& mesh, const float ri,
      const float gi, const float bi, const float ro, const float go,
      const float bo);

  void draw_bg_vertices(const Triangulation& mesh, const float red,
      const float green, const float blue);

  void draw_vertex_faces(Vertex_handle vertex, const Triangulation& mesh,
      const float red, const float green, const float blue,
      const float alpha);

  void draw_vertex_edges(Vertex_handle vertex, const Triangulation& mesh,
      const float ri, const float gi, const float bi, const float ro,
      const float go, const float bo);

  void save_edges(std::ofstream& ofs, const int nb);

  void save_one_edge(std::ofstream& ofs, const Edge& edge);

  void setViewer(GlViewer* vi) { viewer = vi; }
protected:
  GlViewer* viewer;
};
#endif
