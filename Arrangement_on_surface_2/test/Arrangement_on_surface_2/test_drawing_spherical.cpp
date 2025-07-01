#include "CGAL/Constrained_triangulation_2.h"
#include "CGAL/Graphics_scene.h"
#include "CGAL/Triangulation_data_structure_2.h"
#include "CGAL/Triangulation_vertex_base_with_info_2.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>
#include "CGAL/Random.h"
#include "CGAL/IO/Color.h"
#include "CGAL/draw_triangulation_2.h"
#include "CGAL/mark_domain_in_triangulation.h"
#include <CGAL/Graphics_scene_options.h>

template <class PT>
struct Polygon_triangulation_gs_options : public CGAL::Graphics_scene_options<typename PT::Triangulation,
                                                                              typename PT::Vertex_handle,
                                                                              typename PT::Finite_edges_iterator,
                                                                              typename PT::Finite_faces_iterator>
{
  using T2 = typename PT::Triangulation;
  template <class IPM>
  Polygon_triangulation_gs_options(IPM ipm) {
    this->colored_face = [](const T2&, const typename PT::Finite_faces_iterator) -> bool { return true; };

    this->face_color = [ipm](const T2&, const typename PT::Finite_faces_iterator fh) -> CGAL::IO::Color {
      if(!get(ipm, fh)) {
        std::cout << "Face not in domain, returning black color." << std::endl;
        return CGAL::IO::Color(0, 0, 0);
      }
      CGAL::Random random((unsigned int)(std::size_t)(&*fh));
      return CGAL::get_random_color(random);
    };

    this->draw_face = [ipm](const T2&, const typename PT::Finite_faces_iterator fh) -> bool { return true; };

    this->draw_edge = [ipm](const T2& pt, const typename PT::Finite_edges_iterator eh) -> bool {
      typename PT::Face_handle fh1 = eh->first;
      typename PT::Face_handle fh2 = pt.mirror_edge(*eh).first;
      return get(ipm, fh1) || get(ipm, fh2);
    };
  }
};

int main() {
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Vb = CGAL::Triangulation_vertex_base_with_info_2<size_t, K>;
  using Fb = CGAL::Constrained_triangulation_face_base_2<K>;
  using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Cdt = CGAL::Constrained_triangulation_2<K, Tds, CGAL::Exact_predicates_tag>;
  using Point = K::Point_2;

  Cdt cdt;
  std::vector<Point> outer = {{4, 1},       {4, 4},       {2.5, 4},     {2.32026, 3.28104}, {2.32026, 2.64052},
                              {2.5, 3},     {3, 2},       {2.32026, 2}, {2.32026, 3.28104}, {2.5, 4},
                              {2.32026, 4}, {2.32026, 1}, {4, 1}};
  std::vector<Point> inner = {Point(0, 0), Point(0.5, 0), Point(0.5, 0.5), Point(0.5, 0),
                              Point(1, 0), Point(1, 1),   Point(0, 1)};
  std::vector<Point> outer_constraint = {{4, 1},
                                         {4, 4},
                                         {2.5, 4},
                                         {2.32026, 3.28104},
                                         {1.32026, 3.28104},
                                         {1.32026, 2.64052},
                                         {2.32026, 2.64052},
                                         {2.5, 3},
                                         {3, 2},
                                         {2.32026, 2},
                                         {1.32026, 2},
                                         {1.32026, 3.28104},
                                         {2.32026, 3.28104},
                                         {2.5, 4},
                                         {2.32026, 4},
                                         {1.32026, 4},
                                         {1.32026, 1},
                                         {2.32026, 1},
                                         {4, 1}};

  auto add_info = [](const Point& p) { return std::make_pair(p, 1); };

  cdt.insert_with_info<std::pair<Point, size_t>>(boost::make_transform_iterator(outer.begin(), add_info),
                                                 boost::make_transform_iterator(outer.end(), add_info));
  // cdt.insert_with_info<std::pair<Point, size_t>>(boost::make_transform_iterator(inner.begin(), add_info),
  //                                                boost::make_transform_iterator(inner.end(), add_info));
  cdt.insert_constraint(outer_constraint.begin(), outer_constraint.end(), true);
  // cdt.insert_constraint(inner.begin(), inner.end(), true);
  using In_domain_map = CGAL::unordered_flat_map<typename Cdt::Face_handle, bool>;
  In_domain_map in_domain_map;
  boost::associative_property_map<In_domain_map> in_domain(in_domain_map);

  CGAL::mark_domain_in_triangulation(cdt, in_domain);

  std::cout << "Number of faces in triangulation: " << cdt.number_of_faces() << std::endl;

  CGAL::Graphics_scene_options<Cdt::Triangulation, Cdt::Vertex_handle, Cdt::Finite_edges_iterator,
                               Cdt::Finite_faces_iterator>
      gso;
  gso.face_color = [&](const Cdt::Triangulation&, const typename Cdt::Finite_faces_iterator fh) -> CGAL::IO::Color {
    if(!in_domain_map[fh]) {
      return CGAL::IO::Color(255, 255, 255); // black for faces not in domain
    }
    std::array<std::size_t, 3> vertices;
    for(int i = 0; i < 3; ++i) {
      vertices[i] = fh->vertex(i)->info();
    }
    if(std::any_of(vertices.begin(), vertices.end(), [](auto idx) { return idx == 0; })) {
      return CGAL::IO::Color(255, 255, 255);
    }
    CGAL::Random rand((std::size_t)(&*fh));
    return CGAL::get_random_color(rand);
  };
  gso.colored_face = [&](const Cdt::Triangulation&, const typename Cdt::Finite_faces_iterator) -> bool {
    return true; // always color faces
  };

  CGAL::draw(cdt, gso, "Polygon Triangulation Viewer");
}