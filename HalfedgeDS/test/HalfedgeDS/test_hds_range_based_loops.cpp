#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HalfedgeDS_list.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_items_3.h>

using Kernel     = CGAL::Simple_cartesian<double>;
using HDS_list   = CGAL::HalfedgeDS_list<Kernel, CGAL::Polyhedron_items_3>;
using HDS_vector = CGAL::HalfedgeDS_vector<Kernel, CGAL::Polyhedron_items_3>;

// void test_vertex_handles_and_points(
//   Polyhedron& polyhedron) {

//   auto pit = polyhedron.points_begin();
//   auto vit = polyhedron.vertices_begin();
//   for (auto vh : polyhedron.vertex_handles()) {
//     assert(vh == vit);
//     assert(vh->point() == vit->point());
//     assert(vh->point() == *pit);
//     ++vit; ++pit;
//   }
//   assert(pit == polyhedron.points_end());
//   assert(vit == polyhedron.vertices_end());

//   pit = polyhedron.points_begin();
//   vit = polyhedron.vertices_begin();
//   for (auto& point : polyhedron.points()) {
//     assert(*pit == point);
//     assert(vit->point() == point);
//     ++vit; ++pit;
//   }
//   assert(pit == polyhedron.points_end());
//   assert(vit == polyhedron.vertices_end());
// }

// void test_const_vertex_handles_and_points(
//   const Polyhedron& polyhedron) {

//   auto pit = polyhedron.points_begin();
//   auto vit = polyhedron.vertices_begin();
//   for (const auto vh : polyhedron.vertex_handles()) {
//     assert(vh == vit);
//     assert(vh->point() == vit->point());
//     assert(vh->point() == *pit);
//     ++vit; ++pit;
//   }
//   assert(pit == polyhedron.points_end());
//   assert(vit == polyhedron.vertices_end());

//   pit = polyhedron.points_begin();
//   vit = polyhedron.vertices_begin();
//   for (const auto& point : polyhedron.points()) {
//     assert(*pit == point);
//     assert(vit->point() == point);
//     ++vit; ++pit;
//   }
//   assert(pit == polyhedron.points_end());
//   assert(vit == polyhedron.vertices_end());
// }

// void test_facet_handles_and_planes(
//   Polyhedron& polyhedron) {

//   auto pit = polyhedron.planes_begin();
//   auto fit = polyhedron.facets_begin();
//   for (auto fh : polyhedron.facet_handles()) {
//     assert(fh == fit);
//     assert(fh->plane() == fit->plane());
//     assert(fh->plane() == *pit);
//     ++fit; ++pit;
//   }
//   assert(pit == polyhedron.planes_end());
//   assert(fit == polyhedron.facets_end());

//   pit = polyhedron.planes_begin();
//   fit = polyhedron.facets_begin();
//   for (auto& plane : polyhedron.planes()) {
//     assert(*pit == plane);
//     assert(fit->plane() == plane);
//     ++fit; ++pit;
//   }
//   assert(pit == polyhedron.planes_end());
//   assert(fit == polyhedron.facets_end());
// }

// void test_const_facet_handles_and_planes(
//   const Polyhedron& polyhedron) {

//   auto pit = polyhedron.planes_begin();
//   auto fit = polyhedron.facets_begin();
//   for (const auto fh : polyhedron.facet_handles()) {
//     assert(fh == fit);
//     assert(fh->plane() == fit->plane());
//     assert(fh->plane() == *pit);
//     ++fit; ++pit;
//   }
//   assert(pit == polyhedron.planes_end());
//   assert(fit == polyhedron.facets_end());

//   pit = polyhedron.planes_begin();
//   fit = polyhedron.facets_begin();
//   for (const auto& plane : polyhedron.planes()) {
//     assert(*pit == plane);
//     assert(fit->plane() == plane);
//     ++fit; ++pit;
//   }
//   assert(pit == polyhedron.planes_end());
//   assert(fit == polyhedron.facets_end());
// }

// void test_halfedge_handles_and_edges(
//   Polyhedron& polyhedron) {

//   auto hit = polyhedron.halfedges_begin();
//   for (auto hh : polyhedron.halfedge_handles()) {
//     assert(hh == hit);
//     assert(hh->facet() == hit->facet());
//     assert(hh->vertex() == hit->vertex());
//     ++hit;
//   }
//   assert(hit == polyhedron.halfedges_end());

//   auto eit = polyhedron.edges_begin();
//   for (auto& edge : polyhedron.edges()) {
//     assert((*eit).facet() == edge.facet());
//     assert((*eit).vertex() == edge.vertex());
//     ++eit;
//   }
//   assert(eit == polyhedron.edges_end());
// }

// void test_const_halfedge_handles_and_edges(
//   const Polyhedron& polyhedron) {

//   auto hit = polyhedron.halfedges_begin();
//   for (const auto hh : polyhedron.halfedge_handles()) {
//     assert(hh == hit);
//     assert(hh->facet() == hit->facet());
//     assert(hh->vertex() == hit->vertex());
//     ++hit;
//   }
//   assert(hit == polyhedron.halfedges_end());

//   auto eit = polyhedron.edges_begin();
//   for (const auto& edge : polyhedron.edges()) {
//     assert((*eit).facet() == edge.facet());
//     assert((*eit).vertex() == edge.vertex());
//     ++eit;
//   }
//   assert(eit == polyhedron.edges_end());
// }

int main() {

  HDS_list hds_list(1, 2, 2);
  HDS_vector hds_vector(1, 2, 2);

  // test_vertex_handles_and_points(polyhedron);
  // test_const_vertex_handles_and_points(polyhedron);
  // test_facet_handles_and_planes(polyhedron);
  // test_const_facet_handles_and_planes(polyhedron);
  // test_halfedge_handles_and_edges(polyhedron);
  // test_const_halfedge_handles_and_edges(polyhedron);

  std::cout << "test_hds_range_based_loops: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
