#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HalfedgeDS_list.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_items_3.h>

using Kernel     = CGAL::Simple_cartesian<double>;
using HDS_list   = CGAL::HalfedgeDS_list<Kernel, CGAL::Polyhedron_items_3>;
using HDS_vector = CGAL::HalfedgeDS_vector<Kernel, CGAL::Polyhedron_items_3>;

using VL = typename HDS_list::Vertex;
using HL = typename HDS_list::Halfedge;
using FL = typename HDS_list::Face;

using VV = typename HDS_vector::Vertex;
using HV = typename HDS_vector::Halfedge;
using FV = typename HDS_vector::Face;

void test_vertex_handles(
  HDS_list& hds_list,
  HDS_vector& hds_vector) {

  auto lit = hds_list.vertices_begin();
  assert(hds_list.vertex_handles().size() == 1);
  for (auto vh : hds_list.vertex_handles()) {
    assert(vh == lit);
    assert(vh->point() == lit->point());
    assert(vh->halfedge() == lit->halfedge());
    ++lit;
  }
  assert(lit == hds_list.vertices_end());

  auto vit = hds_vector.vertices_begin();
  assert(hds_vector.vertex_handles().size() == 1);
  for (auto vh : hds_vector.vertex_handles()) {
    assert(vh == vit);
    assert(vh->point() == vit->point());
    assert(vh->halfedge() == vit->halfedge());
    ++vit;
  }
  assert(vit == hds_vector.vertices_end());
}

void test_const_vertex_handles(
  const HDS_list& hds_list,
  const HDS_vector& hds_vector) {

  auto lit = hds_list.vertices_begin();
  assert(hds_list.vertex_handles().size() == 1);
  for (const auto& vh : hds_list.vertex_handles()) {
    assert(vh == lit);
    assert(vh->point() == lit->point());
    assert(vh->halfedge() == lit->halfedge());
    ++lit;
  }
  assert(lit == hds_list.vertices_end());

  auto vit = hds_vector.vertices_begin();
  assert(hds_vector.vertex_handles().size() == 1);
  for (const auto& vh : hds_vector.vertex_handles()) {
    assert(vh == vit);
    assert(vh->point() == vit->point());
    assert(vh->halfedge() == vit->halfedge());
    ++vit;
  }
  assert(vit == hds_vector.vertices_end());
}

void test_face_handles(
  HDS_list& hds_list,
  HDS_vector& hds_vector) {

  auto lit = hds_list.faces_begin();
  assert(hds_list.face_handles().size() == 2);
  for (auto fh : hds_list.face_handles()) {
    assert(fh == lit);
    assert(fh->plane() == lit->plane());
    assert(fh->halfedge() == lit->halfedge());
    ++lit;
  }
  assert(lit == hds_list.faces_end());

  auto vit = hds_vector.faces_begin();
  assert(hds_vector.face_handles().size() == 2);
  for (auto fh : hds_vector.face_handles()) {
    assert(fh == vit);
    assert(fh->plane() == vit->plane());
    assert(fh->halfedge() == vit->halfedge());
    ++vit;
  }
  assert(vit == hds_vector.faces_end());
}

void test_const_face_handles(
  const HDS_list& hds_list,
  const HDS_vector& hds_vector) {

  auto lit = hds_list.faces_begin();
  assert(hds_list.face_handles().size() == 2);
  for (const auto& fh : hds_list.face_handles()) {
    assert(fh == lit);
    assert(fh->plane() == lit->plane());
    assert(fh->halfedge() == lit->halfedge());
    ++lit;
  }
  assert(lit == hds_list.faces_end());

  auto vit = hds_vector.faces_begin();
  assert(hds_vector.face_handles().size() == 2);
  for (const auto& fh : hds_vector.face_handles()) {
    assert(fh == vit);
    assert(fh->plane() == vit->plane());
    assert(fh->halfedge() == vit->halfedge());
    ++vit;
  }
  assert(vit == hds_vector.faces_end());
}

void test_halfedge_handles(
  HDS_list& hds_list,
  HDS_vector& hds_vector) {

  auto lit = hds_list.halfedges_begin();
  assert(hds_list.halfedge_handles().size() == 2);
  for (auto hh : hds_list.halfedge_handles()) {
    assert(hh == lit);
    assert(hh->face() == lit->face());
    assert(hh->vertex() == lit->vertex());
    ++lit;
  }
  assert(lit == hds_list.halfedges_end());

  auto vit = hds_vector.halfedges_begin();
  assert(hds_vector.halfedge_handles().size() == 2);
  for (auto hh : hds_vector.halfedge_handles()) {
    assert(hh == vit);
    assert(hh->face() == vit->face());
    assert(hh->vertex() == vit->vertex());
    ++vit;
  }
  assert(vit == hds_vector.halfedges_end());
}

void test_const_halfedge_handles(
  const HDS_list& hds_list,
  const HDS_vector& hds_vector) {

  auto lit = hds_list.halfedges_begin();
  assert(hds_list.halfedge_handles().size() == 2);
  for (const auto& hh : hds_list.halfedge_handles()) {
    assert(hh == lit);
    assert(hh->face() == lit->face());
    assert(hh->vertex() == lit->vertex());
    ++lit;
  }
  assert(lit == hds_list.halfedges_end());

  auto vit = hds_vector.halfedges_begin();
  assert(hds_vector.halfedge_handles().size() == 2);
  for (const auto& hh : hds_vector.halfedge_handles()) {
    assert(hh == vit);
    assert(hh->face() == vit->face());
    assert(hh->vertex() == vit->vertex());
    ++vit;
  }
  assert(vit == hds_vector.halfedges_end());
}

int main() {

  HDS_list hds_list(1, 2, 2);
  hds_list.vertices_push_back(VL());
  hds_list.edges_push_back(HL(), HL());
  hds_list.faces_push_back(FL());
  hds_list.faces_push_back(FL());

  HDS_vector hds_vector(1, 2, 2);
  hds_vector.vertices_push_back(VV());
  hds_vector.edges_push_back(HV(), HV());
  hds_vector.faces_push_back(FV());
  hds_vector.faces_push_back(FV());

  assert(hds_list.size_of_vertices() == 1);
  assert(hds_list.size_of_halfedges() == 2);
  assert(hds_list.size_of_faces() == 2);

  assert(hds_vector.size_of_vertices() == 1);
  assert(hds_vector.size_of_halfedges() == 2);
  assert(hds_vector.size_of_faces() == 2);

  test_vertex_handles(hds_list, hds_vector);
  test_const_vertex_handles(hds_list, hds_vector);
  test_face_handles(hds_list, hds_vector);
  test_const_face_handles(hds_list, hds_vector);
  test_halfedge_handles(hds_list, hds_vector);
  test_const_halfedge_handles(hds_list, hds_vector);

  std::cout << "test_hds_range_based_loops: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
