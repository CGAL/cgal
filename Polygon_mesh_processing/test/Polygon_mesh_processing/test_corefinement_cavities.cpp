#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

template <class Point, class Mesh>
void make_hexa(double x, double y, double z,
               double X, double Y, double Z,
               Mesh& mesh, int t)
{
  CGAL::make_hexahedron(
    Point(x,y,Z),
    Point(X,y,Z),
    Point(X,y,z),
    Point(x,y,z),
    Point(x,Y,z),
    Point(x,Y,Z),
    Point(X,Y,Z),
    Point(X,Y,z),
  mesh);

  using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;

  std::vector<face_descriptor> fcs(faces(mesh).begin(), faces(mesh).end());
  for (face_descriptor f : fcs)
  {
    halfedge_descriptor h = halfedge(f, mesh);
    if (t==1) h=next(h,mesh);
    halfedge_descriptor h2=next(next(h, mesh), mesh);
    CGAL::Euler::split_face(h, h2, mesh);
  }
}

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

namespace PMP = CGAL::Polygon_mesh_processing;


void test_operations(Mesh A, Mesh B,
                     bool reverse_A, bool reverse_B,
                     std::string round,
                     std::size_t union_v, std::size_t inter_v, std::size_t diff1_v, std::size_t diff2_v)
{
#ifndef VERBOSE
  CGAL_USE(round);
#endif
  if (reverse_A) PMP::reverse_face_orientations(A);
  if (reverse_B) PMP::reverse_face_orientations(B);

  Mesh out_union, out_inter, out_diff1, out_diff2;
  std::array<boost::optional<Mesh*>, 4> output;
  output[PMP::Corefinement::UNION] = &out_union;
  output[PMP::Corefinement::INTERSECTION] = &out_inter;
  output[PMP::Corefinement::TM1_MINUS_TM2] = &out_diff1;
  output[PMP::Corefinement::TM2_MINUS_TM1] = &out_diff2;

  Mesh lA=A, lB=B;
  PMP::corefine_and_compute_boolean_operations(lA,lB,output);
#ifdef VERBOSE
  std::ofstream("out_union_"+round+".off") << out_union;
  std::ofstream("out_inter_"+round+".off") << out_inter;
  std::ofstream("out_diff1_"+round+".off") << out_diff1;
  std::ofstream("out_diff2_"+round+".off") << out_diff2;
#endif
  assert(vertices(out_union).size()==union_v);
  assert(vertices(out_inter).size()==inter_v);
  assert(vertices(out_diff1).size()==diff1_v);
  assert(vertices(out_diff2).size()==diff2_v);
}

int main()
{

  Mesh A, mh, B;
  make_hexa<Point_3>(0, 0, 0,
                     4, 4, 4,
                     A, 0);
  make_hexa<Point_3>(1, 1, 1,
                     2, 2, 2,
                     mh, 0);
  make_hexa<Point_3>(1, 1, 1,
                     2, 2, 2,
                     B, 1);

  Mesh A2, mh2, B2;
  make_hexa<Point_3>(5, 0, 0,
                     9, 4, 4,
                     A2, 0);
  make_hexa<Point_3>(6, 1, 1,
                     7, 2, 2,
                     mh2, 0);
  make_hexa<Point_3>(6, 1, 1,
                     7, 2, 2,
                     B2, 1);

  A.join(A2);
  mh.join(mh2);
  PMP::reverse_face_orientations(mh);
  A.join(mh);
  B.join(B2);

#ifdef VERBOSE
  std::ofstream("A.off") << A;
  std::ofstream("B.off") << B;
#endif

  test_operations(A, B, false, false, "r00", 16, 0, 44, 28);
  test_operations(A, B, false, true,  "r01", 28, 44, 0, 16);
  test_operations(A, B, true,  false, "r10", 44, 28, 16, 0);
  test_operations(A, B, true,  true,  "r11", 0, 16, 28, 44);

  test_operations(A, A, false, false, "a00", 32, 32, 0, 0);
  test_operations(A, A, false, true,  "a01", 0, 0, 32, 32);
  test_operations(A, A, true,  false, "a10", 0, 0, 32, 32);
  test_operations(A, A, true,  true,  "a11", 32, 32, 0, 0);

  test_operations(B, B, false, false, "b00", 16, 16, 0, 0);
  test_operations(B, B, false, true,  "b01", 0, 0, 16, 16);
  test_operations(B, B, true,  false, "b10", 0, 0, 16, 16);
  test_operations(B, B, true,  true,  "b11", 16, 16, 0, 0);
}
