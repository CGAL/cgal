#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#ifdef USE_POLYHEDRON
#include <CGAL/Polyhedron_3.h>
#else
#include <CGAL/Surface_mesh.h>
#endif
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>

#include <iostream>
#include <fstream>
#include <sstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
#ifdef USE_POLYHEDRON
typedef CGAL::Polyhedron_3<Kernel> Surface_mesh;
#else
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
#endif


namespace PMP = CGAL::Polygon_mesh_processing;
int main()
{
// testing decimate function
  bool OK = true;
  const int nb_meshes=5;

  for (int i=1; i<=nb_meshes; ++i)
  {
    std::cout << "handling decimation of data/decimation/m" << i << ".off\n";
    std::stringstream ss;
    ss << "data/decimation/m" << i << ".off";
    Surface_mesh sm;
    std::ifstream in(ss.str().c_str());
    in >> sm;

    // call the decimation function
    Surface_mesh out;
    if (!PMP::remesh_planar_patches(sm, out))
    {
      std::cerr << "ERROR: decimate failed to remesh some patches\n";
      OK=false;
    }
    ss=std::stringstream();
    ss << "out" << i << ".off";
    std::ofstream out(ss.str().c_str());
    out << out;
    std::cout << "  output written to out" << i << ".off\n";
    assert(is_valid_polygon_mesh(out));
  }
// testing border non-manifold vertex: not working for now, test kept
/*
// in case we find a solution
  {
    std::cout << "testing handling of non-manifold patches\n";
    Surface_mesh sm;
    std::ifstream("data/decimation/m1.off") >> sm;
    auto f1 = *std::next(faces(sm).begin(), 594);
    auto f2 = *std::next(faces(sm).begin(), 2378);
    CGAL::Euler::remove_face(halfedge(f1, sm), sm);
    CGAL::Euler::remove_face(halfedge(f2, sm), sm);
    if (!PMP::remesh_planar_patches(sm))
    {
      OK=false;
      std::cerr << "ERROR: decimate failed to remesh some patches\n";
    }
    std::ofstream("nm_m1.off") << std::setprecision(17) << sm;
    assert(is_valid_polygon_mesh(sm));
  }
*/
// test duplicated vertex
  {
    std::cout << "testing handling of duplicated non-manifold vertex\n";
    Surface_mesh sm;
    std::ifstream("data/decimation/m1.off") >> sm;
    auto f1 = *std::next(faces(sm).begin(), 594);
    auto f2 = *std::next(faces(sm).begin(), 2378);
    CGAL::Euler::remove_face(halfedge(f1, sm), sm);
    CGAL::Euler::remove_face(halfedge(f2, sm), sm);
    PMP::duplicate_non_manifold_vertices(sm);
    std::size_t nbv_before = vertices(sm).size();
    Surface_mesh out;
    if (!PMP::remesh_planar_patches(sm, out))
      std::cerr << "decimate failed to remesh some patches (expected)\n";
    assert(vertices(out).size()<nbv_before);
    std::ofstream("nmd_m1.off") << std::setprecision(17) << out;
    assert(is_valid_polygon_mesh(out));
  }
  // test duplicated vertex at patch interface
  {
    std::cout << "testing handling of duplicated non-manifold vertex at patch interface\n";
    Surface_mesh sm;
    std::ifstream("data/decimation/m1.off") >> sm;
    auto f1 = *std::next(faces(sm).begin(), 244);
    auto f2 = *std::next(faces(sm).begin(), 2279);
    CGAL::Euler::remove_face(halfedge(f1, sm), sm);
    CGAL::Euler::remove_face(halfedge(f2, sm), sm);
    PMP::duplicate_non_manifold_vertices(sm);
    std::size_t nbv_before = vertices(sm).size();
    Surface_mesh out;
    if (!PMP::remesh_planar_patches(sm, out))
      std::cerr << "decimate failed to remesh some patches (expected)\n";
    assert(vertices(out).size()<nbv_before);
    std::ofstream("nmdi_m1.off") << std::setprecision(17) << out;
    assert(is_valid_polygon_mesh(out));
  }
  assert(OK);

// testing decimate function with almost coplanar/collinear tests
  for (int i=1; i<=nb_meshes; ++i)
  {
    std::cout << "handling decimation of data/decimation/am" << i << ".off (approximate coplanar/collinear)\n";
    std::stringstream ss;
    ss << "data/decimation/am" << i << ".off";
    Surface_mesh sm;
    std::ifstream in(ss.str().c_str());
    in >> sm;

    // call the decimation function
    Surface_mesh out;
    if (!PMP::remesh_planar_patches(sm, out, CGAL::parameters::cosinus_threshold(-0.99)))
    {
      OK=false;
      std::cerr << "ERROR: decimate failed to remesh some patches\n";
    }
    ss=std::stringstream();
    ss << "out_a" << i << ".off";
    std::ofstream out(ss.str().c_str());
    out << out;
    std::cout << "  output written to out_a" << i << ".off\n";
    assert(is_valid_polygon_mesh(out));
  }

//testing decimation of meshes, preserving common interface
  const int nb_meshes_range=9;
  std::vector<Surface_mesh> meshes(nb_meshes_range);
  for (int i=1; i<=nb_meshes_range; ++i)
  {
    std::stringstream ss;
    ss << "data/decimation/range/m" << i << ".off";
    std::ifstream in(ss.str().c_str());
    in >> meshes[i-1];
  }

  std::cout << "decimate a range of meshes with common interfaces\n";
  if (!PMP::decimate_meshes_with_common_interfaces(meshes))
  {
    OK=false;
    std::cerr << "ERROR: decimate failed to remesh some patches\n";
  }
  std::cout << "  output written to";
  for (int i=0; i<nb_meshes_range; ++i)
  {
    std::stringstream ss;
    ss << "out_r" << i+1 << ".off";
    std::ofstream out(ss.str().c_str());
    out << meshes[i];
    std::cout << " out_r" << i+1 << ".off";
  }
  std::cout << "\n";
  for (int i=0; i<nb_meshes_range; ++i)
    assert(is_valid_polygon_mesh(meshes[i]));

  //testing decimation of meshes, preserving common interface and a patch that fails to simplify at the interface
  meshes.clear();
  meshes.resize(nb_meshes_range);
  for (int i=1; i<=nb_meshes_range; ++i)
  {
    std::stringstream ss;
    ss << "data/decimation/range/m" << i << ".off";
    std::ifstream in(ss.str().c_str());
    in >> meshes[i-1];
  }
  auto f1 = *std::next(faces(meshes[4]).begin(), 1);
  auto f2 = *std::next(faces(meshes[4]).begin(), 109);
  CGAL::Euler::remove_face(halfedge(f1, meshes[4]), meshes[4]);
  CGAL::Euler::remove_face(halfedge(f2, meshes[4]), meshes[4]);
  PMP::duplicate_non_manifold_vertices(meshes[4]);
  f1 = *std::next(faces(meshes[0]).begin(), 322);
  f2 = *std::next(faces(meshes[0]).begin(), 963);
  CGAL::Euler::remove_face(halfedge(f1, meshes[0]), meshes[0]);
  CGAL::Euler::remove_face(halfedge(f2, meshes[0]), meshes[0]);
  PMP::duplicate_non_manifold_vertices(meshes[0]);
  f1 = *std::next(faces(meshes[8]).begin(), 23);
  f2 = *std::next(faces(meshes[8]).begin(), 164);
  CGAL::Euler::remove_face(halfedge(f1, meshes[8]), meshes[8]);
  CGAL::Euler::remove_face(halfedge(f2, meshes[8]), meshes[8]);
  PMP::duplicate_non_manifold_vertices(meshes[8]);

  std::cout << "decimate a range of meshes with common interfaces and issue at the interface\n";
  if (!PMP::decimate_meshes_with_common_interfaces(meshes))
    std::cerr << "decimate failed to remesh some patches (expected)\n";
  std::cout << "  output written to";
  for (int i=0; i<nb_meshes_range; ++i)
  {
    std::stringstream ss;
    ss << "out_fi_r" << i+1 << ".off";
    std::ofstream out(ss.str().c_str());
    out << meshes[i];
    std::cout << " out_fi_r" << i+1 << ".off";
  }
  std::cout << "\n";
  for (int i=0; i<nb_meshes_range; ++i)
    assert(is_valid_polygon_mesh(meshes[i]));

  //testing decimation of meshes, preserving common interface with almost coplanar/collinear tests
  meshes.clear();
  meshes.resize(nb_meshes_range);
  for (int i=1; i<=nb_meshes_range; ++i)
  {
    std::stringstream ss;
    ss << "data/decimation/range/am" << i << ".off";
    std::ifstream in(ss.str().c_str());
    in >> meshes[i-1];
  }

  std::cout << "decimate a range of meshes with common interfaces (approximate coplanar/collinear)\n";
  if (!PMP::decimate_meshes_with_common_interfaces(meshes, -0.99))
  {
    OK=false;
    std::cerr << "ERROR: decimate failed to remesh some patches\n";
  }
  std::cout << "  output written to";
  for (int i=0; i<nb_meshes_range; ++i)
  {
    std::stringstream ss;
    ss << "out_ar" << i+1 << ".off";
    std::ofstream out(ss.str().c_str());
    out << meshes[i];
    std::cout << " out_ar" << i+1 << ".off";
  }
  std::cout << "\n";
  for (int i=0; i<nb_meshes_range; ++i)
    assert(is_valid_polygon_mesh(meshes[i]));

#if 0 // tests to be re-enable when using region growing
// testing decimate function with almost coplanar/collinear tests using PCA
  for (int i=1; i<=nb_meshes; ++i)
  {
    std::cout << "handling decimation of data/decimation/am" << i << ".off (approximate coplanar/collinear with PCA)\n";
    std::stringstream ss;
    ss << "data/decimation/am" << i << ".off";
    Surface_mesh sm;
    std::ifstream in(ss.str().c_str());
    in >> sm;

    // call the decimation function

    if (!PMP::decimate_with_pca_for_coplanarity(sm, 1e-5, -0.99))
    {
      std::cerr << "ERROR: decimate failed to remesh some patches\n";
      OK=false;
    }
    ss=std::stringstream();
    ss << "out_a_pca" << i << ".off";
    std::ofstream out(ss.str().c_str());
    out << sm;
    std::cout << "  output written to out_a_pca" << i << ".off\n";
    assert(is_valid_polygon_mesh(sm));
  }

// testing decimation of meshes,  preserving common interface with almost coplanar/collinear tests using PCA
  std::cout << "decimate a range of meshes with common interfaces (approximate coplanar/collinear with PCA)\n";
  if (!PMP::decimate_meshes_with_common_interfaces_and_pca_for_coplanarity(meshes, 0.99, -0.99))
    std::cerr << "ERROR: decimate failed to remesh some patches\n";
  std::cout << "  output written to";
  for (int i=0; i<nb_meshes_range; ++i)
  {
    std::stringstream ss;
    ss << "out_ar_pca" << i+1 << ".off";
    std::ofstream out(ss.str().c_str());
    out << meshes[i];
    std::cout << " out_ar_pca" << i+1 << ".off";
  }
  std::cout << "\n";
  for (int i=0; i<nb_meshes_range; ++i)
    assert(is_valid_polygon_mesh(meshes[i]));

// two examples that fails with approximate but works with PCA
  //PCA first
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere.off using PCA\n";
    std::ifstream in("data/decimation/sphere.off");
    in >> sm;
    if (!PMP::decimate_with_pca_for_coplanarity(sm,1e-5,-0.99))
    {
      std::cerr << "ERROR: decimate failed to remesh some patches\n";
    }
    std::ofstream out("sphere_pca.off");
    out << sm;
    std::cout << "output written to sphere_pca.off\n";
    assert(is_valid_polygon_mesh(sm));
  }
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere_selection.off using PCA\n";
    std::ifstream in("data/decimation/sphere_selection.off");
    in >> sm;
    if (!PMP::decimate_with_pca_for_coplanarity(sm,1e-5,-0.99))
      std::cerr << "decimate failed to remesh some patches\n";
    std::ofstream out("sphere_selection_pca.off");
    out << sm;
    std::cout << "output written to sphere_selection_pca.off\n";
    assert(is_valid_polygon_mesh(sm));
  }
#endif
  // Approximation then
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere.off using approximate predicates\n";
    std::ifstream in("data/decimation/sphere.off");
    in >> sm;
    Surface_mesh out;
    if (!PMP::remesh_planar_patches(sm, out, CGAL::parameters::cosinus_threshold(-0.99)))
      std::cerr << "decimate failed to remesh some patches\n";
  }
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere_selection.off using approximate predicates\n";
    std::ifstream in("data/decimation/sphere_selection.off");
    in >> sm;
    Surface_mesh out;
    if (!PMP::remesh_planar_patches(sm, out, CGAL::parameters::cosinus_threshold(-0.99)))
      std::cout << "decimate failed to remesh some patches (this is the expected behavior)\n";
    std::ofstream out("sphere_selection_app.off");
    out << out;
    std::cout << "output written to sphere_selection_app.off\n";
    assert(is_valid_polygon_mesh(out));
  }

  assert(OK);

  return 0 ;
}
