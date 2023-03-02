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
    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out);
    ss=std::stringstream();
    ss << "out" << i << ".off";
    std::ofstream out(ss.str().c_str());
    out << sm_out;
    std::cout << "  output written to out" << i << ".off\n";
    assert(CGAL::is_valid_polygon_mesh(sm_out));
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
    assert(CGAL::is_valid_polygon_mesh(sm));
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
    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out);
    assert(vertices(sm_out).size()<nbv_before);
    std::ofstream("nmd_m1.off") << std::setprecision(17) << sm_out;
    assert(CGAL::is_valid_polygon_mesh(sm_out));
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
    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out);
    assert(vertices(sm_out).size()<nbv_before);
    std::ofstream("nmdi_m1.off") << std::setprecision(17) << sm_out;
    assert(CGAL::is_valid_polygon_mesh(sm_out));
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
    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out, CGAL::parameters::cosine_of_maxium_angle(-0.99));
    ss=std::stringstream();
    ss << "out_a" << i << ".off";
    std::ofstream out(ss.str().c_str());
    out << sm_out;
    std::cout << "  output written to out_a" << i << ".off\n";
    assert(CGAL::is_valid_polygon_mesh(sm_out));
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
    assert(CGAL::is_valid_polygon_mesh(meshes[i]));

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
    assert(CGAL::is_valid_polygon_mesh(meshes[i]));

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
    assert(CGAL::is_valid_polygon_mesh(meshes[i]));

  // test face/vertex maps
  {
    std::cout << "check face patch ids\n";
    std::cout << "  face map alone\n";
    Surface_mesh in, out;
    std::ifstream(CGAL::data_file_path("meshes/cube-meshed.off")) >> in;
    assert(vertices(in).size()!=0);
    auto fmap_in = in.add_property_map<Surface_mesh::Face_index, std::size_t>("f:pid").first;
    auto fmap_out = out.add_property_map<Surface_mesh::Face_index, std::size_t>("f:pid").first;
    auto vmap_in = in.add_property_map<Surface_mesh::Vertex_index, std::size_t>("v:pid").first;
    auto vmap_out = out.add_property_map<Surface_mesh::Vertex_index, std::size_t>("v:pid").first;
    PMP::remesh_planar_patches(in, out,
                               CGAL::parameters::face_patch_map(fmap_in),
                               CGAL::parameters::face_patch_map(fmap_out));

    auto get_id = [](Surface_mesh::Face_index f, Surface_mesh& sm)
    {
      auto h = halfedge(f, sm);
      Point_3 p=sm.point(source(h, sm)), q=sm.point(target(h, sm)), r=sm.point(target(next(h, sm), sm));
      if (p.x()==-1 && q.x()==-1 && r.x()==-1) return 0;
      if (p.x()== 1 && q.x()== 1 && r.x()== 1) return 1;
      if (p.y()==-1 && q.y()==-1 && r.y()==-1) return 2;
      if (p.y()== 1 && q.y()== 1 && r.y()== 1) return 3;
      if (p.z()==-1 && q.z()==-1 && r.z()==-1) return 4;
      assert(p.z()==1 && q.z()==1 && r.z()==1);
      return 5;
    };

    std::array<std::size_t, 6> pids; // xi, Xi, yi, Yi, zi, Zi;
    for (auto f : faces(out)) pids[get_id(f, out)]=get(fmap_out, f);
    for (auto f : faces(in)) assert(pids[get_id(f, in)]==get(fmap_in, f));
    //------------------------------------------------------
    std::cout << "  face and vertex maps\n";
    out.clear_without_removing_property_maps();
    PMP::remesh_planar_patches(in, out,
                               CGAL::parameters::face_patch_map(fmap_in).vertex_corner_map(vmap_in),
                               CGAL::parameters::face_patch_map(fmap_out).vertex_corner_map(vmap_out));

    auto check_corner_ids = [&]()
    {
      std::vector<Point_3> id2pt(vertices(out).size());
      for (auto v : vertices(out))
        id2pt[get(vmap_out, v)]=out.point(v);
      for (auto v : vertices(in))
      {

        if (!( get(vmap_in, v)==std::size_t(-1) || id2pt[get(vmap_in, v)]==in.point(v) ))
          std::cout << get(vmap_in, v) << " vs " << std::size_t(-1) << " vs " << std::size_t(-2) << "\n";
        assert( get(vmap_in, v)==std::size_t(-1) || id2pt[get(vmap_in, v)]==in.point(v) );
      }

    };
    check_corner_ids();
    for (auto f : faces(out)) pids[get_id(f, out)]=get(fmap_out, f);
    for (auto f : faces(in)) assert(pids[get_id(f, in)]==get(fmap_in, f));
    //------------------------------------------------------
    std::cout << "  vertex map alone\n";
    out.clear_without_removing_property_maps();
    PMP::remesh_planar_patches(in, out,
                               CGAL::parameters::vertex_corner_map(vmap_in),
                               CGAL::parameters::vertex_corner_map(vmap_out));
    check_corner_ids();
    //------------------------------------------------------
    std::cout << "  no simplification face+vertex maps\n";
    out.clear_without_removing_property_maps();
    in.clear_without_removing_property_maps();
    std::ifstream(CGAL::data_file_path("meshes/sphere.off")) >> in;
    assert(vertices(in).size()!=0);
    PMP::remesh_planar_patches(in, out,
                               CGAL::parameters::face_patch_map(fmap_in).vertex_corner_map(vmap_in),
                               CGAL::parameters::face_patch_map(fmap_out).vertex_corner_map(vmap_out));
    check_corner_ids();
    std::map<std::array<Point_3, 3>, std::size_t> pids_map;
    auto get_pt_sorted_array = [](Surface_mesh::Face_index f, Surface_mesh& sm)
    {
      auto h = halfedge(f, sm);
      std::array<Point_3, 3> pts = CGAL::make_array(sm.point(source(h, sm)),
                                                    sm.point(target(h, sm)),
                                                    sm.point(target(next(h, sm), sm)));
      std::sort(pts.begin(), pts.end());
      return pts;
    };
    for (auto f : faces(out)) pids_map[get_pt_sorted_array(f, out)]=get(fmap_out, f);
    for (auto f : faces(in)) assert(pids_map[get_pt_sorted_array(f, in)]==get(fmap_in, f));
  }

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
    assert(CGAL::is_valid_polygon_mesh(sm));
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
    assert(CGAL::is_valid_polygon_mesh(meshes[i]));

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
    assert(CGAL::is_valid_polygon_mesh(sm));
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
    assert(CGAL::is_valid_polygon_mesh(sm));
  }
#endif
  // Approximation then
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere.off using approximate predicates\n";
    std::ifstream in("data/decimation/sphere.off");
    in >> sm;
    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out, CGAL::parameters::cosine_of_maxium_angle(-0.99));
  }
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere_selection.off using approximate predicates\n";
    std::ifstream in("data/decimation/sphere_selection.off");
    in >> sm;
    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out, CGAL::parameters::cosine_of_maxium_angle(-0.99));
    std::ofstream out("sphere_selection_app.off");
    out << sm_out;
    std::cout << "output written to sphere_selection_app.off\n";
    assert(CGAL::is_valid_polygon_mesh(sm_out));
  }

  assert(OK);

  return 0 ;
}
