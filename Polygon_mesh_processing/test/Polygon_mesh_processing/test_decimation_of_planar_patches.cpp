#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#ifdef USE_POLYHEDRON
#include <CGAL/Polyhedron_3.h>
#else
#include <CGAL/Surface_mesh.h>
#endif
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/subdivision_method_3.h>

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

void approximate_remeshing(Surface_mesh& sm, double cos_th, double frechet)
{
  std::vector<std::size_t> region_ids(num_faces(sm));
  std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices
  std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions
  boost::vector_property_map<CGAL::Epick::Vector_3> normal_map; // normal of the supporting planes of the regions detected

  // detect planar regions in the mesh
  std::size_t nb_regions =
    PMP::region_growing_of_planes_on_faces(sm,
                                           CGAL::make_random_access_property_map(region_ids),
                                           CGAL::parameters::cosine_of_maximum_angle(cos_th).
                                                             region_primitive_map(normal_map).
                                                             maximum_distance(frechet));

  // detect corner vertices on the boundary of planar regions
  std::size_t nb_corners =
    PMP::detect_corners_of_regions(sm,
                                   CGAL::make_random_access_property_map(region_ids),
                                   nb_regions,
                                   CGAL::make_random_access_property_map(corner_id_map),
                                   CGAL::parameters::cosine_of_maximum_angle(cos_th).
                                                     maximum_distance(frechet).
                                                     edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

  // run the remeshing algorithm using filled properties
  Surface_mesh out;
  PMP::remesh_almost_planar_patches(sm,
                                    sm,
                                    nb_regions, nb_corners,
                                    CGAL::make_random_access_property_map(region_ids),
                                    CGAL::make_random_access_property_map(corner_id_map),
                                    CGAL::make_random_access_property_map(ecm),
                                    CGAL::parameters::patch_normal_map(normal_map),
                                    CGAL::parameters::visitor([](Surface_mesh& sm){sm.clear_without_removing_property_maps ();}));

}

int main()
{
  CGAL::Aff_transformation_3<Kernel> rot( 0.914409, 0.0707007,-0.39857 , 0,
                                         -0.34672 , 0.644949 ,-0.681048, 0,
                                          0.208907, 0.760948 , 0.61426 , 0);

// testing decimate function
  bool OK = true;
  const int nb_meshes=4;

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
// test on cheese
{
    std::cout << "handling decimation of data/cheese.off\n";
    Surface_mesh sm;
    std::ifstream(CGAL::data_file_path("meshes/cheese.off")) >> sm;
    auto ecm = sm.add_property_map<Surface_mesh::Edge_index, bool>("ecm", false).first;
    PMP::detect_sharp_edges(sm, 60, ecm);
    PMP::isotropic_remeshing(faces(sm), 0.004, sm, CGAL::parameters::edge_is_constrained_map(ecm));

    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out);
    std::ofstream("cheese_out.off") << sm_out;
    assert(CGAL::is_valid_polygon_mesh(sm_out));

    PMP::transform(rot, sm);
    approximate_remeshing(sm, 0.98, 1e-2);
    std::ofstream("cheese_out_rg.off") << sm;
    assert(CGAL::is_valid_polygon_mesh(sm));
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
    std::cout << "handling decimation of transformed data/decimation/m" << i << ".off (approximate coplanar/collinear)\n";
    std::stringstream ss;
    ss << "data/decimation/m" << i << ".off";
    Surface_mesh sm;
    std::ifstream in(ss.str().c_str());
    in >> sm;
    PMP::transform(rot, sm);

    // call the decimation function
    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out, CGAL::parameters::cosine_of_maximum_angle(0.99));
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

#if 0 // tests to be re-enable when using region growing with common interface
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
#endif
// testing decimate function with almost coplanar/collinear tests using RG
  for (int i=1; i<=nb_meshes; ++i)
  {
    std::cout << "handling decimation of transformed data/decimation/m" << i << ".off (approximate coplanar/collinear with RG)\n";
    std::stringstream ss;
    ss << "data/decimation/m" << i << ".off";
    Surface_mesh sm;
    std::ifstream in(ss.str().c_str());
    in >> sm;
    PMP::transform(rot, sm);

    // call the decimation function

    approximate_remeshing(sm, 0.98, 1e-2);
    ss=std::stringstream();
    ss << "out_a_rg" << i << ".off";
    std::ofstream out(ss.str().c_str());
    out << sm;
    std::cout << "  output written to out_a_rg" << i << ".off\n";
    assert(CGAL::is_valid_polygon_mesh(sm));
  }

// two examples that fails with approximate but works with RG
  //PCA first
  {
    Surface_mesh sm;
    std::cout << "decimate of refined data/decimation/sphere.off using RG\n";
    std::ifstream(CGAL::data_file_path("meshes/sphere.off")) >> sm;
    CGAL::Subdivision_method_3::Sqrt3_subdivision(sm,
                                                  CGAL::parameters::number_of_iterations(3));

    approximate_remeshing(sm, 0.98, 1e-2);
    std::ofstream out("sphere_rg.off");
    out << sm;
    std::cout << "output written to sphere_rg.off\n";
    assert(CGAL::is_valid_polygon_mesh(sm));
  }
  // Approximation then
  {
    Surface_mesh sm;
    std::ifstream(CGAL::data_file_path("meshes/sphere.off")) >> sm;
    CGAL::Subdivision_method_3::Sqrt3_subdivision(sm,
                                                  CGAL::parameters::number_of_iterations(3));
    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out, CGAL::parameters::cosine_of_maximum_angle(0.99));
  }
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere_selection.off using approximate predicates\n";
    std::ifstream in("data/decimation/sphere_selection.off");
    in >> sm;
    Surface_mesh sm_out;
    PMP::remesh_planar_patches(sm, sm_out, CGAL::parameters::cosine_of_maximum_angle(0.99));
    std::ofstream out("sphere_selection_app.off");
    out << sm_out;
    std::cout << "output written to sphere_selection_app.off\n";
    assert(CGAL::is_valid_polygon_mesh(sm_out));
  }

  assert(OK);

  return 0 ;
}
