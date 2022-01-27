#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#ifdef USE_POLYHEDRON
#include <CGAL/Polyhedron_3.h>
#else
#include <CGAL/Surface_mesh.h>
#endif
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>

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
    if (!PMP::remesh_planar_patches(sm))
    {
      std::cerr << "ERROR: decimate cannot be done correctly\n";
      continue;
    }
    ss=std::stringstream();
    ss << "out" << i << ".off";
    std::ofstream out(ss.str().c_str());
    out << sm;
    std::cout << "  output written to out" << i << ".off\n";
  }

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
    if (!PMP::remesh_planar_patches(sm, 0.9801))
    {
      std::cerr << "ERROR: decimate cannot be done correctly\n";
      continue;
    }
    ss=std::stringstream();
    ss << "out_a" << i << ".off";
    std::ofstream out(ss.str().c_str());
    out << sm;
    std::cout << "  output written to out_a" << i << ".off\n";
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
    std::cerr << "ERROR: decimate cannot be done correctly\n";
  else
  {
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
  }

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
  if (!PMP::decimate_meshes_with_common_interfaces(meshes, 0.9801))
    std::cerr << "ERROR: decimate cannot be done correctly\n";
  else
  {
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
  }

#ifndef CGAL_DO_NOT_USE_PCA
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

    if (!PMP::decimate_with_pca_for_coplanarity(sm, 1e-5, 0.9801))
    {
      std::cerr << "ERROR: decimate cannot be done correctly\n";
      continue;
    }
    ss=std::stringstream();
    ss << "out_a_pca" << i << ".off";
    std::ofstream out(ss.str().c_str());
    out << sm;
    std::cout << "  output written to out_a_pca" << i << ".off\n";
  }

// testing decimation of meshes,  preserving common interface with almost coplanar/collinear tests using PCA
  std::cout << "decimate a range of meshes with common interfaces (approximate coplanar/collinear with PCA)\n";
  if (!PMP::decimate_meshes_with_common_interfaces_and_pca_for_coplanarity(meshes, 0.99, 0.9801))
    std::cerr << "ERROR: decimate cannot be done correctly\n";
  else
  {
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
  }

// two examples that fails with approximate but works with PCA
  //PCA first
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere.off using PCA\n";
    std::ifstream in("data/decimation/sphere.off");
    in >> sm;
    if (!PMP::decimate_with_pca_for_coplanarity(sm,1e-5,0.9801))
      std::cout << "ERROR: decimate cannot be done correctly\n";
    else
    {
      std::ofstream out("sphere_pca.off");
      out << sm;
      std::cout << "output written to sphere_pca.off\n";
    }
  }
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere_selection.off using PCA\n";
    std::ifstream in("data/decimation/sphere_selection.off");
    in >> sm;
    if (!PMP::decimate_with_pca_for_coplanarity(sm,1e-5,0.9801))
      std::cout << "ERROR: decimate cannot be done correctly\n";
    else
    {
      std::ofstream out("sphere_selection_pca.off");
      out << sm;
      std::cout << "output written to sphere_selection_pca.off\n";
    }
  }
#endif
  // Approximation then
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere.off using approximate predicates\n";
    std::ifstream in("data/decimation/sphere.off");
    in >> sm;
    if (!PMP::remesh_planar_patches(sm,0.9801))
      std::cout << "ERROR: decimate cannot be done correctly (this is the expected behavior)\n";
    else
    {
      std::ofstream out("sphere_app.off");
      out << sm;
      std::cout << "output written to sphere_app.off\n";
    }
  }
  {
    Surface_mesh sm;
    std::cout << "decimate of data/decimation/sphere_selection.off using approximate predicates\n";
    std::ifstream in("data/decimation/sphere_selection.off");
    in >> sm;
    if (!PMP::remesh_planar_patches(sm,0.9801))
      std::cout << "ERROR: decimate cannot be done correctly (this is the expected behavior)\n";
    else
    {
      std::ofstream out("sphere_selection_app.off");
      out << sm;
      std::cout << "output written to sphere_selection_app.off\n";
    }
  }

  return 0 ;
}
