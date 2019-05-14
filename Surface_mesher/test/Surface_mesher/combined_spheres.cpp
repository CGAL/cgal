#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Multi_surface_3.h>

#include <CGAL/Surface_mesher/Standard_criteria.h>
#include <CGAL/Surface_mesher/Vertices_on_the_same_surface_criterion.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <CGAL/Surface_mesher/Sphere_oracle_3.h>
#include <CGAL/Surface_mesher/Point_surface_indices_oracle_visitor.h>

#include <CGAL/Point_traits.h>
#include <CGAL/Point_with_surface_index_geom_traits.h>

#include <iostream>
#include <fstream>
#include <string>

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Point_with_surface_index_geom_traits<K> My_traits;
// Multi_surface_traits<Regular_traits> ?
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<My_traits> Vb;
typedef CGAL::Surface_mesh_cell_base_3<My_traits> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<My_traits, Tds> Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef My_traits::Point_3 Point_3;
typedef My_traits::Sphere_3 Sphere_3;
typedef My_traits::FT FT;

using CGAL::Surface_mesher::Refine_criterion;
using CGAL::Surface_mesher::Standard_criteria;
using CGAL::Surface_mesher::Point_surface_indices_visitor;
typedef Refine_criterion<Tr> Criterion;
// changer le nom du Refine_criterion en Criterion_virtual_base
typedef Standard_criteria <Criterion> Multi_criterion;
// changer son nom en Multi_criterion, ou Combined_criterion

typedef Point_surface_indices_visitor Set_indices;

struct Point_with_surface_index_creator
{
  typedef FT    argument_type;
  typedef FT    argument1_type;
  typedef FT    argument2_type;
  typedef FT    argument3_type;
  typedef My_traits::Point_3 result_type;

  result_type operator()(const FT& x, const FT& y, const FT& z) const
  {
    return result_type(K::Point_3(x, y, z));
  }
};

typedef CGAL::Surface_mesher::Sphere_oracle_3<
  My_traits,
  Point_with_surface_index_creator,
  Set_indices> Single_oracle;

using CGAL::Surface_mesher::Combining_oracle;

typedef Combining_oracle<Single_oracle, Single_oracle> Oracle_2;
typedef Combining_oracle<Oracle_2, Single_oracle> Oracle_3;
typedef Combining_oracle<Oracle_3, Single_oracle> Oracle_4;
typedef Combining_oracle<Oracle_4, Single_oracle> Oracle_5;
typedef Oracle_5 Oracle;

#include <vector>

int main(int, char**)
{
  /*** Spheres radiuss ***/
  FT r1; // 93 milimeters
  FT r2;
  FT r3;
  FT r4;
  FT r5;
  std::vector<double> size_bounds(5);
  std::vector<double> radii(5);
  
  std::cout << "Input r1, r2, r3, r4, r5:" << std::endl;
  std::cin >> r1 >> r2 >> r3 >> r4 >> r5;
  std::cout << "Input the corresponding 5 size bounds:" << std::endl;
  std::cin >> size_bounds[0]
           >> size_bounds[1]
           >> size_bounds[2]
           >> size_bounds[3]
           >> size_bounds[4];
  if(!std::cin)
    return EXIT_FAILURE;

  radii[0] = CGAL::to_double(r1);
  radii[1] = CGAL::to_double(r2);
  radii[2] = CGAL::to_double(r3);
  radii[3] = CGAL::to_double(r4);
  radii[4] = CGAL::to_double(r5);

  const int number_of_initial_points = 20;
  
  const double facets_uniform_size_bound = 0.5; // mm
  const double facets_aspect_ratio_bound = 30; // degres

  Sphere_3 sphere1(CGAL::ORIGIN, r1*r1);
  Sphere_3 sphere2(CGAL::ORIGIN, r2*r2);
  Sphere_3 sphere3(CGAL::ORIGIN, r3*r3);
  Sphere_3 sphere4(CGAL::ORIGIN, r4*r4);
  Sphere_3 sphere5(CGAL::ORIGIN, r5*r5);

//   const Sphere_3 bounding_sphere(CGAL::ORIGIN, 
//                                  bounding_sphere_radius*bounding_sphere_radius);

//   Implicit_sphere sphere1(Sphere(r1), bounding_sphere, precision);
//   Implicit_sphere sphere2(Sphere(r2), bounding_sphere, precision);
//   Implicit_sphere sphere3(Sphere(r3), bounding_sphere, precision);
//   Implicit_sphere sphere4(Sphere(r4), bounding_sphere, precision);
//   Implicit_sphere sphere5(Sphere(r5), bounding_sphere, precision);

  typedef CGAL::Multi_surface_3<Sphere_3, Sphere_3> Surface_2;
  typedef CGAL::Multi_surface_3<Surface_2, Sphere_3> Surface_3;
  typedef CGAL::Multi_surface_3<Surface_3, Sphere_3> Surface_4;
  typedef CGAL::Multi_surface_3<Surface_4, Sphere_3> Surface_5;
  typedef Surface_5 Surface;

//   typedef CGAL::Multi_surface_3<Implicit_sphere, Implicit_sphere> Surface_2;
//   typedef CGAL::Multi_surface_3<Surface_2, Implicit_sphere> Surface_3;
//   typedef CGAL::Multi_surface_3<Surface_3, Implicit_sphere> Surface_4;
//   typedef CGAL::Multi_surface_3<Surface_4, Implicit_sphere> Surface_5;
//   typedef Surface_5 Surface;

  Surface_2 surface_2(sphere1, sphere2);
  Surface_3 surface_3(surface_2, sphere3);
  Surface_4 surface_4(surface_3, sphere4);
  Surface surface(surface_4, sphere5);

  Tr tr;
  C2t3 c2t3(tr);

//   Single_oracle 
//     single_oracle_1 (sphere1,
//                      K::Point_3(CGAL::ORIGIN), // center of the bounding sphere
//                      bounding_sphere_radius,   // its radius (in mm)
//                      precision,
//                      use_bipolar_oracle,       // bipolar oracle
//                      false,                    // debug off
//                      Set_indices(1));
  Single_oracle single_oracle_1(Set_indices(1));
  Single_oracle single_oracle_2(Set_indices(2));
  Single_oracle single_oracle_3(Set_indices(3));
  Single_oracle single_oracle_4(Set_indices(4));
  Single_oracle single_oracle_5(Set_indices(5));

  Oracle_2 oracle_2(single_oracle_1, single_oracle_2);
  Oracle_3 oracle_3(oracle_2, single_oracle_3);
  Oracle_4 oracle_4(oracle_3, single_oracle_4);
  Oracle_5 oracle(oracle_4, single_oracle_5);

  CGAL::Surface_mesher::Uniform_size_criterion<Tr>
    uniform_size_criterion (facets_uniform_size_bound); 
  CGAL::Surface_mesher::Aspect_ratio_criterion<Tr>
    aspect_ratio_criterion (facets_aspect_ratio_bound);
  CGAL::Surface_mesher::Vertices_on_the_same_surface_criterion<Tr>
    vertices_on_the_same_surface_criterion;

  std::vector<Criterion*> criterion_vector;
  criterion_vector.push_back(&uniform_size_criterion);
  criterion_vector.push_back(&vertices_on_the_same_surface_criterion);
  criterion_vector.push_back(&aspect_ratio_criterion);
  Multi_criterion multi_criterion (criterion_vector);

  oracle.construct_initial_points_object()(surface,
                                           CGAL::inserter(tr),
                                           number_of_initial_points);


  CGAL::make_surface_mesh(c2t3, surface, oracle, multi_criterion,
                          CGAL::Non_manifold_tag());

  std::string filename;
  std::cout << "Ouput file name (without extension):" << std::endl;
  std::cin >> filename;

  std::ofstream out_cgal((filename+".off").c_str());

  CGAL::output_surface_facets_to_off (out_cgal, c2t3);
}
