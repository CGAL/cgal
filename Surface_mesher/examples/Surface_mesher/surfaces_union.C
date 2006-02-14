#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>
#include <CGAL/Complex_2_in_triangulation_surface_mesh_cell_base_3.h>
#include <CGAL/Mesh_3/Complex_2_in_triangulation_cell_base_3.h>
#include <CGAL/Implicit_surfaces_mesher_3.h>

#include <CGAL/Surface_mesher/Surface_mesher.h>

#include <CGAL/Surface_mesher/Criteria/Standard_criteria.h>
#include <CGAL/Surface_mesher/Criteria/Vertices_on_the_same_surface_criterion.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Surface_mesher/Oracles/Implicit_oracle.h>
#include <CGAL/Surface_mesher/Oracles/Polyhedral.h>
#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Surface_mesher/Oracles/Combining_oracle.h>

#include <CGAL/Point_with_surface_index_geom_traits.h>
#include <CGAL/Surface_mesher/Oracles/Point_surface_indices_visitor.h>

#include <iostream>
#include <fstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/File_medit.h>

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Robust_circumcenter_traits_3<K>  K2;
typedef CGAL::Point_with_surface_index_geom_traits<K2> My_traits;
typedef CGAL::Triangulation_vertex_base_3<My_traits> Vb1;
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<My_traits, Vb1> Vb;
typedef CGAL::Triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Complex_2_in_triangulation_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Complex_2_in_triangulation_surface_mesh_cell_base_3<My_traits, Cb2> Cb3;
typedef CGAL::Mesh_3::Complex_2_in_triangulation_cell_base_3<My_traits, Cb3> Cb;

typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<My_traits, Tds> Tr;

typedef CGAL::Complex_2_in_triangulation_3_surface_mesh<Tr> C2t3;

typedef K::FT FT;

using CGAL::Surface_mesher::Implicit_oracle;
using CGAL::Surface_mesher::Polyhedral;
using CGAL::Surface_mesher::Combining_oracle;
using CGAL::Surface_mesher::Refine_criterion;
using CGAL::Surface_mesher::Standard_criteria;
using CGAL::Surface_mesher::Point_surface_indices_visitor;
typedef Refine_criterion<Tr> Criterion;
// changer le nom du Refine_criterion en Criterion_virtual_base
typedef Standard_criteria <Criterion> Multi_criterion;
// changer son nom en Multi_criterion, ou Combined_criterion

typedef Point_surface_indices_visitor<Tr> Set_indices;

class Sphere {
private:
  const FT r; 
public:
  Sphere(FT rayon) : r(rayon) {}

  FT operator()(const FT x, const FT y, const FT z) const
  {
    return x*x+y*y+z*z-r*r;
  }
};

typedef CGAL::Creator_uniform_3<
  FT,
  K::Point_3> Point_creator_for_implicit_oracle;

typedef Implicit_oracle<My_traits,
                        Sphere,
                        Set_indices, // visitor that sets indices of points
                        Point_creator_for_implicit_oracle // to create
                                                          // points from
                                                          // three FT
                        > Implicite_sphere;
typedef Polyhedral<Tr, Set_indices> Polyhedron;
typedef Combining_oracle<Implicite_sphere, Polyhedron> Union_oracle;

using CGAL::Surface_mesher::Surface_mesher;
typedef Surface_mesher<Tr, Union_oracle, Multi_criterion> Surface_mesher_t;

int main(int, char**)
{
  /*** Sphere radius ***/
  const FT r = 0.6;
  const FT precision = 0.1; // mm
  const int number_of_initial_points = 50;
  const FT bounding_sphere_radius = 5.;
  const bool use_bipolar_oracle = true;

  const double facets_uniform_size_bound = 0.1;
  const double facets_aspect_ratio_bound = 0; // degres

  Sphere sphere(r);
  Tr tr;
  Implicite_sphere implicite_sphere(sphere,                                    
                                    K::Point_3(CGAL::ORIGIN),
                                    // center of the bounding sphere
                                    bounding_sphere_radius,
                                    // its radius (in mm)
                                    precision,
                                    use_bipolar_oracle,
                                    false, // debug off
                                    Set_indices(1));

  std::ifstream ifs("inputs/geosphere.off");
  Polyhedron polyhedron(ifs, Set_indices(2));
  ifs.close();

  Union_oracle union_oracle(implicite_sphere, polyhedron);

  Union_oracle::Points vector_of_initial_points =
    union_oracle.random_points(number_of_initial_points);

  for(Union_oracle::Points::iterator pit = vector_of_initial_points.begin();
      pit != vector_of_initial_points.end();
      ++pit)
  {
    tr.insert(*pit);
  }
  std::cerr << "Number of initial points, before refinement: "
            << tr.number_of_vertices() << std::endl;

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

  C2t3 c2t3(tr);

  {
    std::ofstream medit_before("sphere_union-before.mesh");
    CGAL::output_pslg_to_medit(medit_before, c2t3);
  }

  Surface_mesher_t mesher(tr, c2t3, union_oracle, multi_criterion);

  mesher.refine_mesh(true);

  std::ofstream ofs("sphere_union.off");
  CGAL::output_surface_facets_to_off(ofs, tr);

  std::ofstream medit("sphere_union.mesh");
  CGAL::output_pslg_to_medit(medit, c2t3);
}
