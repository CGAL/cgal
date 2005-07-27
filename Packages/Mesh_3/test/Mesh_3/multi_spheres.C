#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>
#include <CGAL/Complex_2_in_triangulation_surface_mesh_cell_base_3.h>
#include <CGAL/Mesh_3/Complex_2_in_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Implicit_surfaces_mesher_3.h>

#include <CGAL/Surface_mesher/Criteria/Standard_criteria.h>
#include <CGAL/Surface_mesher/Criteria/Vertices_on_the_same_surface_criterion.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Surface_mesher/Oracles/Multi_implicit_oracle.h>

#include <CGAL/Mesh_3/Slivers_exuder.h>

#include <CGAL/Point_traits.h>
#include <CGAL/Point_with_surface_index.h>

#include <iostream>
#include <fstream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/File_medit.h>

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Regular_traits;
typedef CGAL::Point_with_surface_index_geom_traits<Regular_traits> My_traits;
// Multi_surface_traits<Regular_traits> ?
typedef CGAL::Triangulation_vertex_base_with_info_3<bool, My_traits> Vb1;
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<My_traits, Vb1> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Complex_2_in_triangulation_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Complex_2_in_triangulation_surface_mesh_cell_base_3<My_traits, Cb2> Cb3;
typedef CGAL::Mesh_3::Complex_2_in_triangulation_cell_base_3<My_traits, Cb3> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<My_traits, Tds> Tr;

typedef K::FT FT;

// typedef FT surface_function(const FT, const FT, const FT); 
// // surface_function is a typedef of the function (FT,FT,FT)->FT
// // defining a surface.

using CGAL::Surface_mesher::Multi_implicit_oracle;
using CGAL::Surface_mesher::Refine_criterion;
using CGAL::Surface_mesher::Standard_criteria;

typedef Refine_criterion<Tr> Criterion;
// changer le nom du Refine_criterion en Criterion_virtual_base
typedef Standard_criteria <Criterion> Multi_criterion;
// changer son nom en Multi_criterion, ou Combined_criterion

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

typedef Multi_implicit_oracle<My_traits, Sphere> Oracle;

/*** criteria for Mesh_3 ***/
typedef CGAL::Mesh_criteria_3<Tr> Tets_criteria;

typedef CGAL::Implicit_surfaces_mesher_3<Tr,
					 Oracle,
                                         Multi_criterion,
                                         Tets_criteria> Mesher;

#include <vector>

int main(int argc, char** argv)
{
  /*** Spheres radiuss ***/
  const FT r1 = 93.; // 93 milimeters
  const FT r2 = 94.;
  const FT r3 = 97.;
  const FT r4 = 100.;
  const FT r5 = 267.;

  const FT precision = 0.1; // mm
  const int number_of_initial_points = 10;
  const FT bounding_sphere_radius = 500.;
  const bool use_bipolar_oracle = true;

  const int number_of_oracles = 5;

  const double facets_uniform_size_bound = 30.; // mm
  const double facets_aspect_ratio_bound = 0; // degres
  const double tets_radius_radius_ratio_bound = 2.5;
  const double tets_squared_size_bound = 0;



  Sphere sphere1(r1);
  Sphere sphere2(r2);
  Sphere sphere3(r3);
  Sphere sphere4(r4);
  Sphere sphere5(r5);

  std::vector<Sphere*> functions;
  functions.push_back(&sphere1);
  functions.push_back(&sphere2);
  functions.push_back(&sphere3);
  functions.push_back(&sphere4);
  functions.push_back(&sphere5);
  
  Tr tr;

  Oracle oracle1 (sphere1,
		  K::Point_3(CGAL::ORIGIN), // center of the bounding sphere
		  bounding_sphere_radius, // its radius (in mm)
		  precision,
		  use_bipolar_oracle); // bipolar oracle
  Oracle oracle2 (sphere2,
		  K::Point_3(CGAL::ORIGIN), // center of the bounding sphere
		  bounding_sphere_radius, // its radius (in mm)
		  precision,
		  use_bipolar_oracle);
  Oracle oracle3 (sphere3,
		  K::Point_3(CGAL::ORIGIN), // center of the bounding sphere
		  bounding_sphere_radius, // its radius (in mm)
		  precision,
		  use_bipolar_oracle);
  Oracle oracle4 (sphere4,
		  K::Point_3(CGAL::ORIGIN), // center of the bounding sphere
		  bounding_sphere_radius, // its radius (in mm)
		  precision,
		  use_bipolar_oracle);
  Oracle oracle5 (sphere5,
		  K::Point_3(CGAL::ORIGIN), // center of the bounding sphere
		  bounding_sphere_radius, // its radius (in mm)
		  precision,
		  use_bipolar_oracle);

  Oracle* oracles[number_of_oracles] = { &oracle1, 
					 &oracle2,
					 &oracle3,
					 &oracle4, 
					 &oracle5 };

  for(int i = 0; i < number_of_oracles; ++i)
  {
    Oracle::Points vector_of_initial_points =
      oracles[i]->random_points(number_of_initial_points);

    for(Oracle::Points::iterator pit = vector_of_initial_points.begin();
	pit != vector_of_initial_points.end();
	++pit)
    {
      pit->set_surface_index(i+1); // temporary
      tr.insert(*pit);
    }
  }

  Oracle multi_oracle (functions.begin(),
		       functions.end(),
		       K::Point_3(CGAL::ORIGIN), // center of bounding sphere
		       bounding_sphere_radius, // bounding radius
		       use_bipolar_oracle); // bipolar oracle

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

  Tets_criteria tets_criteria(tets_radius_radius_ratio_bound,
			      tets_squared_size_bound);

  Mesher mesher (tr, multi_oracle, multi_criterion, tets_criteria);
  mesher.refine_mesh();

  std::ofstream out("multi_spheres.mesh");
  CGAL::output_pslg_to_medit(out, tr);
}
