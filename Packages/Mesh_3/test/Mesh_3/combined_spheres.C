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

#include <CGAL/Surface_mesher/Oracles/Implicit_oracle.h>
#include <CGAL/Surface_mesher/Oracles/Combining_oracle.h>

#include <CGAL/Mesh_3/Slivers_exuder.h>

#include <CGAL/Point_traits.h>
#include <CGAL/Weighted_point_with_surface_index_geom_traits.h>
#include <CGAL/Surface_mesher/Oracles/Point_surface_indices_visitor.h>

#include <iostream>
#include <fstream>
#include <string>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/File_medit.h>

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Regular_traits;
typedef CGAL::Weighted_point_with_surface_index_geom_traits<Regular_traits> My_traits;
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

using CGAL::Surface_mesher::Implicit_oracle;
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

struct Point_with_surface_index_creator
{
  typedef FT    argument_type;
  typedef FT    argument1_type;
  typedef FT    argument2_type;
  typedef FT    argument3_type;
  typedef My_traits::Point result_type;
  typedef CGAL::Arity_tag<3> Arity;

  result_type operator()(const FT& x, const FT& y, const FT& z) const
  {
    return result_type(result_type::Point(x, y, z));
  }
};

typedef Implicit_oracle<My_traits,
                        Sphere,
                        Set_indices,
                        Point_with_surface_index_creator> Single_oracle;
typedef Combining_oracle<Single_oracle, Single_oracle> Oracle_2;
typedef Combining_oracle<Oracle_2, Single_oracle> Oracle_3;
typedef Combining_oracle<Oracle_3, Single_oracle> Oracle_4;
typedef Combining_oracle<Oracle_4, Single_oracle> Oracle_5;
typedef Oracle_5 Oracle;

/*** criteria for Mesh_3 ***/
typedef CGAL::Mesh_criteria_3<Tr> Tets_criteria;

class Special_tets_criteria
{
  double size_bound;
  double shape_bound;
  double special_size;
  double r;

public:
  Special_tets_criteria(const double special_size, const double r1,
                        const double radius_edge_bound = 2,
                        const double squared_radius_bound = 0)
    : size_bound(squared_radius_bound),
      shape_bound(radius_edge_bound),
      special_size(special_size),
      r(r1) {}

  typedef Tr::Cell_handle Cell_handle;
  typedef Tets_criteria::Quality Quality;

  class Is_bad
  {
  protected:
    const double shape_bound;
    const double size_bound;
    const double r1;
    const double squared_special_size;
  public:
    typedef Tr::Point Point_3;
      
    Is_bad(const double radius_edge_bound, 
	   const double squared_radius_bound,
           const double r1,
           const double size)
      : shape_bound(radius_edge_bound),
	size_bound(squared_radius_bound),
        r1(r1),
        squared_special_size(size*size){};
      
    bool operator()(const Cell_handle& c,
                    Quality& qual) const
    {
      const Point_3& p = c->vertex(0)->point();
      const Point_3& q = c->vertex(1)->point();
      const Point_3& r = c->vertex(2)->point();
      const Point_3& s = c->vertex(3)->point();

      typedef Tr::Geom_traits Geom_traits;
      typedef Geom_traits::Compute_squared_radius_3 Radius;
      typedef Geom_traits::Compute_squared_distance_3 Distance;
      typedef Geom_traits::Construct_circumcenter_3 Circumcenter;
      typedef Geom_traits::FT FT;

      Radius radius = Geom_traits().compute_squared_radius_3_object();
      Distance distance = Geom_traits().compute_squared_distance_3_object();
      Circumcenter circumcenter = 
        Geom_traits().construct_circumcenter_3_object();

      const double size = CGAL::to_double(radius(p, q, r, s));
      const double sq_distance_from_origin = 
        CGAL::to_double(distance(CGAL::ORIGIN, circumcenter(p, q, r, s)));

      double min_size_bound = size_bound;
      if( squared_special_size != 0 && 
          sq_distance_from_origin < r1*r1 )
      {
        if( min_size_bound == 0 )
          min_size_bound = squared_special_size;
        else
          if( squared_special_size < min_size_bound )
            min_size_bound = squared_special_size;
      }

      if( min_size_bound != 0)
      {
        qual.second = size / min_size_bound;
        // normalized by size bound to deal
        // with size field
        if( qual.sq_size() > 1 )
        {
          qual.first = 1; // (do not compute aspect)
          return true;
        }
      }
      if( shape_bound == 0 )
      {
        qual = Quality(0,1);
        return false;
      }

      double min_sq_length = CGAL::to_double(distance(p, q));
      min_sq_length = CGAL::min(min_sq_length,
                                CGAL::to_double(distance(p, r)));
      min_sq_length = CGAL::min(min_sq_length,
                                CGAL::to_double(distance(p, s)));
      min_sq_length = CGAL::min(min_sq_length,
                                CGAL::to_double(distance(q, r)));
      min_sq_length = CGAL::min(min_sq_length,
                                CGAL::to_double(distance(q, s)));
      min_sq_length = CGAL::min(min_sq_length,
                                CGAL::to_double(distance(r, s)));

      qual.first = size / min_sq_length;

      return (qual.first > shape_bound);
    }

  }; // end Is_bad


  Is_bad is_bad_object() const
  { return Is_bad(shape_bound, size_bound, r, special_size); }

}; // end Special_tets_criteria

typedef CGAL::Implicit_surfaces_mesher_3<Tr,
					 Oracle,
                                         Multi_criterion,
                                         Special_tets_criteria> Mesher;

#include <vector>

int main(int, char**)
{
  /*** Spheres radiuss ***/
  FT r1; // 93 milimeters
  FT r2;
  FT r3;
  FT r4;
  FT r5;
  double special_size_bound;

  std::cout << "Input r1, r2, r3, r4, r5:" << std::endl;
  std::cin >> r1 >> r2 >> r3 >> r4 >> r5;
  std::cout << "Input size bound for sphere <= r1:" << std::endl;
  std::cin >> special_size_bound;
  if(!cin)
    return EXIT_FAILURE;

  const FT precision = 0.1; // mm
  const int number_of_initial_points = 10;
  const FT bounding_sphere_radius = 500.;
  const bool use_bipolar_oracle = true;

  const double facets_uniform_size_bound = 30.; // mm
  const double facets_aspect_ratio_bound = 0; // degres
  const double tets_radius_radius_ratio_bound = 2.5;
  const double tets_squared_size_bound = 0;

  Sphere sphere1(r1);
  Sphere sphere2(r2);
  Sphere sphere3(r3);
  Sphere sphere4(r4);
  Sphere sphere5(r5);

  Tr tr;

  Single_oracle 
    single_oracle_1 (sphere1,
                     K::Point_3(CGAL::ORIGIN), // center of the bounding sphere
                     bounding_sphere_radius,   // its radius (in mm)
                     precision,
                     use_bipolar_oracle,       // bipolar oracle
                     false,                    // debug off
                     Set_indices(1));
  Single_oracle
    single_oracle_2 (sphere2,
                     K::Point_3(CGAL::ORIGIN),
                     bounding_sphere_radius,
                     precision,
                     use_bipolar_oracle,
                     false,
                     Set_indices(2));
  Single_oracle 
    single_oracle_3 (sphere3,
                     K::Point_3(CGAL::ORIGIN),
                     bounding_sphere_radius,
                     precision,
                     use_bipolar_oracle,
                     false,
                     Set_indices(3));
  Single_oracle 
    single_oracle_4 (sphere4,
                     K::Point_3(CGAL::ORIGIN),
                     bounding_sphere_radius,
                     precision,
                     use_bipolar_oracle,
                     false,
                     Set_indices(4));
  Single_oracle 
    single_oracle_5 (sphere5,
                     K::Point_3(CGAL::ORIGIN),
                     bounding_sphere_radius,
                     precision,
                     use_bipolar_oracle,
                     false,
                     Set_indices(5));

  Oracle_2 oracle_2(single_oracle_1, single_oracle_2);
  Oracle_3 oracle_3(oracle_2, single_oracle_3);
  Oracle_4 oracle_4(oracle_3, single_oracle_4);
  Oracle_5 oracle(oracle_4, single_oracle_5);

  Oracle::Points vector_of_initial_points =
    oracle.random_points(number_of_initial_points);

    for(Oracle::Points::iterator pit = vector_of_initial_points.begin();
	pit != vector_of_initial_points.end();
	++pit)
    {
      tr.insert(*pit);
    }

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

  Special_tets_criteria tets_criteria(tets_radius_radius_ratio_bound,
                                      tets_squared_size_bound,
                                      r1,
                                      special_size_bound);

  Mesher mesher (tr, oracle, multi_criterion, tets_criteria);
  mesher.refine_mesh();

  std::string filename;
  std::cout << "Input filename (default: combined_spheres.mesh):" << std::endl;
  std::cin >> filename;

  if(filename.empty())
    filename = "combined_spheres.mesh";

  std::ofstream out(filename.c_str());
  CGAL::output_pslg_to_medit(out, mesher.complex_2_in_triangulation_3());
}
