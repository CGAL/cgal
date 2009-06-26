#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>

// IO.h must be included before vertex and cell bases.
#include <CGAL/Mesh_3/IO.h>

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Multi_surface_3.h>

#include <CGAL/Volume_mesher_cell_base_3.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Implicit_surfaces_mesher_3.h>

#include <CGAL/Surface_mesher/Standard_criteria.h>
#include <CGAL/Surface_mesher/Vertices_on_the_same_surface_criterion.h>

#include <CGAL/Mesh_3/Slivers_exuder.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/File_medit.h>

// #include <CGAL/Surface_mesher/Sphere_oracle_3.h>
#include <CGAL/Surface_mesher/Point_surface_indices_oracle_visitor.h>

#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Point_traits.h>
#include <CGAL/Weighted_point_with_surface_index_geom_traits.h>

#include <iostream>
#include <fstream>
#include <string>

struct K2 : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Robust_circumcenter_traits_3<K2>  K;
// struct K : public CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt {};
typedef CGAL::Regular_triangulation_filtered_traits_3<K> Regular_traits;
typedef CGAL::Weighted_point_with_surface_index_geom_traits<Regular_traits> My_traits;
// Multi_surface_traits<Regular_traits> ?
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<My_traits> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Surface_mesh_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Volume_mesher_cell_base_3<My_traits, Cb2> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<My_traits, Tds> Tr;
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

class Sphere {
private:
  const FT r; 
public:
  Sphere(FT rayon) : r(rayon) {}

  FT operator()(const Point_3& p) const
  {
    const FT& x = p.x();
    const FT& y = p.y();
    const FT& z = p.z();
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
  //typedef CGAL::Arity_tag<3> Arity;

  result_type operator()(const FT& x, const FT& y, const FT& z) const
  {
    return result_type(result_type::Point(x, y, z));
  }
};

typedef CGAL::Implicit_surface_3<My_traits, Sphere> Implicit_sphere;

typedef CGAL::Surface_mesher::Implicit_surface_oracle_3<
  My_traits,
  Implicit_sphere,
  Point_with_surface_index_creator,
  Set_indices> Single_oracle;

// typedef CGAL::Surface_mesher::Sphere_oracle_3<
//   My_traits,
//   Point_with_surface_index_creator,
//   Set_indices> Single_oracle;

using CGAL::Surface_mesher::Combining_oracle;

typedef Combining_oracle<Single_oracle, Single_oracle> Oracle_2;
typedef Combining_oracle<Oracle_2, Single_oracle> Oracle_3;
typedef Combining_oracle<Oracle_3, Single_oracle> Oracle_4;
typedef Combining_oracle<Oracle_4, Single_oracle> Oracle_5;
typedef Oracle_5 Oracle;

/*** A size criterion for multi surfaces ***/
class Uniform_size_multi_surface_criterion : public Criterion {
public:
  typedef Criterion::Quality Quality;

private:
  std::vector<double> squared_sphere_radii;

  typedef Tr::Facet Facet;
  typedef Tr::Vertex_handle Vertex_handle;
  typedef Tr::Cell_handle Cell_handle;
  
public:
  Uniform_size_multi_surface_criterion(const std::vector<double> sphere_radii)
  {
    std::vector<double>::size_type size = sphere_radii.size();

    squared_sphere_radii.resize(size);

    for(std::vector<double>::size_type i = 0; i < size; ++i)
    {
      squared_sphere_radii[i] = sphere_radii[i] * sphere_radii[i];
    }
  }

  bool is_bad (const Facet& f, Quality& q) const {
    const Cell_handle& ch = f.first;
    const int i = f.second;
    const Vertex_handle& v1 = ch->vertex((i+1)&3);
    const Vertex_handle& v2 = ch->vertex((i+2)&3);
    const Vertex_handle& v3 = ch->vertex((i+3)&3);

    const int& number = v1->point().surface_index();
    if ( number == 0 ||
         (v2->point().surface_index() != number) ||
         (v3->point().surface_index() != number ) )
    {
      q = Quality(0);
      return true;
    }
    else
      if(squared_sphere_radii[number] == 0.)
      {
        q = Quality(0);
        return false;
      }
      else 
      {
        const Point_3& p1 = v1->point();

        q =  squared_sphere_radii[number] / 
          CGAL::squared_distance(p1,
                                 f.first->get_facet_surface_center(f.second));
        return q < 1;
      }
  }
}; // end Uniform_size_multi_surface_criterion

/*** criteria for Mesh_3 ***/
typedef CGAL::Mesh_criteria_3<Tr> Tets_criteria;

class Special_tets_criteria
{
  std::vector<double> squared_sphere_radii;
  std::vector<double> squared_radius_bounds;
  double squared_radius_edge_bound;

public:
  Special_tets_criteria(const std::vector<double> sphere_radii,
			const std::vector<double> radius_bounds,
                        const double radius_edge_bound = 2)
    : squared_sphere_radii(5),
      squared_radius_bounds(5),
      squared_radius_edge_bound(radius_edge_bound * radius_edge_bound)
  {
    for(int i = 0; i < 5; ++i)
    {
      squared_sphere_radii[i] = sphere_radii[i] * sphere_radii[i];
      squared_radius_bounds[i] = radius_bounds[i] * radius_bounds[i];
    }
  }

  typedef Tr::Cell_handle Cell_handle;
  typedef Tets_criteria::Quality Quality;

  class Is_bad {
  protected:
    std::vector<double> squared_sphere_radii;
    std::vector<double> squared_radius_bounds;
    double squared_radius_edge_bound;
  public:
    typedef Tr::Point Point_3;
      
    Is_bad(const std::vector<double> squared_sphere_radii_,
           const std::vector<double> squared_radius_bounds_,
           const double squared_radius_edge_bound_)
      : squared_sphere_radii(5),
        squared_radius_bounds(5),
        squared_radius_edge_bound(squared_radius_edge_bound_)
    {
      for(int i = 0; i < 5; ++i)
      {
        squared_sphere_radii[i] = squared_sphere_radii_[i];
        squared_radius_bounds[i] = squared_radius_bounds_[i];
      }
    }
      
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

      const double sq_distance_from_origin = 
        CGAL::to_double(distance(CGAL::ORIGIN, circumcenter(p, q, r, s)));

      double current_squared_radius_bound = 0;
      
      if(sq_distance_from_origin < squared_sphere_radii[0] )
      {
        current_squared_radius_bound = squared_radius_bounds[0];
      }
      else if(sq_distance_from_origin < squared_sphere_radii[1] )
      {
        current_squared_radius_bound = squared_radius_bounds[1];
      } 
      else if(sq_distance_from_origin < squared_sphere_radii[2] )
      {
        current_squared_radius_bound = squared_radius_bounds[2];
      } 
      else if(sq_distance_from_origin < squared_sphere_radii[3] )
      {
        current_squared_radius_bound = squared_radius_bounds[3];
      } 
      else 
      {
        current_squared_radius_bound = squared_radius_bounds[4];
      }

      const double sq_radius = CGAL::to_double(radius(p, q, r, s));

      if( current_squared_radius_bound != 0)
      {
        qual.second = sq_radius / current_squared_radius_bound;
        // normalized by size bound to deal
        // with size field
        if( qual.sq_size() > 1 )
        {
          qual.first = 1; // (do not compute aspect)
#ifdef CGAL_MESH_3_DEBUG_CRITERIA
          std::cerr << "bad tetrahedron: squared radius=" << sq_radius << "\n";
#endif
          return true;
        }
      }
      else
	qual.second = 1;
      if( squared_radius_edge_bound == 0 )
      {
        qual = Quality(0,1);
        return false;
      }

      double min_sq_length = CGAL::to_double(distance(p, q));
      min_sq_length = (CGAL::min)(min_sq_length,
				  CGAL::to_double(distance(p, r)));
      min_sq_length = (CGAL::min)(min_sq_length,
				  CGAL::to_double(distance(p, s)));
      min_sq_length = (CGAL::min)(min_sq_length,
				  CGAL::to_double(distance(q, r)));
      min_sq_length = (CGAL::min)(min_sq_length,
				  CGAL::to_double(distance(q, s)));
      min_sq_length = (CGAL::min)(min_sq_length,
				  CGAL::to_double(distance(r, s)));

      qual.first = sq_radius / min_sq_length;
#ifdef CGAL_MESH_3_DEBUG_CRITERIA
      if( qual.first > squared_radius_edge_bound )
        std::cerr << "bad tetrahedron: radius-edge ratio = " 
                  << CGAL::sqrt(qual.first)
                  << "\n";
#endif
      return (qual.first > squared_radius_edge_bound);
    }

  }; // end Is_bad


  Is_bad is_bad_object() const
  { return Is_bad(squared_sphere_radii, squared_radius_bounds,
		  squared_radius_edge_bound); }

}; // end Special_tets_criteria

#include <vector>

int main(int, char**)
{
  /*** Spheres radiuss ***/
  FT r1; // 93 milimeters
  FT r2;
  FT r3;
  FT r4;
  FT r5;
  std::vector<double> size_bounds(6);
  std::vector<double> radii(5);
  std::vector<double> facets_size_bounds(6);
  
  std::cout << "Input r1, r2, r3, r4, r5:" << std::endl;
  std::cin >> r1 >> r2 >> r3 >> r4 >> r5;
  std::cout << "Input the corresponding 5 size bounds:" << std::endl;
  std::cin >> size_bounds[0]
           >> size_bounds[1]
           >> size_bounds[2]
           >> size_bounds[3]
           >> size_bounds[4];

  size_bounds[5] = std::numeric_limits<double>::infinity(); // trick

  if(!std::cin)
    return EXIT_FAILURE;

  radii[0] = CGAL::to_double(r1);
  radii[1] = CGAL::to_double(r2);
  radii[2] = CGAL::to_double(r3);
  radii[3] = CGAL::to_double(r4);
  radii[4] = CGAL::to_double(r5);

  for(int i = 1; i < 6; ++i)
    facets_size_bounds[i] = (std::min)(size_bounds[i-1], size_bounds[i]);
  facets_size_bounds[0] = 0;

  const FT precision = 0.001;
  const int number_of_initial_points = 20;
  const FT bounding_sphere_radius = 2*r5;

  const double facets_aspect_ratio_bound = 0; // degres
  const double tets_radius_edge_bound = 2.5;

//   Sphere_3 sphere1(CGAL::ORIGIN, r1*r1);
//   Sphere_3 sphere2(CGAL::ORIGIN, r2*r2);
//   Sphere_3 sphere3(CGAL::ORIGIN, r3*r3);
//   Sphere_3 sphere4(CGAL::ORIGIN, r4*r4);
//   Sphere_3 sphere5(CGAL::ORIGIN, r5*r5);

  const Sphere_3 bounding_sphere(CGAL::ORIGIN, 
                                 bounding_sphere_radius*bounding_sphere_radius);

  Implicit_sphere sphere1(Sphere(r1), bounding_sphere, precision);
  Implicit_sphere sphere2(Sphere(r2), bounding_sphere, precision);
  Implicit_sphere sphere3(Sphere(r3), bounding_sphere, precision);
  Implicit_sphere sphere4(Sphere(r4), bounding_sphere, precision);
  Implicit_sphere sphere5(Sphere(r5), bounding_sphere, precision);

  typedef CGAL::Multi_surface_3<Implicit_sphere, Implicit_sphere> Surface_2;
  typedef CGAL::Multi_surface_3<Surface_2, Implicit_sphere> Surface_3;
  typedef CGAL::Multi_surface_3<Surface_3, Implicit_sphere> Surface_4;
  typedef CGAL::Multi_surface_3<Surface_4, Implicit_sphere> Surface_5;
  typedef Surface_5 Surface;

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

  CGAL::Surface_mesher::Aspect_ratio_criterion<Tr>
    aspect_ratio_criterion (facets_aspect_ratio_bound);
  Uniform_size_multi_surface_criterion 
    uniform_size_multi_surface_criterion(facets_size_bounds);

  std::vector<Criterion*> criterion_vector;
  criterion_vector.push_back(&uniform_size_multi_surface_criterion);
  criterion_vector.push_back(&aspect_ratio_criterion);
  Multi_criterion multi_criterion (criterion_vector);

  Special_tets_criteria tets_criteria(radii,
                                      size_bounds,
				      tets_radius_edge_bound);

  typedef CGAL::Implicit_surfaces_mesher_3<C2t3,
    Surface,
    Multi_criterion,
    Special_tets_criteria,
    Oracle> Mesher;

  oracle.construct_initial_points_object()(surface,
                                           CGAL::inserter(tr),
                                           number_of_initial_points);

  Mesher mesher (c2t3, surface, multi_criterion, tets_criteria, oracle);
  mesher.refine_mesh();

  CGAL::Mesh_3::Slivers_exuder<C2t3> exuder(tr);
  exuder.pump_vertices(0.2);

  std::string filename;
  std::cout << "Ouput file name (without extension):" << std::endl;
  std::cin >> filename;

  std::ofstream out((filename+".mesh").c_str());
  CGAL::output_to_medit(out, mesher.complex_2_in_triangulation_3());
  out.close();

  std::ofstream out_cgal((filename+".cgal").c_str());

  CGAL::Mesh_3::output_mesh(out_cgal,
                            mesher.complex_2_in_triangulation_3());
  std::cerr << "ok\n";
}
