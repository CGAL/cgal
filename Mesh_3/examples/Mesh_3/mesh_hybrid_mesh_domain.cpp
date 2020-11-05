// #define CGAL_MESH_3_PROTECTION_DEBUG 1
// #define CGAL_MESH_3_VERBOSE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

// Ouput
#include <CGAL/Mesh_3/Dump_c3t3.h>

// Read 1D features from input file
#include "read_polylines.h"

// Sphere Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Labeled_mesh_domain_3<K> Implicit_domain;

// Polyhedral Domain
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Polyhedron_domain;

class Hybrid_domain
{
  const Implicit_domain& implicit_domain;
  const Polyhedron_domain& polyhedron_domain;

public:
  Hybrid_domain(const Implicit_domain& implicit_domain,
                const Polyhedron_domain& polyhedron_domain)
    : implicit_domain(implicit_domain)
    , polyhedron_domain(polyhedron_domain)
  {}

  // types required by the `MeshDomain_3` concept
  typedef int Surface_patch_index;
  typedef int Subdomain_index;
  typedef int Index;

  typedef K R;
  typedef K::Point_3 Point_3;
  typedef std::tuple<Point_3, Index, int> Intersection;

  CGAL::Bbox_3 bbox() const {
    return implicit_domain.bbox() + polyhedron_domain.bbox();
  }

  struct Construct_initial_points
  {
    Construct_initial_points(const Hybrid_domain& domain)
      : r_domain_(domain) {}

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 20) const
    {
      //construct initial points on implicit domain
      typedef Implicit_domain::Index Implicit_Index;
      std::vector<std::pair<Point_3,
                            Implicit_Index> > implicit_points_vector;
      Implicit_domain::Construct_initial_points cstr_implicit_initial_points =
        r_domain_.implicit_domain.construct_initial_points_object();
      cstr_implicit_initial_points(std::back_inserter(implicit_points_vector),
                                   n/2);
      for(std::size_t i = 0, end = implicit_points_vector.size(); i<end; ++i) {
        *pts++ = std::make_pair(implicit_points_vector[i].first, 2);
      }

      //construct initial points on polyhedral domain
      typedef Polyhedron_domain::Index Polyhedron_Index;
      std::vector<std::pair<Point_3,
                            Polyhedron_Index> > polyhedron_points_vector;
      Polyhedron_domain::Construct_initial_points cstr_polyhedron_initial_points =
        r_domain_.polyhedron_domain.construct_initial_points_object();
      cstr_polyhedron_initial_points(std::back_inserter(polyhedron_points_vector),
                                     n/2);
      for(std::size_t i = 0, end = polyhedron_points_vector.size(); i<end; ++i) {
        *pts++ = std::make_pair(polyhedron_points_vector[i].first, 1);
      }
      return pts;
    }

  private:
    const Hybrid_domain& r_domain_;
  }; // end Construct_initial_points_object

  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  struct Is_in_domain
  {
    Is_in_domain(const Hybrid_domain& domain) : r_domain_(domain) {}

    boost::optional<Subdomain_index> operator()(const K::Point_3& p) const
    {
      boost::optional<Subdomain_index> subdomain_index =
        r_domain_.implicit_domain.is_in_domain_object()(p);
      if(subdomain_index) return 2;
      else return r_domain_.polyhedron_domain.is_in_domain_object()(p);
    }
  private:
    const Hybrid_domain& r_domain_;
  };

  Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }

  struct Construct_intersection
  {
    Construct_intersection(const Hybrid_domain& domain)
      : r_domain_(domain) {}

    template <typename Query>
    Intersection operator()(const Query& query) const
    {
      using boost::get;
      //intersection with implicit domain
      Implicit_domain::Intersection implicit_inter =
        r_domain_.implicit_domain.construct_intersection_object()(query);
      //if found, return it
      if(get<2>(implicit_inter) != 0) {
        return Intersection(get<0>(implicit_inter), 2, get<2>(implicit_inter));
      }

      //intersection with polyhedral domain
      Polyhedron_domain::Intersection polyhedron_inter =
        r_domain_.polyhedron_domain.construct_intersection_object()(query);
      //if found, return it
      if(get<2>(polyhedron_inter) != 0) {
        const Point_3 inter_point = get<0>(polyhedron_inter);
        if(!r_domain_.implicit_domain.is_in_domain_object()(inter_point)) {
          return Intersection(inter_point, 1, get<2>(polyhedron_inter));
        }
      }
      //no intersection found
      return Intersection();
    }
  private:
    const Hybrid_domain& r_domain_;
  }; // end Construct_intersection

  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }

  //Index types converters
  Index index_from_surface_patch_index(const Surface_patch_index& index) const
  { return index; }

  Index index_from_subdomain_index(const Subdomain_index& index) const
  { return index; }

  Surface_patch_index surface_patch_index(const Index& index) const
  { return index; }

  Subdomain_index subdomain_index(const Index& index) const
  { return index; }
}; // end class Hybrid_domain

typedef CGAL::Mesh_domain_with_polyline_features_3<Hybrid_domain> Domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Domain, K>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>     Mesh_criteria;
typedef Mesh_criteria::Edge_criteria  Edge_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria  Cell_criteria;

// Function
typedef K::FT FT;
typedef K::Point_3 Point;

FT sphere_centered_at_111 (const Point& p)
{
  const FT dx=p.x()-1;
  const FT dy=p.y()-1;
  const FT dz=p.z()-1;

  return dx*dx+dy*dy+dz*dz-1;
}

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  const char* fname = "data/cube.off";
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(fname);
  input >> polyhedron;
  if(input.bad()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  // Polyhedra domain
  Polyhedron_domain polyhedron_domain(polyhedron);

  // Implicit domain
  // - the first argument is the function pointer
  // - the second argument is a bounding sphere of the domain
  // (Warning: Sphere_3 constructor uses square radius !)
  Implicit_domain sphere_domain =
    Implicit_domain::create_implicit_mesh_domain(sphere_centered_at_111,
                                                 K::Sphere_3(K::Point_3(1, 1, 1),
                                                             2.));

  Domain domain(sphere_domain, polyhedron_domain);

  // Read polyline features
  const char* lines_fname = "data/hybrid_example.polylines.txt";
  std::vector<std::vector<Point> > featured_curves;
  if (!read_polylines(lines_fname, featured_curves))
  { // see file "read_polylines.h"
    std::cerr << "Error: Cannot read file " << lines_fname << std::endl;
    return EXIT_FAILURE;
  }

  // Add features for protection
  domain.add_features(featured_curves.begin(), featured_curves.end());

  // Criteria
  Edge_criteria edge_criteria(0.08);
  Facet_criteria facet_criteria(30, 0.08, 0.025); // angle, size, approximation
  Cell_criteria cell_criteria(2, 0.1); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

  // Mesh generation (without optimization)
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      no_perturb(), no_exude());

  // Output
  dump_c3t3(c3t3, "out");

  return EXIT_SUCCESS;
}
