#include "test_meshing_utilities.h"

#include <CGAL/Labeled_mesh_domain_3.h>

#include <CGAL/disable_warnings.h>

template <typename K, typename Concurrency_tag = CGAL::Sequential_tag>
struct Implicit_tester : public Tester<K>
{
  typedef typename K::Point_3 Point;
  typedef typename K::FT FT;
  static FT sphere_function (const Point& p)
  {
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
    return x2+y2+z2-1;
  }

  void implicit() const
  {

    typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

    typedef typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      typename CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;

    typedef typename K::Sphere_3 Sphere_3;

    typedef typename Mesh_domain::Surface_patch_index Surface_patch_index;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;

    namespace p = CGAL::parameters;
    Mesh_domain domain =
      Mesh_domain::create_implicit_mesh_domain
      (p::function = Implicit_tester<K>::sphere_function,
       p::bounding_object = Sphere_3(CGAL::ORIGIN,2.),
       p::p_rng = &CGAL::get_default_random(),
       p::relative_error_bound = 1e-3);

    // Set mesh criteria
    Facet_criteria facet_criteria(0, 0, 0.3);
    Cell_criteria cell_criteria(0, 0.5);
    Mesh_criteria criteria(facet_criteria, cell_criteria);

    std::vector<Point> initial_points;
    initial_points.push_back(Point(1,0,0));
    initial_points.push_back(Point(0,1,0));
    initial_points.push_back(Point(0,0,1));
    initial_points.push_back(Point(-1,0,0));
    initial_points.push_back(Point(0,-1,0));
    initial_points.push_back(Point(0,0,-1));

    // Mesh generation
    C3t3 c3t3;
    c3t3.insert_surface_points(initial_points.begin(),
                               initial_points.end(),
                               domain.index_from_surface_patch_index(Surface_patch_index(0,1)));

    CGAL::refine_mesh_3(c3t3, domain, criteria,
                        CGAL::parameters::no_exude(),
                        CGAL::parameters::no_perturb());

    CGAL::remove_far_points_in_mesh_3(c3t3);

#ifdef CGAL_LINKED_WITH_TBB
    // Parallel
    if (std::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value)
    {
      this->verify(c3t3, domain, criteria, Bissection_tag(), 40, 65, 60, 110);
    }
    else
#endif //CGAL_LINKED_WITH_TBB
    {
      // Verify
      this->verify(c3t3, domain, criteria, Bissection_tag(), 50, 58, 80, 90);
    }
  }
};


int main()
{
  Implicit_tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from an implicit function:\n";
  test_epic.implicit();

#ifdef CGAL_LINKED_WITH_TBB
  Implicit_tester<K_e_i, CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation from an implicit function:\n";
  test_epic_p.implicit();
#endif
  return EXIT_SUCCESS;
}
