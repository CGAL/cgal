#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

#include "test_meshing_utilities.h"

#include <CGAL/Image_3.h>
#include <CGAL/Gray_image_mesh_domain_3.h>
#include <CGAL/use.h>

#include <functional>

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

template<typename T>
struct Greater_than {
  typedef T argument_type;
  Greater_than(const T& second) : second(second) {}
  bool operator()(const T& first) const {
    return std::greater<T>()(first, second);
  }
  T second;
};

template <typename Concurrency_tag = CGAL::Sequential_tag>
struct Image_tester : public Tester<K_e_i>
{
public:
  void image() const
  {
    typedef float                                       Image_word_type;
    typedef CGAL::Image_3                               Image;
    typedef CGAL::Gray_image_mesh_domain_3<
      Image,
      K_e_i,
      Image_word_type,
      Greater_than<double> >                            Mesh_domain;

    typedef typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type                            Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    typedef CGAL::Mesh_criteria_3<Tr>                   Mesh_criteria;

    CGAL_USE_TYPE(typename Mesh_domain::Surface_patch_index);

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Image image;
    if (!image.read("data/skull_2.9.inr"))
    {
      std::cout << "Image reading error. Exit test.\n";
      return;
    }

    std::cout << "\tSeed is\t"
              << CGAL::get_default_random().get_seed() << std::endl;

    // Domain
    Mesh_domain domain(image,
      2.9f, //isovalue
      0.f,  //value_outside
      1e-3, //error_bound
      &CGAL::get_default_random());//random generator for determinism

    // Mesh criteria
    Mesh_criteria criteria(facet_angle = 30,
                           facet_size = 6,
                           facet_distance = 2,
                           facet_topology = CGAL::MANIFOLD,
                           cell_radius_edge_ratio = 3,
                           cell_size = 8);

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        no_perturb(),
                                        no_exude(),
      mesh_3_options(number_of_initial_points = 30),
      non_manifold()
      );

    // Verify
    this->verify_c3t3_volume(c3t3, 1236086 * 0.95, 1236086 * 1.05);
    this->verify(c3t3, domain, criteria, Bissection_tag());
  }
};


int main()
{
  Image_tester<> test_epic;
  std::cerr << "Mesh generation from a 3D image:\n";
  test_epic.image();

#ifdef CGAL_LINKED_WITH_TBB
  Image_tester<CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation from a 3D image:\n";
  test_epic_p.image();
#endif

  return EXIT_SUCCESS;
}
