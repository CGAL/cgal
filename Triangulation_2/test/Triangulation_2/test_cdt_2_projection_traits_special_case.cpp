#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>

#include <string>
#include <fstream>


typedef const double vec[3];
const vec input[10] =
{ { -0.37503900000000001, 0.5, 0.47710999999999998                     },
  { -0.37503900000000001, 0.5, 0.31415300000000002                     },
  { -0.37503900000000001, 0.5, 0.151196                                },
  { -0.37503900000000001, 0.10979194793565403, 0.151196                },
  { -0.37503899999999996, -0.049577392484619259, 0.151196              },
  { -0.37503900000000001, -0.18584700000000001, 0.151196               },
  { -0.37503900000000001, -0.062128898913781927, 0.20998675245268067   },
  { -0.37503900000000001, 0.11124288530418507, 0.29237289933619037     },
  { -0.37503900000000001, 0.20409758712826165, 0.33649738670770629     },
  { -0.37503900000000001, 0.29695228895233827, 0.38062187407922232     } };

const vec input_normal = { 1., 0., 0. };

template <typename Kernel, typename CDT_2_traits>
bool test(std::string test_name)
{
  std::cerr << "Testing " << test_name << std::endl;
  const unsigned nb_input_points = sizeof(input)/sizeof(vec);

  typedef CGAL::Constrained_triangulation_face_base_2<CDT_2_traits>     CDT_2_fb;
  typedef CGAL::Triangulation_vertex_base_2<CDT_2_traits>               CDT_2_vb;
  typedef CGAL::Triangulation_data_structure_2<CDT_2_vb, CDT_2_fb>      CDT_2_tds;
  typedef CGAL::No_constraint_intersection_requiring_constructions_tag  CDT_2_itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<CDT_2_traits,
                                                     CDT_2_tds,
                                                     CDT_2_itag>        CDT_2;

  typename Kernel::Vector_3 normal(input_normal[0],
                                   input_normal[1],
                                   input_normal[2]);
  CDT_2_traits traits(normal);
  CDT_2 cdt(traits);
  typedef typename CDT_2::Vertex_handle Vh;

  Vh first, previous;
  for(unsigned i = 0; i < nb_input_points; ++i)
  {
    typename Kernel::Point_3 p(input[i][0], input[i][1], input[i][2]);
    Vh current = cdt.insert(p);
    assert(current != previous);
    if(previous != Vh()) cdt.insert_constraint(previous, current);
    if(first == Vh()) first = current;
    previous = current;
  }
  assert(first != Vh());
  assert(previous != first);
  cdt.insert_constraint(previous, first);
  return cdt.is_valid();
}

int main()
{
  std::cerr.precision(17);
  typedef CGAL::Exact_predicates_exact_constructions_kernel   Epeck;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;

  bool ok = true;
  ok = ok &&
      test<CGAL::Epick,
           CGAL::Triangulation_2_projection_traits_3<Epick> >
      ("CDT_2 in a 3D plane, with Epick");
  ok = ok &&
      test<CGAL::Epeck,
           CGAL::Triangulation_2_projection_traits_3<Epeck> >
      ("CDT_2 in a 3D plane, with Epeck");
  return ok ? 0 : 1;
}
