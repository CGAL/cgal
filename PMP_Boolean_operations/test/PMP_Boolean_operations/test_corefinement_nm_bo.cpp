#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        SM;
typedef CGAL::Polyhedron_3<K>                                   Poly;

namespace PMP = CGAL::Polygon_mesh_processing;

template<class Mesh, class ID>
void test_one(std::string filename1, std::string filename2)
{
  std::cout << "Running test with " << filename1 << " " << filename2 << "\n";
  Mesh mesh1, mesh2;
  if(!PMP::IO::read_polygon_mesh(filename1, mesh1) || !PMP::IO::read_polygon_mesh(filename2, mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    exit(1);
  }


  PMP::Corefinement::Non_manifold_output_visitor<Mesh> visitor(mesh1, mesh2);

  std::array<Mesh,4> out_meshes;
  std::array<std::optional<Mesh*>, 4> output = {&out_meshes[0], &out_meshes[1], &out_meshes[2], &out_meshes[3]};

  std::array<bool,4> res = PMP::corefine_and_compute_boolean_operations(mesh1, mesh2, output, CGAL::parameters::visitor(visitor));

  std::array<std::string, 4> out_name_prefixes = {"union", "intersection", "tm1_minus_tm2", "tm2_minus_tm1"};

  for (int k=0; k<4; ++k)
  {
    std::cout << "  checking " << out_name_prefixes[k] << "\n";
    std::vector<K::Point_3> points;
    std::vector< std::array<ID, 3> > polygons;
    switch(k)
    {
      case 0: visitor.extract_union(points, polygons); break;
      case 1: visitor.extract_intersection(points, polygons); break;
      case 2: visitor.extract_tm1_minus_tm2(points, polygons); break;
      case 3: visitor.extract_tm2_minus_tm1(points, polygons); break;
    }
    CGAL::IO::write_polygon_soup(out_name_prefixes[k]+".off", points, polygons, CGAL::parameters::stream_precision(17));
    assert(res[k] == PMP::is_polygon_soup_a_polygon_mesh(polygons));
    if (res[k])
    {
      Mesh tmp;
      PMP::polygon_soup_to_polygon_mesh(points, polygons, tmp);
      std::size_t i=0;
      PMP::match_faces(tmp, *(output[k].value()), CGAL::Emptyset_iterator(),
                       CGAL::Counting_output_iterator(&i),
                       CGAL::Counting_output_iterator(&i));
      assert(i==0); // might be breaking with coplanar patches
    }
  }
}

void test(std::string filename1, std::string filename2)
{
  test_one<SM, std::size_t>(filename1, filename2);
  test_one<Poly, std::size_t>(filename1, filename2);
  test_one<SM, int>(filename1, filename2);
  test_one<Poly, int>(filename1, filename2);
  test_one<SM, unsigned int>(filename1, filename2);
  test_one<Poly, unsigned int>(filename1, filename2);
}

int main(int argc, char* argv[])
{
  if (argc<3)
  {
    test("data-coref/nm/tm1.off", "data-coref/nm/tm_union.off");
    test("data-coref/nm/tm1.off", "data-coref/nm/tm_inter.off");
    test("data-coref/nm/tm1.off", "data-coref/nm/tm_minus.off");
    test("data-coref/nm/tm_minus.off", "data-coref/nm/tm1.off");
    test("data-coref/nm/tm1.off", "data-coref/nm/tm_multi.off");
    test("data-coref/nm/tm_multi.off", "data-coref/nm/tm1.off");
    test("data-coref/sphere.off", "data-coref/elephant.off");
  }
  else
  {
    const std::string filename1 = argv[1];
    const std::string filename2 = argv[2];

    test(filename1, filename2);
  }

  return 0;
}
