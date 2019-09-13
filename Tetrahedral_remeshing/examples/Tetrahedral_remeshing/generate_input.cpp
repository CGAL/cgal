#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/IO/File_binary_mesh_3.h> 

#include <CGAL/Random.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K, int> T3;
typedef int Corner_index;
typedef int Curve_segment_index;

typedef CGAL::Mesh_complex_3_in_triangulation_3<T3, Corner_index, Curve_segment_index> C3t3;


int main(int argc, char* argv[])
{
  const std::size_t nbv = 1000;

  int input_id = (argc > 1) ? atoi(argv[1]) : 1;
  char* filename;

  T3 tr;
  C3t3 c3t3;
  c3t3.triangulation() = tr;

  CGAL::Random rng;

  if (input_id == 1) //sphere and only one subdomain
  {
    filename = "data/triangulation_one_subdomain.binary.cgal";

    while (tr.number_of_vertices() < nbv)
      tr.insert(T3::Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));

    for (T3::Finite_cells_iterator cit = tr.finite_cells_begin();
      cit != tr.finite_cells_end(); ++cit)
    {
      c3t3.add_to_complex(cit, 1);
    }
  }
  else if (input_id == 2) //sphere separated in 2 subdomains by a plane
  {
    filename = "data/triangulation_two_subdomains.binary.cgal";

    while (c3t3.triangulation().number_of_vertices() < nbv)
      c3t3.triangulation().insert(
        T3::Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));

    const K::Plane_3 plane(K::Point_3(0,0,0), K::Point_3(0,1,0), K::Point_3(0,0,1));

    for (T3::Finite_cells_iterator cit = c3t3.triangulation().finite_cells_begin();
         cit != c3t3.triangulation().finite_cells_end(); ++cit)
    {
      int index;
      if(plane.has_on_positive_side(
        CGAL::centroid(cit->vertex(0)->point(), cit->vertex(1)->point(),
                       cit->vertex(2)->point(), cit->vertex(3)->point())))
        index = 1;
      else
        index = 2;

      c3t3.add_to_complex(cit, index);
    }
  }

  std::ofstream out(filename, std::ios_base::out | std::ios_base::binary);
  CGAL::Mesh_3::save_binary_file(out, c3t3);

  std::string file_in(filename);
  std::string file_out = file_in.substr(0, file_in.find_first_of("."));
  file_out.append(".mesh");
  std::ofstream medit_out(file_out.c_str(), std::ios_base::out);
  c3t3.output_to_medit(medit_out);

  return (!out.bad());
}
