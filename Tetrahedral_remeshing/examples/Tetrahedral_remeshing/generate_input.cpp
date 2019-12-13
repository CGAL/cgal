#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>

#include <CGAL/Random.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K, int> T3;


bool load_binary_triangulation(std::istream& is, T3& t3)
{
  std::string s;
  if (!(is >> s)) return false;
  bool binary = (s == "binary");
  if (binary) {
    if (!(is >> s)) return false;
  }
  if (s != "CGAL" || !(is >> s) || s != "c3t3")
    return false;

  std::getline(is, s);
  if (binary) CGAL::set_binary_mode(is);
  is >> t3;
  return bool(is);
}

bool save_binary_triangulation(std::ostream& os, const T3& t3)
{
  typedef T3::Geom_traits::FT FT;
  os << "binary CGAL c3t3\n";
  CGAL::set_binary_mode(os);
  return !!(os << t3);
}

int main(int argc, char* argv[])
{
  const std::size_t nbv = 1000;

  int input_id = (argc > 1) ? atoi(argv[1]) : 1;
  char* filename;

  T3 tr;

  CGAL::Random rng;

  if (input_id == 1) //sphere and only one subdomain
  {
    filename = "data/triangulation_one_subdomain.binary.cgal";

    while (tr.number_of_vertices() < nbv)
      tr.insert(T3::Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));

    for (T3::Finite_cells_iterator cit = tr.finite_cells_begin();
      cit != tr.finite_cells_end(); ++cit)
    {
      cit->set_subdomain_index(1);
    }
  }
  else if (input_id == 2) //sphere separated in 2 subdomains by a plane
  {
    filename = "data/triangulation_two_subdomains.binary.cgal";

    while (tr.number_of_vertices() < nbv)
      tr.insert(
        T3::Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));

    const K::Plane_3 plane(K::Point_3(0,0,0), K::Point_3(0,1,0), K::Point_3(0,0,1));

    for (T3::Finite_cells_iterator cit = tr.finite_cells_begin();
         cit != tr.finite_cells_end(); ++cit)
    {
      if(plane.has_on_positive_side(
        CGAL::centroid(cit->vertex(0)->point(), cit->vertex(1)->point(),
                       cit->vertex(2)->point(), cit->vertex(3)->point())))
        cit->set_subdomain_index(1);
      else
        cit->set_subdomain_index(2);
    }
  }

  std::ofstream out(filename, std::ios_base::out | std::ios_base::binary);
  save_binary_triangulation(out, tr);

//  std::string file_in(filename);
//  std::string file_out = file_in.substr(0, file_in.find_first_of("."));
//  file_out.append(".mesh");
//  std::ofstream medit_out(file_out.c_str(), std::ios_base::out);
//  c3t3.output_to_medit(medit_out);

  return (!out.bad());
}
