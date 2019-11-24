#include <CGAL/Polygonal_schema.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>

typedef CGAL::Surface_mesh_topology::Polygonal_schema_with_combinatorial_map<> PS;
typedef typename PS::Dart_handle  Dart_handle;

using namespace CGAL::Surface_mesh_topology;

unsigned int starting_seed;

///////////////////////////////////////////////////////////////////////////////
int myrandom (int i) { return std::rand()%i;}
//-----------------------------------------------------------------------------
void generate_random_map(PS& ps, std::size_t nb_labels,
                         std::size_t max_dart_per_face=30,
                         bool closed=true,
                         std::size_t percentage_of_perforated=20) // 0 => no perforated faces
{
  std::vector<std::string> all_labels(nb_labels*2);
  for (int i=0; i<nb_labels; ++i)
  {
    all_labels[2*i]=std::to_string(i+1);
    all_labels[(2*i)+1]=std::to_string(-(i+1));
  }

  std::random_shuffle (all_labels.begin(), all_labels.end(), myrandom);

  std::size_t endlabel=all_labels.size();
  if (!closed)
  { endlabel-=(all_labels.size()/10); } // We remove 10% of labels.

  for (std::size_t i=0; i<endlabel; )
  {
    ps.begin_facet();
    for (std::size_t j=0, nb=1+rand()%(max_dart_per_face);
         i<endlabel && j<=nb; ++i, ++j)
    {
      ps.add_edges_to_facet(all_labels[i]);
    }
    Dart_handle dh=ps.end_facet();

    if (rand()%101<percentage_of_perforated)
    { ps.perforate_facet(dh); }
  }

  ps.keep_biggest_connected_component();
}
///////////////////////////////////////////////////////////////////////////////
bool test_two_random_paths(const PS& ps,
                           const Curves_on_surface_topology<PS>& cst,
                           std::size_t nbtests,
                           std::size_t lmin=5,
                           std::size_t lmax=20)
{
  CGAL_assertion(lmin>0 && lmin<lmax);
  CGAL::Random random(starting_seed+nbtests);
  Path_on_surface<PS> p1(ps), p2(ps);
  internal::generate_random_closed_path(p1, random.get_int(lmin, lmax), random);

  p2=p1;
  p2.update_path_randomly(100, random);

  bool res=true;
  if (!cst.are_freely_homotopic(p1, p2))
  {
    res=false;
    std::cout<<"ERROR for path number "<<nbtests<<" (seed="<<random.get_seed()
             <<")"<<std::endl; res=false;
  }
  
  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool run_n_random_paths_tests(const PS& ps,
                              std::size_t n,
                              const Curves_on_surface_topology<PS>& cst,
                              std::size_t lmin=5,
                              std::size_t lmax=20)
{
  bool res=true;
  for (std::size_t i=0; i<n; ++i)
  {
    std::cout<<"."<<std::flush;
    if (!test_two_random_paths(ps, cst, i, lmin, lmax))
    { res=false; }
  }
  std::cout<<std::endl;
  return res;
}
///////////////////////////////////////////////////////////////////////////////
int main()
{
  starting_seed=unsigned(std::time(0));
  std::cout<<"Start tests (seed="<<starting_seed<<")"<<std::endl;
  srand(starting_seed);
  CGAL::Random r(starting_seed);

  PS ps1, ps2, ps3, ps4;
  generate_random_map(ps1, 20000, 30, true, 0);   // Closed, no perforated faces
  generate_random_map(ps2, 20000, 30, true, 15);  // Closed, 15% of perforated faces
  generate_random_map(ps3, 20000, 30, false, 0);  // Open, no perforated faces
  generate_random_map(ps4, 20000, 30, false, 15); // Open, 15% of perforated faces

  Curves_on_surface_topology<PS> cst1(ps1);
  Curves_on_surface_topology<PS> cst2(ps2);
  Curves_on_surface_topology<PS> cst3(ps3);
  Curves_on_surface_topology<PS> cst4(ps4);

  bool res=run_n_random_paths_tests(ps1, 20, cst1) &&
      run_n_random_paths_tests(ps2, 20, cst2)  &&
      run_n_random_paths_tests(ps3, 20, cst3)  &&
      run_n_random_paths_tests(ps4, 20, cst4);

  if (res)
  { std::cout<<"test_homotopy_with_polygonal_schema ALL TESTS OK."<<std::endl; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
