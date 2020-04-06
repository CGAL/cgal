#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <fstream>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>

using namespace CGAL::Surface_mesh_topology;

typedef CGAL::Surface_mesh<CGAL::Simple_cartesian<double>::Point_3> SM;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<2, 3> LCC_CM;
typedef CGAL::Linear_cell_complex_for_generalized_map<2, 3> LCC_GM;
typedef CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> Poly;

template<typename LCC>
void load_lcc(LCC& lcc, const std::string& filename)
{
  std::ifstream inp(filename);
  if (inp.fail())
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    exit(EXIT_FAILURE);
  }
  CGAL::load_off(lcc, inp);
}

template<typename SM>
void load_sm(SM& sm, const std::string& filename)
{
  std::ifstream inp(filename);
  if (inp.fail())
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    exit(EXIT_FAILURE);
  }
  inp>>sm;
}

template<typename Mesh>
bool test_one_data_structure(const Mesh& mesh, std::size_t nbedges, double length)
{
  bool res=true;
  
  CGAL::Surface_mesh_topology::Euclidean_length_weight_functor<Mesh> wf(mesh);
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<Mesh>      cst(mesh);

  std::cout<<"."<<std::flush;
  CGAL::Surface_mesh_topology::Path_on_surface<Mesh> cycle=cst.compute_edgewidth();
  if (cycle.length()!=nbedges)
  {
    std::cout<<"[ERROR]: number of edges for the cycle is not minimal ("<<cycle.length()
             <<">"<<nbedges<<")."<<std::endl;
    res=false;
  }
  
  std::cout<<"."<<std::flush;
  cycle=cst.compute_edgewidth(wf);
  double l=0.;
  for (int i=0; i<cycle.length(); ++i) { l+=wf(cycle[i]); }
  if (l>length)
  {
    std::cout<<"[ERROR] length of the cycle is not minimal. ("<<l
             <<">"<<length<<")."<<std::endl;
    res=false;
  }
  return res;
}


/// Compute edge width for the given mesh (loaded in filename), and test if
/// the size of the cycle is nbedges for unit weight and the length of the cycle
/// is smaller than length for Euclidean length weight.
bool test(const char* filename, std::size_t nbedges, double length)
{
  bool res=true;
  
  LCC_CM lcc_cm;
  load_lcc(lcc_cm, filename);
  if (!test_one_data_structure(lcc_cm, nbedges, length))
  {
    std::cout<<"[ERROR] for Linear_cell_complex_for_combinatorial_map."<<std::endl;
    res=false;
  }
  
  LCC_GM lcc_gm;
  load_lcc(lcc_gm, filename);
  if (!test_one_data_structure(lcc_cm, nbedges, length))
  {
    std::cout<<"[ERROR] for Linear_cell_complex_for_generalized_map."<<std::endl;
    res=false;
  }

  SM sm;
  load_sm(sm, filename);
  if (!test_one_data_structure(lcc_cm, nbedges, length))
  {
    std::cout<<"[ERROR] for Surface_mesh."<<std::endl;
    res=false;
  }

  Poly poly;
  load_sm(sm, filename);
  if (!test_one_data_structure(lcc_cm, nbedges, length))
  {
    std::cout<<"[ERROR] for Polyhedron_3."<<std::endl;
    res=false;
  }

  return res;
}

int main()
{
  bool res=true;
  std::cout<<"[BEGIN] Test shortest cycle non contractible "<<std::flush;

  if (!test("data/elephant-with-holes.off", 5, 0.0515047))
  {
    std::cout<<"[ERROR] for data/elephant-with-holes.off."<<std::endl;
    res=false;
  }
    
  if (!test("data/obj1.off", 4, 16.2814))
  { std::cout<<"[ERROR] for data/obj1.off."<<std::endl;
    res=false;
  }

  if (!test("data/obj1-with-holes.off", 4, 16.6082))
  {
    std::cout<<"[ERROR] for data/elephant-with-holes.off."<<std::endl;
    res=false;
  }

  if (!res)
  { return EXIT_FAILURE; }
  else
  { std::cout<<" All tests OK."<<std::endl; }
  
  return EXIT_SUCCESS;
}
