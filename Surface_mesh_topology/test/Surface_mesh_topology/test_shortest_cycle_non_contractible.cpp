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

template<typename Map>
struct Marked_weight_functor
{
  using Weight_t=unsigned int;
  using Dart_const_descriptor=typename Map::Dart_const_descriptor;

  Marked_weight_functor(const Map& m, typename Map::size_type mark) : m_map(m), m_mark(mark)
  {}

  Weight_t operator() (Dart_const_descriptor dh) const
  {
    if (m_map.is_marked(dh, m_mark))
    { return 1; }

    return 100;
  }

protected:
  const Map& m_map;
  typename Map::size_type m_mark;
};

template<typename Mesh>
bool test_weighted(const Mesh& , const CGAL::Surface_mesh_topology::Path_on_surface<Mesh>& )
{ return true; }

template<>
bool test_weighted<LCC_CM>(const LCC_CM& map,
                           const CGAL::Surface_mesh_topology::Path_on_surface<LCC_CM>& cycle)
{
  bool res=true;

  // Cycle is the smallest non contractible cycle with unary weight on double-torus-2-b.off.
  // Its length is 12.
  // 1) We create a cycle parallel to this first one.
  typename LCC_CM::Dart_const_descriptor dh=nullptr;

  // We iterate through the cycle until a junction
  for (std::size_t i=0; dh==nullptr && i<cycle.length(); ++i)
  {
    if (map.opposite(map.next(map.opposite(map.next(cycle.get_ith_real_dart(i)))))!=
        cycle.get_ith_real_dart(i))
    { dh=map.next(cycle.get_ith_real_dart(i)); }
  }

  // We iterate through the border of the face until the next junction
  while(map.opposite(map.next(map.opposite(map.next(dh))))==dh)
  { dh=map.next(dh); }

  dh=map.next(dh);
  // 2) Here dh is  on the parallel of the first cycle. We mark darts of the cycle parallel
  //    to the first one. Its length is 24.
  auto mark=map.get_new_mark();
  std::size_t nbedges=0;
  typename LCC_CM::Dart_const_descriptor dh2=dh;
  do
  {
    map.mark(dh2, mark);
    ++nbedges;
    if (map.opposite(map.next(map.opposite(map.next(dh2))))==dh2)
    { dh2=map.next(dh2); }
    else
    { dh2=map.next(map.opposite(map.next(dh2))); }
  }
  while(dh2!=dh);

  if (nbedges!=24)
  {
    std::cout<<"[ERROR] in test_weighted: the length of the second cycle is wrong: "
             <<nbedges<<"!=24"<<std::endl;
    res=false;
  }

  Marked_weight_functor<LCC_CM> wf(map, mark);
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_CM> cst(map);

  auto cycle2=cst.compute_shortest_non_contractible_cycle(wf);
  if (cycle2.length()!=nbedges)
  {
    std::cout<<"[ERROR] in test_weighted for double-torus-2-b.off: the length"
             <<" of the weighted smallest non contractible cycle is wrong: "
             <<cycle2.length()<<"!="<<nbedges<<std::endl;
    res=false;
  }

  auto cycle3=cycle2;
  cycle3.reverse();
  if (!cst.are_freely_homotopic(cycle, cycle2) &&
      !cst.are_freely_homotopic(cycle, cycle3))
  {
    std::cout<<"[ERROR] in test_weighted for double-torus-2-b.off: "
             <<"cycle and cycle2 are not homotopic."<<std::endl;
    res=false;
  }

  return res;
}


template<typename Mesh>
bool test_one_data_structure(const Mesh& mesh, std::size_t nbedges, double length,
                             std::size_t nbfaces, bool is_double_torus_2_b)
{
  bool res=true;

  CGAL::Surface_mesh_topology::Euclidean_length_weight_functor<Mesh> wf(mesh);
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<Mesh>      cst(mesh);

  std::cout<<"."<<std::flush;
  CGAL::Surface_mesh_topology::Path_on_surface<Mesh> cycle=cst.compute_edge_width();
  if (cycle.length()!=nbedges)
  {
    std::cout<<"[ERROR]: number of edges for the cycle is not correct ("
             <<cycle.length()<<"!="<<nbedges<<")."<<std::endl;
    res=false;
  }

  if (is_double_torus_2_b)
  {
    if (!test_weighted(mesh, cycle))
    { res=false; }
  }

  std::cout<<"."<<std::flush;
  cycle=cst.compute_shortest_non_contractible_cycle(wf);
  double l=0.;
  for (std::size_t i=0; i<cycle.length(); ++i) { l+=wf(cycle[i]); }
  if (l>length)
  {
    std::cout<<"[ERROR] length of the cycle is not minimal ("
             <<l<<">"<<length<<")."<<std::endl;
    res=false;
  }

  auto facewidth=cst.compute_face_width();
  if (facewidth.size()!=nbfaces)
  {
    std::cout<<"[ERROR] number of faces for the face width is not correct ("
             <<facewidth.size()<<"!="<<nbfaces<<")."<<std::endl;
    res=false;
  }

  return res;
}


/// Compute edge width for the given mesh (loaded in filename), and test if
/// the size of the cycle is nbedges for unit weight and the length of the cycle
/// is smaller than length for Euclidean length weight.
bool test(const std::string filename, std::size_t nbedges, double length, std::size_t nbfaces)
{
  bool res=true;

  {
  LCC_CM lcc_cm;
  load_lcc(lcc_cm, filename);
  if (!test_one_data_structure(lcc_cm, nbedges, length, nbfaces,
                               std::string(filename)=="data/double-torus-2-b.off"))
  {
    std::cout<<"[ERROR] for Linear_cell_complex_for_combinatorial_map."<<std::endl;
    res=false;
  }
  }

  {
  LCC_GM lcc_gm;
  load_lcc(lcc_gm, filename);
  if (!test_one_data_structure(lcc_gm, nbedges, length, nbfaces,
                               std::string(filename)=="data/double-torus-2-b.off"))
  {
    std::cout<<"[ERROR] for Linear_cell_complex_for_generalized_map."<<std::endl;
    res=false;
  }
  }

  {
  SM sm;
  load_sm(sm, filename);
  if (!test_one_data_structure(sm, nbedges, length, nbfaces,
                               std::string(filename)=="data/double-torus-2-b.off"))
  {
    std::cout<<"[ERROR] for Surface_mesh."<<std::endl;
    res=false;
  }
  }

  {
  Poly poly;
  load_sm(poly, filename);
  if (!test_one_data_structure(poly, nbedges, length, nbfaces,
                               std::string(filename)=="data/double-torus-2-b.off"))
  {
    std::cout<<"[ERROR] for Polyhedron_3."<<std::endl;
    res=false;
  }
  }

  return res;
}

int main()
{
  bool res=true;
  std::cout<<"[BEGIN] Test shortest cycle non contractible "<<std::flush;

  if (!test("data/double-torus-2-a.off", 6, 5.41404, 6))
  { std::cout<<"[ERROR] for data/double-torus-2-a.off."<<std::endl; res=false; }

  if (!test("data/double-torus-2-b.off", 12, 5.41404, 6))
  { std::cout<<"[ERROR] for data/double-torus-2.off."<<std::endl; res=false; }

  // if (!test("data/double-torus-2-c.off", 24, 5.41404, 6))
  // { std::cout<<"[ERROR] for data/double-torus-2.off."<<std::endl; res=false; }

  // if (!test("data/double-torus-2-d.off", 48, 5.41404, 6))
  // { std::cout<<"[ERROR] for data/double-torus-2.off."<<std::endl; res=false; }

  // if (!test("data/elephant-with-holes.off", 5, 0.0515047, 7))
  // { std::cout<<"[ERROR] for data/elephant-with-holes.off."<<std::endl; res=false; }

  // if (!test("data/obj1.off", 4, 16.2814, 4))
  // { std::cout<<"[ERROR] for data/obj1.off."<<std::endl; res=false; } */

  // if (!test("data/obj1-with-holes.off", 4, 16.6082, 4))
  // { std::cout<<"[ERROR] for data/obj1-with-holes.off."<<std::endl; res=false; }

  if (!res)
  { return EXIT_FAILURE; }
  else
  { std::cout<<" All tests OK."<<std::endl; }

  return EXIT_SUCCESS;
}
