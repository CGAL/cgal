#include <CGAL/Surface_mesh.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                         SM;

using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
bool test_constructions()
{
  bool res=true;
  SM sm;
  std::ifstream in("data/cylinder-with-two-borders.off");
  if (!in.is_open())
  {
    std::cout<<"ERROR reading file data/cylinder-with-two-borders.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  in>>sm;

  Path_on_surface<SM> p1(sm), p2(sm);

  p1.push_back_by_index(24);
  p1.extend_positive_turn(1);
  p1.extend_negative_turn(1);
  p1.extend_straight_positive();
  p1.extend_straight_negative();
  p1.push_back_by_index(93, true);
  p1.extend_positive_turn(3);
  p1.extend_positive_turn(3);
  p1.extend_negative_turn(1);
  p1.extend_straight_positive(5);
  p1.extend_straight_negative(5);
  if (!p1.is_valid())
  { std::cout<<"ERROR: Path p1 NOT is valid, test failed"<<std::endl; res=false; }

  p2.push_back_by_index(2);
  p2.extend_positive_turn(1);
  p2.extend_negative_turn(1);
  p2.extend_straight_positive();
  p2.extend_straight_negative();
  p2.push_back_by_index(93, true);
  p2.extend_positive_turn(3);
  p2.extend_positive_turn(3);
  p2.extend_negative_turn(1);
  p2.extend_straight_positive(5);
  p2.extend_straight_negative(5);
  if (p2.is_valid())
  { std::cout<<"ERROR: Path p2 is valid, test failed"<<std::endl; res=false; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool test_random_path()
{
  SM sm;
  bool res=true;
  std::ifstream in("data/rond_point_saucisse.off");
  if (!in.is_open())
  {
    std::cout<<"ERROR reading file data/rond_point_saucisse.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  in>>sm;

  Path_on_surface<SM> p1(sm), p2(sm), p3(sm);
  CGAL::Random random;

  p1.generate_random_path(100, random, true);
  if (!p1.is_valid())
  { std::cout<<"ERROR: Path p1 is NOT valid, test failed."<<std::endl; res=false; }

  p3.generate_random_path(100, random, false);
  if (!p3.is_valid())
  { std::cout<<"ERROR: Path p3 is NOT valid, test failed."<<std::endl; res=false; }

  p2.generate_random_closed_path(100, random);
  p2.update_is_closed();
  if (!p2.is_valid())
  { std::cout<<"ERROR: Path p2 is NOT valid, test failed."<<std::endl; res=false; }
  if (!p2.is_closed())
  { std::cout<<"ERROR: Path p2 is NOT closed, test failed."<<std::endl; res=false; }

  if (!res)
  { std::cout<<"Wrong seed="<<random.get_seed()<<std::endl; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
int main()
{
  bool success=true;
  success=test_constructions() && success;
  for (unsigned i=0; i<1000; ++i)
  { success=test_random_path() && success; }

  std::cout<<"tests_path_on_surface "<<(success?"OK":"failed") <<std::endl;

  if (!success)
  { return EXIT_FAILURE; }

  return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
