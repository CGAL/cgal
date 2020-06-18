/*
FUNDAMENTAL GROUP OF THE CIRCLE

We know that the fundamental group of the torus is Z
Hence we can choose a path on a circle by choosing an integer

The test generates all pairs of integer (i, j) between -10 and 10 and the associated paths pi and pj with two different basepoints
Then it verify that
-> pi is contractible iff i==0
-> pj is contractible iff j==0
-> pi is homotopic to himself
-> pj is homotopic to himself
-> pi and pj are homotopic iff pj and pi are
-> pi and pj are homotopic iff i==j
*/
#include <CGAL/Surface_mesh.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>

// If you want to use a viewer, you can use qglviewer.
#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/draw_face_graph_with_paths.h>
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                         SM;

using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
void create_positive_loop_88(Path_on_surface<SM>& p, unsigned int n)
{
  p.clear();
  p.push_back_by_index(88, false, false);
  if (n==0)
  {
    p.push_back_by_index(88, true, true);
    CGAL_assertion(p.is_closed());
  }
  else
  {
    p.extend_straight_positive(10*n-1);
    CGAL_assertion(p.is_closed());
  }
}
///////////////////////////////////////////////////////////////////////////////
void create_negative_loop_88(Path_on_surface<SM>& p, unsigned int n)
{
  p.clear();
  p.push_back_by_index(88, true, false);
  if (n==0)
  {
    p.push_back_by_index(88, false, true);
    CGAL_assertion(p.is_closed());
  }
  else
  {
    p.extend_straight_positive(10*n-1);
    CGAL_assertion(p.is_closed());
  }
}
///////////////////////////////////////////////////////////////////////////////
void create_positive_loop_24(Path_on_surface<SM>& p, unsigned int n)
{
  p.clear();
  p.push_back_by_index(24, false, false);
  if (n==0)
  {
    p.push_back_by_index(24, true, true);
    CGAL_assertion(p.is_closed());
  }
  else
  {
    p.extend_straight_positive((10*n)-1);

#ifdef CGAL_USE_BASIC_VIEWER
  /* std::vector<Path_on_surface<SM> > v;
  v.push_back(p);
  CGAL::draw(p.get_mesh(), v, "Title"); */
#endif // CGAL_USE_BASIC_VIEWER

    CGAL_assertion(p.is_closed());
  }
}
///////////////////////////////////////////////////////////////////////////////
void create_negative_loop_24(Path_on_surface<SM>& p, unsigned int n)
{
  p.clear();
  p.push_back_by_index(24, true, false);
  if (n==0)
  {
    p.push_back_by_index(24, false, true);
    CGAL_assertion(p.is_closed());
  }
  else
  {
    p.extend_straight_negative((10*n)-1);
    CGAL_assertion(p.is_closed());
  }
}
///////////////////////////////////////////////////////////////////////////////
void create_loop_88(Path_on_surface<SM>& p, int n)
{
  if (n>0)
  { create_positive_loop_88(p, n); }
  else
  { create_negative_loop_88(p, -n); }
}
///////////////////////////////////////////////////////////////////////////////
void create_loop_24(Path_on_surface<SM>& p, int n)
{
  if (n>0)
  { create_positive_loop_24(p, n); }
  else
  { create_negative_loop_24(p, -n); }
}
///////////////////////////////////////////////////////////////////////////////
int main()
{
  SM sm;
  std::ifstream in("data/cylinder-with-two-borders.off");
  if (!in.is_open())
  {
    std::cout<<"ERROR reading file data/cylinder-with-two-borders.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  in>>sm;

  Curves_on_surface_topology<SM> cst(sm);
  Path_on_surface<SM> p1(sm), p2(sm);
  bool c1, c2, h11, h12, h21, h22;
  bool test_valid=true;

  for (int i=-10; i<=10; ++i)
  {
    for (int j=-10; j<=10; ++j)
    {
      create_loop_88(p1, i);
      create_loop_24(p2, j);
      c1=cst.is_contractible(p1);
      c2=cst.is_contractible(p2);
      h11=cst.are_freely_homotopic(p1, p1);
      h12=cst.are_freely_homotopic(p1, p2);
      h21=cst.are_freely_homotopic(p2, p1);
      h22=cst.are_freely_homotopic(p2, p2);

      if (i==0 && !c1)
      { std::cout<<"FAILURE : a path associated with int "<<i<<" is not contractible"<<std::endl; test_valid=false; }
      else if (i!=0 && c1)
      { std::cout<<"FAILURE : a path associated with int "<<i<<" is contractible"<<std::endl; test_valid=false; }
      if (j==0 && !c2)
      { std::cout<<"FAILURE : a path associated with int "<<j<<" is not contractible"<<std::endl; test_valid=false; }
      else if (j!=0 && c2)
      { std::cout<<"FAILURE : a path associated with int "<<j<<" is contractible"<<std::endl; test_valid=false; }

      if (!h11)
      { std::cout<<"FAILURE : a path associated with int "<<i<<" is not homotopic to himself"<<std::endl; test_valid=false; }
      if (!h22)
      { std::cout<<"FAILURE : a path associated with int "<<j<<" is not homotopic to himself"<<std::endl; test_valid=false; }
      if (h12!=h21)
      { std::cout<<"FAILURE : the homotopy relation is not symetric on paths associated with ints "<<i<<", "<<j<<std::endl; test_valid=false; }

      if (i==j && (!h12 || !h21))
      { std::cout<<"FAILURE : paths both associated with int "<<i<<" are not homotopic"<<std::endl; test_valid=false; }
      if (i!=j && (h12 || h21))
      { std::cout<<"FAILURE : paths associated with ints "<<i<<", "<<j<<" are homotopic"<<std::endl; test_valid=false; }
    }
  }

  if (test_valid)
  {
    std::cout<<"All tests are OK"<<std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
///////////////////////////////////////////////////////////////////////////////
