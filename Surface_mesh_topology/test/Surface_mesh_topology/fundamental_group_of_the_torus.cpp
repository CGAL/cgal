/*
FUNDAMENTAL GROUP OF THE TORUS

We know that the fundamental group of the torus is the direct product ZxZ
Hence we can choose a path on a circle by choosing an pair of integers

The test generates all groups of four integers (i, j, k, l) between -5 and 5 and the associated paths pij and pkl with two different basepoints
Then it verify that
-> pij is contractible iff i==0 and j==0
-> pkl is contractible iff k==0 and l==0
-> pik is homotopic to himself
-> pkl is homotopic to himself
-> pij and pkl are homotopic iff pkl and pij are
-> pij and pkl are homotopic iff i==j and k==l
*/
#include <CGAL/Surface_mesh.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_face_graph_with_paths.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                         SM;

using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
void extend_loop_positive_i(Path_on_surface<SM>& pij, int i)
{// pre i>0
  pij.push_back_by_index(23);
  pij.extend_straight_positive(5*i-1);
}
///////////////////////////////////////////////////////////////////////////////
void extend_loop_positive_j(Path_on_surface<SM>& pij, int j)
{// pre j>0
  pij.push_back_by_index(29);
  pij.extend_straight_positive(5*j-1);
}
///////////////////////////////////////////////////////////////////////////////
void extend_loop_positive_k(Path_on_surface<SM>& pkl, int k)
{// pre i>0
  pkl.push_back_by_index(33);
  pkl.extend_straight_positive(5*k-1);
}
///////////////////////////////////////////////////////////////////////////////
void extend_loop_positive_l(Path_on_surface<SM>& pkl, int l)
{// pre j>0
  pkl.push_back_by_index(8);
  pkl.extend_straight_positive(5*l-1);
}
///////////////////////////////////////////////////////////////////////////////
void extend_loop_negative_i(Path_on_surface<SM>& pij, int i)
{// pre i<0
  pij.push_back_by_index(98);
  pij.extend_straight_positive((-5)*i-1);
}
///////////////////////////////////////////////////////////////////////////////
void extend_loop_negative_j(Path_on_surface<SM>& pij, int j)
{// pre j<0
  pij.push_back_by_index(24);
  pij.extend_straight_positive((-5)*j-1);
}
///////////////////////////////////////////////////////////////////////////////
void extend_loop_negative_k(Path_on_surface<SM>& pkl, int k)
{// pre i<0
  pkl.push_back_by_index(2);
  pkl.extend_straight_positive((-5)*k-1);
}
///////////////////////////////////////////////////////////////////////////////
void extend_loop_negative_l(Path_on_surface<SM>& pkl, int l)
{// pre j<0
  pkl.push_back_by_index(1);
  pkl.extend_straight_positive((-5)*l-1);
}
///////////////////////////////////////////////////////////////////////////////
void create_loop_ij(Path_on_surface<SM>& pij, int i, int j)
{
  pij.clear();

  if (i>0)
  { extend_loop_positive_i(pij, i); }
  else if (i<0)
  { extend_loop_negative_i(pij, i); }

  if (j>0)
  { extend_loop_positive_j(pij, j); }
  else if (j<0)
  { extend_loop_negative_j(pij, j); }
}
///////////////////////////////////////////////////////////////////////////////
void create_loop_kl(Path_on_surface<SM>& pkl, int k, int l)
{
  pkl.clear();

  if (k>0)
  { extend_loop_positive_k(pkl, k); }
  else if (k<0)
  { extend_loop_negative_k(pkl, k); }

  if (l>0)
  { extend_loop_positive_l(pkl, l); }
  else if (l<0)
  { extend_loop_negative_l(pkl, l); }
}
///////////////////////////////////////////////////////////////////////////////
int main()
{
  SM sm;
  std::ifstream in(CGAL::data_file_path("meshes/torus_quad.off"));
  if (!in.is_open())
  {
    std::cout<<"ERROR reading file data/torus_quad.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  in>>sm;

  Curves_on_surface_topology<SM> cst(sm);
  Path_on_surface<SM> pij(sm), pkl(sm);
  bool cij, ckl, hij_ij, hij_kl, hkl_ij, hkl_kl;
  bool test_valid=true;
/*
#ifdef CGAL_USE_BASIC_VIEWER
  std::vector<Path_on_surface<SM> > paths={pij, pkl};
  CGAL::draw(sm, paths); // Enable only if CGAL was compiled with Qt6
#endif // CGAL_USE_BASIC_VIEWER
*/
  for (int i=-4; i<=4; ++i)
  {
  for (int j=-4; j<=4; ++j)
  {
    for (int k=-4; k<=4; ++k)
    {
    for (int l=-4; l<=4; ++l)
    {
      create_loop_ij(pij, i, j);
      create_loop_kl(pkl, k, l);
      cij=cst.is_contractible(pij);
      ckl=cst.is_contractible(pkl);
      hij_ij=cst.are_freely_homotopic(pij, pij);
      hij_kl=cst.are_freely_homotopic(pij, pkl);
      hkl_ij=cst.are_freely_homotopic(pkl, pij);
      hkl_kl=cst.are_freely_homotopic(pkl, pkl);

      if (i==0 && j==0 && !cij)
      { std::cout<<"FAILURE : a path associated with ints "<<i<<", "<<j<<" is not contractible"<<std::endl; test_valid=false; }
      else if ((i!=0 || j!=0) && cij)
      { std::cout<<"FAILURE : a path associated with ints "<<i<<", "<<j<<" is contractible"<<std::endl; test_valid=false; }
      if (k==0 && l==0 && !ckl)
      { std::cout<<"FAILURE : a path associated with ints "<<k<<", "<<l<<" is not contractible"<<std::endl; test_valid=false; }
      else if ((k!=0 || l!=0) && ckl)
      { std::cout<<"FAILURE : a path associated with ints "<<k<<", "<<l<<" is contractible"<<std::endl; test_valid=false; }

      if (!hij_ij)
      { std::cout<<"FAILURE : a path associated with ints "<<i<<", "<<j<<" is not homotopic to himself"<<std::endl; test_valid=false; }
      if (!hkl_kl)
      { std::cout<<"FAILURE : a path associated with ints "<<k<<", "<<l<<" is not homotopic to himself"<<std::endl; test_valid=false; }
      if (hij_kl!=hkl_ij)
      { std::cout<<"FAILURE : the homotopy relation is not symmetric on paths associated with ints "<<i<<", "<<j<<" and "<<k<<", "<<l<<std::endl; test_valid=false; }

      if (i==k && j==l && (!hij_kl || !hkl_ij))
      { std::cout<<"FAILURE : paths both associated with ints "<<i<<", "<<j<<" are not homotopic"<<std::endl; test_valid=false; }
      if ((i!=k || j!=l) && (hij_kl || hkl_ij))
      { std::cout<<"FAILURE : paths associated with ints "<<i<<", "<<j<<" and "<<k<<", "<<l<<" are homotopic"<<std::endl; test_valid=false; }
    }
    }
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
