#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/draw_face_graph_with_paths.h>

using Kernel         =CGAL::Simple_cartesian<double>;
using Point          =Kernel::Point_3;
using Mesh           =CGAL::Surface_mesh<Point>;
using Path_on_surface=CGAL::Surface_mesh_topology::Path_on_surface<Mesh>;
using CST            =CGAL::Surface_mesh_topology::Curves_on_surface_topology<Mesh>;

struct Weight_functor
{
  using Weight_t=double;
  Weight_functor(const Mesh& mesh) : m_mesh(mesh) {}
  double operator()(Mesh::Halfedge_index he) const
  {
    const Point& A=m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 0));
    const Point& B=m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 1));
    return CGAL::sqrt(CGAL::squared_distance(A, B));
  }
private:
  const Mesh& m_mesh;
};

int main(int argc, char* argv[])
{
  std::cout<<"Program edgewidth_surface_mesh started."<<std::endl;
  std::string filename(argc==1?"data/3torus.off":argv[1]);
  bool draw=(argc<3?false:(std::string(argv[2])=="-draw"));
  std::ifstream inp(filename);
  if (inp.fail())
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    return EXIT_FAILURE;
  }
  Mesh sm;
  inp>>sm;
  std::cout<<"File '"<<filename<<"' loaded. Finding edge-width of the mesh..."<<std::endl;

  Weight_functor wf(sm);
  CST            cst(sm, true);

  Path_on_surface cycle=cst.compute_edgewidth(wf, true);
  if (cycle.length()==0)
  { std::cout<<"  Cannot find edge-width. Stop."<<std::endl; }
  else
  {
    double cycle_length=0;
    for (int i=0; i<cycle.length(); ++i)
    { cycle_length+=wf(cycle[i]); }

    std::cout<<"  Number of edges in cycle: "<<cycle.length()<<std::endl;
    std::cout<<"  Cycle length: "<<cycle_length<<std::endl;
    std::cout<<"  Root: "<<sm.point(sm.vertex(sm.edge(cycle[0]), 0))<<std::endl;
    if (draw) { CGAL::draw(sm, cycle); }
  }

  return EXIT_SUCCESS;
}
