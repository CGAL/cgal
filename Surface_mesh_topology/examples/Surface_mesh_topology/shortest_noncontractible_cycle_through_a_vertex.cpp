#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/draw_face_graph_with_paths.h>
#include <CGAL/Random.h>

///////////////////////////////////////////////////////////////////////////////
using LCC_3            =CGAL::Linear_cell_complex_for_generalized_map<2, 3>;
using Dart_const_handle=LCC_3::Dart_const_handle;
using Dart_handle=LCC_3::Dart_handle; // TODO REMOVE
using Path_on_surface  =CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;
using CST              =CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;
///////////////////////////////////////////////////////////////////////////////
struct Weight_functor
{
  Weight_functor(const LCC_3& lcc) : m_lcc(lcc) { }
  using Weight_t=double;
  Weight_t operator()(Dart_const_handle dh) const
  {
    const LCC_3::Point& x=m_lcc.point(dh);
    const LCC_3::Point& y=m_lcc.point(m_lcc.alpha<0>(dh));
    return CGAL::sqrt(CGAL::squared_distance(x, y));
  }
private:
  const LCC_3& m_lcc;
};
///////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void usage(int /*argc*/, char* argv[])
{
  std::cout<<"Usage "<<argv[0]<<" [-draw] [-dist] [-time] [-seed S] file.off"<<std::endl
           <<"  Compute the shortest non contractible cycle passing through a vertex on the given mesh."
           <<"  Vertex is choosen randomly."<<std::endl
           <<"   -draw: draw the mesh and the cycle."<<std::endl
           <<"   -dist: use Euclidean distance (instead of bfs by default)."<<std::endl
           <<"   -time: show computation times."<<std::endl
           <<"   -seed S: uses S as seed of random generator. Otherwise use a different seed at each run (based on time)."<<std::endl
           <<std::endl;
  exit(EXIT_FAILURE);
}
///////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void error_command_line(int argc, char** argv, const char* msg)
{
  std::cout<<"[ERROR] "<<msg<<std::endl;
  usage(argc, argv);
}
///////////////////////////////////////////////////////////////////////////////
void process_command_line(int argc, char** argv,
                          std::string& file,
                          bool& draw,
                          bool& dist,
                          bool& time,
                          CGAL::Random& random)
{
  std::string arg;
  bool usage_required=false;
  for (int i=1; i<argc; ++i)
  {
    arg=argv[i];
    if (arg=="-draw")
    { draw=true; }
    else if (arg=="-dist")
    { dist=true; }
    else if (arg=="-time")
    { time=true; }
    else if (arg=="-seed")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "no number after -seed option."); }
      random=CGAL::Random(static_cast<unsigned int>
                          (std::stoi(std::string(argv[++i]))));
      // initialize the random generator with the given seed
    }
    else if (arg=="-h" || arg=="--help" || arg=="-?")
    { usage_required=true; }
    else if (arg[0]=='-')
    { std::cout<<"[WARNING] Unknown option "<<arg<<", ignored."<<std::endl; }
    else
    { file=arg; }
  }

  if (usage_required) {usage(argc, argv); }
}
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  std::cout<<"Program shortest_noncontractible_cycle_through_a_vertex started."
           <<std::endl;
  std::string filename("data/3torus.off");
  bool draw=false;
  bool dist=false;
  bool time=false;
  CGAL::Random random; // Used when user do not provide its own seed.
  process_command_line(argc, argv, filename, draw, dist, time, random);
  
  std::ifstream inp(filename);
  if (inp.fail())
  {
    std::cout<<"[ERROR] Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    usage(argc, argv);
  }
  
  LCC_3 lcc;
  CGAL::load_off(lcc, inp);
  std::cout<<"File '"<<filename<<"' loaded. Running the main program (seed="
           <<random.get_seed()<<")"<<std::endl;

  Weight_functor wf(lcc);
  CST            cst(lcc, time);

  /// Change the value of `root` to test the algorithm at another vertex
  Dart_handle root=lcc.dart_handle(random.get_int(0, lcc.number_of_darts()));
  std::cout<<"Finding the shortest noncontractible cycle..."<<std::endl;
  Path_on_surface cycle(lcc);
  if (dist)
  { cycle=cst.compute_shortest_noncontractible_cycle_with_basepoint(root, wf, time); }
  else
  { cycle=cst.compute_shortest_noncontractible_cycle_with_basepoint(root, time); }
  
  if (cycle.length()==0)
  { std::cout<<"  Cannot find such cycle. Stop."<<std::endl; }
  else
  {
    double cycle_length=0;
    for (int i=0; i<cycle.length(); ++i)
    { cycle_length+=wf(cycle[i]); }

    std::cout<<"  Number of edges in cycle: "<<cycle.length()<<std::endl;
    std::cout<<"  Cycle length: "<<cycle_length<<std::endl; 
    std::cout<<"  Root: "<<lcc.point(root)<<std::endl;
    if (draw) { CGAL::draw(lcc, cycle); }
  }

  return EXIT_SUCCESS;
}
