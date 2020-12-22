#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>

/* If you want to use a viewer, you can use qglviewer. */
#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/draw_face_graph_with_paths.h>
#endif

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void usage(int /*argc*/, char** argv)
{
  std::cout<<"usage: "<<argv[0]<<" file [-draw] [-L l1 l2] [-D d1 d2] "
           <<" [-N n] [-seed S] [-time]"<<std::endl
           <<"   Load the given off file, compute one random path, deform it "
           <<"into a second path and test that the two paths are homotopic."
           <<std::endl
           <<"   -draw: draw mesh and the two paths." <<std::endl
           <<"   -L l1 l2: create a path of length >= l: a random number between [l1, l2] (by default [10, 100])."<<std::endl
           <<"   -D d1 d2: use d deformations to generate the second path: d is a random number between [d1, d2] (by default [10, 100]."<<std::endl
           <<"   -N n: do n tests of homotopy (using 2*n random paths) (by default 1)."<<std::endl
           <<"   -seed S: uses S as seed of random generator. Otherwise use a different seed at each run (based on time)."<<std::endl
           <<"   -time: display computation times."<<std::endl
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
                          int& l1,
                          int& l2,
                          int& d1,
                          int& d2,
                          unsigned int& N,
                          CGAL::Random& random,
                          bool& time)
{
  std::string arg;
  bool usage_required=false;
  for (int i=1; i<argc; ++i)
  {
    arg=argv[i];
    if (arg=="-draw")
    { draw=true; }
    else if (arg=="-D")
    {
     if (i+2>=argc)
     { error_command_line(argc, argv, "not enough number after -D option."); }
     d1=std::stoi(std::string(argv[++i]));
     d2=std::stoi(std::string(argv[++i]));
    }
    else if (arg=="-L")
    {
     if (i+2>=argc)
     { error_command_line(argc, argv, "not enough number after -L option."); }
     l1=std::stoi(std::string(argv[++i]));
     l2=std::stoi(std::string(argv[++i]));
    }
    else if (arg=="-N")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "no number after -nbtests option."); }
      N=static_cast<unsigned int>(std::stoi(std::string(argv[++i])));
    }
    else if (arg=="-seed")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "no number after -seed option."); }
      random=CGAL::Random(static_cast<unsigned int>
                          (std::stoi(std::string(argv[++i]))));
      // initialize the random generator with the given seed
    }
    else if (arg=="-time")
    { time=true; }
    else if (arg=="-h" || arg=="--help" || arg=="-?")
    { usage_required=true; }
    else if (arg[0]=='-')
    { std::cout<<"[WARNING] Unknown option "<<arg<<", ignored."<<std::endl; }
    else
    { file=arg; }
  }

  if (N==0) { N=1; }
  if (l2<l1) l2=l1;
  if (d2<d1) d2=d1;

  if (usage_required) {usage(argc, argv); }
}
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  std::string file="data/elephant.off";
  bool draw=false;
  int l1=10, l2=100;
  int d1=10, d2=100;
  unsigned int N=1;
  CGAL::Random random; // Used when user do not provide its own seed.
  bool time=false;

  process_command_line(argc, argv, file, draw, l1, l2, d1, d2, N, random, time);

  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, file.c_str()))
  {
    std::cout<<"PROBLEM reading file "<<file<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout<<"Initial map: ";
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid() << std::endl;

  Curves_on_surface_topology<LCC_3_cmap> cst(lcc, time);

  cst.compute_minimal_quadrangulation(time);

  std::cout<<"Reduced map: ";
  cst.get_minimal_quadrangulation().display_characteristics(std::cout)
      <<", valid="<< cst.get_minimal_quadrangulation().is_valid()<<std::endl;

  unsigned int nbcontractible=0;
  std::vector<std::size_t> errors_seeds;
  unsigned int length, defo;

  for (unsigned int i=0; i<N; ++i)
  {
    if (i!=0)
    {
      random=CGAL::Random(static_cast<unsigned int>
                          (random.get_int(0, (std::numeric_limits<int>::max)())));
    }
    std::cout<<"Random seed: "<<random.get_seed()<<": ";

    length=static_cast<unsigned int>(random.get_int(l1, l2+1));
    defo=static_cast<unsigned int>(random.get_int(d1, d2+1));

    std::vector<Path_on_surface<LCC_3_cmap> > paths;
    std::vector<Path_on_surface<LCC_3_cmap> > transformed_paths;

    Path_on_surface<LCC_3_cmap> path1(lcc);
    path1.generate_random_closed_path(length, random);

    std::cout<<"Path1 size: "<<path1.length()<<" (from "<<length<<" darts); ";
    paths.push_back(path1);

    Path_on_surface<LCC_3_cmap> path2(path1);
    path2.update_path_randomly(defo, random);
    std::cout<<"Path2 size: "<<path2.length()<<" (from "<<defo<<" deformations).";
    paths.push_back(path2);
    std::cout<<std::endl;

    if (cst.is_contractible(path1, time))
    { ++nbcontractible; }

    bool res=cst.are_freely_homotopic(path1, path2, time);
    if (!res)
    { errors_seeds.push_back(random.get_seed()); }

#ifdef CGAL_USE_BASIC_VIEWER
    if (draw)
    { CGAL::draw(lcc, paths); }
#endif
  }

  if (errors_seeds.empty())
  {
    if (N==1) { std::cout<<"Test OK: both paths are homotopic."<<std::endl; }
    else
    { std::cout<<"All the "<<N
               <<" tests OK: each pair of paths were homotopic."<<std::endl; }
  }
  else
  {
    std::cout<<"ERRORS for "<<errors_seeds.size()<<" tests among "<<N
            <<" (i.e. "<<static_cast<double>(errors_seeds.size()*100)/
              static_cast<double>(N)<<"%)."<<std::endl;
    std::cout<<"Errors for seeds: ";
    for (std::size_t i=0; i<errors_seeds.size(); ++i)
    { std::cout<<errors_seeds[i]<<"  "; }
    std::cout<<std::endl;
  }

  std::cout<<"Number of contractible paths: "<<nbcontractible<<" among "<<N
           <<" (i.e. "<<static_cast<double>(nbcontractible*100)/
             static_cast<double>(N)<<"%)."<<std::endl;

  return EXIT_SUCCESS;
}
