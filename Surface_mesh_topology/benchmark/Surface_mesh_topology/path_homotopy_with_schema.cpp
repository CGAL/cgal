#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Polygonal_schema.h>

///////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void usage(int /*argc*/, char** argv)
{
  std::cout<<"usage: "<<argv[0]<<" scheme [-nbdefo D] [-nbedges E] "
           <<" [-nbtests N] [-seed S] [-time]"<<std::endl
           <<"   Load the given polygonal scheme (by default \"a b -a -b\"),"
           << "compute one random path, deform it "
           <<"into a second path and test that the two paths are homotope."
           <<std::endl
           <<"   -nbdefo D: use D deformations to generate the second path (by default a random number between 10 and 100)."<<std::endl
           <<"   -nbedges E: generate paths of length E (by default a random number between 10 and 100)."<<std::endl
           <<"   -nbtests N: do N tests of homotopy (using 2*N random paths) (by default 1)."<<std::endl
           <<"   -seed S: uses S as seed of random generator. Otherwise use a different seed at each run (based on time)."<<std::endl
           <<"   -time: display computation times."<<std::endl
           <<std::endl;
  exit(EXIT_FAILURE);
}
///////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void error_command_line(int argc, char** argv, const char* msg)
{
  std::cout<<"ERROR: "<<msg<<std::endl;
  usage(argc, argv);
}
///////////////////////////////////////////////////////////////////////////////
void process_command_line(int argc, char** argv,
                          std::string& file,
                          bool& withD,
                          unsigned int& D,
                          bool& withE,
                          unsigned int& E,
                          unsigned int& N,
                          CGAL::Random& random,
                          bool& time)
{
  std::string arg;
  for (int i=1; i<argc; ++i)
  {
    arg=argv[i];
    if (arg=="-nbdefo")
    {
     if (i==argc-1)
     { error_command_line(argc, argv, "Error: no number after -nbdefo option."); }
     withD=true;
     D=static_cast<unsigned int>(std::stoi(std::string(argv[++i])));
    }
    else if (arg=="-nbedges")
    {
     if (i==argc-1)
     { error_command_line(argc, argv, "Error: no number after -nbedges option."); }
     withE=true;
     E=static_cast<unsigned int>(std::stoi(std::string(argv[++i])));
    }
    else if (arg=="-nbtests")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "Error: no number after -nbtests option."); }
      N=static_cast<unsigned int>(std::stoi(std::string(argv[++i])));
    }
    else if (arg=="-seed")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "Error: no number after -seed option."); }
      random=CGAL::Random(static_cast<unsigned int>(std::stoi(std::string(argv[++i])))); // initialize the random generator with the given seed
    }
    else if (arg=="-time")
    { time=true; }
    else if (arg=="-h" || arg=="--help" || arg=="-?")
    { usage(argc, argv); }
    else if (arg[0]=='-')
    { std::cout<<"Unknown option "<<arg<<", ignored."<<std::endl; }
    else
    { file=arg; }
  }

  if (N==0) { N=1; }
}
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  using namespace CGAL::Surface_mesh_topology;
  std::string scheme="a b -a -b";
  bool withD=false;
  unsigned int D=0;
  bool withE=false;
  unsigned int E=0;
  unsigned int N=1;
  CGAL::Random random; // Used when user do not provide its own seed.
  bool time=false;

  process_command_line(argc, argv, scheme, withD, D, withE, E, N, random, time);

  Polygonal_schema_with_combinatorial_map<> cm;
  cm.add_facet(scheme);
  std::cout<<"Initial map: ";
  cm.display_characteristics(std::cout) << ", valid="
                                         << cm.is_valid() << std::endl;

  Curves_on_surface_topology<Polygonal_schema_with_combinatorial_map<> > smct(cm, time);
  smct.compute_minimal_quadrangulation(time);
  std::cout<<"Reduced map: ";
  smct.get_minimal_quadrangulation().display_characteristics(std::cout)
    << ", valid="<< smct.get_minimal_quadrangulation().is_valid() << std::endl;

  unsigned int nbcontractible=0;
  std::vector<std::size_t> errors_seeds;

  for (unsigned int i=0; i<N; ++i)
  {
    if (i!=0)
    {
      random=CGAL::Random(random.get_int(0, std::numeric_limits<int>::max()));
    }
    std::cout<<"Random seed: "<<random.get_seed()<<": ";

    if (!withE)
    { E=static_cast<unsigned int>(random.get_int
                                  (10, std::max(std::size_t(11),
                                                cm.number_of_darts()/10))); }

    if (!withD)
    { D=static_cast<unsigned int>(random.get_int
                                  (10, std::max(std::size_t(11),
                                                cm.number_of_darts()/10))); }



    Path_on_surface<Polygonal_schema_with_combinatorial_map<>> path1(cm);
    path1.generate_random_closed_path(E, random);

    std::cout<<"Path1 size: "<<path1.length()<<" (from "<<E<<" darts); ";
    Path_on_surface<Polygonal_schema_with_combinatorial_map<>> deformed_path1(path1);
    deformed_path1.update_path_randomly(D, random);
    std::cout<<"Path2 size: "<<deformed_path1.length()<<" (from "<<D<<" deformations).";
    std::cout<<std::endl;

    if (smct.is_contractible(path1, time))
    { ++nbcontractible; }

    bool res=smct.are_freely_homotopic(path1, deformed_path1, time);
    if (!res)
    { errors_seeds.push_back(random.get_seed()); }
  }

  if (errors_seeds.empty())
  {
    if (N==1) { std::cout<<"Test OK: both paths are homotopic."<<std::endl; }
    else { std::cout<<"All the "<<N<<" tests OK: each pair of paths were homotopic."<<std::endl; }
  }
  else
  {
    std::cout<<"ERRORS for "<<errors_seeds.size()<<" tests among "<<N
            <<" (i.e. "<<(double)(errors_seeds.size()*100)/double(N)<<"%)."<<std::endl;
    std::cout<<"Errors for seeds: ";
    for (std::size_t i=0; i<errors_seeds.size(); ++i)
    { std::cout<<errors_seeds[i]<<"  "; }
    std::cout<<std::endl;
  }

  std::cout<<"Number of contractible paths: "<<nbcontractible<<" among "<<N
           <<" (i.e. "<<(double)(nbcontractible*100)/double(N)<<"%)."<< std::endl
           << "==============" <<std::endl;

  return EXIT_SUCCESS;
}
