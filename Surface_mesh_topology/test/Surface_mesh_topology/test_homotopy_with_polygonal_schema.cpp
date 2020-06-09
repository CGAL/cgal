#include <CGAL/Polygonal_schema.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <iostream>
#include <ctime>
#include <cstdlib>

typedef CGAL::Surface_mesh_topology::Polygonal_schema_with_combinatorial_map<> PS;
typedef typename PS::Dart_handle  Dart_handle;

using namespace CGAL::Surface_mesh_topology;

static unsigned int starting_seed;

///////////////////////////////////////////////////////////////////////////////
bool test_two_random_paths(const PS& ps,
                           const Curves_on_surface_topology<PS>& cst,
                           unsigned int nbtests,
                           unsigned int& seed,
                           int lmin=5,
                           int lmax=20)
{
  CGAL_assertion(lmin>0 && lmin<lmax);
  CGAL::Random random(seed++);
  Path_on_surface<PS> p1(ps), p2(ps);
  internal::generate_random_closed_path
    (p1, static_cast<std::size_t>(random.get_int(lmin, lmax)), random);

  p2=p1;
  p2.update_path_randomly(100, random);

  // std::cout<<"("<<p1.length()<<", "<<p2.length()<<")"<<std::flush;
  bool res=true;
  if (!cst.are_freely_homotopic(p1, p2))
  {
    res=false;
    std::cout<<std::endl<<"ERROR for path number "<<nbtests<<" (seed="
             <<random.get_seed()<<")"<<std::endl;
  }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool run_n_random_paths_tests(const PS& ps,
                              unsigned int n,
                              const Curves_on_surface_topology<PS>& cst,
                              unsigned int& seed,
                              int lmin=5,
                              int lmax=20)
{
  static std::size_t nbtest=0;
  std::cout<<"***** Test "<<nbtest++<<" "<<std::flush;
  ps.display_characteristics(std::cout)<<"   ";
  bool res=true;
  for (unsigned int i=0; i<n; ++i)
  {
    std::cout<<"."<<std::flush;
    if (!test_two_random_paths(ps, cst, i, seed, lmin, lmax))
    { res=false; }
  }
  std::cout<<(res?" OK.":" FAILED.")<<std::endl;
  return res;
}
///////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void usage(int /*argc*/, char** argv)
{
  std::cout<<"usage: "<<argv[0]<<" [-seed S]"<<std::endl
           <<"   Run different homotopy on 2-maps built using polygonal schema."
           <<std::endl
           <<"   -seed S: uses S as seed of random generator. Otherwise use a"
           <<" different seed at each run (based on time)."<<std::endl
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
                          unsigned int& seed)
{
  seed=static_cast<unsigned int>(std::time(nullptr));
  std::string arg;
  for (int i=1; i<argc; ++i)
  {
    arg=argv[i];
    if (arg=="-seed")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "Error: no number after -seed option."); }
      seed=static_cast<unsigned int>(std::stoi(std::string(argv[++i])));
      // initialize the random generator with the given seed
    }
    else if (arg=="-h" || arg=="--help" || arg=="-?")
    { usage(argc, argv); }
    else if (arg[0]=='-')
    { std::cout<<"Unknown option "<<arg<<", ignored."<<std::endl; }
  }
}
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  process_command_line(argc, argv, starting_seed);
  std::cout<<"Start tests="<<starting_seed<<"."<<std::endl;

  bool closed=true;
  std::size_t percent=0;
  bool res=true;
  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<5; ++j)
    { // 5 tests for the same "type" of polygonal schema
      //    (open/closed; with  or without perforated)
      PS ps;
      Curves_on_surface_topology<PS> cst(ps);
      generate_random_polygonal_schema(ps, 1000, starting_seed++,  30, closed,  percent);   // Closed, no perforated faces
      res&=run_n_random_paths_tests(ps, 10, cst, starting_seed);
    }
    switch(i)
    {
      case 0: closed=true; percent=15; break;
      case 1: closed=false; percent=0; break;
      case 2: closed=false; percent=15; break;
    }
  }

  if (res)
  {
    std::cout<<"test_homotopy_with_polygonal_schema ALL TESTS OK."<<std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
///////////////////////////////////////////////////////////////////////////////
