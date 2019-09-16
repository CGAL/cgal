#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3;
typedef LCC_3::Dart_handle             Dart_handle;
typedef LCC_3::Point                   Point;
typedef LCC_3::FT                      FT;

void load_and_simplify_off(LCC_3& lcc, const std::string& filename,
                           bool updateattribs, int percent)
{
  std::ifstream ifile(filename.c_str());
  if (ifile)
  {
    CGAL::load_off(lcc, ifile);
    CGAL::Timer timer;
    Dart_handle dh;
    std::size_t nb=(lcc.number_of_darts()*percent)/200;
    timer.start();

    if (!updateattribs) lcc.set_automatic_attributes_management(false);
    for (LCC_3::Dart_range::iterator it=lcc.darts().begin(),
           itend=lcc.darts().end(); it!=itend && nb>0; )
    {
      dh=it++;
      if ( it!=itend && it==lcc.beta<2>(dh) ) ++it;
      lcc.remove_cell<1>(dh);
      --nb;
    }
    if ( !updateattribs ) lcc.set_automatic_attributes_management(true);

    timer.stop();
    lcc.display_characteristics(std::cout);
    std::cout<<", valid="<< lcc.is_valid()
             <<" time: "<<timer.time()<<" seconds." << std::endl;
  }
}

int main(int narg, char** argv)
{
  if (narg>1 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"-?")) )
  {
    std::cout<<"Usage: a.out file.off [percentage]"<<std::endl;
    return EXIT_FAILURE;
  }

  std::string filename;
  if ( narg==1 )
  {
    filename=std::string("data/armadillo.off");
    std::cout<<"No filename given: use data/armadillo.off by default."<<std::endl;
  }
  else filename=std::string(argv[1]);

  int percent = 30; // remove 30 percent of edges
  if ( narg>2 ) { percent = atoi(argv[2]); }
  std::cout<<percent<<"% edges to remove."<<std::endl;

  LCC_3 lcc;
  std::cout<<"Update attribute DURING operations: ";
  load_and_simplify_off(lcc, filename, true, percent);

  LCC_3 lcc2;
  std::cout<<"Update attribute AFTER operations: ";
  load_and_simplify_off(lcc2, filename, false, percent);

  return EXIT_SUCCESS;
}
