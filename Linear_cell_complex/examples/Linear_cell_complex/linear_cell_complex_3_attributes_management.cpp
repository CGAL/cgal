#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>

typedef CGAL::Linear_cell_complex<2,3> LCC_3;
typedef LCC_3::Dart_handle             Dart_handle;
typedef LCC_3::Point                   Point;
typedef LCC_3::FT                      FT;

void load_and_simplify_off(LCC_3& lcc, char* filename,
                           bool updateattribs, int percent)
{
  std::ifstream ifile(filename);
  int nb=0;
  if (ifile)
  {
    CGAL::load_off(lcc, ifile);
    CGAL::Timer timer;
    Dart_handle dh;
    unsigned int nb=(lcc.number_of_darts()*percent)/200;
    timer.start();  
    for (LCC_3::Dart_range::iterator it=lcc.darts().begin(),
           itend=lcc.darts().end(); it!=itend && nb>0; )
    {
      dh=it++;
      // if ( dh < lcc.beta<2>(dh) )
      // {
        if ( it!=itend && it==lcc.beta<2>(dh) ) ++it;      
        CGAL::remove_cell<LCC_3, 1>(lcc, dh, updateattribs);
        --nb;
      // }
    }
    if ( !updateattribs ) lcc.correct_invalid_attributes();
    
    timer.stop();
    lcc.display_characteristics(std::cout);
    std::cout<<", valid="<< lcc.is_valid()
             <<" time: "<<timer.time()<<" seconds." << std::endl;
  }
}

int main(int narg, char** argv)
{
  if ( narg==1 )
  {
    std::cout<<"Usage: a.out file.off [percentage]"<<std::endl;
    return EXIT_FAILURE;
  }
  
  int percent = 30; // remove 30 percent of edges 
  if ( narg>2 ) { percent = atoi(argv[2]); }
  std::cout<<percent<<"% edges to remove."<<std::endl;
  {
    LCC_3 lcc;
    std::cout<<"Update attribute DURING operations: ";
    load_and_simplify_off(lcc, argv[1], true, percent);
  }
  {
    LCC_3 lcc2;
    std::cout<<"Update attribute AFTER operations: ";
    load_and_simplify_off(lcc2, argv[1], false, percent);
  }
  return EXIT_SUCCESS;
}
