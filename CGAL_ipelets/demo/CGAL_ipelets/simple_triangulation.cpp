#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/CGAL_Ipelet_base.h>

namespace my_triangulation{

typedef CGAL::Exact_predicates_inexact_constructions_kernel               Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel>                            Delaunay;

//Function names of the ipelet
const std::string labels[] = {  "Delaunay","Help" };
//Help message associated to the first function
const std::string hmsg[] = {
  "Draw a Delaunay triangulation of a set of points"
};

class Triangulation_ipelet
  : public CGAL::Ipelet_base<Kernel,2>{
public:
  //declare an ipelet called CGAL Delaunay, with 2 functions (including help message).
  Triangulation_ipelet()
    :CGAL::Ipelet_base<Kernel,2>("CGAL Delaunay",labels,hmsg){}
  void protected_run(int);
};

//function called when using the ipelet.
void Triangulation_ipelet::protected_run(int fn)
{
  switch (fn){
    case 1:
      show_help();//print an help message
      return;
    default:
    std::list<Point_2> pt_lst;

    //Recovering points using output iterator of type
    //Dispatch_or_drop_output_iterator
    read_active_objects(
      CGAL::dispatch_or_drop_output<Point_2>(std::back_inserter(pt_lst))
    );

    if (pt_lst.empty()) {
      print_error_message("No mark selected");
      return;
    }

    Delaunay dt;
    dt.insert(pt_lst.begin(),pt_lst.end());

    //draw the triangulation.
    draw_in_ipe(dt);
  };
}

}//namespace my_triangulation

//register the ipelet in Ipe
CGAL_IPELET(my_triangulation::Triangulation_ipelet)

