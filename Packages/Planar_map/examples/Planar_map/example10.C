// file: examples/Planar_map/example10.C


#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/write_pm.h>
#include <iostream>
#include <string>

template <class Pt>
class Pm_my_vertex : public CGAL::Pm_vertex_base<Pt>
{
public:
  Pm_my_vertex() : CGAL::Pm_vertex_base<Pt>() { }

  void set_color(const std::string & c) { color = c; }
  std::string get_color() const { return color;}
  
private:
  std::string color;
};

// building new dcel with my vertex base.
template <class Traits>
class Pm_my_dcel : 
  public CGAL::Pm_dcel<Pm_my_vertex<typename Traits::Point_2>,
               CGAL::Pm_halfedge_base<typename Traits::X_monotone_curve_2>, 
                       CGAL::Pm_face_base> 
{
public:  // Creation
  Pm_my_dcel() { }
};

// extend the drawer to print the color as well. 
template <class PM>
class Pm_my_file_writer : public CGAL::Pm_file_writer<PM>
{
public:
  typedef typename PM::Vertex_handle             Vertex_handle;
  typedef typename PM::Vertex_const_handle       Vertex_const_handle;
  typedef typename PM::Vertex_iterator           Vertex_iterator;
  typedef typename PM::Vertex_const_iterator     Vertex_const_iterator;

  Pm_my_file_writer(std::ostream & o, const PM & pm, bool verbose = false) : 
    CGAL::Pm_file_writer<PM>(o, pm, verbose) { }
  
  void write_vertex(Vertex_const_handle v) const
  {
    this->out() << v->point() <<"  ";
    this->out() << v->get_color()<< std::endl;
  }
};

typedef CGAL::Quotient<int>               NT;
typedef CGAL::Cartesian<NT>               Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel> Traits;
typedef Pm_my_dcel<Traits>                Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>   Planar_map;
typedef Planar_map::Vertex_iterator       Vertex_iterator;

int main()
{
  Planar_map pm;
  std::cin >> pm;
 
  std::cout << "* * * Demonstrating definition of user attributes for "
            << "Planar map components" << std::endl << std::endl
            << std::endl;
  
  // Update the colors for halfedge and vertex:
  for (Vertex_iterator v_iter = pm.vertices_begin(); 
       v_iter != pm.vertices_end(); 
       ++v_iter)
    v_iter->set_color("BLUE");

 // Print the map to output stream with the user attributes:
  std::cout << "* * * Printing the Planar map" << std::endl;
  std::cout << std::endl;
  
  Pm_my_file_writer<Planar_map>  writer(std::cout, pm); 
  CGAL::write_pm(pm, writer, std::cout);

  return 0;
}
