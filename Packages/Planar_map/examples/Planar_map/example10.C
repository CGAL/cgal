// example10

#include <CGAL/basic.h> //CGAL definitions that need to come before anything

#include <CGAL/Quotient.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_segment_exact_traits.h>

#include <CGAL/IO/Pm_iostream.h>
#include <iostream>
#include<string>

CGAL_BEGIN_NAMESPACE

template <class Pt>
class Pm_my_vertex : public Pm_vertex_base<Pt>
{
public:
  Pm_my_vertex() : Pm_vertex_base<Pt>(){}

  void  set_color(const std::string& c) { color = c;}
  
  std::string get_color() const { return color;}
  
private:
  std::string color;
};

// building new dcel with my vertex base.
template <class Traits>
class Pm_my_dcel : 
  public Pm_dcel<Pm_my_vertex<typename Traits::Point>,
                 Pm_halfedge_base<typename Traits::X_curve>, 
                 Pm_face_base> 
{
public:  // Creation
  Pm_my_dcel() {}
};

// extend the drawer to print the color as well. 
template <class PM>
class Pm_my_file_writer:  public Pm_file_writer<PM> {
public:
  
  typedef typename PM::Vertex_handle             Vertex_handle;
  typedef typename PM::Vertex_const_handle       Vertex_const_handle;
  typedef typename PM::Vertex_iterator           Vertex_iterator;
  typedef typename PM::Vertex_const_iterator     Vertex_const_iterator;

  Pm_my_file_writer(std::ostream& o, const PM& pm, bool verbose = false) : 
    Pm_file_writer<PM>(o, pm, verbose) {}
  
  void write_vertex(Vertex_handle v) {
    out() << v->point() <<"  ";
    out() << v->get_color()<< std::endl;
  }
  
  void write_vertex(Vertex_const_handle v) {
    out() << v->point() <<"  ";
    out() << v->get_color()<< std::endl;
  }
   
  void write_vertices(Vertex_iterator Vertices_begin, 
		      Vertex_iterator Vertices_end) {
    for (Vertex_iterator v_iter = Vertices_begin; 
	 v_iter != Vertices_end;
	 v_iter++)
      write_vertex(v_iter);
  }
  
  void write_vertices(Vertex_const_iterator Vertices_begin, 
		      Vertex_const_iterator Vertices_end) {
    for (Vertex_const_iterator v_iter = Vertices_begin;
	 v_iter !=  Vertices_end; 
	 v_iter++)
      write_vertex(v_iter);
  }
};

CGAL_END_NAMESPACE

typedef CGAL::Quotient<int>                NT;
typedef CGAL::Cartesian<NT>                R;
typedef CGAL::Pm_segment_exact_traits<R>   Traits;

typedef CGAL::Pm_my_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>    PM;

typedef PM::Vertex_iterator                Vertex_iterator;
//typedef PM::Vertex_const_iterator          Vertex_const_iterator;

int main(int argc, char* argv[])
{
  PM pm;
  cin >> pm;
 
  // updating the colors for halfedge and vertex.
  for (Vertex_iterator v_iter = pm.vertices_begin(); 
       v_iter != pm.vertices_end(); 
       v_iter++)
    v_iter->set_color("BLUE");

 //printing pm to output stream with the user attributes.
  CGAL::Pm_my_file_writer<PM>  writer(cout, pm); 
  write_pm(pm, writer, cout);

  return 0;
}
