#include <CGAL/basic.h>
#include <cassert>
#include <fstream>
#include <list>
#include <stack>
#include <set>

#if !defined(CGAL_USE_LEDA) && !defined(CGAL_USE_CGAL_WINDOW)
int main(int argc, char* argv[])
{
  std::cout <<"Sorry, this demo needs LEDA or CGAL::WINDOW for visualisation.";
  std::cout << std::endl;
  return 0;
}

#else
#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_yellow CGAL::yellow
#endif //CGAL_USE_CGAL_WINDOW


#include <CGAL/Cartesian.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/IO/Window_stream.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Gt;

typedef Gt::Point_2    Point;
typedef Gt::Segment_2  Segment;
typedef Gt::Triangle_2 Triangle;

typedef CGAL::Constrained_Delaunay_triangulation_2<Gt>  CDT;

typedef CDT::Constraint     Constraint;
typedef CDT::Face_handle    Face_handle;
typedef CDT::Vertex_handle  Vertex_handle;

typedef CGAL::Window_stream Window_stream;

void
any_button(Window_stream &win)
{
    double x, y;
    std::cerr << "Press any button to continue" << std::endl;
    win.read_mouse(x,y);
}

void
draw_constraints(Window_stream &win, std::list<Constraint> & lc)
{
  win << CGAL::RED;
  std::list<Constraint>::iterator cit=lc.begin();
  for( ; cit != lc.end(); ++cit) {
    win << Segment((*cit).first,(*cit).second);
  }
  win << CGAL::BLUE;    
} 

void
input_constraints_from_file(std::list<Constraint> & constraints,
			    std::ifstream& is)
{
  int n;
  is >> n;
    std::cerr << "Reading " << n << " constraints" << std::endl;

  Point p,q;
  for(; n > 0; n--) {
    is >> p >> q;
    constraints.push_back(std::make_pair(p,q));
  }
}

void
draw_connected_component(const Point&  p, 
			 const CDT& ct,
			 Window_stream& win)
{
  Face_handle fh = ct.locate(p);
  std::set<Face_handle> component; 
  std::stack<Face_handle, std::list<Face_handle> > st; 
  // component includes the faces of the connected_component
  // stack includes the faces in component whose neighbors
  // have not yet been looked at

  st.push(fh);
  component.insert(fh);
  while (! st.empty()){
    fh = st.top();
    st.pop();
    for(int i = 0 ; i < 3 ; ++i){
      if ( (! fh->is_constrained(i)) && 
	   component.find(fh->neighbor(i)) == component.end() ) {
	component.insert(fh->neighbor(i));
	st.push(fh->neighbor(i));
      }
    }
  }

  // draw
  //win << CGAL::GREEN;
  win.set_fill_color(leda_yellow);
  std::set<Face_handle>::iterator it;
  for ( it = component.begin(); it != component.end(); it++) {
    if (! ct.is_infinite( *it)) win << ct.triangle( *it);
    else win << ct.segment(*it, (*it)->index(ct.infinite_vertex()));
  }
  return;
}

int
main( )
{
  CGAL::Window_stream win(400,400); // physical window size
  win.init(-1.1, 1.1, -1.1);   // logical window size
  win.display_help_text("help/constrained");
  win << CGAL::BLUE;
  CGAL::cgalize( win);
  win.display();
  win.set_frame_label("Cells of a Constrained Triangulation");

  std::ifstream is("data/poisson");
  std::list<Constraint> lc;
  input_constraints_from_file(lc,is);
  CDT ct(lc);
  assert(ct.is_valid());

  win << ct;
  draw_constraints(win,lc);
  
    std::cerr << "Enter points with the left button" << std::endl;
  std::cerr << "Terminate with right button " << std::endl;
  Point p;
    while(1) {
    double x, y;
    int b = win.get_mouse(x,y);
    bool button_pressed = (b == MOUSE_BUTTON(1)) ||
                          (b == MOUSE_BUTTON(2)) ||
                          (b == MOUSE_BUTTON(3));
    if (button_pressed && b == MOUSE_BUTTON(1)) {
      p = Point(Coord_type(x), Coord_type(y));
      win.clear();
      win << CGAL::BLUE <<ct;
      win << p ;
      draw_connected_component(p, ct, win);
      draw_constraints(win,lc);
    }
    if (button_pressed && b == MOUSE_BUTTON(3)) {
      break;
    }
 
  }
  return 0;
}
  
#endif // CGAL_USE_LEDA
