#ifdef CGAL_USE_LEDA
#include <CGAL/basic.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Delaunay_d.h>
#include <CGAL/IO/Delaunay_d_window_stream.h>
#include <LEDA/d3_window.h>
#include <iostream>

typedef leda_integer RT;
typedef CGAL::Homogeneous_d<RT> Kernel;
typedef CGAL::Delaunay_d<Kernel> Delaunay_d;
typedef Delaunay_d::Point_d Point_d;
typedef Delaunay_d::Simplex_handle Simplex_handle;
typedef Delaunay_d::Vertex_handle Vertex_handle;

static int wsize = 500;

static leda_rat_vector transform(const Point_d& p,int i)
{ return leda_rat_vector(p.homogeneous(0),
                         p.homogeneous(1),
                         p.homogeneous(2)/i,
                         p.homogeneous(p.dimension())); }

enum { ChangeWin = 30, Nearest, Furthest, Exit };

static Delaunay_d* pT;
static leda_window* pW;
static Delaunay_d::Delaunay_voronoi_kind kind = Delaunay_d::NEAREST;
static void redraw_diagram()
{ CGAL::d2_show(*pT,*pW,kind); }

int main()
{
  CGAL::set_pretty_mode ( std::cerr );
  SETDTHREAD(191);
  leda_string startmess = "input points with left mouse button and ";
  startmess += "exit program with right mouse button!";
  leda_window W1(wsize,wsize,"Delaunay Diagram (projected)"),
              W2(wsize,wsize,"Convex Hull (lifted)");
  pW = &W1;
  W1.set_redraw(redraw_diagram);
  int but = 0;
  int scale = 30;
  W1.init(-50,50,-50);
  W1.set_show_coordinates(true);
  W1.set_grid_mode(5);
  leda_list<leda_string> L;
  W1.button("Nearest Site Delaunay",Nearest);
  W1.button("Furthest Site Delaunay",Furthest);
  W1.button("Change to 3D Window",ChangeWin);
  W1.button("Exit",Exit);
  W1.int_item("z-axis scale",scale,1,100,
  "determines the extend of the z-axis");
  W2.init(-50,50,-50);
  W1.display(leda_window::min,leda_window::min); 
  W2.display(wsize+20,0);
  W1.message(startmess);
  double a,b;  // coordinates of a point in the window
  but = W1.read_mouse(a,b);
  W1.del_messages();

  Delaunay_d T(2);  
  pT = &T;
  // we are working in space
  // see dd_delaunay_traits.h for the adaptations
  std::ofstream To("delddemo.log");

  GRAPH< Point_d, int > G;
  leda_node_array<leda_rat_vector> pos(G);
  leda_d3_window W3(W2,G,pos);
  leda_node v;
 
  while (but != MOUSE_BUTTON(3)) {
    // while mouse click is not the right button

    // read the window coordinates into a and b
    W1.clear();       
    RT ia(a), ib(b);

    To << a << "," << b << std::endl; To.flush(); 
    Point_d x(ia,ib); 
    switch (but) {
    case Nearest:
      kind = Delaunay_d::NEAREST;
      break;
    case Furthest:
      kind = Delaunay_d::FURTHEST;
      break;
    case ChangeWin:
    case MOUSE_BUTTON(2):
      forall_nodes(v,G) pos[v] = transform(G[v],scale);
      W3.init(pos);
      CGAL::d2_show(T,W1,kind); 
      while (but != MOUSE_BUTTON(3))
        but = W3.move();
      break;
    case Exit:
      but = MOUSE_BUTTON(3);
      continue;
      break;
    default:
      T.insert(x);
      T.is_valid();
      CGAL::d3_surface_map(T,G);
      pos.init(G);
      forall_nodes(v,G) pos[v] = transform(G[v],scale);
      W3.init(pos);
      W3.draw();
      break;
    }
    CGAL::d2_show(T,W1,kind); 
    but = W1.read_mouse(a,b);  
  }
  return 0;
}

#else
#include <iostream>

int main()
{ 
  std::cout << "this program requires LEDA" << std::endl;
  return 0;
}

#endif // CGAL_USE_LEDA


