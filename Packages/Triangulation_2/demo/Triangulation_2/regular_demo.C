#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) && !defined(CGAL_USE_CGAL_WINDOW )
int main()
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else
#include <CGAL/Cartesian.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/IO/Window_stream.h>


typedef CGAL::Cartesian<double> Rp;
typedef double W;
typedef CGAL::Regular_triangulation_euclidean_traits_2<Rp,W>  Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Regular_triangulation_face_base_2<Gt> Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb > Tds;
typedef CGAL::Regular_triangulation_2<Gt, Tds> Regular_triangulation;
typedef CGAL::Window_stream  Window_stream;

void
any_button(Window_stream &W)
{
    double x, y;
    std::cerr << "Press any button to continue" << std::endl;
    W.read_mouse(x,y);
}


int main()
{
  Regular_triangulation rt;
  //std::ifstream in("regular.cin.plante");

  Window_stream W(400, 400); // physical window size
  W.init(-1, 10, -1);   // logical window size
  CGAL::cgalize( W);
  W.display();

  Gt::Weighted_point wp;
  std::list< Gt::Weighted_point> lp;
  std::list< Gt::Weighted_point>::iterator lpit;
  std::cout << " input weighted points with left mouse button" 
	    <<    std::endl;
  std::cout << " end with Ctrl C " << std::endl;
  while(W >> wp){
     std::cout << wp << std::endl;
     rt.insert(wp);
     lp.push_back(wp);
     rt.is_valid();
     W.clear();
     W << CGAL::GREEN;
     for ( lpit=lp.begin(); lpit != lp.end(); lpit++){
       W << *lpit;
     }
     W << CGAL::BLUE << rt;
     W << CGAL::RED ;
     rt.draw_dual(W);
     }
   return 0;	
}

#endif // CGAL_USE_LEDA || CGAL_USE_CGAL_WINDOW
