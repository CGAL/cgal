// examples/Planar_map/example3.C
// ------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Pm_segment_epsilon_traits.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

typedef CGAL::Cartesian<double>                   Coord_t;
typedef CGAL::Pm_segment_epsilon_traits<Coord_t>  Pmtraits;
//typedef CGAL::Pm_segment_epsilon_traits<Coord_t,1,4>  Pmtraits;
typedef Pmtraits::Point                           Point;
typedef Pmtraits::X_curve                         Curve;
typedef CGAL::Pm_default_dcel<Pmtraits>           Pmdcel;

typedef CGAL::Planar_map_2<Pmdcel,Pmtraits>       Planar_map;
typedef Planar_map::Halfedge_handle	          Halfedge_handle;
typedef Planar_map::Halfedge_const_handle	  Halfedge_const_handle;
typedef Planar_map::Locate_type                   Locate_type;
typedef Planar_map::Ccb_halfedge_circulator       Ccb_halfedge_circulator;
typedef Planar_map::Ccb_halfedge_const_circulator Ccb_halfedge_const_circulator;

// use only for this program
void draw_point_locate(const Point &, const Planar_map&);
// end part of use only for this program

int main()
{
  // creating an instance of Planar_map
  Planar_map pm;

  Curve cv[18];
  int i;

  CGAL::set_ascii_mode(std::cout);

  Point a1(6, 1), a2(1, 3), a3(4, 3), a4(8, 3), a5(11,3) ;
  Point a6(3, 5), a7(9, 5), a8(1, 7), a9(4, 7), a10(8,7) ;
  Point a11(11, 7), a12(6, 9);

	// those curves are about to enter to pm
  cv[0] = Curve(a1, a3);
  cv[1] = Curve(a1, a4);
  cv[2] = Curve(a2, a3);
  cv[3] = Curve(a2, a6);
  cv[4] = Curve(a3, a6);
  cv[5] = Curve(a3, a4);
  cv[6] = Curve(a4, a5);
  cv[7] = Curve(a4, a7);
  cv[8] = Curve(a5, a7);
  cv[9] = Curve(a6, a8);
  cv[10] = Curve(a6, a9);
  cv[11] = Curve(a7, a10);
  cv[12] = Curve(a7, a11);
  cv[13] = Curve(a8, a9);
  cv[14] = Curve(a9, a10);
  cv[15] = Curve(a9, a12);
  cv[16] = Curve(a10, a11);
  cv[17] = Curve(a10, a12);
   
  std::cout << "the curves of the map :" << std::endl; 
  for (i = 0; i < 18; i++)
    std::cout << cv[i] << std::endl;

  std::cout << std::endl;

  // insert the curves to the map
  std::cout << "inserting the curves to the map..." << std::endl;
  for (i = 0; i < 18; i++)
    pm.insert(cv[i]);

  // check the validity of the map
  std::cout << "check map validity... ";
  if (pm.is_valid())
    std::cout << "map valid!" << std::endl;
  else{
    std::cout << "map invalid!" << std::endl;
    return 1;
  }
  std::cout << std::endl;
  
  // draw the map 
  std::cout << "    1  3 4  6  8 9 11		 " << std::endl;
  std::cout << "                                 "
               "                                " << std::endl;
  std::cout << "           a12            9             ";
  std::cout << "a1(" << a1 << ")  a2(" << a2 << ")  a3(" << a3 << ")";
  std::cout << std::endl;        
  std::cout << "           *  *           8             ";
  std::cout << "a4(" << a4 << ")  a5(" << a5 << ")  a6(" << a6 << ")";
  std::cout << std::endl;
  std::cout << "    a8===a9===a10==a11    7             ";
  std::cout << "a7(" << a7 << ")  a8(" << a8 << ")  a9(" << a9 << ")";
  std::cout << std::endl;
  std::cout << "     *  *       *  *      6             ";
  std::cout<< "a10(" << a10 << ")  a11(" << a11 << ")  a12(" << a12 << ")";
  std::cout << std::endl;
  std::cout << "      a6         a7       5" << std::endl;
  std::cout << "     *  *       *  *      4" << std::endl;
  std::cout << "    a2===a3===a4==a5      3" << std::endl;
  std::cout << "          *   *           2" << std::endl;
  std::cout << "           a1             1" << std::endl;

  // try to find some point locations
  draw_point_locate( Point(6,1) , pm );
  draw_point_locate( Point(2,3) , pm );
  draw_point_locate( Point(6,8) , pm );
  std::cout << std::endl ;

  std::cout << "Enter x coordinate for point location:" << std::endl;
  double x;
  std::cin >> x;
  std::cout << "Enter y coordinate for point location:" << std::endl;
  double y;
  std::cin >> y;
        
  draw_point_locate( Point(x,y) ,pm );

  //        std::cout << std::endl << pm;

  return 0;

}

void draw_point_locate(const Point & p, const Planar_map & pm){
  Locate_type lt;	
  Halfedge_const_handle edge  = pm.locate( p,  lt);
  Ccb_halfedge_const_circulator curr,first;
  std::cout << "The location of point " << p << " is of type " ;
  switch ( lt ) {
  case Planar_map::VERTEX	:	
    std::cout << "VERTEX" << std::endl;
    std::cout << "The vertex is : (" << edge->target()->point() << ")" 
	      << std::endl;
    break;
  case Planar_map::UNBOUNDED_VERTEX	:	
    std::cout << "UNBOUNDED_VERTEX" << std::endl;
    std::cout << "The vertex is : (" << edge->target()->point() << ")" 
	      << std::endl;
    break;	case Planar_map::EDGE	:
      std::cout << "EDGE" << std::endl;
    std::cout << "The edge is : {(" << edge->source()->point() ;
    std::cout<< ")->(" << edge->target()->point() << ")}" 
	     << std::endl;
    break;
  case Planar_map::UNBOUNDED_EDGE	:
    std::cout << "UNBOUNDED_EDGE" << std::endl;
    std::cout << "The edge is : {(" << edge->source()->point() ;
    std::cout<< ")->(" << edge->target()->point() << ")}" 
	     << std::endl;
    break;
  case Planar_map::FACE	:
    first = Ccb_halfedge_const_circulator(edge);
    curr=first;
    std::cout << "FACE" << std::endl;
    std::cout << "The face is : " ;
    std::cout << "[" ;
    do {
      std::cout << "(" << curr->target()->point() << ")-" ;
      ++curr;
    }while (curr!=first) ;
    std::cout << "(" << curr->target()->point() << ")]" << std::endl;
    break;
  case Planar_map::UNBOUNDED_FACE	:
    std::cout << "UNBOUNDED_FACE" << std::endl;
    break;
  }

}

