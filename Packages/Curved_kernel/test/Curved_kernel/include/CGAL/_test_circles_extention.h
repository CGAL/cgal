#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circular_kernel.h>
#include <CGAL/Circular_arc_traits.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Curved_kernel/function_objects_polynomial_circular.h>
#include <CGAL/Curved_kernel/Circular_arc_2.h>

#include <CGAL/NT_extensions_Root_of/CGAL_Quotient.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Gmpq.h>


#include <CGAL/Random.h>

template <class CK>
void _test_circle_bbox(CK ck)
{
  typedef typename CK::Circle_2                    Circle_2;
  typedef typename CK::Circular_arc_2              Circular_arc_2;
  typedef typename CK::Point_2                     Point_2;
  typedef typename CK::Line_2                      Line_2;
  typedef typename CK::Circular_arc_point_2     Circular_arc_point_2;
  typedef typename CK::Construct_circle_2          Construct_circle_2;
  typedef typename CK::Intersect_2   Intersect_2;
  typedef typename CK::Make_x_monotone_2           Make_x_monotone_2;
  typedef typename CK::Split_2                     Split_2;  
  typedef typename CK::Get_equation                Get_equation;
  typedef typename CK::Polynomial_for_circles_2_2  Polynomial_for_circles_2_2;
  typedef typename CK::Compare_xy_2                Compare_xy_2;

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 5;
  int random_min = -5;
  
  int x;
  int y;
  int r;

  for (int i = 0; i < 50 ; i++){
    x = theRandom.get_int(random_min,random_max);
    y = theRandom.get_int(random_min,random_max);
    r = theRandom.get_int(1,random_max);
    Point_2 center_circ1(x,y);
    Circle_2 circ1(center_circ1, r);
    x = theRandom.get_int(random_min,random_max);
    y = theRandom.get_int(random_min,random_max);
    r = theRandom.get_int(1,random_max);
    Point_2 center_circ2(x,y);
    Circle_2 circ2(center_circ2, r);
    Circular_arc_2 arc1(circ1);
    Circular_arc_2 arc2(circ2);
    CGAL::Bbox_2 box1 = arc1.bbox();
    CGAL::Bbox_2 box2 = arc2.bbox();
    bool box_overlap = do_overlap(box1, box2);
    
    Intersect_2 theConstruct_intersect_2 
      = ck.intersect_2_object();
    std::vector< CGAL::Object > 
      vector_for_intersection_1;
    theConstruct_intersect_2(arc1, 
			     arc2,
			     std::back_inserter(vector_for_intersection_1));
    if(vector_for_intersection_1.size() > 0){
      std::cout << " intersection " << std::endl;
      assert(box_overlap);
    }
    if(!box_overlap){
      std::cout << " box_overlap " << std::endl;
      assert(vector_for_intersection_1.size() == 0);
    }
  }
}
  template <class CK>
    void _test_circular_arc_bbox(CK ck)
  {
    typedef typename CK::Circle_2                    Circle_2;
    typedef typename CK::Circular_arc_2              Circular_arc_2;
    typedef typename CK::Point_2                     Point_2;
    typedef typename CK::Line_2                      Line_2;
    typedef typename CK::Circular_arc_point_2     Circular_arc_point_2;
    typedef typename CK::Construct_circle_2          Construct_circle_2;
    typedef typename CK::Intersect_2   Intersect_2;
    typedef typename CK::Make_x_monotone_2           Make_x_monotone_2;
    typedef typename CK::Split_2                     Split_2;  
    typedef typename CK::Get_equation                Get_equation;
    typedef typename CK::Polynomial_for_circles_2_2  Polynomial_for_circles_2_2;
    typedef typename CK::Compare_xy_2                Compare_xy_2;

    CGAL::Random generatorOfgenerator;
    int random_seed = generatorOfgenerator.get_int(0, 123456);
    std::cout << "random_seed = " << random_seed << std::endl;
    CGAL::Random theRandom(random_seed);
    int random_max = 127;
    int random_min = -127;

    for (int i = 0; i < 50 ; i++){
      Point_2 center_random(theRandom.get_int(random_min,random_max),
			   theRandom.get_int(random_min,random_max));
      Circle_2 circle_random(center_random,
		     theRandom.get_int(1,random_max));
      int x_random1, y_random1;
      int x_random2, y_random2;
      do{
	x_random1 = theRandom.get_int(random_min, random_max);
	y_random1 = theRandom.get_int(random_min, random_max);
      }while(x_random1 == 0 && y_random1 ==0);
    
      do{
	x_random2 = theRandom.get_int(random_min, random_max);
	y_random2 = theRandom.get_int(random_min, random_max);
      }while(x_random2 == 0 && y_random2 ==0);
    
      Line_2 line_random_1(center_random,
			   Point_2(center_random.x() +
				   x_random1,
				   center_random.y() + 
				   y_random1));
      Line_2 line_random_2(center_random,
			   Point_2(center_random.x() +
				   x_random2,
				   center_random.y() + 
				   y_random2));
      Circular_arc_2 arc_random(circle_random,
				line_random_1, theRandom.get_bool(),
				line_random_2, theRandom.get_bool());
      
      CGAL::Bbox_2 box1 = arc_random.bbox();
      assert(typename CK::FT(box1.xmin()) <= arc_random.source().x());
      assert(typename CK::FT(box1.xmin()) <= arc_random.target().x());
      assert(typename CK::FT(box1.xmax()) >= arc_random.source().x());
      assert(typename CK::FT(box1.xmax()) >= arc_random.target().x());
      assert(typename CK::FT(box1.ymin()) <= arc_random.source().y());
      assert(typename CK::FT(box1.ymin()) <= arc_random.target().y());
      assert(typename CK::FT(box1.ymax()) >= arc_random.source().y());
      assert(typename CK::FT(box1.ymax()) >= arc_random.target().y());
//      assert(((typename CK::FT(box1.xmin()) - arc_random.center().x())
//	     *(typename CK::FT(box1.xmin()) - arc_random.center().x()))
//	     <= arc_random.supporting_circle().squared_radius());
      
    }
  }
    
template <class CK>
    void _test_has_on(CK ck)
  {
    typedef typename CK::Circle_2                    Circle_2;
    typedef typename CK::Circular_arc_2              Circular_arc_2;
    typedef typename CK::Point_2                     Point_2;
    typedef typename CK::Circular_arc_point_2     Circular_arc_point_2;
    typedef typename CK::Intersect_2   Intersect_2;
    typedef typename CK::Make_x_monotone_2           Make_x_monotone_2;
    typedef typename CK::Line_arc_2              Line_arc_2;
    
    Point_2 center_circ(0,0);
    Circle_2 circ(center_circ, 100);
    Circular_arc_2 arc(circ);
    Line_arc_2 line_vertical(Point_2(0, 15),
			     Point_2(0, -15));
    Intersect_2 theConstruct_intersect_2 
      = ck.intersect_2_object();
    std::vector< CGAL::Object > 
      vector_for_intersection_1;
    theConstruct_intersect_2(arc, 
			     line_vertical,
			     std::back_inserter(vector_for_intersection_1));
    Circular_arc_point_2 point_top;
    Circular_arc_point_2 point_down;
    std::pair< Circular_arc_point_2, std::size_t > aux;
    assign(aux, vector_for_intersection_1[0]);
    point_down = aux.first;
    assign(aux, vector_for_intersection_1[1]);
    point_top = aux.first;
    Make_x_monotone_2 theMake_x_monotone = ck.make_x_monotone_2_object();
    std::vector< CGAL::Object > outputIterator1;
    theMake_x_monotone(arc,
		       std::back_inserter(outputIterator1));
    Circular_arc_2 arc_top;
    Circular_arc_2 arc_down;
    assign(arc_top,outputIterator1[1]);
    assign(arc_down, outputIterator1[0]);
    assert(!CGAL::CircularFunctors::has_on<CK>(arc_top,
					    line_vertical.source()));
    assert(CGAL::CircularFunctors::has_on<CK>(arc_top,
					      arc_top.source()));
    assert(CGAL::CircularFunctors::has_on<CK>(arc_top,
					      arc_top.target()));
    assert(CGAL::CircularFunctors::has_on<CK>(arc_top,
					      point_top));
    assert(!CGAL::CircularFunctors::has_on<CK>(arc_top,
					      point_down));
    assert(CGAL::CircularFunctors::has_on<CK>(arc_down,
					      point_down));
    assert(!CGAL::CircularFunctors::has_on<CK>(arc_down,
					      point_top));
  }
  
