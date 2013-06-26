// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ 
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#include <CGAL/Random.h>
#include <cassert>

template <class CK>
void _test_circle_bbox(CK ck)
{
  typedef typename CK::Circle_2                    Circle_2;
  typedef typename CK::Circular_arc_2              Circular_arc_2;
  typedef typename CK::Point_2                     Point_2;
  typedef typename CK::Intersect_2   Intersect_2;

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
    void _test_circular_arc_bbox(CK)
  {
    typedef typename CK::Circle_2                    Circle_2;
    typedef typename CK::Circular_arc_2              Circular_arc_2;
    typedef typename CK::Point_2                     Point_2;
    typedef typename CK::Line_2                      Line_2;

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
    std::pair< Circular_arc_point_2, unsigned int> aux;
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
    assert(!ck.has_on_2_object()(arc_top,
				line_vertical.source()));
    assert(ck.has_on_2_object()(arc_top,
				arc_top.source()));
    assert(ck.has_on_2_object()(arc_top,
				arc_top.target()));
    assert(ck.has_on_2_object()(arc_top,
				point_top));
    assert(!ck.has_on_2_object()(arc_top,
				point_down));
    assert(ck.has_on_2_object()(arc_down,
				point_down));
    assert(!ck.has_on_2_object()(arc_down,
				point_top));
  }

template <class CK>
    void _test_circular_arc_point_bbox(CK ck)
  {
    typedef typename CK::Circle_2                    Circle_2;
    typedef typename CK::Circular_arc_2              Circular_arc_2;
    typedef typename CK::Point_2                     Point_2;
    typedef typename CK::Circular_arc_point_2     Circular_arc_point_2;
    typedef typename CK::Intersect_2   Intersect_2;


    
    int test_suite_cases[3][2][3] = {{{-7,-8,113},{9,9,162}},
                                     {{0,0,1},{1,2,4}},
                                     {{0,0,1},{1,0,1}}};
    int n_cases = 3;
    
    for(int i=0; i<n_cases; i++) {
      Point_2 center_circ(test_suite_cases[i][0][0],test_suite_cases[i][0][1]);
      Circle_2 circ(center_circ, test_suite_cases[i][0][2]);
      Circular_arc_2 arc(circ);
      Point_2 center_circ2(test_suite_cases[i][1][0],test_suite_cases[i][1][1]);
      Circle_2 circ2(center_circ2, test_suite_cases[i][1][2]);
      Circular_arc_2 arc2(circ2);

      Intersect_2 theConstruct_intersect_2 
        = ck.intersect_2_object();
      std::vector< CGAL::Object > 
        vector_for_intersection_1;
      theConstruct_intersect_2(arc, 
	  		     arc2,
			     std::back_inserter(vector_for_intersection_1));
      std::pair< Circular_arc_point_2, unsigned int> aux;
      assign(aux, vector_for_intersection_1[0]);
      CGAL::Bbox_2 box1 = aux.first.bbox(); 
      assert(typename CK::FT(box1.xmin()) <= aux.first.x());
      assert(typename CK::FT(box1.xmax()) >= aux.first.x());
      assert(typename CK::FT(box1.ymin()) <= aux.first.y());
      assert(typename CK::FT(box1.ymax()) >= aux.first.y());
      std::cout << "Ok" << std::endl;
      assign(aux, vector_for_intersection_1[1]);
      CGAL::Bbox_2 box2 = aux.first.bbox(); 
      assert(typename CK::FT(box2.xmin()) <= aux.first.x());
      assert(typename CK::FT(box2.xmax()) >= aux.first.x());
      assert(typename CK::FT(box2.ymin()) <= aux.first.y());
      assert(typename CK::FT(box2.ymax()) >= aux.first.y());
      std::cout << "Ok" << std::endl;
    }
  }
  
  
