#include <CGAL/Random.h>

template <class AK>
void _test_solve(AK ak)
{
  typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
  typename AK::Solve theSolve =
    ak.solve_object();

  //Polynomial_for_circles_2_2
  typename AK::Construct_polynomial_circle_2_2 theConstruct_2_2 =
    ak.construct_polynomial_circle_2_2_object();

  std::vector< std::pair<Root_for_circles_2_2, size_t> > res1;
  theSolve(theConstruct_2_2(5, 5, 25),
	   theConstruct_2_2(0, 5, 100),
	   std::back_inserter(res1));
  assert(res1.size() == 1);
  assert(res1[0].second == 2u);
  assert(res1[0].first == Root_for_circles_2_2(10, 5));

  std::vector< std::pair<Root_for_circles_2_2, size_t> > res2;
  theSolve(theConstruct_2_2(-5, 5, 25),
	   theConstruct_2_2(0, 5, 100),
	   std::back_inserter(res2));
  assert(res2.size() == 1);
  assert(res2[0].second == 2u);
  assert(res2[0].first == Root_for_circles_2_2(-10, 5));

  std::vector< std::pair<Root_for_circles_2_2, size_t> > res3;
  theSolve(theConstruct_2_2(0, 5, 25),
	   theConstruct_2_2(0, 0, 100),
	   std::back_inserter(res3));
  assert(res3.size() == 1);
  assert(res3[0].second == 2u);
  assert(res3[0].first == Root_for_circles_2_2(0, 10));
 
  std::vector< std::pair<Root_for_circles_2_2, size_t> > res4;
  theSolve(theConstruct_2_2(0, -5, 25),
	   theConstruct_2_2(0, 0, 100),
	   std::back_inserter(res4));
  assert(res4.size() == 1);
  assert(res4[0].second == 2u);
  assert(res4[0].first == Root_for_circles_2_2(0, -10));
 
  std::vector< std::pair<Root_for_circles_2_2, size_t> > res5;
  theSolve(theConstruct_2_2(-5, 5, 25),
	   theConstruct_2_2(0, 0, 25),
	   std::back_inserter(res5));
  assert(res5.size() == 2);
  assert(res5[0].second == 1u);
  assert(res5[0].first == Root_for_circles_2_2(-5, 0));
  assert(res5[1].second == 1u);
  assert(res5[1].first == Root_for_circles_2_2(0, 5));
 
  //Polynomial_1_2 Polynomial_for_circles_2_2
  typename AK::Construct_polynomial_1_2 theConstruct_1_2 =
    ak.construct_polynomial_1_2_object();

  std::vector< std::pair<Root_for_circles_2_2, size_t> > res6;
  theSolve(theConstruct_1_2(0, 1, -5),
	   theConstruct_2_2(0, 5, 100),
	   std::back_inserter(res6));
  assert(res6.size() == 2);
  assert(res6[0].second == 1u);
  assert(res6[0].first == Root_for_circles_2_2(-10, 5));
  assert(res6[1].second == 1u);
  assert(res6[1].first == Root_for_circles_2_2(10, 5));

  std::vector< std::pair<Root_for_circles_2_2, size_t> > res7;
  theSolve(theConstruct_1_2(1, 0, -5),
	   theConstruct_2_2(5, 5, 100),
	   std::back_inserter(res7));
  assert(res7.size() == 2);
  assert(res7[0].second == 1u);
  assert(res7[0].first == Root_for_circles_2_2(5, -5));
  assert(res7[1].second == 1u);
  assert(res7[1].first == Root_for_circles_2_2(5, 15));

  std::vector< std::pair<Root_for_circles_2_2, size_t> > res8;
  theSolve(theConstruct_1_2(1, 0, 5),
	   theConstruct_2_2(5, 5, 100),
	   std::back_inserter(res8));
  assert(res8.size() == 1);
  assert(res8[0].second == 2u);
  assert(res8[0].first == Root_for_circles_2_2(-5, 5));

  std::vector< std::pair<Root_for_circles_2_2, size_t> > res9;
  theSolve(theConstruct_1_2(1, 0, -15),
	   theConstruct_2_2(5, 5, 100),
	   std::back_inserter(res9));
  assert(res9.size() == 1);
  assert(res9[0].second == 2u);
  assert(res9[0].first == Root_for_circles_2_2(15, 5));

  std::vector< std::pair<Root_for_circles_2_2, size_t> > res10;
  theSolve(theConstruct_1_2(0, 1, -15),
	   theConstruct_2_2(5, 5, 100),
	   std::back_inserter(res10));
  assert(res10.size() == 1);
  assert(res10[0].second == 2u);
  assert(res10[0].first == Root_for_circles_2_2(5, 15));
  
  std::vector< std::pair<Root_for_circles_2_2, size_t> > res11;
  theSolve(theConstruct_1_2(0, 1, 5),
	   theConstruct_2_2(5, 5, 100),
	   std::back_inserter(res11));
  assert(res11.size() == 1);
  assert(res11[0].second == 2u);
  assert(res11[0].first == Root_for_circles_2_2(5, -5));
  
}

template <class AK>
void _test_sign_at(AK ak)
{
  typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
  typename AK::Sign_at theSigh_at =
    ak.sign_at_object();
  
  //Polynomial_for_circles_2_2
  typename AK::Construct_polynomial_circle_2_2 theConstruct_2_2 =
    ak.construct_polynomial_circle_2_2_object();

  assert(theSigh_at(theConstruct_2_2(5, 5, 100),
		    Root_for_circles_2_2(-5,5)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_2(5, 5, 100),
		    Root_for_circles_2_2(15,5)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_2(5, 5, 100),
		    Root_for_circles_2_2(5,15)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_2(5, 5, 100),
		    Root_for_circles_2_2(5,-5)) == CGAL::ZERO);
  
  assert(theSigh_at(theConstruct_2_2(5, 5, 100),
		    Root_for_circles_2_2(5,5)) != CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_2(5, 5, 100),
		    Root_for_circles_2_2(5,16)) != CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_2(5, 5, 100),
		    Root_for_circles_2_2(5,-6)) != CGAL::ZERO);

  //Polynomial_1_2
  typename AK::Construct_polynomial_1_2 theConstruct_1_2 =
    ak.construct_polynomial_1_2_object();

  assert(theSigh_at(theConstruct_1_2(1, 0, -5),
		    Root_for_circles_2_2(5,-6)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_1_2(1, 0, -5),
		    Root_for_circles_2_2(6,-6)) != CGAL::ZERO);
  assert(theSigh_at(theConstruct_1_2(0, 1, -5),
		    Root_for_circles_2_2(5,-6)) != CGAL::ZERO);
  assert(theSigh_at(theConstruct_1_2(0, 1, -5),
		    Root_for_circles_2_2(5, 5)) == CGAL::ZERO);
  
}

template <class AK>
void _test_critical_points(AK ak)
{
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;

  typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
  typename AK::Construct_polynomial_circle_2_2 theConstruct_2_2 =
    ak.construct_polynomial_circle_2_2_object();
  typename AK::X_critical_points theX_critical_points =
    ak.x_critical_points_object();
  
  typename AK::Y_critical_points theY_critical_points =
    ak.y_critical_points_object();
  
  for(int i = 0; i < 20; i++){
    int x = theRandom.get_int(random_min,random_max);
    int y = theRandom.get_int(random_min,random_max);
    int r = theRandom.get_int(1,random_max);
  
    assert(theX_critical_points(theConstruct_2_2(x,y,r*r),true)
	   == Root_for_circles_2_2(x - r, y));
    assert(theX_critical_points(theConstruct_2_2(x,y,r*r),false)
	   == Root_for_circles_2_2(x + r, y));
    assert(theY_critical_points(theConstruct_2_2(x,y,r*r),true)
	   == Root_for_circles_2_2(x, y - r));
    assert(theY_critical_points(theConstruct_2_2(x,y,r*r),false)
	   == Root_for_circles_2_2(x, y + r));
  }
  
}

template <class AK>
void _test_compare_Root_for_circles(AK ak)
{
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;

  typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;

  typename AK::Compare_x theCompare_x =
    ak.compare_x_object();
  
  typename AK::Compare_y theCompare_y =
    ak.compare_y_object();
  
  typename AK::Compare_xy theCompare_xy =
    ak.compare_xy_object();

  for (int i = 0; i < 20; i++){
    Root_for_circles_2_2 r1(theRandom.get_int(random_min,random_max),
			    theRandom.get_int(random_min,random_max));
    Root_for_circles_2_2 r2(theRandom.get_int(random_min,random_max),
			    theRandom.get_int(random_min,random_max));
    if(r1.x() > r2.x()){
      assert(theCompare_x(r1, r2) == CGAL::LARGER);
      assert(theCompare_xy(r1, r2) == CGAL::LARGER);
    }
    else if(r1.x() == r2.x()){
      assert(theCompare_x(r1, r2) == CGAL::EQUAL);
      if(r1.y() < r2.y()){
	assert(theCompare_y(r1, r2) == CGAL::SMALLER);
	assert(theCompare_xy(r1, r2) == CGAL::SMALLER);
      }
      else if(r1.y() > r2.y()){
	assert(theCompare_y(r1, r2) == CGAL::LARGER);
	assert(theCompare_xy(r1, r2) == CGAL::LARGER);
      }
      else {
	assert(theCompare_y(r1, r2) == CGAL::EQUAL);
	assert(theCompare_xy(r1, r2) == CGAL::EQUAL);
      }
    }
    else {
      assert(theCompare_x(r1, r2) == CGAL::SMALLER);
      assert(theCompare_xy(r1, r2) == CGAL::SMALLER);
    }
    if(r1.y() > r2.y())
      assert(theCompare_y(r1, r2) == CGAL::LARGER);
    else if(r1.y() < r2.y())
      assert(theCompare_y(r1, r2) == CGAL::SMALLER);
    else
      assert(theCompare_y(r1, r2) == CGAL::EQUAL);
  }
}
