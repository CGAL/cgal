#include <CGAL/basic.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>
#include <vector>

template <typename T>
void _test_rc_random_1()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	typedef typename D3::Finite_vertices_iterator  Finite_vertices_iterator;	
	typedef typename D3::Geom_traits::Sphere_3     Sphere_3;
	typedef Point                                  Point_3;

	std::cout << "_test_rc_random_1" << std::endl;
		
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	CGAL::Random_points_in_cube_3<Point> cube_insert(1.0);
	
  const int number_of_points = (1<<14); 

	for(int i=0; i<number_of_points; i++) to_insert.push_back(*cube_insert++);
	dt.insert(to_insert.begin(), to_insert.end());

	for(Finite_vertices_iterator fit = dt.finite_vertices_begin(); 
	fit != dt.finite_vertices_end(); fit++) {
		Sphere_3 s1 = Sphere_3(Point_3(0.25, 0.25, 0.25), 0.01);
		Sphere_3 s2 = Sphere_3(Point_3(0.75, 0.75, 0.75), 0.01);
		if(s1.has_on_unbounded_side(fit->point()) &&
		   s2.has_on_unbounded_side(fit->point())) to_remove.push_back(fit);
  }

	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_random_2()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	typedef typename D3::Finite_vertices_iterator  Finite_vertices_iterator;	
	typedef typename D3::Geom_traits::Sphere_3     Sphere_3;
	typedef Point                                  Point_3;
		
	std::cout << "_test_rc_random_2" << std::endl;		
		
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	CGAL::Random_points_in_cube_3<Point> cube_insert(1.0);
	
  const int number_of_points = (1<<14); 

	for(int i=0; i<number_of_points; i++) to_insert.push_back(*cube_insert++);
	dt.insert(to_insert.begin(), to_insert.end());

	for(Finite_vertices_iterator fit = dt.finite_vertices_begin(); 
	fit != dt.finite_vertices_end(); fit++) {
		Sphere_3 s1 = Sphere_3(Point_3(0.25, 0.25, 0.25), 0.01);
		Sphere_3 s2 = Sphere_3(Point_3(0.75, 0.75, 0.75), 0.01);
		if(s1.has_on_bounded_side(fit->point())) to_remove.push_back(fit);
		if(s2.has_on_bounded_side(fit->point())) to_remove.push_back(fit);
  }

	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_random_3()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	typedef typename D3::Finite_vertices_iterator  Finite_vertices_iterator;	
	
	std::cout << "_test_rc_random_3" << std::endl;	
	
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	CGAL::Random random;
	CGAL::Random_points_in_cube_3<Point> cube_insert(1.0);
	
  const int number_of_points = (1<<14); 

	for(int i=0; i<number_of_points; i++) to_insert.push_back(*cube_insert++);
	dt.insert(to_insert.begin(), to_insert.end());

	for(Finite_vertices_iterator fit = dt.finite_vertices_begin(); 
	fit != dt.finite_vertices_end(); fit++) {
		if(random.get_double() < 0.05) to_remove.push_back(fit);
  }

	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_random_4()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	typedef typename D3::Finite_vertices_iterator  Finite_vertices_iterator;	
	
	std::cout << "_test_rc_random_4" << std::endl;	
	
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	CGAL::Random random;
	CGAL::Random_points_in_cube_3<Point> cube_insert(1.0);
	
  const int number_of_points = (1<<14); 

	for(int i=0; i<number_of_points; i++) to_insert.push_back(*cube_insert++);
	dt.insert(to_insert.begin(), to_insert.end());

	for(Finite_vertices_iterator fit = dt.finite_vertices_begin(); 
	fit != dt.finite_vertices_end(); fit++) {
		if(random.get_double() < 0.9) to_remove.push_back(fit);
  }

	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_trivial()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	
	std::cout << "_test_rc_trivial" << std::endl;	
	
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	Vertex_handle a1, a2, a3, a4, a5;
		
	a1 = dt.insert(Point(0,0,0));
	a2 = dt.insert(Point(0,0,1));
	a3 = dt.insert(Point(0,1,1));
	a4 = dt.insert(Point(0,1,0));
	a5 = dt.insert(Point(1,0,0));
		
	to_remove.push_back(a5);
				
	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_specific_cases_1()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	
	std::cout << "_test_rc_specific_cases_1" << std::endl;		
	
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	Vertex_handle a1, a2, a3, a4, a5, a6, a7, a8;
	Vertex_handle b1, b2, b3, b4, b5, b6, b7, b8;
		
	a1 = dt.insert(Point(0,0,0));
	a2 = dt.insert(Point(0,0,1));
	a3 = dt.insert(Point(0,1,1));
	a4 = dt.insert(Point(0,1,0));
	a5 = dt.insert(Point(1,0,0));
	a6 = dt.insert(Point(1,0,1));
	a7 = dt.insert(Point(1,1,1));
	a8 = dt.insert(Point(1,1,0));
	
	b1 = dt.insert(Point(0.25,0.25,0.25));
	b2 = dt.insert(Point(0.25,0.25,0.75));
	b3 = dt.insert(Point(0.25,0.75,0.75));
	b4 = dt.insert(Point(0.25,0.75,0.25));
	b5 = dt.insert(Point(0.75,0.25,0.25));
	b6 = dt.insert(Point(0.75,0.25,0.75));
	b7 = dt.insert(Point(0.75,0.75,0.75));
	b8 = dt.insert(Point(0.75,0.75,0.25));
	
	to_remove.push_back(a1);
	to_remove.push_back(a2);
	to_remove.push_back(a3);
	to_remove.push_back(a4);
	to_remove.push_back(a5);
	to_remove.push_back(a6);
	to_remove.push_back(a7);
	to_remove.push_back(a8);
				
	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_specific_cases_2()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	
	std::cout << "_test_rc_specific_cases_2" << std::endl;	
	
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	Vertex_handle a1, a2, a3, a4, a5, a6, a7, a8;
	Vertex_handle b1, b2, b3, b4, b5, b6, b7, b8;
		
	a1 = dt.insert(Point(0,0,0));
	a2 = dt.insert(Point(0,0,1));
	a3 = dt.insert(Point(0,1,1));
	a4 = dt.insert(Point(0,1,0));
	a5 = dt.insert(Point(1,0,0));
	a6 = dt.insert(Point(1,0,1));
	a7 = dt.insert(Point(1,1,1));
	a8 = dt.insert(Point(1,1,0));
	
	b1 = dt.insert(Point(0.25,0.25,0.25));
	b2 = dt.insert(Point(0.25,0.25,0.75));
	b3 = dt.insert(Point(0.25,0.75,0.75));
	b4 = dt.insert(Point(0.25,0.75,0.25));
	b5 = dt.insert(Point(0.75,0.25,0.25));
	b6 = dt.insert(Point(0.75,0.25,0.75));
	b7 = dt.insert(Point(0.75,0.75,0.75));
	b8 = dt.insert(Point(0.75,0.75,0.25));
	
	to_remove.push_back(a1);
	to_remove.push_back(a2);
	to_remove.push_back(a3);
	to_remove.push_back(a4);
	to_remove.push_back(a5);
	to_remove.push_back(a6);
				
	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_specific_cases_3()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;

	std::cout << "_test_rc_specific_cases_3" << std::endl;
	
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	Vertex_handle a1, a2, a3, a4, a5, a6, a7, a8;
	Vertex_handle b1, b2, b3, b4, b5, b6, b7, b8;
		
	a1 = dt.insert(Point(0,0,0));
	a2 = dt.insert(Point(0,0,1));
	a3 = dt.insert(Point(0,1,1));
	a4 = dt.insert(Point(0,1,0));
	a5 = dt.insert(Point(1,0,0));
	a6 = dt.insert(Point(1,0,1));
	a7 = dt.insert(Point(1,1,1));
	a8 = dt.insert(Point(1,1,0));
	
	b1 = dt.insert(Point(0.25,0.25,0.25));
	b2 = dt.insert(Point(0.25,0.25,0.75));
	b3 = dt.insert(Point(0.25,0.75,0.75));
	b4 = dt.insert(Point(0.25,0.75,0.25));
	b5 = dt.insert(Point(0.75,0.25,0.25));
	b6 = dt.insert(Point(0.75,0.25,0.75));
	b7 = dt.insert(Point(0.75,0.75,0.75));
	b8 = dt.insert(Point(0.75,0.75,0.25));
	
	to_remove.push_back(b1);
	to_remove.push_back(b2);
	to_remove.push_back(b3);
	to_remove.push_back(a4);
	to_remove.push_back(a5);
	to_remove.push_back(a6);
				
	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_specific_cases_4()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	
	std::cout << "_test_rc_specific_cases_4" << std::endl;	
	
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	Vertex_handle a1, a2, a3, a4, a5, a6, a7, a8;
	Vertex_handle b1, b2, b3, b4, b5, b6, b7, b8;
		
	a1 = dt.insert(Point(0,0,0));
	a2 = dt.insert(Point(0,0,1));
	a3 = dt.insert(Point(0,1,1));
	a4 = dt.insert(Point(0,1,0));
	a5 = dt.insert(Point(1,0,0));
	a6 = dt.insert(Point(1,0,1));
	a7 = dt.insert(Point(1,1,1));
	a8 = dt.insert(Point(1,1,0));
	
	b1 = dt.insert(Point(0.25,0.25,0.25));
	b2 = dt.insert(Point(0.25,0.25,0.75));
	b3 = dt.insert(Point(0.25,0.75,0.75));
	b4 = dt.insert(Point(0.25,0.75,0.25));
	b5 = dt.insert(Point(0.75,0.25,0.25));
	b6 = dt.insert(Point(0.75,0.25,0.75));
	b7 = dt.insert(Point(0.75,0.75,0.75));
	b8 = dt.insert(Point(0.75,0.75,0.25));
	
	to_remove.push_back(b1);
				
	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_specific_cases_5()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	
	std::cout << "_test_rc_specific_cases_5" << std::endl;	
	
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	dt.insert(Point(-1000,-1000,-1000));
	dt.insert(Point(-1000,-1000,1000));
	dt.insert(Point(-1000,1000,1000));
	dt.insert(Point(-1000,1000,-1000));
	dt.insert(Point(1000,-1000,-1000));
	dt.insert(Point(1000,-1000,1000));
	dt.insert(Point(1000,1000,1000));
	dt.insert(Point(1000,1000,-1000));
	
	Vertex_handle a1, a2, a3, a4, a5, a6, a7, a8;
	Vertex_handle b1, b2, b3, b4, b5, b6, b7, b8;
	Vertex_handle c1;
		
	a1 = dt.insert(Point(0,0,0));
	a2 = dt.insert(Point(0,0,1));
	a3 = dt.insert(Point(0,1,1));
	a4 = dt.insert(Point(0,1,0));
	a5 = dt.insert(Point(1,0,0));
	a6 = dt.insert(Point(1,0,1));
	a7 = dt.insert(Point(1,1,1));
	a8 = dt.insert(Point(1,1,0));
	
	b1 = dt.insert(Point(0.25,0.25,0.25));
	b2 = dt.insert(Point(0.25,0.25,0.75));
	b3 = dt.insert(Point(0.25,0.75,0.75));
	b4 = dt.insert(Point(0.25,0.75,0.25));
	b5 = dt.insert(Point(0.75,0.25,0.25));
	b6 = dt.insert(Point(0.75,0.25,0.75));
	b7 = dt.insert(Point(0.75,0.75,0.75));
	b8 = dt.insert(Point(0.75,0.75,0.25));
	
	c1 = dt.insert(Point(0.5,0.51,0.52));
	
	to_remove.push_back(b1);
	to_remove.push_back(b2);
	to_remove.push_back(b3);
	to_remove.push_back(b4);
	to_remove.push_back(b5);
	to_remove.push_back(b6);
	to_remove.push_back(b7);
	to_remove.push_back(b8);
	to_remove.push_back(a2);
	to_remove.push_back(a4);
	to_remove.push_back(a6);
	to_remove.push_back(a8);	
		
	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_specific_cases_6()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;

	
	std::cout << "_test_rc_specific_cases_6" << std::endl;	
	
	std::size_t s1, s2, ns;
	
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	Vertex_handle a1, a2, a3, a4, a5, a6, a7, a8;
  Vertex_handle v1, v2;
	
	dt.insert(Point(-1000,-1000,-1000));
	dt.insert(Point(-1000,-1000,1000));
	dt.insert(Point(-1000,1000,1000));
	dt.insert(Point(-1000,1000,-1000));
	dt.insert(Point(1000,-1000,-1000));
	dt.insert(Point(1000,-1000,1000));
	dt.insert(Point(1000,1000,1000));
	dt.insert(Point(1000,1000,-1000));
	
	a1 = dt.insert(Point(0,0,0));
	a2 = dt.insert(Point(0,0,1));
	a3 = dt.insert(Point(0,1,1));
	a4 = dt.insert(Point(0,1,0));
	a5 = dt.insert(Point(1,0,0));
	a6 = dt.insert(Point(1,0,1));
	a7 = dt.insert(Point(2,2,2));
	a8 = dt.insert(Point(1,1,0));
	
	v1 = dt.insert(Point(0.25,0.25,0.25));
	v2 = dt.insert(Point(0.75,0.75,0.75));	
	
	//////////////////////////
	
	to_remove.push_back(v1);
	s1 = dt.number_of_vertices();
	s2 = to_remove.size();
	ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
	
	///////////////////////
	
	v1 = dt.insert(Point(0.25,0.25,0.25));
	to_remove.clear();
	to_remove.push_back(v1); 
	to_remove.push_back(v2);
	s1 = dt.number_of_vertices();
	s2 = to_remove.size();
	ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
	
	///////////////////////
	
	v1 = dt.insert(Point(0.25,0.25,0.25));
	v2 = dt.insert(Point(0.75,0.75,0.75));	
	
	to_remove.clear();
	to_remove.push_back(a1);
	s1 = dt.number_of_vertices();
	s2 = to_remove.size();
	ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
	
	///////////////////////
	
	a1 = dt.insert(Point(0,0,0));	
	
	to_remove.clear();
	to_remove.push_back(a1);
	to_remove.push_back(a2);
	s1 = dt.number_of_vertices();
	s2 = to_remove.size();
	ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
	
	///////////////////////
	
	a1 = dt.insert(Point(0,0,0));	
	a2 = dt.insert(Point(0,0,1));

	to_remove.clear();
	to_remove.push_back(a1);
	to_remove.push_back(a2);
	to_remove.push_back(a3);
	to_remove.push_back(a4);
	to_remove.push_back(a7);
	to_remove.push_back(a8);
	
	s1 = dt.number_of_vertices();
	s2 = to_remove.size();
	ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_grid_1()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	typedef typename D3::Finite_vertices_iterator  Finite_vertices_iterator;	
	typedef Point                                  Point_3;

	std::cout << "_test_rc_grid_1" << std::endl;
	
	CGAL::Random random;	
		
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	for(int i=1; i<=(1<<10); i<<=1)
	  for(int j=1; j<=(1<<10); j<<=1)
		  for(int k=1; k<=(1<<10); k<<=1)
	    to_insert.push_back(Point_3(i,j,k));
	
	dt.insert(to_insert.begin(), to_insert.end());
	
	for(Finite_vertices_iterator fit = dt.finite_vertices_begin(); 
	fit != dt.finite_vertices_end(); fit++) {
		if(random.get_double() < 0.9) to_remove.push_back(fit);
  }

	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_grid_2()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	typedef typename D3::Finite_vertices_iterator  Finite_vertices_iterator;	
  typedef Point                                  Point_3;
	std::cout << "_test_rc_grid_2" << std::endl;
	
	CGAL::Random random;	
		
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	for(int i=1; i<=(1<<10); i<<=1)
	  for(int j=1; j<=(1<<10); j<<=1)
		  for(int k=1; k<=(1<<10); k<<=1)
	    to_insert.push_back(Point_3(i,j,k));
	
	dt.insert(to_insert.begin(), to_insert.end());
	
	for(Finite_vertices_iterator fit = dt.finite_vertices_begin(); 
	fit != dt.finite_vertices_end(); fit++) {
		if(random.get_double() < 0.05) to_remove.push_back(fit);
  }

	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}

template <typename T>
void _test_rc_grid_3()
{
	typedef T                                      D3;
	typedef typename D3::Point                     Point;
	typedef typename D3::Vertex_handle             Vertex_handle;
	typedef typename D3::Finite_vertices_iterator  Finite_vertices_iterator;	
  typedef Point                                  Point_3;

	std::cout << "_test_rc_grid_3" << std::endl;
	
	CGAL::Random random;	
		
	std::vector<Point> to_insert;
	std::vector<Vertex_handle> to_remove;
	D3 dt;
	
	for(int i=1; i<=(1<<10); i<<=1)
	  for(int j=1; j<=(1<<10); j<<=1)
		  for(int k=1; k<=(1<<10); k<<=1)
	    to_insert.push_back(Point_3(i,j,k));
	
	dt.insert(to_insert.begin(), to_insert.end());
	
	for(Finite_vertices_iterator fit = dt.finite_vertices_begin(); 
	fit != dt.finite_vertices_end(); fit++) {
		if(fit->point().x() == (1<<5)) to_remove.push_back(fit);
  }

	std::size_t s1 = dt.number_of_vertices();
	std::size_t s2 = to_remove.size();
	std::size_t ns = dt.remove_cluster(to_remove.begin(), to_remove.end());
	assert(dt.is_valid());
	assert(ns == s2);
	assert(dt.number_of_vertices() + s2 == s1);
}


template <typename T>
void _test_remove_cluster() {
	std::cout << "Removing clusters of points from triangulations..." << std::endl;
	
	// small cases
	{
	  _test_rc_trivial<T>();
	  _test_rc_specific_cases_1<T>();
	  _test_rc_specific_cases_2<T>();
	  _test_rc_specific_cases_3<T>();
	  _test_rc_specific_cases_4<T>();
	  _test_rc_specific_cases_5<T>();			
	  _test_rc_specific_cases_6<T>();	
  }

  // big cases
  {
	  // random
    _test_rc_random_1<T>();
	  _test_rc_random_2<T>();
	  _test_rc_random_3<T>();
	  _test_rc_random_4<T>();
		
		// grid
		_test_rc_grid_1<T>();
		_test_rc_grid_2<T>();
		_test_rc_grid_3<T>();
  }
}


// Local Variables:
// tab-width: 2
// End:
