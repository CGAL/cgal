#include "MyMink.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <Graphics.h>
#include <iostream>
#include <fstream>
#include <boost\timer.hpp>

#include <CGAL/minkowski_sum_2.h>
#include <SweepCollisionDetection.h>
//#include <CGAL/Filtered_kernel.h>
using namespace std;
using namespace CGAL;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//typedef CGAL::Filtered_kernel<Kernel_> Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon;
typedef CGAL::Arr_segment_traits_2<Kernel> Segment_traits_2;
typedef Segment_traits_2::Point_2 Point;

typedef CGAL::Polygon_with_holes_2<Kernel>            PolygonWithHoles;   
std::istream &operator>>(std::istream & is, Polygon& pl)
{
    pl.clear();
    int size;
    is >> size;
    for (int i=0; i<size; i++)
    {
	Point p ;
        is >> p;
        pl.push_back(p);
    }


    return is;
}

// Global counters
int _my_global_counter =0;	
double rc_time =0;
double arr_build_time =0;
double final_stage_time =0;
double degenerate_stage =0;

int main(int argc,char* argv[]){
	Polygon a;
	Polygon b;
	PolygonWithHoles sum;
	
	
	if (argc < 3)
	{
		cout<< "Usage: mink.exe file1 file2" <<endl;
		return 0;
	}
	//fstream data("data.txt");
	//data >> a >> b;
	fstream data1(argv[1]);
	data1 >> a;
	fstream data2(argv[2]);
	data2 >> b;
	if (a.is_clockwise_oriented())
		a.reverse_orientation();
	if (b.is_clockwise_oriented())
		b.reverse_orientation();
	if (!a.is_simple ())
	{
		 cout<< "a is not simple" <<endl;
		 return -1;
	}
	if (!b.is_simple ())
	{
		 cout<< "b is not simple" <<endl;
		 return -1;
	}
	if (!a.is_counterclockwise_oriented ())
	{
		 cout<< "a is not counter clockwise oriented." <<endl;
		 return -1;
	}
	if (!b.is_counterclockwise_oriented ())
	{
		 cout<< "b is not counter clockwise oriented." <<endl;
		 return -1;
	}
	//if (a.is_convex())
	//{
	//	cout << "a is convex" <<endl;
	//}
	//if (b.is_convex())
	//{
	//	cout << "b is convex" <<endl;
	//}
	//double bottomx = (CGAL::to_double(a.bottom_vertex()->x())+CGAL::to_double(b.bottom_vertex()->x()));
	double bottomy = (CGAL::to_double(a.bottom_vertex()->y())+CGAL::to_double(b.bottom_vertex()->y()));
	double leftx = (CGAL::to_double(a.left_vertex()->x())+CGAL::to_double(b.left_vertex()->x()));
	//double lefty = (CGAL::to_double(a.bottom_vertex()->x())+CGAL::to_double(b.bottom_vertex()->x()));
	//double topx = (CGAL::to_double(a.bottom_vertex()->x())+CGAL::to_double(b.bottom_vertex()->x()));
	double topy = (CGAL::to_double(a.top_vertex()->y())+CGAL::to_double(b.top_vertex()->y()));
	double rightx = (CGAL::to_double(a.right_vertex()->x())+CGAL::to_double(b.right_vertex()->x()));
	//double righty = (CGAL::to_double(a.bottom_vertex()->x())+CGAL::to_double(b.bottom_vertex()->x()));
	global_graphics = new Graphics(0,0,leftx,rightx,bottomy,topy); 
	QColor c(0,255,0);

if (SHOW_STAGES)
{
	

	global_graphics->draw_polygon(a,c);
	global_graphics->display();
	global_graphics->clear();
	global_graphics->draw_polygon(b,c);
	global_graphics->display();
	global_graphics->clear();
}
//#endif
	

	/*a.push_back(Point(0.197005,0.504386));
	a.push_back(Point(0.199309,0.384503));
	a.push_back(Point(0.800691,0.375731));
	a.push_back(Point(0.796083,0.500000));*/

/*	a.push_back(Point(0.5,0.7));
	a.push_back(Point(0.2,0.2));
	a.push_back(Point(0.5,0.4));
	a.push_back(Point(0.8,0.2));
	*/
/*
	b.push_back(Point(0,0));
	b.push_back(Point(1,0));
  //b.push_back(Point(2,0));
	b.push_back(Point(1,1));

	*/

	/*b.push_back(Point(0.906682,0.890351));
	b.push_back(Point(0.261521,0.896199));
	b.push_back(Point(0.256912,0.703216));
	b.push_back(Point(0.839862,0.711988));
	b.push_back(Point(0.849078,0.583333));
	b.push_back(Point(0.236175,0.606725));
	b.push_back(Point(0.240783,0.472222));
	b.push_back(Point(0.879032,0.472222));
	b.push_back(Point(0.892857,0.311404));
	b.push_back(Point(0.164747,0.270468));
	b.push_back(Point(0.169355,0.156433));
	b.push_back(Point(0.975806,0.179825));*/


	/*b.push_back(Point(0.307604,0.899123));
	b.push_back(Point(0.100230,0.878655));
	b.push_back(Point(0.102535,0.097953));
	b.push_back(Point(0.899770,0.118421));
	b.push_back(Point(0.899770,0.896199));
	b.push_back(Point(0.809908,0.896199));
	b.push_back(Point(0.802995,0.796784));
	b.push_back(Point(0.853687,0.796784));
	b.push_back(Point(0.849078,0.200292));
	b.push_back(Point(0.148618,0.200292));
	b.push_back(Point(0.148618,0.711988));
	b.push_back(Point(0.300000,0.703216));*/

	//// test collision detection.
	//SweepCollisionDetector<Kernel,Polygon::Container> coll;
	//bool collide = coll.checkCollision(a,b);

	//cout << "collide: " << collide << "\n";

	double timer_end;

	//Polygon sum_bound;
	boost::timer t;
	//cout<< "before";
	try{
	sum = minkowski_sum_2_(a,b);
/*	global_graphics->draw_polygon(sum,c);
	global_graphics->display();
	global_graphics->clear();*/
	} 
	catch(exception ex)
	{
		cout <<ex.what();
	}
	//cout << _my_global_counter << std::endl;
//	sum = minkowski_sum_2(a,b);
	//cout << "after";

	timer_end = t.elapsed();
/*	cout << "Conv cycles time stage : " <<rc_time << std::endl;
	cout << "Arrangement time stage : " <<arr_build_time << std::endl;
	cout << "degenerate time stage : " << degenerate_stage <<std::endl;
	cout << "Final stage : " << final_stage_time << std::endl;*/
	cout << "timer: " << timer_end << std::endl;
//	cout << "Num of convolution segments : " <<_my_global_counter << std::endl;

	/// WEIN RUN
///**** commented out For Mem test
	boost::timer t2;
	sum = minkowski_sum_2(a,b);	
	timer_end = t2.elapsed();

	cout << "timer: " << timer_end << std::endl;
	
//*******/


	/*	cout << "Conv cycles time stage : " <<global_time_counter_first_stage << std::endl;
	cout << "Arrangement time stage : " <<global_time_counter_second_stage << std::endl;
	cout << "Final stage : " << global_time_counter_third_stage << std::endl;*/
	//cout << "Num of convolution segments : " << _ron_global_counter << std::endl;
	//char e;
	//cin >> e;
/*	global_graphics->draw_polygon(sum.outer_boundary(),c);
	PolygonWithHoles::Hole_iterator itr = sum.holes_begin();
	for (;itr != sum.holes_end();++itr)
	{
		global_graphics->draw_polygon(*itr,c);
	}
	global_graphics->display();
	global_graphics->clear();*/
	//std::list<Polygon>  sum_holes;
	//Minkowski_sum_by_convolution_lien_2<Kernel,
	//Minkowski_sum_by_convolution_lien_2(a,b,sum_bound,std::back_inserter(sum_holes));
}
