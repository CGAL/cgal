// Descriptions of the file format can be found at
// http://www.autodesk.com/techpubs/autocad/acad2000/dxf/
// http://www.tnt.uni-hannover.de/soft/compgraph/fileformats/docs/DXF.ascii

#ifndef CGAL_IO_DXF_VARIANT_READER_H
#define CGAL_IO_DXF_VARIANT_READER_H

#include <CGAL/IO/Dxf_reader.h>
#include <iostream>
#include <string>
#include <list>
#include <boost/variant.hpp>

//Comment define JUJU_CAST to don't use static_cast<int>
#define COEF_MULT 1000000
#define CGAL_CAST_INT
namespace CGAL {


template<class CK,class Circular_arc_2, class Line_arc_2, class OutputIterator>
  OutputIterator variant_load(std::istream& is, OutputIterator res)
{
  typedef typename CK::Circular_arc_point_2 Circular_arc_point_2;
  typedef typename CK::Line_2  Line_2;
  typedef typename CK::Point_2 Point_2;
  typedef typename CK::Circle_2 Circle_2;
  typedef typename boost::variant< Circular_arc_2, Line_arc_2 >        Arc;
  typedef std::list<std::pair<Point_2, double> > Polygon;
  typedef std::list<Polygon> Polygons;
  typedef std::list<Circle_2> Circles;
  typedef typename CK::Intersect_2   Intersect_2;
  const double PI = 3.14159265359;

  Polygons polygons;
  Circles circles;
  CGAL::Dxf_reader<CK> reader;

  reader(is, polygons, circles);

  std::cout << "Read " << polygons.size() << " polygons, and " 
	    << circles.size() << " circles" << std::endl;
  for(typename Circles::iterator it = circles.begin(); it != circles.end(); it++){
#ifndef CGAL_CAST_INT
    Arc arc = *it;
#else
    Arc arc = Circle_2(Point_2(it->center().x()*COEF_MULT,
			       it->center().y()*COEF_MULT),
		       it->squared_radius()*COEF_MULT*COEF_MULT);
#endif
    *res++ = arc;
  }

  Point_2 first_point;
  Point_2 begin_point;
  Point_2 end_arc ;
  Point_2 center;
  double angle;
  for(typename Polygons::iterator it = polygons.begin(); it != polygons.end(); it++){
    typename Polygon::iterator pit = it->begin();
    first_point = pit->first;
    begin_point = pit->first;
    angle = pit->second;
    pit++;
    do{
      if(!angle){
#ifndef CGAL_CAST_INT
	Arc arc = Line_arc_2(begin_point, pit->first);
#else
	//For testing
	 Arc arc = Line_arc_2(Point_2(static_cast<int>(to_double(begin_point.x())*COEF_MULT),
				     static_cast<int>(to_double(begin_point.y())*COEF_MULT)),
			     Point_2(static_cast<int>(to_double(pit->first.x())*COEF_MULT),
				     static_cast<int>(to_double(pit->first.y())*COEF_MULT))); 
#endif
	if(begin_point != pit->first)
	  *res++ = arc;
	begin_point = pit->first;
	angle = pit->second;
      }
      else {
	end_arc = pit->first;
	double diff_x = to_double(begin_point.x() - end_arc.x());
	double diff_y = to_double(begin_point.y() - end_arc.y());
	double dist = sqrt(diff_x*diff_x  + diff_y*diff_y);
	double radius_sqt;
	if(angle == 1){
	  radius_sqt = (dist/sqrt(2.0));
	  radius_sqt *= radius_sqt;  
	}
	else{
	  double half_dist = dist/2.0;
	  double rayon = half_dist / sin(angle*PI/2.0);
	  double haut = rayon - (half_dist/tan(angle*PI/2.0));
	  radius_sqt = half_dist*half_dist + haut*haut;
	} 
	Intersect_2 theConstruct_intersect_2;
	std::vector<CGAL::Object> vector_objects;
	theConstruct_intersect_2(Circle_2(begin_point, radius_sqt),
				 Circle_2(end_arc, radius_sqt),
				 std::back_inserter(vector_objects));
	std::pair<Circular_arc_point_2, std::size_t> the_pair;
	//there are always 2 intersection point
	//we must choose the good
	if(begin_point.y() < end_arc.y())
	  assert(assign(the_pair, vector_objects[1]));
	else if((begin_point.y() == end_arc.y()) && (begin_point.x() > end_arc.x()))
	  assert(assign(the_pair, vector_objects[1]));
	else
	  assert(assign(the_pair, vector_objects[0]));
	center = Point_2(to_double(the_pair.first.x()),to_double(the_pair.first.y()));
#ifndef CGAL_CAST_INT
	Circular_arc_2 circular_arc(begin_point,
				    center,
				    end_arc);
#else
	//For testing
	Circular_arc_2 circular_arc(Point_2(static_cast<int>(to_double(begin_point.x())*COEF_MULT),
					    static_cast<int>(to_double(begin_point.y())*COEF_MULT)),
				    Point_2(static_cast<int>(to_double(center.x())*COEF_MULT),
					    static_cast<int>(to_double(center.y())*COEF_MULT)),
				    Point_2(static_cast<int>(to_double(end_arc.x())*COEF_MULT),
					    static_cast<int>(to_double(end_arc.y())*COEF_MULT)));

#endif
	Arc arc = circular_arc;
	*res++ = arc;
	begin_point = end_arc;
	angle = pit->second;
      }
      pit++;
    }while(pit != it->end());
    if(!angle){
#ifndef CGAL_CAST_INT
      Arc arc = Line_arc_2(begin_point, first_point);
#else
      //for testing
      Arc arc = Line_arc_2(Point_2(static_cast<int>(to_double(begin_point.x())*COEF_MULT),
				   static_cast<int>(to_double(begin_point.y())*COEF_MULT)),
			   Point_2(static_cast<int>(to_double(first_point.x())*COEF_MULT),
				   static_cast<int>(to_double(first_point.y())*COEF_MULT)));
#endif
      if(begin_point != first_point)
	*res++ = arc;
    }
    else {
      end_arc = first_point;
      double diff_x = to_double(begin_point.x() - end_arc.x());
      double diff_y = to_double(begin_point.y() - end_arc.y());
      double dist = sqrt(diff_x*diff_x  + diff_y*diff_y);
      double radius_sqt;
      if(angle == 1){
	radius_sqt = (dist/sqrt(2.0));
	radius_sqt *= radius_sqt; 
      }
      else{
	double half_dist = dist/2.0;
	double rayon = half_dist / sin(angle*PI/2.0);
	double haut = rayon - (half_dist/tan(angle*PI/2.0));
	radius_sqt = half_dist*half_dist + haut*haut;

      }
      Intersect_2 theConstruct_intersect_2;
      std::vector<CGAL::Object> vector_objects;
      theConstruct_intersect_2(Circular_arc_2(Circle_2(begin_point, radius_sqt)),
			       Circular_arc_2(Circle_2(end_arc, radius_sqt)),
			       std::back_inserter(vector_objects));
      std::pair<Circular_arc_point_2, std::size_t> the_pair;
      //there are always 2 intersection point
      //we must choose the good
      if(begin_point.y() < end_arc.y())
	assert(assign(the_pair, vector_objects[1]));
      else if((begin_point.y() == end_arc.y()) && (begin_point.x() > end_arc.x()))
	assert(assign(the_pair, vector_objects[1]));
      else
	assert(assign(the_pair, vector_objects[0]));
      center = Point_2(to_double(the_pair.first.x()),to_double(the_pair.first.y()));
#ifndef CGAL_CAST_INT
      Circular_arc_2 circular_arc(begin_point,
				  center,
				  end_arc);
#else
      	//For testing
      Circular_arc_2 circular_arc(Point_2(static_cast<int>(to_double(begin_point.x())*COEF_MULT),
					  static_cast<int>(to_double(begin_point.y())*COEF_MULT)),
				  Point_2(static_cast<int>(to_double(center.x())*COEF_MULT),
					  static_cast<int>(to_double(center.y())*COEF_MULT)),
				  Point_2(static_cast<int>(to_double(end_arc.x())*COEF_MULT),
					  static_cast<int>(to_double(end_arc.y())*COEF_MULT)));
#endif
      Arc arc = circular_arc;
      *res++ = arc;
    }
  }
  std::cout << " Loaded" << std::endl;
  return res;
}

}// namespace CGAL

#endif // CGAL_IO_DXF_VARIANT_READER_H
