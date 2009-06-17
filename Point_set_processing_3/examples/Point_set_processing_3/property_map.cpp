// property_map.cpp

//----------------------------------------------------------
// Test various types of containers and property maps
// compatible with Point_set_processing_3 package.
//----------------------------------------------------------
// Usage: no parameters

#include <CGAL/Simple_cartesian.h>
#include <CGAL/property_map.h>
#include <algorithm>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point_3, Vector_3> PointVectorPair; 

// Data type = index, followed by the point, followed by a boolean
// that tells us whether the normal is oriented or not, followed by the normal vector.
typedef boost::tuple<int, Point_3, bool, Vector_3> IndexedPointWithOrientableNormalTuple;



// This is an implementation detail of the process_point_set function.
// We need this function because in process_point_set() we use std::sort.
// We sort arbitrary objects of type T, and the property map will allow
// us to access the Point_3 associated to it.

template <typename T, typename PointPMap>
struct MyLess {

  PointPMap pmap;

  MyLess(const PointPMap& p)
    : pmap(p)
  {}

  bool operator()(const T& t0, const T& t1) const
  {
    return get(pmap, &t0) < get(pmap, &t1);
  }
};


// In this example we have a function that only operates on the point part.
// It sorts them lexicographically.

template <typename Iterator, typename PointPMap >
void process_point_set(Iterator beg, Iterator end, PointPMap pmap)
{
  MyLess<typename std::iterator_traits<Iterator>::value_type,PointPMap> less(pmap);
  std::sort(beg,end,less);
}


// We can call it just with points. Then interally we use a property map
// that maps point iterators on points.

template <typename Iterator>
void process_point_set(Iterator beg, Iterator end)
{
  process_point_set(beg,end, CGAL::make_dereference_property_map(beg));
}



// Here comes a function that changes the orientation and the normal

template <typename Iterator, typename PointPMap, typename OrientationPMap, typename NormalPMap >
void orient_normals(Iterator beg, Iterator end, PointPMap point_pmap, OrientationPMap orient_pmap, NormalPMap normal_pmap)
{
  for(;beg!= end;++beg){
    Vector_3& v = get(normal_pmap, beg);
    
    boost::put(orient_pmap, beg, (v == CGAL::NULL_VECTOR));

    if(v.x() < 0){
      v = -v;
      boost::put(normal_pmap, beg, v);
    }
  }
}



int main()
{
  CGAL::set_pretty_mode(std::cout);

  // Here we run it on plain points. No need for a property map
  {
    std::vector<Point_3> points;
    
    process_point_set(points.begin(), points.end());
  }


  // Here we run it on points with normal vectors stored in a std::pair.
  // We use a property map that accesses pair::first.
  {
    std::vector<PointVectorPair> points;

    for(int i = 0; i < 10; i++){
      points.push_back(std::make_pair(Point_3(9-i,0,0), Vector_3(i,0,0)));
    }
    
    process_point_set(points.begin(),
                      points.end(),
                      CGAL::First_of_pair_property_map<PointVectorPair>());
    
    for(int i = 0; i < 10; i++){
      std::cout << points[i].first << "\t" << points[i].second << std::endl;
    }
  }


  // Here we run it on tuples. To see the interest I made up my own
  // data type = index, followed by the point, followed by a boolean
  // that tells us whether the normal is oriented or not, followed by the normal vector.
  // As the point is the second element of the tuple (that is with index 1)
  // we use a property map that accesses the 1st element of the tuple.
  {
    std::vector<IndexedPointWithOrientableNormalTuple> points;

    for(int i = 0; i < 10; i++){
      double x = (i%2)?i:-i;
      points.push_back(boost::make_tuple(i,Point_3(9-i,0,0), false, Vector_3(x,0,0)));
    }
    
    process_point_set(points.begin(),
                      points.end(),
                      CGAL::Nth_of_tuple_property_map<1,IndexedPointWithOrientableNormalTuple>());
    
    std::cout << boost::tuples::set_open('[') << boost::tuples::set_close(']') << boost::tuples::set_delimiter(','); 

    for(int i = 0; i < 10; i++){
      std::cout << points[i]  << std::endl;
    }

    //We keep the sequence in order, but determine the normal and if it is different from zero set the Boolean to true 
    orient_normals(points.begin(),
                   points.end(),
                   CGAL::make_nth_of_tuple_property_map<1>(points.begin()),
                   CGAL::make_nth_of_tuple_property_map<2>(points.begin()),
                   CGAL::make_nth_of_tuple_property_map<3>(points.begin()));
    
    std::cout << "\nAfter orient_normals\n";
    for(int i = 0; i < 10; i++){
      std::cout << points[i]  << std::endl;
    }     
  }
  
  return 0;
}
