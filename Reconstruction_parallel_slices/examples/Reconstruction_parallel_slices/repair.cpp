#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/property_map.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include <boost/lexical_cast.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Direction_2 Direction_2;

typedef std::pair<Point_2,std::string> Point;
typedef std::list<Point>::iterator iterator;

typedef std::pair<iterator,iterator> Segment;
typedef std::vector<Segment>::iterator segment_iterator;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,segment_iterator> Box;
/*
// callback function that reports all truly intersecting triangles
void report_inters( const Box& a, const Box& b) 
{
  std::cout << "Box " << (a.handle() - triangles.begin()) << " and "
            << (b.handle() - triangles.begin()) << " intersect\n";
}
*/

struct Less {
  bool operator()(const std::pair<Direction_2, int> d1,
                  const std::pair<Direction_2, int> d2) const
  {
    return CGAL::compare_angle_with_x_axis(d1.first,d2.first) == CGAL::SMALLER;
  }
};

int
crossing(const Point_2& p, const Point_2& q, const Point_2& q2,
         const Point_2& r, const Point_2& s, const Point_2& s2)
{
  //std::cerr << "crossing test:\n" << p << std::endl;
  //std::cerr << q2 << std::endl;
  //std::cerr << r << std::endl;
  //std::cerr << s2 << std::endl;
  std::vector<std::pair<Direction_2, int> > directions;
  directions.push_back(std::make_pair(Direction_2(p-q),0));
  directions.push_back(std::make_pair(Direction_2(q2-q),0));
  directions.push_back(std::make_pair(Direction_2(r-q),1));
  directions.push_back(std::make_pair(Direction_2(s2-q),1));

  Less less;
  sort(directions.begin(), directions.end(),less);
  if((directions[0].first == directions[1].first) ||
     (directions[1].first == directions[2].first) ||
     (directions[2].first == directions[3].first)){
    return 0; // collinear segments 
  }
  if((directions[0].second == directions[1].second) ||
     (directions[1].second == directions[2].second) ||
     (directions[2].second == directions[3].second)){
    return -1; // no crossing 
  }
  return 1; // a crossing
}

template <class KeyType, class PointPmap>
bool
find_safe_start(std::list< KeyType >& points, PointPmap ppmap)
{
  iterator pit = points.begin();
  iterator qit = pit; ++qit;
  iterator rit = qit; ++rit;
  std::size_t count = points.size();

  points.pop_back();
  while(std::find(qit, points.end(),*pit) != points.end()){
   std::cerr << "change the start point (degree > 2)" << std::endl;
    points.splice(points.end(), points, pit);
    if(--count == 0){
      points.push_back(points.front());
      return false;
    }
    pit = qit;
    ++qit;
    ++rit;
  }

  while(collinear( get(ppmap, *pit),
                   get(ppmap, *qit),
                   get(ppmap, *rit) ) &&
        collinear_are_strictly_ordered_along_line(get(ppmap, *pit),
                                                  get(ppmap, *qit),
                                                  get(ppmap, *rit))
  ){
    std::cerr << "change the start point (collinear overlapping)" << std::endl;
    points.splice(points.end(), points, pit);
    if(--count == 0){
      points.push_back(points.front());
      return false;
    }
    pit = qit;
    ++qit;
    ++rit;
  }
  points.push_back(points.front());
  return true;
}

template <class KeyType, class PointPmap>
void fix_intersections(std::list< KeyType >& points, PointPmap ppmap)
{
  iterator pit = points.begin();
  iterator end = points.end();

  while(true){
    bool swapped = false;
    const Point_2& p = get(ppmap, *pit);
    iterator qit = pit;
    ++qit;
    if(qit == end){
      break;
    }

    const Point_2& q = get(ppmap, *qit);
    //std::cerr << "\nouter loop: " << p << "  " << q << std::endl;  
    iterator rit = qit;
    ++rit;
    while(true){
      if(rit == end){
        break;
      }
      const Point_2& r = get(ppmap, *rit);
      iterator sit = rit;
      ++sit;
      if(sit == end){
        break;
      }

      const Point_2& s = get(ppmap, *sit);
      //std::cerr << "  inner loop: " << r << "  " << s << std::endl;  
      Segment_2 segA(p,q), segB(r,s); 
      if(CGAL::do_intersect(segA,segB)){
        //std::cerr << "intersection" << std::endl;
        if(p!=r && p!=s && q!=r && q!=s){
          //std::cerr << "real intersection" << std::endl;
          //std::cerr << segA << std::endl;
          //std::cerr << segB << std::endl;
          std::list<std::pair<Point_2,std::string> > tmp;
          tmp.splice(tmp.begin(),
                     points,
                     qit,
                     sit);
          tmp.reverse();
          points.splice(sit,tmp);
          swapped = true;
          break;
        } 
        else if (q == s){
          iterator q2it = qit; ++q2it;
          const Point_2& q2 = get(ppmap, *q2it);
          iterator s2it = sit; ++s2it;
          const Point_2& s2 = get(ppmap, *s2it);
          int sign = crossing(p,q,q2,
                              r,s,s2);
          if(sign == 1){
            //std::cerr << "crossing--------------------" << std::endl;
            std::list<std::pair<Point_2,std::string> > tmp;
            tmp.splice(tmp.begin(),
                       points,
                       q2it,
                       sit);
            tmp.reverse();
            points.splice(s2it,tmp);
            if((! CGAL::do_intersect(Segment_2(p,r),
                                     Segment_2(s,s2))) && 
               (! CGAL::do_intersect(Segment_2(p,r),
                                     Segment_2(q,q2)))){
              points.erase(qit);
            } else if((! CGAL::do_intersect(Segment_2(q2,s2),
                                           Segment_2(p,q))) &&
                      (! CGAL::do_intersect(Segment_2(q2,s2),
                                            Segment_2(q,r)))
                      ) {
              points.erase(sit);

            } else {
              //std::cerr << "shortcut of pqq2 or rss2 would introduce an intersection" << std::endl;
            }
            points.splice(qit,tmp);
          } else if(sign != 1){
            //std::cerr << "touching--------------------" << std::endl;
            //std::cerr << "2\n" << p << std::endl << q2 << std::endl;
            //std::cerr << "2\n" << r << std::endl << s2 << std::endl;
            if((! CGAL::do_intersect(Segment_2(p,q2),
                                     Segment_2(r,s))) && 
               (! CGAL::do_intersect(Segment_2(p,q2),
                                     Segment_2(s2,s)))){
              points.erase(qit);
            } else if((! CGAL::do_intersect(Segment_2(r,s2),
                                           Segment_2(p,q))) &&
                      (! CGAL::do_intersect(Segment_2(r,s2),
                                            Segment_2(q2,q)))
                      ) {
              points.erase(sit);

            } else {
              //std::cerr << "shortcut of pqq2 or rss2 would introduce an intersection" << std::endl;
            }
          } else {
            //std::cerr << "todo: treat overlapping segments" << std::endl;
          }
             
        }          
      }
      ++rit;
        }
    if(! swapped){
      ++pit;
    }
  }
}



template <class KeyType, class PointPmap>
void fix_consecutive_overlapping_segments(
  std::list< KeyType >& points, PointPmap ppmap)
{
 iterator pit = points.begin();
  while(true){
    if(points.size() == 4){
      break;
    }
    iterator qit = pit;
    ++qit;
    if(qit == points.end()){
      break;
    }
    const Point_2& p = get(ppmap, *pit);
    const Point_2& q = get(ppmap, *qit);
    
    iterator rit = qit;
    ++rit;
    if(rit == points.end()){
      const Point_2& r = get(ppmap, *(++points.begin()));
      if(CGAL::collinear(p,q,r) && 
         (! CGAL::collinear_are_strictly_ordered_along_line(p,q,r))){
        points.pop_back();
        points.pop_front();
        points.push_back(points.front()); // duplicate the first point
      }
      break;
    }
    const Point_2& r = get(ppmap, *rit);
    
    if(CGAL::collinear(p,q,r)){
      if(! CGAL::collinear_are_strictly_ordered_along_line(p,q,r)){
        points.erase(qit);

        if(get(ppmap, *pit) == get(ppmap, *rit)){
          std::cerr << "identical consecutive points " << get(ppmap, *pit) << std::endl;
          points.erase(rit);
          
        }
        if(pit != points.begin()){
            --pit;
          }
        continue;  // and do not advance pit, because its two successors are different
      }
    }
    ++pit; 
  }
}

struct Get_point_from_pair{
  //classical typedefs
  typedef const std::pair<Point_2, std::string>& key_type;
  typedef Point_2 value_type;
  typedef const value_type& reference;
  typedef boost::readable_property_map_tag category;
  friend reference get(const Get_point_from_pair&, key_type p) {
    return p.first;
  }
};

void contour(int count, int n, int dim)
{
  std::cerr << "\ncontour  " << count << " with " << n << " vertices" << std::endl;
 
  std::cin.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
  std::list<std::pair<Point_2,std::string> > points;

  for(int i = 0; i <n; i++){
    std::string line;
    std::getline(std::cin, line);
    std::istringstream iss(line);

    double x,y,z;
    iss >> x >> y >> z;
    
    Point_2 p = (dim==0) ? Point_2(y,z) : (dim==1) ? Point_2(x,z) : Point_2(x,y);
    if(! points.empty()){
      if(p != points.back().first){
      points.push_back(std::make_pair(p, line));
      } else {
        std::cerr << "identical consecutive points " << p << std::endl;
      }
    } else {
     points.push_back(std::make_pair(p, line)); 
    }
  }
  if(points.size() <= 3){
    std::cerr << "ignore segment" <<  std::endl;
    std::cerr << points.front().first << " -- " << (++points.begin())->first << std::endl;
    return;
  }

  Get_point_from_pair ppmap;
  
  if(find_safe_start(points, ppmap)){ 
    fix_consecutive_overlapping_segments(points, ppmap);
    
    if(points.size() <= 3){
      std::cerr << "ignore segment" <<  std::endl;
      std::cerr << points.front().first << " -- " << get(ppmap, *(++points.begin())) << std::endl;
      return;
    }
  
    fix_intersections(points, ppmap);
  } else {
    std::cerr << "No safe start point found" << std::endl;
  }

  std::cout << points.size() <<std::endl;
  for(iterator it = points.begin(); it!= points.end(); ++it){
    std::cout << it->second << std::endl;
  }
}

int main(int argc, char** argv)
{
  if (argc!=2)
  {
    std::cerr << "Please provide the constant coordinate index (0,1,2)\n";
    return 1;
  }
  int constant_coordinate = boost::lexical_cast<int>(argv[1]);
  int n;
  int count = 0;
  while(std::cin >> n){
    contour(count, n, constant_coordinate);
    count++;
  }
  std::cerr << "done" << std::endl;
  return 0;
}
