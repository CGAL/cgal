#ifndef INTERSECTION_LINE_2_CIRCLE_2_MAP_H
#define INTERSECTION_LINE_2_CIRCLE_2_MAP_H

#include <map>
#include <vector>
#include <CGAL/Object.h>

namespace CGAL {
namespace CGALi {

class Intersection_line_2_circle_2_map {

typedef struct inter_map_pair {
  int x, y;
  inter_map_pair(int xx=0, int yy=0) : x(xx), y(yy) {}
  inter_map_pair(const inter_map_pair &i) : x(i.x), y(i.y) {}
  bool operator<(const inter_map_pair &i) const {
    if(x < i.x) return true;
    if(x > i.x) return false; 
    if(y < i.y) return true;
    return false;
  }
} inter_map_pair;
typedef std::map< inter_map_pair , CGAL::Object > Table;

private:
  Table intersection_map;
  unsigned int id_gen;

public:
  Intersection_line_2_circle_2_map() : id_gen(0) { intersection_map.clear(); }
  ~Intersection_line_2_circle_2_map() { intersection_map.clear(); }
  
  unsigned int get_new_id() {
    return ++id_gen;
  } 

  template < class T >
  bool find(int id1, int id2, T& res) const {
    Table::const_iterator p = intersection_map.find(
      inter_map_pair(id1,id2));
    if(p == intersection_map.end()) return false;
    assign(res, p->second);
    return true;
  }

  template < class T >
  void put(const int id1, const int id2, const T& res) {
    intersection_map[inter_map_pair(id1,id2)] = CGAL::make_object(res);
  }
};

} // endof internal cgal namespace
} //endof cgal namespace
#endif // INTERSECTION_LINE_2_CIRCLE_2_MAP_H
