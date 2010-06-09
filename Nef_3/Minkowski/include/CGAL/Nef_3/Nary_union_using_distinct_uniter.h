#ifndef CGAL_NEF_NARY_UNION_USING_DISTINCT_UNITER_H
#define CGAL_NEF_NARY_UNION_USING_DISTINCT_UNITER_H

#include <CGAL/Nef_3/Nary_union_strategy.h>
#include <CGAL/Box_intersection_d/Box_d.h>
#include <CGAL/Nef_3/Distinct_polyhedron_uniter.h>
#include <list>

namespace CGAL {

template<class Polyhedron>
class Nary_union_using_distinct_uniter {

  typedef typename Polyhedron::Vertex_const_iterator 
    Vertex_const_iterator;
  typedef Box_intersection_d::Box_d<double, 3> Box;
  typedef typename std::pair<Polyhedron, Box> PolyBox;
  std::list<PolyBox> queue;
  typedef typename std::list<PolyBox>::iterator pit;
  typedef CGAL::Distinct_polyhedron_uniter<Polyhedron> 
    Distinct_polyhedron_uniter;

  void print_box(const Box& b) {
    std::cerr << b.min_coord(0) << ","
	      << b.min_coord(1) << ","
	      << b.min_coord(2) << "->"
	      << b.max_coord(0) << ","
	      << b.max_coord(1) << ","
	      << b.max_coord(2) << std::endl;
  }

  Box compute_box(const Polyhedron& P) {
    bool first = true;
    Box b;
    Vertex_const_iterator vi;
    CGAL_forall_vertices(vi,P) {
      std::pair<double, double> q[3];
      q[0] = CGAL::to_interval( vi->point().x() );
      q[1] = CGAL::to_interval( vi->point().y() );
      q[2] = CGAL::to_interval( vi->point().z() );
      if(first) {
	double min[3], max[3];
	min[0] = q[0].first;
	min[1] = q[1].first;
	min[2] = q[2].first;
	max[0] = q[0].second;
	max[1] = q[1].second;
	max[2] = q[2].second;
	b = Box(min, max);
	first = false;
      } else
	b.extend(q);
    }
    return b;
  }

  Box unite_boxes(const Box& b1, const Box& b2) {
    Box b = b1;
    double q[3];
    q[0] = b2.min_coord(0);
    q[1] = b2.min_coord(1);
    q[2] = b2.min_coord(2);
    b.extend(q);
    q[0] = b2.max_coord(0);
    q[1] = b2.max_coord(1);
    q[2] = b2.max_coord(2);
    b.extend(q);
    return b;
  }

  bool boxes_are_distinct(const Box& b1, const Box& b2) {
    if(((b1.min_coord(0) >= b2.min_coord(0) &&
	 b1.min_coord(0) <= b2.max_coord(0)) ||
	(b2.min_coord(0) >= b1.min_coord(0) &&
	 b2.min_coord(0) <= b1.max_coord(0))) &&
       ((b1.min_coord(1) >= b2.min_coord(1) &&
	 b1.min_coord(1) <= b2.max_coord(1)) ||
	(b2.min_coord(1) >= b1.min_coord(1) &&
	 b2.min_coord(1) <= b1.max_coord(1))) &&
       ((b1.min_coord(2) >= b2.min_coord(2) &&
	 b1.min_coord(2) <= b2.max_coord(2)) ||
	(b2.min_coord(2) >= b1.min_coord(2) &&
	 b2.min_coord(2) <= b1.max_coord(2))))
      return false;
    return true;
  }

 public:
  void add_polyhedron(const Polyhedron& P) {
    Box b = compute_box(P);
    queue.push_back(PolyBox(P,b));
    if(queue.size()%100 == 0)
      std::cerr << "queue size: " << queue.size() << std::endl;
  }
  
  Polyhedron get_union() {
    
    Polyhedron empty;
    while(queue.size() > 1) {
      if(queue.size()%100 == 0 || queue.size()<30)     
	std::cerr << queue.size() << " polyhedra in the queue" << std::endl;
      
      pit i1(queue.begin()), i2(i1++);
      if(boxes_are_distinct(i1->second, i2->second)) {
	std::list<Polyhedron> N_list;
	N_list.push_back(i1->first);
	N_list.push_back(i2->first);
	Polyhedron tmp(empty);
	Distinct_polyhedron_uniter dpu(N_list);
	tmp.delegate(dpu,false,false);
	PolyBox pb(tmp,
		   unite_boxes(i1->second, i2->second));
	queue.push_back(pb);
      } else {
	PolyBox pb(i1->first + i2->first, 
		   unite_boxes(i1->second, i2->second));
	queue.push_back(pb);
      }
      queue.pop_front();
      queue.pop_front();
    }
    return queue.front().first;
  }
};

} //namespace CGAL
#endif // CGAL_NEF_NARY_UNION_USING_DISTINCT_UNITER_H
