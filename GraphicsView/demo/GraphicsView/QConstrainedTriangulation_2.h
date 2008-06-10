#ifndef CGAL_Q_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_Q_CONSTRAINED_TRIANGULATION_2_H

#include "QDatastructure.h"
#include <list>

namespace CGAL {

template <typename CT>
class QConstrainedTriangulation_2 : public QDatastructure
{

public:

  QConstrainedTriangulation_2(CT* ct_)
    : ct(ct_)
  {}

  void clear()
  {
    ct->clear();
    emit(changed());
  }

  template <typename Iterator> 
  void insert(Iterator b, Iterator e)
  {
    ct->insert(b,e);
    emit(changed());
  }

  template <typename Iterator> 
  void insert_constraints(Iterator b, Iterator e)
  {
    for(; b!= e; ++b){
      ct->insert_constraint(b->first, b->second);
    }
    emit(changed());
  }

 template <typename Iterator> 
  void insert_polyline(Iterator b, Iterator e)
  {
    typename CT::Point p, q;
    typename CT::Vertex_handle vh, wh;
    Iterator it = b;
    vh = ct->insert(*it);
    p = *it;
    ++it;
    for(; it != e; ++it){
      q = *it;
      if(p != q){
	wh = ct->insert(*it);
	ct->insert_constraint(vh,wh);
	vh = wh;
	p = q;
      } else {
	std::cout << "duplicate point: " << p << std::endl; 
      }
    }
    emit(changed());
  }
  

  virtual void insert(CGAL::Object o)
  {
    typedef typename CT::Point Point;
    std::list<Point> points;
    if(CGAL::assign(points, o)){
      if(points.size() == 1){
	ct->insert(points.front());
      } else if(points.size() == 2){
	ct->insert_constraint(points.front(), points.back());
      } else {
	std::list<Point>::iterator it = points.begin();
	++it;
	typename CT::Vertex_handle vh, wh;
	vh = ct->insert(points.front());
	for(; it != points.end(); ++it){
	  wh = ct->insert(*it);
	  ct->insert_constraint(vh,wh);
	  vh = wh;
	}
      }
    }
    emit(changed());
  }

private:

  CT * ct;
};

} // namespace CGAL

#endif // CGAL_Q_CONSTRAINED_TRIANGULATION_2_H
