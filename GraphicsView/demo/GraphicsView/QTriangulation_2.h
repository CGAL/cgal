#ifndef CGAL_Q_TRIANGULATION_2_H
#define CGAL_Q_TRIANGULATION_2_H

#include "QDatastructure.h"
#include <list>

namespace CGAL {

template <typename T>
class QTriangulation_2 : public QDatastructure
{

public:

  QTriangulation_2(T* t_)
    : t(t_)
  {}

  void clear()
  {
    t->clear();
    emit(changed());
  }

  template <typename Iterator> 
  void insert(Iterator b, Iterator e)
  {
    t->insert(b,e);
    emit(changed());
  }

  void insert(CGAL::Object o)
  {
    typedef typename T::Point Point;
    std::list<Point> points;
    if(CGAL::assign(points, o)){
      t->insert(points.begin(), points.end());
    }
    emit(changed());
  }

private:

  T * t;
};

} // namespace CGAL

#endif // CGAL_Q_TRIANGULATION_2_H
