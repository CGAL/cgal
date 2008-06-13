
#ifndef CGAL_Q_POLYLINE_INPUT_2_H
#define CGAL_Q_POLYLINE_INPUT_2_H

#include "QConverter.h"

#include "QPolylineInput_2_non_templated_base.h"

namespace CGAL {

template <typename K>
class QPolylineInput_2 : public QPolylineInput_2_non_templated_base
{
public:
  QPolylineInput_2(QGraphicsScene* s, int n = 0, bool closed = true)
    : QPolylineInput_2_non_templated_base(s, n, closed)
  {
  }

protected:
  void generate_polygon() {
    std::list<typename K::Point_2> points;
    QConverter<K> convert;
    convert(points, this->polygon); 
    emit(generate(CGAL::make_object(points)));
  }
};

} // namespace CGAL

#endif // CGAL_Q_POLYLINE_INPUT_2_H
