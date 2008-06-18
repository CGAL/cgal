
#ifndef CGAL_Q_POLYLINE_INPUT_2_H
#define CGAL_Q_POLYLINE_INPUT_2_H

#include <CGAL/IO/QtConverter.h>

#include <CGAL/IO/QtPolylineInput_non_templated_base.h>

namespace CGAL {

template <typename K>
class QtPolylineInput : public QtPolylineInput_non_templated_base
{
public:
  QtPolylineInput(QGraphicsScene* s, int n = 0, bool closed = true)
    : QtPolylineInput_non_templated_base(s, n, closed)
  {
  }

protected:
  void generate_polygon() {
    std::list<typename K::Point_2> points;
    QtConverter<K> convert;
    convert(points, this->polygon); 
    emit(generate(CGAL::make_object(points)));
  }
};

} // namespace CGAL

#endif // CGAL_Q_POLYLINE_INPUT_2_H
