#ifndef CGAL_INPUT_H
#define CGAL_INPUT_H

#include <CGAL/Object.h>
#include <QObject>

namespace CGAL {

class Input  : public QObject
{
  Q_OBJECT

signals:

  void produce(CGAL::Object o);

};


} // namespace CGAL

#endif // CGAL_INPUT_H
