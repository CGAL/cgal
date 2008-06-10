#ifndef CGAL_Q_INPUT_H
#define CGAL_Q_INPUT_H

#include <CGAL/Object.h>
#include <QObject>

namespace CGAL {

class QInput  : public QObject
{
  Q_OBJECT

signals:

  void generate(CGAL::Object o);

};


} // namespace CGAL

#endif // CGAL_Q_INPUT_H
