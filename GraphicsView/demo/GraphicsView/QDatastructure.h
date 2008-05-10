#ifndef CGAL_QDATASTRUCTURE_H
#define CGAL_QDATASTRUCTURE_H

#include <QObject>
#include <CGAL/Object.h>
#include <list>

namespace CGAL {

class QDatastructure : public QObject
{
  Q_OBJECT

protected:

  virtual void insert(CGAL::Object o) = 0;


public slots:

  void consume(CGAL::Object o)
  {
    insert(o);
  }


  signals:

  void changed();
};




} // namespace CGAL

#endif // CGAL_QDATASTRUCTURE_H
