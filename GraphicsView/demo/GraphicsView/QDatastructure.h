#ifndef CGAL_QDATASTRUCTURE_H
#define CGAL_QDATASTRUCTURE_H

#include <QObject>
#include <CGAL/Object.h>
#include <list>

namespace CGAL {

class QDatastructure : public QObject
{
  Q_OBJECT

public slots:
  virtual void insert(CGAL::Object ) = 0;

  signals:

  void changed();
};




} // namespace CGAL

#endif // CGAL_QDATASTRUCTURE_H
