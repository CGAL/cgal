#ifndef CGAL_Q_GRAPHICS_ITEM_2_H
#define CGAL_Q_GRAPHICS_ITEM_2_H

#include <CGAL/Object.h>


#include <QAbstractGraphicsShapeItem>

namespace CGAL {

class QGraphicsItem_2 : public QObject, public QAbstractGraphicsShapeItem  {

  Q_OBJECT

public slots:

virtual void modelChanged() = 0;
};


} // namespace CGAL

#endif // CGAL_Q_GRAPHICS_ITEM_2_H
