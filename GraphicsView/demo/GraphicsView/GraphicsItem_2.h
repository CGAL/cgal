#ifndef CGAL_GRAPHICS_ITEM_2_H
#define CGAL_GRAPHICS_ITEM_2_H

#include <CGAL/Object.h>


#include <QAbstractGraphicsShapeItem>

namespace CGAL {

class GraphicsItem_2 : public QObject, public QAbstractGraphicsShapeItem  {

  Q_OBJECT

protected:

  virtual void vModelChanged() = 0;

public slots:

  void modelChanged()
  {
    vModelChanged();
  }

};


} // namespace CGAL

#endif // CGAL_GRAPHICS_ITEM_2_H
