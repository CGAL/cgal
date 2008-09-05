#ifndef CGAL_QT_GRAPHICS_ITEM_H
#define CGAL_QT_GRAPHICS_ITEM_H

#include <CGAL/auto_link/Qt4.h>
#include <QObject>
#include <QGraphicsItem>
#include <CGAL/Object.h>



namespace CGAL {
namespace Qt {

class GraphicsItem : public QObject, public QGraphicsItem {

  Q_OBJECT

public slots:

  virtual void modelChanged() = 0;
};


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_ITEM_H
