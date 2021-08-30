
#ifndef CGAL_QT_DELAUNAY_MESH_INSERT_SEEDS_H
#define CGAL_QT_DELAUNAY_MESH_INSERT_SEEDS_H

#include <CGAL/Qt/GraphicsViewInput.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>
#include <CGAL/Qt/Converter.h>

namespace CGAL {
namespace Qt {

template <typename CDT>
class DelaunayMeshInsertSeeds
  : public GraphicsViewInput
{
public:
  typedef typename CDT::Point Point;

  DelaunayMeshInsertSeeds(QGraphicsScene* s,
                          CDT* cdt,
                          QObject* parent)
    : GraphicsViewInput(parent)
    , cdt_(cdt)
    , scene_(s)
  {}

protected:

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

private:

  CDT* cdt_;
  QGraphicsScene* scene_;
};


template <typename T>
void
DelaunayMeshInsertSeeds<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  if((event->modifiers() & ::Qt::ShiftModifier)
    && event->button() == ::Qt::LeftButton)
  {
    Converter<typename T::Geom_traits> convert;

    typename T::Point seed = convert(event->scenePos());
    Q_EMIT(generate(CGAL::make_object(seed)));
  }
}



template <typename T>
bool
DelaunayMeshInsertSeeds<T>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMousePress)
  {
    QGraphicsSceneMouseEvent *mouseEvent
      = static_cast<QGraphicsSceneMouseEvent *>(event);
    mousePressEvent(mouseEvent);
    return false;
  }
  else
  {
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_DELAUNAY_MESH_INSERT_SEEDS_H
