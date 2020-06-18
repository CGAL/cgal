#ifndef CGAL_ARRANGEMENTS_DEMO_GRAPHICS_SCENE_MIXIN_H
#define CGAL_ARRANGEMENTS_DEMO_GRAPHICS_SCENE_MIXIN_H

#include <QRectF>
#include <QPoint>
#include <QPointF>

class QGraphicsScene;
class QGraphicsView;

class QGraphicsSceneMixin
{
public:
  /*! Costructor */
  QGraphicsSceneMixin(QGraphicsScene* scene_ = nullptr);

  /*! Destructor (virtual) */
  virtual ~QGraphicsSceneMixin();
  virtual void setScene(QGraphicsScene* scene_);
  QGraphicsScene* getScene() const;
  QRectF viewportRect() const;
  QPoint fromScene(QPointF p);
  QPointF toScene(QPoint p) const;
  QGraphicsView* getView() const;

protected: // fields
  QGraphicsScene* scene;
};

#endif
