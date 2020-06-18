#include "GraphicsSceneMixin.h"
#include "ArrangementDemoGraphicsView.h"
#include <QGraphicsScene>
#include <QGraphicsView>

QGraphicsSceneMixin::QGraphicsSceneMixin(QGraphicsScene* scene_) : scene{scene_}
{
}

QGraphicsSceneMixin::~QGraphicsSceneMixin() { }

void QGraphicsSceneMixin::setScene(QGraphicsScene* scene_)
{
  this->scene = scene_;
}

QGraphicsScene* QGraphicsSceneMixin::getScene() const { return this->scene; }

QRectF QGraphicsSceneMixin::viewportRect() const
{
  QGraphicsView* viewport = this->getView();
  return viewport ? ArrangementDemoGraphicsView::viewportRect(viewport)
                  : QRectF{};
}

QPoint QGraphicsSceneMixin::fromScene(QPointF p)
{
  QGraphicsView* viewport = this->getView();
  return viewport ? viewport->mapFromScene(p) : QPoint{};
}

QPointF QGraphicsSceneMixin::toScene(QPoint p) const
{
  QGraphicsView* viewport = this->getView();
  return viewport ? viewport->mapToScene(p) : QPointF{};
}

QGraphicsView* QGraphicsSceneMixin::getView() const
{
  if (!this->scene) return nullptr;

  QList<QGraphicsView*> views = this->scene->views();
  if (views.size() == 0) return nullptr;
  // assumes the first view is the right one
  return views.first();
}
