#ifndef ARRANGEMENT_DEMO_POINT_SNAPPER_H
#define ARRANGEMENT_DEMO_POINT_SNAPPER_H

#include "Utils.h"

template <typename ArrTraits>
class SnapStrategy : public QGraphicsSceneMixin
{
public:
  // typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;

  virtual Point_2 snapPoint(QGraphicsSceneMouseEvent* event) = 0;
  Point_2 snapQPoint(QGraphicsSceneMouseEvent* event);

protected:
  SnapStrategy(QGraphicsScene* scene_);
}; // class SnapStrategy

template <typename ArrTraits>
SnapStrategy<ArrTraits>::SnapStrategy(QGraphicsScene* scene_) : 
  QGraphicsSceneMixin(scene_)
{
}

template <typename ArrTraits>
class SnapToGridStrategy : public SnapStrategy<ArrTraits>
{
public:
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::Point_2 Kernel_point_2;
  typedef SnapStrategy<ArrTraits> Superclass;

  /*! Constructors */
  SnapToGridStrategy() : Superclass(NULL), gridSize(50) { }

  SnapToGridStrategy(QGraphicsScene* scene) :
      Superclass(scene), gridSize(50) { }

  /*! Destructors (virtual) */
  ~SnapToGridStrategy() { }

  Point_2 snapPoint(QGraphicsSceneMouseEvent* event)
  {
    return this->snapPoint(event, ArrTraits());
  }

  template <typename TTraits>
  Point_2 snapPoint(QGraphicsSceneMouseEvent* event, TTraits /* traits */)
  {
    QPointF clickedPoint = event->scenePos();
    QRectF viewportRect = this->viewportRect();
    if (viewportRect == QRectF())
    { // fallback case; we usually shouldn't end up here
      Kernel_point_2 res = this->convert(event->scenePos());
      return Point_2(CGAL::to_double(res.x()), CGAL::to_double(res.y()));
    }

    qreal d(this->gridSize / 2.0);
    int left =
      int(viewportRect.left()) - (int(viewportRect.left()) % this->gridSize);
    int right = int(viewportRect.right()) +
                (this->gridSize - int(viewportRect.right()) % this->gridSize);
    int x = int(clickedPoint.x());
    int y = int(clickedPoint.y());
    for (int i = left - this->gridSize; i <= right; i += this->gridSize)
    {
      if (i - d <= clickedPoint.x() && clickedPoint.x() <= i + d)
      {
        x = i;
        break;
      }
    }
    int top =
      int(viewportRect.top()) - (int(viewportRect.top()) % this->gridSize);
    int bottom = int(viewportRect.bottom()) +
                 (this->gridSize - int(viewportRect.bottom()) % this->gridSize);
    for (int i = top - this->gridSize; i <= bottom; i += this->gridSize)
    {
      if (i - d <= clickedPoint.y() && clickedPoint.y() <= i + d)
      {
        y = i;
        break;
      }
    }
    // return this->convert( QPointF( x, y ) );
    Point_2 res(x, y);
    return res;
  }

  void setGridSize(int size) { this->gridSize = size; }

protected:
  int gridSize;
  CGAL::Qt::Converter<Kernel> convert;
}; // class SnapToGridStrategy

template <typename Arr_>
class SnapToArrangementVertexStrategy :
    public SnapStrategy<typename Arr_::Geometry_traits_2>
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename ArrTraitsAdaptor<Traits>::Kernel Kernel;
  typedef SnapStrategy<Traits> Superclass;
  typedef typename Arrangement::Vertex_iterator Vertex_iterator;
  typedef
    typename Kernel::Compute_squared_distance_2 Compute_squared_distance_2;
  typedef typename Kernel::FT FT;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Kernel::Point_2 Kernel_point_2;

  SnapToArrangementVertexStrategy() : Superclass(NULL), arrangement(NULL) { }

  SnapToArrangementVertexStrategy(Arrangement* arr, QGraphicsScene* scene_) :
      Superclass(scene_), arrangement(arr)
  {
  }

  Point_2 snapPoint(QGraphicsSceneMouseEvent* event)
  {
    Kernel_point_2 clickedPoint = this->convert(event->scenePos());
    return this->snapPoint(clickedPoint, Traits());
  }

  template <typename TTraits>
  Point_2 snapPoint(const Kernel_point_2& clickedPoint, TTraits /* traits */)
  {
    Point_2 initialPoint(
      CGAL::to_double(clickedPoint.x()), CGAL::to_double(clickedPoint.y()));
    Point_2 closestPoint(
      CGAL::to_double(clickedPoint.x()), CGAL::to_double(clickedPoint.y()));
    bool first = true;
    FT minDist(0);
    QRectF viewportRect = this->viewportRect();
    if (viewportRect == QRectF()) { return initialPoint; }

    FT maxDist((viewportRect.right() - viewportRect.left()) / 4.0);
    for (Vertex_iterator vit = this->arrangement->vertices_begin();
         vit != this->arrangement->vertices_end(); ++vit)
    {
      Point_2 point = vit->point();
      Kernel_point_2 thisPoint(
        CGAL::to_double(point.x()), CGAL::to_double(point.y()));
      FT dist = this->compute_squared_distance_2(clickedPoint, thisPoint);
      if (first || (dist < minDist))
      {
        first = false;
        minDist = dist;
        closestPoint = point;
      }
    }
    if (!first && minDist < maxDist) { return closestPoint; }
    else
    {
      return initialPoint;
    }
  }

  void setArrangement(Arrangement* arr) { this->arrangement = arr; }

protected:
  Arrangement* arrangement;
  Compute_squared_distance_2 compute_squared_distance_2;
  CGAL::Qt::Converter<Kernel> convert;
}; // class SnapToArrangementVertexStrategy

#endif
