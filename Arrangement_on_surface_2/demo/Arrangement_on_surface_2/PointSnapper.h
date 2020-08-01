#ifndef ARRANGEMENT_DEMO_POINT_SNAPPER_H
#define ARRANGEMENT_DEMO_POINT_SNAPPER_H

#include "GridGraphicsItem.h"
#include "GraphicsSceneMixin.h"

#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/ipower.h>
#include <boost/optional.hpp>

class GridGraphicsItem;
class QGraphicsScene;

class PointSnapperBase : public QGraphicsSceneMixin
{
public:
  using Nt_traits = CGAL::CORE_algebraic_number_traits;
  using Integer = Nt_traits::Integer;
  using Rational = Nt_traits::Rational;
  using RatKernel = CGAL::Cartesian<Rational>;
  using Point_2 = CGAL::Point_2<RatKernel>;
  using Compute_squared_distance_2 = RatKernel::Compute_squared_distance_2;

public:
  PointSnapperBase(QGraphicsScene* scene, GridGraphicsItem* grid);

  Point_2 snapPoint(const QPointF& qpt);
  Point_2 snapToGrid(const QPointF& qpt);;
  virtual boost::optional<Point_2> snapToArrangement(const QPointF& qpt) = 0;
  void setSnapToGrid(bool val);
  void setSnapToArrangement(bool val);
  bool isSnapToGridEnabled();
  bool isSnapToArrangementEnabled();

protected:
  GridGraphicsItem* gridGraphicsItem;
  bool snapToGridEnabled;
  bool snapToArrangementEnabled;
  Compute_squared_distance_2 compute_squared_distance_2;
};

template <typename Arr_>
class PointSnapper : public PointSnapperBase
{
  using Arrangement = Arr_;
  using Compute_squared_distance_2 = RatKernel::Compute_squared_distance_2;

public:
  PointSnapper(QGraphicsScene*, GridGraphicsItem*, Arrangement*);
  boost::optional<Point_2> snapToArrangement(const QPointF& qpt) override;

private:
  Arrangement* arr;
};

#endif
