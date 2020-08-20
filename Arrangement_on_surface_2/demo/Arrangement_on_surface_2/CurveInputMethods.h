#ifndef ARRANGEMENT_DEMO_CURVE_INPUT_METHODS_H
#define ARRANGEMENT_DEMO_CURVE_INPUT_METHODS_H

#include <QGraphicsLineItem>
#include <vector>

#include "GraphicsSceneMixin.h"
#include "PointSnapper.h"
#include "PointsGraphicsItem.h"

class QEvent;
class QGraphicsSceneMouseEvent;

namespace CGAL
{

namespace Qt
{
enum class CurveType
{
  Segment,
  Polyline,
  Ray,
  Line,
  Circle,
  Ellipse,
  ThreePointCircularArc,
  FivePointConicArc,
  AlgebraicEquation,
  Bezier,
  None,
};

class CurveInputMethodCallback
{
public:
  using Point_2 = PointSnapperBase::Point_2;

  virtual void curveInputDoneEvent(
    const std::vector<Point_2>& clickedPoints, CurveType type) = 0;
};

class CurveInputMethod : public GraphicsSceneMixin
{
public:
  using Point_2 = PointSnapperBase::Point_2;

  CurveInputMethod(CurveType, int numPoints_ = -1);
  virtual ~CurveInputMethod() { }

  void setCallback(CurveInputMethodCallback*);
  void setPointSnapper(PointSnapperBase*);
  CurveType curveType() const;

  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
  virtual void mousePressEvent(QGraphicsSceneMouseEvent* event);
  void reset();
  void beginInput_();

  void setColor(QColor);
  QColor getColor() const;

  Point_2 snapPoint(QGraphicsSceneMouseEvent* event);
  QPointF snapQPoint(QGraphicsSceneMouseEvent* event);

protected:
  virtual void beginInput();
  virtual void resetInput();
  virtual void updateVisualGuideNewPoint(const std::vector<QPointF>&);
  virtual void
  updateVisualGuideMouseMoved(const std::vector<QPointF>&, const QPointF&);
  void appendGraphicsItem(QGraphicsItem* item);

protected:
  QColor color;

private:
  const int numPoints;
  std::vector<QPointF> clickedPoints;
  std::vector<Point_2> clickedBigPoints;
  PointsGraphicsItem pointsGraphicsItem;
  CurveInputMethodCallback* callback;
  const CurveType type;
  std::vector<QGraphicsItem*> items;
  bool itemsAdded;
  PointSnapperBase* snapper;
};

class SegmentInputMethod : public CurveInputMethod
{
public:
  SegmentInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsLineItem segmentGuide;
};

class RayInputMethod : public CurveInputMethod
{
public:
  RayInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsLineItem rayGuide;
};

class LineInputMethod : public CurveInputMethod
{
public:
  LineInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsLineItem lineGuide;
};

class PolylineInputMethod : public CurveInputMethod
{
public:
  PolylineInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

  void updateVisualGuideNewPoint(const std::vector<QPointF>&) override;

private:
  QGraphicsPathItem polylineGuide;
  QGraphicsLineItem lastLine;
  QPainterPath painterPath;
};

class CircleInputMethod : public CurveInputMethod
{
public:
  CircleInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsEllipseItem circleGuide;
};

class EllipseInputMethod : public CurveInputMethod
{
public:
  EllipseInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsEllipseItem ellipseGuide;
};

class ThreePointCircularInputMethod : public CurveInputMethod
{
public:
  ThreePointCircularInputMethod();
};

class FivePointConicInputMethod : public CurveInputMethod
{
public:
  FivePointConicInputMethod();
};

class BezierInputMethod : public CurveInputMethod
{
public:
  BezierInputMethod();
  void beginInput() override;
  void updateVisualGuideNewPoint(const std::vector<QPointF>&) override;
  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsPathItem bezierOldGuide;
  QGraphicsPathItem bezierGuide;
  QPainterPath painterOldPath;
  QPainterPath painterPath;
  std::vector<QPointF> controlPoints;
  std::vector<QPointF> cache;
};

} // namespace Qt
} // namespace CGAL

#endif
