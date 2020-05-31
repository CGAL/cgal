#include "CurveGraphicsItem.h"
#include "ArrangementTypes.h"

namespace CGAL
{
namespace Qt
{

template <class ArrTraits>
CurveGraphicsItem<ArrTraits>::CurveGraphicsItem() :
	painterOstream(0), boundingBox(), m_edgeColor(::Qt::red), m_edgeWidth(0),
	m_vertexColor(::Qt::red), m_vertexRadius(1)
{
	this->setZValue(4);
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::paint(
	QPainter* painter, const QStyleOptionGraphicsItem* /* option */,
	QWidget* /* widget */)
{
	// draw the curves
	QPen edgesPen(this->m_edgeColor, this->m_edgeWidth);
	painter->setPen(edgesPen);
	QRectF clippingRectangle = this->viewportRect();

	this->painterOstream =
		ArrangementPainterOstream<Traits>(painter, clippingRectangle);
	this->painterOstream.setScene(this->getScene());

	for (auto& curve : this->curves)
	{
		this->painterOstream << curve;
	}

	// draw the points
	QPen verticesPen(this->m_vertexColor, this->m_vertexRadius);
	painter->setPen(verticesPen);
	for (auto point : this->points)
	{
		this->painterOstream << point;
	}
}

template <class ArrTraits>
QRectF CurveGraphicsItem<ArrTraits>::boundingRect() const
{
	QRectF boundingRectangle = this->convert(this->boundingBox);
	return boundingRectangle;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::insert(const X_monotone_curve_2& curve)
{
	this->curves.push_back(curve);

	this->updateBoundingBox();
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::insert(const Point_2& point)
{
	this->points.push_back(point);

	this->updateBoundingBox();
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::clear()
{
	this->curves.clear();
	this->points.clear();

	this->updateBoundingBox();
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::modelChanged()
{
	this->updateBoundingBox();
	this->update();
}

template <class ArrTraits>
const QColor& CurveGraphicsItem<ArrTraits>::edgeColor() const
{
	return this->m_edgeColor;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::setEdgeColor(const QColor& color)
{
	this->m_edgeColor = color;
}

template <class ArrTraits>
int CurveGraphicsItem<ArrTraits>::edgeWidth() const
{
	return this->m_edgeWidth;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::setEdgeWidth(int width)
{
	this->m_edgeWidth = width;
}

template <class ArrTraits>
const QColor& CurveGraphicsItem<ArrTraits>::vertexColor() const
{
	return this->m_vertexColor;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::setVertexColor(const QColor& color)
{
	this->m_vertexColor = color;
}

template <class ArrTraits>
int CurveGraphicsItem<ArrTraits>::vertexRadius() const
{
	return this->m_vertexRadius;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::setVertexRadius(int radius)
{
	this->m_vertexRadius = radius;
}

static Bbox_2 remove_infs(const Bbox_2& box, const QRectF& rect) {
  double xmin = std::min(rect.left(), box.xmax());
  double ymin = std::min(rect.bottom(), box.ymax());
  double xmax = std::max(rect.right(), box.xmin());
  double ymax = std::max(rect.top(), box.ymin());
  if (!std::isinf(box.xmin()))
    xmin = box.xmin();
  if (!std::isinf(box.ymin()))
    ymin = box.ymin();
  if (!std::isinf(box.xmax()))
    xmax = box.xmax();
  if (!std::isinf(box.ymax()))
    ymax = box.ymax();
  return CGAL::Bbox_2{xmin, ymin, xmax, ymax};
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::updateBoundingBox()
{
	this->prepareGeometryChange();

	this->boundingBox = {};
	for (auto& curve : this->curves)
	{
		// some algebraic curves throw exceptions when asking about their bbox
		try
		{
			this->boundingBox += curve.bbox();
		}
		catch (...) { }
	}

	for (auto& pt : this->points)
	{
		double x = CGAL::to_double(pt.x());
		double y = CGAL::to_double(pt.y());
		this->boundingBox += CGAL::Bbox_2{x, y, x, y};
	}
  this->boundingBox = remove_infs(this->boundingBox, this->viewportRect());
}

template class CurveGraphicsItem<Seg_traits>;
template class CurveGraphicsItem<Pol_traits>;
template class CurveGraphicsItem<Conic_traits>;
template class CurveGraphicsItem<Lin_traits>;
template class CurveGraphicsItem<Alg_seg_traits>;

} // namespace Qt
} // namespace CGAL
