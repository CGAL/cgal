// Copyright (c) 2008, 2012, 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#include "ArrangementGraphicsItem.h"
#include "ArrangementGraphicsItemTyped.h"
#include "ArrangementPainterOstream.h"
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "PointsGraphicsItem.h"
#include "Utils/ConstructBoundingBox.h"
#include "Utils/PointLocationFunctions.h"

#include <CGAL/Qt/Converter.h>
#include <QGraphicsScene>
#include <QPainter>
#include <limits>

namespace CGAL {
namespace Qt {

static constexpr int margin = 2;

ArrangementGraphicsItemBase::ArrangementGraphicsItemBase() :
    verticesPen(QPen(::Qt::blue, 3.)), edgesPen(QPen(::Qt::blue, 1.))
{
  this->verticesPen.setCosmetic(true);
  this->edgesPen.setCosmetic(true);
  this->edgesPen.setWidth(2);
  this->facesPen.setCosmetic(true);
  this->facesPen.setColor(::Qt::transparent);
}

// msvc2015 doesn't play well with polymorphic lambdas
namespace
{
struct ExplicitLambda
{
  template <typename Arrangement>
  void operator()(demo_types::TypeHolder<Arrangement>)
  {
    Arrangement* arr = nullptr;
    CGAL::assign(arr, arr_obj);
    agi = new ArrangementGraphicsItem<Arrangement>(arr);
  }

  ArrangementGraphicsItemBase*& agi;
  CGAL::Object& arr_obj;
};
} // anonymous namespace

ArrangementGraphicsItemBase* ArrangementGraphicsItemBase::create(
  demo_types::TraitsType tt, CGAL::Object arr_obj)
{
  ArrangementGraphicsItemBase* agi;
  ExplicitLambda explicit_lambda{agi, arr_obj};
  demo_types::visitArrangementType(tt, explicit_lambda);
  return agi;
}

const QPen& ArrangementGraphicsItemBase::getVerticesPen() const
{
  return this->verticesPen;
}

const QPen& ArrangementGraphicsItemBase::getEdgesPen() const
{
  return this->edgesPen;
}

void ArrangementGraphicsItemBase::setVerticesPen(const QPen& pen)
{
  this->verticesPen = pen;
}

void ArrangementGraphicsItemBase::setEdgesPen(const QPen& pen)
{
  this->edgesPen = pen;
}

template <typename Arr_>
ArrangementGraphicsItem<Arr_>::ArrangementGraphicsItem(Arrangement* arr_) :
    arr(arr_), pointsGraphicsItem(new PointsGraphicsItem(this))
{
  this->updatePointsItem();
  this->updateBoundingBox();
}

template < typename Arr_ >
QRectF
ArrangementGraphicsItem< Arr_ >::
boundingRect( ) const
{
  qreal xmin = -(std::numeric_limits<qreal>::max)() / 4;
  qreal ymin = -(std::numeric_limits<qreal>::max)() / 4;
  qreal xmax = (std::numeric_limits<qreal>::max)() / 4;
  qreal ymax = (std::numeric_limits<qreal>::max)() / 4;
  if (this->bb.xmin() > xmin) xmin = this->bb.xmin();
  if (this->bb.ymin() > ymin) ymin = this->bb.ymin();
  if (this->bb.xmax() < xmax) xmax = this->bb.xmax();
  if (this->bb.ymax() < ymax) ymax = this->bb.ymax();
  if (xmin > xmax || ymin > ymax)
  {
    xmin = 0;
    xmax = 0;
    ymin = 0;
    ymax = 0;
  }
  return {QPointF{xmin, ymin}, QPointF{xmax, ymax}};
}

template < typename Arr_ >
void
ArrangementGraphicsItem< Arr_ >::
paint(QPainter* painter,
      const QStyleOptionGraphicsItem* /* option */,
      QWidget*  /*widget*/)
{
  this->paint(painter, *(arr->traits()));
}

template <typename Arr_>
template <typename TTraits>
void ArrangementGraphicsItem<Arr_>::paint(
  QPainter* painter, const TTraits& traits)
{
  this->paintFaces(painter);
  this->paintEdges(painter, traits);
}

template <typename Arr_>
template <typename Coefficient_>
void ArrangementGraphicsItem<Arr_>::paint(
  QPainter* painter,
  const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>& traits)
{
  this->paintWithFloodFill(painter, traits);
}

template <typename Arr_>
template <typename AlgebraicKernel_d_1>
void ArrangementGraphicsItem<Arr_>::paint(
  QPainter* painter,
  const CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>& traits)
{
  this->paintWithFloodFill(painter, traits);
}

template <typename Arr_>
template <typename TTraits>
void ArrangementGraphicsItem<Arr_>::paintWithFloodFill(
  QPainter* painter, const TTraits& traits)
{
  auto windowRect = painter->window();
  auto width = windowRect.width();
  auto height = windowRect.height();
  if (tempImage.width() != width || tempImage.height() != height)
    tempImage = {width, height, QImage::Format_ARGB32};

  QRgb* st = reinterpret_cast<QRgb*>(tempImage.bits());
  // Transparent (0) is used as a flag for an invalid color
  std::fill(st, st + width * height, 0);
  // draw margins with white
  // prevents the flood algorithm from going over the border
  // also useful with algebraic curves as sometimes pixels on the borders
  // aren't painted, and the flood algorithm 'leaks'
  static constexpr QRgb white = 0xFFFFFFFF;
  for (uint16_t i = 0; i < margin; i++)
  {
    for (uint16_t j = 0; j < width; j++)
    {
      st[i * width + j] = white;
      st[(height - 1 - i) * width + j] = white;
    }
    for (uint16_t j = 0; j < height; j++)
    {
      st[j * width + i] = white;
      st[j * width + (width - 1 - i)] = white;
    }
  }

  QPainter painter2{&tempImage};
  painter2.setTransform(painter->transform());
  painter2.setPen(this->edgesPen);

  // paint bounded faces normally?
  // by experimenting it's faster to just paint all using the flood algo
  // specially with algebraic faces since currenlty all edges have to
  // be recalculated/rendered again for faces
  // this->paintFaces(&painter2);
  this->paintEdges(&painter2, traits);
  // paint remaining faces using floodfill
  this->paintFacesFloodFill(&painter2, tempImage);

  auto& painterTransform = painter->transform();
  painter->resetTransform();
  painter->drawImage(QPoint{0, 0}, tempImage);
  painter->setTransform(painterTransform);
}

template <typename Arr_>
template <typename TTraits>
void ArrangementGraphicsItem<Arr_>::paintEdges(
  QPainter* painter, const TTraits&)
{
  auto painterOstream =
    ArrangementPainterOstream<TTraits>(painter, this->boundingRect());

  painterOstream.setScene(this->getScene());
  painter->setPen(this->edgesPen);
  painter->setBrush(::Qt::transparent);

  for (auto it = this->arr->edges_begin(); it != this->arr->edges_end(); ++it)
  {
    X_monotone_curve_2 curve = it->curve();
    painterOstream << curve;
  }
}

template <typename Arr_>
template <typename Coefficient_>
void ArrangementGraphicsItem<Arr_>::paintEdges(
  QPainter* painter, const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>&)
{
  using TTraits = CGAL::Arr_algebraic_segment_traits_2<Coefficient_>;

  auto painterOstream =
    ArrangementPainterOstream<TTraits>(painter, this->boundingRect());

  painterOstream.setScene(this->getScene());
  painter->setPen(this->edgesPen);
  painter->setBrush(::Qt::transparent);

  painterOstream.paintEdges(arr->edges_begin(), arr->edges_end());
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::paintFacesFloodFill(
  QPainter* painter, QImage& image)
{
  static constexpr QRgb invalid_rgb = 0;

  QRgb* raw_img = reinterpret_cast<QRgb*>(image.bits());
  uint16_t width = image.width();
  uint16_t height = image.height();

  QTransform deviceToScene = painter->deviceTransform().inverted();

  PointLocationFunctions<Arrangement> point_location_functions;
  auto get_face = [&](auto x, auto y) {
    QPointF point =
      deviceToScene.map(QPointF{static_cast<qreal>(x), static_cast<qreal>(y)});
    return point_location_functions.getFace(this->arr, point);
  };

  auto paint_face = [this, raw_img, width, get_face](auto x, auto y) {
    auto face = get_face(x, y);
    QRgb color = face->color().rgb();
    this->flood_fill(raw_img, width, x, y, color);
  };

  static constexpr int tot_margin = margin + 2;
  auto cur_img_line = raw_img + tot_margin * width;
  for (uint16_t y = tot_margin; y + 1 + tot_margin < height; y++)
  {
    for (uint16_t x = tot_margin; x + 1 + tot_margin < width; x++)
    {
      // just to account for rendering errors
      // make sure the pixel falls in the right face
      if (
        cur_img_line[x] == invalid_rgb &&
        cur_img_line[x + 1] == invalid_rgb &&
        cur_img_line[x - 1] == invalid_rgb &&
        cur_img_line[x + 2] == invalid_rgb &&
        cur_img_line[x - 2] == invalid_rgb &&
        cur_img_line[width + x] == invalid_rgb &&
        cur_img_line[-width + x] == invalid_rgb &&
        cur_img_line[2 * width + x] == invalid_rgb &&
        cur_img_line[- 2 * width + x] == invalid_rgb &&
        cur_img_line[width + x + 1] == invalid_rgb &&
        cur_img_line[width + x - 1] == invalid_rgb &&
        cur_img_line[-width + x + 1] == invalid_rgb &&
        cur_img_line[-width + x - 1] == invalid_rgb &&
        cur_img_line[width + x + 2] == invalid_rgb &&
        cur_img_line[width + x - 2] == invalid_rgb &&
        cur_img_line[-width + x + 2] == invalid_rgb &&
        cur_img_line[-width + x - 2] == invalid_rgb &&
        cur_img_line[2 * width + x + 1] == invalid_rgb &&
        cur_img_line[2 * width + x - 1] == invalid_rgb &&
        cur_img_line[-2 * width + x + 1] == invalid_rgb &&
        cur_img_line[-2 * width + x - 1] == invalid_rgb &&
        cur_img_line[2 * width + x + 2] == invalid_rgb &&
        cur_img_line[2 * width + x - 2] == invalid_rgb &&
        cur_img_line[-2 * width + x + 2] == invalid_rgb &&
        cur_img_line[-2 * width + x - 2] == invalid_rgb)
        paint_face(x, y);
    }
    cur_img_line += width;
  }
}

template < typename Arr_ >
void ArrangementGraphicsItem< Arr_ >::updateBoundingBox( )
{
  this->prepareGeometryChange();

  this->bb = {};
  ConstructBoundingBox<Traits> construct_bounding_box;
  for (auto it = this->arr->edges_begin(); it != this->arr->edges_end(); ++it)
    this->bb += construct_bounding_box(it->curve());

  for (auto it = this->arr->vertices_begin(); it != this->arr->vertices_end();
       ++it)
    this->bb += construct_bounding_box(it->point());
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::updatePointsItem()
{
  this->pointsGraphicsItem->clear();
  for (auto it = this->arr->vertices_begin(); it != this->arr->vertices_end();
       ++it)
  {
    Point_2 p = it->point();
    this->pointsGraphicsItem->insert(p);
  }
}

template < typename Arr_ >
void ArrangementGraphicsItem< Arr_ >::modelChanged( )
{
  this->updatePointsItem();
  this->updateBoundingBox( );
  this->update( );
}

template < typename Arr_ >
void
ArrangementGraphicsItem< Arr_ >::
paintFace( Face_handle f, QPainter* painter )
{
  if (f->visited()) return;

  Hole_iterator hit; // holes iterator
  this->paintFace(f, painter, *(arr->traits()));
  f->set_visited(true);

  for (hit = f->holes_begin(); hit != f->holes_end(); ++hit)
  {
    // Traverse in clockwise order
    Ccb_halfedge_circulator cc = *hit;
    do
    {
      Halfedge_handle he = cc;
      Halfedge_handle he2 = he->twin();
      Face_handle inner_face = he2->face();
      if (this->antenna(he)) continue;

      this->visit_ccb_faces(inner_face, painter);
    } while (++cc != *hit);
  }
}

template <typename Arr_>
template <typename Kernel_>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter, const CGAL::Arr_segment_traits_2<Kernel_>&)
{
  if (!f->is_unbounded())  // f is not the unbounded face
  {
    QVector< QPointF > pts; // holds the points of the polygon

    /* running with around the outer of the face and generate from it
     * polygon
     */
    Ccb_halfedge_circulator cc=f->outer_ccb();
    do {
      double x = CGAL::to_double(cc->source()->point().x());
      double y = CGAL::to_double(cc->source()->point().y());
      QPointF coord_source(x , y);
      pts.push_back(coord_source );
      //created from the outer boundary of the face
    } while (++cc != f->outer_ccb());

    // make polygon from the outer ccb of the face 'f'
    QPolygonF pgn (pts);

    QColor color = f->color();
    painter->setBrush(color);
    painter->drawPolygon( pgn );
  }
  else
  {
    QRectF rect = this->viewportRect( );
    QColor color = f->color();
    painter->fillRect(rect, color);
  }
}

template <typename Arr_>
template <typename Kernel_>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter, const CGAL::Arr_polyline_traits_2<Kernel_>&)
{
  if (!f->is_unbounded())
  {
    QVector< QPointF > pts; // holds the points of the polygon

    CGAL::Qt::Converter<Kernel_> convert;
    /* running with around the outer of the face and generate from it
     * polygon
     */
    Ccb_halfedge_circulator cc = f->outer_ccb();
    do {
      if (this->antenna(cc)) continue;

      QPointF src = convert(cc->source()->point());
      QPointF tgt = convert(cc->target()->point());

      Halfedge he = *cc;
      if (he.direction() == ARR_LEFT_TO_RIGHT && src.x() > tgt.x())
        std::swap(src, tgt);
      else if (he.direction() == ARR_RIGHT_TO_LEFT && src.x() < tgt.x())
        std::swap(src, tgt);

      auto&& curve = cc->curve();
      auto first_subcurve = curve.subcurves_begin();
      QPointF src_first = convert(first_subcurve->source());

      if (src_first == src)
      {
        for (auto it = curve.subcurves_begin(); it != curve.subcurves_end();
             ++it)
        {
          pts.push_back(convert(it->source()));
          pts.push_back(convert(it->target()));
        }
      }
      else
      {
        QVector<QPointF> pts_tmp;
        for (auto it = curve.subcurves_begin(); it != curve.subcurves_end();
             ++it)
        {
          pts_tmp.push_front(convert(it->source()));
          pts_tmp.push_front(convert(it->target()));
        }
        pts += pts_tmp;
      }
    } while (++cc != f->outer_ccb());

    // make polygon from the outer ccb of the face 'f'
    QPolygonF pgn( pts );

    painter->setBrush(f->color());
    painter->drawPolygon( pgn );
  }
  else
  {
    QRectF rect = this->viewportRect( );
    QColor color = f->color();
    painter->fillRect(rect, color);
  }
}

template <typename Arr_>
template <typename Coefficient_>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>& /* traits */)
{
  if (f->is_unbounded()) return;

  ArrangementPainterOstream<Traits> painterOstream{
    painter, this->boundingRect()};
  painterOstream.setScene(this->getScene());

  painter->save();
  painter->setTransform(painterOstream.getPointsListMapping());

  QVector<QPointF> pts;
  /* running with around the outer of the face and generate from it
   * polygon
   */
  Ccb_halfedge_circulator cc = f->outer_ccb();
  do
  {
    if (this->antenna(cc)) continue;

    auto points = painterOstream.getPointsList(cc->curve());
    if (points.empty()) continue;

    double src_x = CGAL::to_double(cc->source()->point().x());
    double tgt_x = CGAL::to_double(cc->target()->point().x());

    if (src_x < tgt_x)
    {
      for (auto& vec : points)
          for(auto& vit : vec)
            pts.push_back({vit.first, vit.second});
    }
    else
    {
      for (auto vecit = points.rbegin(); vecit != points.rend(); ++vecit)
        for (auto vit = vecit->rbegin(); vit != vecit->rend(); ++vit)
          pts.push_back({vit->first, vit->second});
    }
  } while (++cc != f->outer_ccb());

  // make polygon from the outer ccb of the face 'f'
  QPolygonF pgn(pts);
  painter->setBrush(f->color());
  painter->drawPolygon(pgn);

  painter->restore();
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::paintFaces(QPainter* painter)
{
  QPen pen = painter->pen();
  painter->setPen(this->facesPen);

  // Prepare all faces for painting
  for (auto fi = this->arr->faces_begin(); fi != this->arr->faces_end(); ++fi)
    fi->set_visited(false);

  for (auto fi = this->arr->faces_begin(); fi != this->arr->faces_end(); ++fi)
    this->paintFace(fi, painter);

  painter->setPen(pen);
}

template <typename Arr_>
void ArrangementGraphicsItem<Arr_>::visit_ccb_faces(
  Face_handle& fh, QPainter* painter)
{
  this->paintFace(fh, painter);

  Ccb_halfedge_circulator cc = fh->outer_ccb();
  do
  {
    Halfedge he = *cc;
    if (!he.twin()->face()->visited())
    {
      Face_handle nei = he.twin()->face();
      this->visit_ccb_faces(nei, painter);
    }
    // created from the outer boundary of the face
  } while (++cc != fh->outer_ccb());
}

template <typename Arr_>
bool ArrangementGraphicsItem<Arr_>::antenna(Halfedge_handle h)
{
  Halfedge_handle twin = h->twin();
  return (twin->face() == h->face());
}

template <typename Arr_>
template <typename ArrTraits>
void ArrangementGraphicsItem<Arr_>::paintFace(Face_handle, QPainter*, ArrTraits)
{
}

template <typename Arr_>
template <typename RatKernel, typename AlgKernel, typename NtTraits>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  const CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>&)
{
  if (!f->is_unbounded()) // f is not the unbounded face
  {
    QVector<QPointF> pts; // holds the points of the polygon
    /* running with around the outer of the face and generate from it
     * polygon
     */
    Ccb_halfedge_circulator cc = f->outer_ccb();
    do
    {
      if (this->antenna(cc)) { continue; }

      Halfedge_handle he = cc;
      X_monotone_curve_2 c = he->curve();
      // Get the co-ordinates of the curve's source and target.
      double sx = CGAL::to_double(he->source()->point().x()),
             sy = CGAL::to_double(he->source()->point().y()),
             tx = CGAL::to_double(he->target()->point().x()),
             ty = CGAL::to_double(he->target()->point().y());

      QPointF coord_source(sx, sy);
      QPointF coord_target(tx, ty);

      // Transform the point coordinates from general coordinate system to
      // Qt scene coordinate system
      QPoint coord_source_viewport = this->fromScene(coord_source);
      QPoint coord_target_viewport = this->fromScene(coord_target);

      if (c.orientation() == CGAL::COLLINEAR) { pts.push_back(coord_source); }
      else
      {
        // If the curve is monotone, than its source and its target has the
        // extreme x co-ordinates on this curve.
        bool is_source_left = (sx < tx);
        int x_min = is_source_left ? coord_source_viewport.x()
                                   : coord_target_viewport.x();
        int x_max = is_source_left ? coord_target_viewport.x()
                                   : coord_source_viewport.x();

        pts.push_back(coord_source);

        // Draw the curve as pieces of small segments
        const int DRAW_FACTOR = 5;
        int start;
        int end;
        int step;
        if (is_source_left)
        {
          start = x_min + DRAW_FACTOR;
          end = x_max;
          step = DRAW_FACTOR;
        }
        else
        {
          start = x_max;
          end = x_min;
          step = -DRAW_FACTOR;
        }
        for (int x = start; x < end; x += step)
        {
          double curr_x = this->toScene(QPoint{x, 0}).x();
          AlgKernel ker;
          Point_2 curr_p(curr_x, 0);

          // If curr_x > x_max or curr_x < x_min
          if (!(ker.compare_x_2_object()(curr_p, c.left()) != CGAL::SMALLER &&
                ker.compare_x_2_object()(curr_p, c.right()) != CGAL::LARGER))
          { continue; }

          auto px = c.point_at_x(curr_p);
          double curr_y = CGAL::to_double(px.y());
          QPointF curr(curr_x, curr_y);
          pts.push_back(curr);
        }
        pts.push_back(coord_target);
      }
    } while (++cc != f->outer_ccb());

    QPolygonF pgn(pts);
    painter->setBrush(f->color());
    painter->drawPolygon(pgn);
  }
  else
  {
    QRectF rect = this->viewportRect();
    QColor color = f->color();
    painter->fillRect(rect, color);
  }
}

template <typename Arr_>
template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  const CGAL::Arr_Bezier_curve_traits_2<
    RatKernel, AlgKernel, NtTraits, BoundingTraits>&)
{
  if (!f->is_unbounded())
  {
    ArrangementPainterOstream<Traits> painterOstream{
      painter, this->boundingRect()};
    painterOstream.setScene(this->getScene());

    QVector<QPointF> pts;
    Ccb_halfedge_circulator cc = f->outer_ccb();
    do
    {
      if (this->antenna(cc)) continue;

      const X_monotone_curve_2& curve = cc->curve();
      auto points = painterOstream.getPoints(curve);

      if (cc->source()->point().is_same(curve.source()))
        for (auto& vit : points)
          pts.push_back({vit.first, vit.second});
      else
        for (auto it = points.rbegin(); it != points.rend(); ++it)
          pts.push_back({it->first, it->second});

    } while(++cc != f->outer_ccb());

    QPolygonF pgn(pts);
    QColor color = f->color();
    painter->setBrush(color);
    painter->drawPolygon(pgn);
  }
  else
  {
    QRectF rect = this->viewportRect();
    QColor color = f->color();
    painter->fillRect(rect, color);
  }
}

template <typename Arr_>
template <typename Kernel_>
void ArrangementGraphicsItem<Arr_>::paintFace(
  Face_handle f, QPainter* painter,
  const CGAL::Arr_linear_traits_2<Kernel_>& /* traits */)
{
  QVector<QPointF> pts; // holds the points of the polygon
  QColor color = f->color();
  painter->setBrush(color);
  if (!f->is_unbounded())
  {
    /* running with around the outer of the face and generate from it
     * polygon
     */
    Ccb_halfedge_circulator cc = f->outer_ccb();
    do
    {
      Halfedge_handle hh = cc;
      if (this->antenna(hh)) continue;

      double x = CGAL::to_double(cc->source()->point().x());
      double y = CGAL::to_double(cc->source()->point().y());
      QPointF coord_source(x, y);
      pts.push_back(coord_source);
    } while (++cc != f->outer_ccb());

    painter->drawPolygon(QPolygonF{pts});
  }
  else
  {
    std::vector<Halfedge_handle> halfedges;
    Ccb_halfedge_circulator cc = f->outer_ccb();
    do
    {
      Halfedge_handle hh = cc;
      if (cc->is_fictitious() || this->antenna(hh))
        continue;

      halfedges.push_back(hh);
    } while (++cc != f->outer_ccb());

    QRectF rect = this->viewportRect();

    switch (halfedges.size())
    {
    // one unbounded face
    case 0: {
      painter->fillRect(rect, color);
      break;
    }
    // two unbounded faces
    case 1: {
      auto&& curve = halfedges.front()->curve();
      auto&& line = curve.line();

      // horizontal or mostly horizontal
      if (CGAL::abs(line.b()) > CGAL::abs(line.a()))
      {
        qreal left_y = CGAL::to_double(line.y_at_x(rect.left()));
        qreal right_y = CGAL::to_double(line.y_at_x(rect.right()));

        pts.push_back({rect.left(), left_y});
        pts.push_back({rect.right(), right_y});
        if (cc->direction() == ARR_LEFT_TO_RIGHT)
        {
          pts.push_back(rect.bottomRight());
          pts.push_back(rect.bottomLeft());
        }
        else
        {
          pts.push_back(rect.topRight());
          pts.push_back(rect.topLeft());
        }
      }
      // vertical or mostly vertical
      else
      {
        qreal top_x = CGAL::to_double(line.x_at_y(rect.top()));
        qreal bottom_x = CGAL::to_double(line.x_at_y(rect.bottom()));

        pts.push_back({top_x, rect.top()});
        pts.push_back({bottom_x, rect.bottom()});
        if (cc->direction() == ARR_LEFT_TO_RIGHT)
        {
          pts.push_back(rect.bottomLeft());
          pts.push_back(rect.topLeft());
        }
        else
        {
          pts.push_back(rect.bottomRight());
          pts.push_back(rect.topRight());
        }
      }

      painter->drawPolygon(QPolygonF{pts});
      break;
    }
    // general case
    default: {
      using Vertex_handle = typename Arrangement::Vertex_handle;

      auto handle_vertex = [&](Halfedge_handle halfedge, Vertex_handle vertex)
      {
        // dx * y = dy * x + c
        typename Kernel_::RT dx, dy, c;

        auto&& curve = halfedge->curve();
        if (curve.is_segment())
        {
          auto&& src = halfedge->source()->point();
          auto&& tgt = halfedge->target()->point();
          dy = tgt.y() - src.y();
          dx = tgt.x() - src.x();
          c = dx * src.y() - dy * src.x();
        }
        else if(curve.is_ray())
        {
          auto&& ray = curve.ray();
          auto&& dir = ray.direction();
          auto&& src = ray.source();
          dy = dir.dy();
          dx = dir.dx();
          c = dx * src.y() - dy * src.x();
        }
        // halfedge can be a line in case of a face between two parallel lines
        else
        {
          auto&& line = curve.line();
          auto param_in_x = vertex->parameter_space_in_x();
          auto param_in_y = vertex->parameter_space_in_y();
          dy = line.a();
          dx = -line.b();
          c = line.c();

          if (
            (param_in_x == LEFT_BOUNDARY && CGAL::is_positive(dx)) ||
            (param_in_x == RIGHT_BOUNDARY && CGAL::is_negative(dx)) ||
            (param_in_y == BOTTOM_BOUNDARY && CGAL::is_positive(dy)) ||
            (param_in_y == TOP_BOUNDARY && CGAL::is_negative(dy)))
          {
            dx = -dx;
            dy = -dy;
            c = -c;
          }
        }
        // horizontal or mostly horizontal
        if (CGAL::abs(dx) > CGAL::abs(dy))
        {
          qreal x = CGAL::is_positive(dx) ? rect.right() + rect.width()
                                          : rect.left() - rect.width();
          // in case vertex is a finite point
          if (!vertex->is_at_open_boundary())
          {
            auto p = vertex->point();
            // check not far away
            if (
              (CGAL::is_positive(dx) && p.x() < x) ||
              (CGAL::is_negative(dx) && p.x() > x))
            {
              pts.push_back({CGAL::to_double(p.x()), CGAL::to_double(p.y())});
              return;
            }
          }
          // in case vertex is far away (finite or at infinity), get a nearer
          // point on the line
          // i.e. in case farther than twice the width of the viewport
          qreal y = CGAL::to_double((dy / dx) * x + (c / dx));
          pts.push_back({x, y});
        }
        // vertical or mostly vertical
        else
        {
          qreal y = CGAL::is_positive(dy) ? rect.bottom() + rect.height()
                                          : rect.top() - rect.height();
          // in case vertex is a finite point
          if (!vertex->is_at_open_boundary())
          {
            auto p = vertex->point();
            // check not far away
            if (
              (CGAL::is_positive(dy) && p.y() < y) ||
              (CGAL::is_negative(dy) && p.y() > y))
            {
              pts.push_back({CGAL::to_double(p.x()), CGAL::to_double(p.y())});
              return;
            }
          }
          // in case vertex is far away (finite or at infinity), get a nearer
          // point on the line
          // i.e. in case farther than twice the height of the viewport
          qreal x = CGAL::to_double((dx / dy) * y - (c / dy));
          pts.push_back({x, y});
        }
      };

      for (auto hh : halfedges)
      {
        handle_vertex(hh, hh->source());
        handle_vertex(hh, hh->target());
      }

      painter->drawPolygon(QPolygonF{pts});
      break;
    }
    }
  }
}

// TODO: clean all this portion up (move it somewhere?), it's ugly!
static CGAL::Bbox_2 reject_not_in_allowable_range(
  const CGAL::Bbox_2& box, const CGAL::Bbox_2& allowable_range)
{
  double xmin = std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();

  if (box.xmin() > allowable_range.xmin()) xmin = box.xmin();
  if (box.ymin() > allowable_range.ymin()) ymin = box.ymin();
  if (box.xmax() < allowable_range.xmax()) xmax = box.xmax();
  if (box.ymax() < allowable_range.ymax()) ymax = box.ymax();
  return {xmin, ymin, xmax, ymax};
}

static bool isFinite(const CGAL::Bbox_2& box)
{
  return !std::isinf(box.xmin()) && !std::isinf(box.xmax()) &&
         !std::isinf(box.ymin()) && !std::isinf(box.ymax());
}

static CGAL::Bbox_2 addMargins(const CGAL::Bbox_2& box)
{
  // add margin to bounding box
  double x_margin;
  double y_margin;
  if (box.xmin() == box.xmax() || box.ymin() == box.ymax())
  {
    static constexpr float const_margin = 50;
    x_margin = const_margin;
    y_margin = const_margin;
  }
  else
  {
    static constexpr double prop_margin = 0.10;
    x_margin = (box.xmax() - box.xmin()) * prop_margin;
    y_margin = (box.ymax() - box.ymin()) * prop_margin;
  }
  return {
    box.xmin() - x_margin, box.ymin() - y_margin, box.xmax() + x_margin,
    box.ymax() + y_margin};
}

template <class Arr_>
CGAL::Bbox_2
findOtherInterestingPoints(const Arr_*, const CGAL::Bbox_2&)
{
  return {};
}

#ifdef CGAL_USE_CORE
template <typename Coefficient_>
static const auto&
getXyCurves(const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>* traits)
{
  // the traits object is only needed the first time
  // this assumes that X_monotone_curves created from the first traits object
  // will work with arrangements with a different object
  using Traits = CGAL::Arr_algebraic_segment_traits_2<Coefficient_>;
  static std::vector<typename Traits::X_monotone_curve_2> xy_curves;
  if (xy_curves.empty())
  {
    typedef typename Traits::Polynomial_2 Polynomial_2;
    auto construct_curve = traits->construct_curve_2_object();
    auto make_x_monotone = traits->make_x_monotone_2_object();

    Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
    Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);
    auto x_cv = construct_curve(x);
    auto y_cv = construct_curve(y);

    std::vector<CGAL::Object> arcs;
    make_x_monotone(x_cv, std::back_inserter(arcs));
    make_x_monotone(y_cv, std::back_inserter(arcs));
    for (auto& arc_obj : arcs)
    {
      typename Traits::X_monotone_curve_2 arc;
      CGAL::assign(arc, arc_obj);
      xy_curves.push_back(arc);
    }
  }
  return xy_curves;
}

template <>
CGAL::Bbox_2 findOtherInterestingPoints<demo_types::DemoTypes::Alg_seg_arr>(
  const demo_types::DemoTypes::Alg_seg_arr* arr,
  const CGAL::Bbox_2& allowable_range)
{
  using Traits = demo_types::DemoTypes::Alg_seg_traits;
  CGAL::Bbox_2 bb = {};
  std::vector<CGAL::Object> intersections;
  for (auto it = arr->edges_begin(); it != arr->edges_end(); ++it)
  {
    for (auto& arc : getXyCurves(arr->traits()))
      if (arc.is_vertical() != it->curve().is_vertical())
        it->curve().intersections(arc, std::back_inserter(intersections));
  }
  for (auto it = intersections.begin(); it != intersections.end(); it++)
  {
    std::pair<typename Traits::Point_2, unsigned int> point_multiplicity;
    CGAL::assign(point_multiplicity, *it);
    auto& point = point_multiplicity.first;
    if (point.location() == CGAL::ARR_INTERIOR)
    {
      auto xy = point.to_double();
      bb += reject_not_in_allowable_range(
        {xy.first, xy.second, xy.first, xy.second}, allowable_range);
    }
  }
  return bb;
}
#endif // CGAL_USE_CORE


template <typename Arr_>
QRectF ArrangementGraphicsItem<Arr_>::getInterestingViewport() const
{
  QRectF scene_rect = this->getScene()->sceneRect();
  CGAL::Bbox_2 scene_bbox = {
    scene_rect.left(), scene_rect.top(), scene_rect.right(),
    scene_rect.bottom()};

  ConstructBoundingBox<Traits> construct_bounding_box;
  CGAL::Bbox_2 bounding_box = {};
  for (auto it = this->arr->edges_begin();
       it != this->arr->edges_end(); ++it)
  {
    bounding_box += reject_not_in_allowable_range(
      construct_bounding_box(it->curve()), scene_bbox);
  }
  for (auto it = this->arr->vertices_begin();
       it != this->arr->vertices_end(); ++it)
  {
    bounding_box += reject_not_in_allowable_range(
      construct_bounding_box(it->point()), scene_bbox);
  }

  if (!isFinite(bounding_box))
    bounding_box += findOtherInterestingPoints(this->arr, scene_bbox);

  // ideally this should happen only if the arrangement is empty
  // (if findOtherInterestingPoints worked for all arrangement types)
  if (!isFinite(bounding_box)) bounding_box += {0, 0, 0, 0};

  bounding_box = addMargins(bounding_box);

  return QRectF(
    bounding_box.xmin(), bounding_box.ymin(),
    bounding_box.xmax() - bounding_box.xmin(),
    bounding_box.ymax() - bounding_box.ymin());
}

} // namespace QT
} // namespace CGAL
