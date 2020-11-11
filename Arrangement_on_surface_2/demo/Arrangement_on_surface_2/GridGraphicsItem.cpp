// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#include "GridGraphicsItem.h"
#include <cmath>
#include <limits>

#include <QGraphicsView>
#include <QPainter>
#include <QtGlobal>

GridGraphicsItem::GridGraphicsItem() :
    gridColor{::Qt::gray}, axesColor{::Qt::black},
    labelsColor{::Qt::black}, spacing{75}
{
}

int GridGraphicsItem::getXPower5() const { return x_power5; }

int GridGraphicsItem::getXPower2() const { return x_power2; }

int GridGraphicsItem::getYPower5() const { return y_power5; }

int GridGraphicsItem::getYPower2() const { return y_power2; }

float GridGraphicsItem::getXUnit() const { return x_unit; }

float GridGraphicsItem::getYUnit() const { return y_unit; }

QColor GridGraphicsItem::getGridColor() const { return this->gridColor; }

QColor GridGraphicsItem::getAxesColor() const { return this->axesColor; }

QColor GridGraphicsItem::getLabelsColor() const { return this->labelsColor; }

int GridGraphicsItem::getSpacing() const { return this->spacing; }

void GridGraphicsItem::setGridColor(const QColor& color)
{
  this->gridColor = color;
}

void GridGraphicsItem::setAxesColor(const QColor& color)
{
  this->axesColor = color;
}

void GridGraphicsItem::setLabelsColor(const QColor& color)
{
  this->labelsColor = color;
}

void GridGraphicsItem::setSpacing(int spacing_)
{
  this->spacing = spacing_;
  this->update();
}

static inline qreal
horizontalAdvance(const QFontMetrics& fm, const QString& text)
{
#if (QT_VERSION >= QT_VERSION_CHECK(5, 11, 0))
  return fm.horizontalAdvance(text);
#else
  return fm.boundingRect(text).width();
#endif
}

void GridGraphicsItem::paint(
  QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*)
{
  // equivalent to m11() and m22() when there are no rotations
  QTransform worldTransform = painter->transform();
  QLineF ux_line = worldTransform.map(QLineF{0, 0, 1, 0});
  QLineF uy_line = worldTransform.map(QLineF{0, 0, 0, 1});
  float scaleX = ux_line.length();
  float scaleY = uy_line.length();

  float gridSceneMaxXLen = spacing / scaleX;
  float gridSceneMaxYLen = spacing / scaleY;

  // 2^power2 * 5^power5 = gridSceneMaxSize
  x_power5 = std::log(gridSceneMaxXLen) / std::log(5.);
  float pow5l = std::pow(5, x_power5);
  x_power2 = std::log(gridSceneMaxXLen / pow5l) / std::log(2.);
  float pow2l = std::pow(2, x_power2);

  x_unit = pow2l * pow5l;

  y_power5 = std::log(gridSceneMaxYLen) / std::log(5.);
  y_power2 = std::log(gridSceneMaxYLen / pow5l) / std::log(2.);
  pow5l = std::pow(5, y_power5);
  pow2l = std::pow(2, y_power2);

  y_unit = pow2l * pow5l;

  QRectF sceneViewport = this->viewportRect(painter);

  QVarLengthArray<QLineF, 100> linesX;
  QVarLengthArray<QLineF, 100> linesY;

  qreal left = int(sceneViewport.left() / x_unit) * x_unit;
  qreal top = int(sceneViewport.top() / y_unit) * y_unit;

  for (qreal x = left; x < sceneViewport.right(); x += x_unit)
    linesX.append(QLineF(x, sceneViewport.top(), x, sceneViewport.bottom()));

  for (qreal y = top; y < sceneViewport.bottom(); y += y_unit)
    linesY.append(QLineF(sceneViewport.left(), y, sceneViewport.right(), y));

  // set up the painter
  QPen gridPen;
  gridPen.setColor(this->gridColor);
  gridPen.setCosmetic(true);
  painter->setPen(gridPen);

  // draw the grid
  painter->drawLines(linesX.data(), linesX.size());
  painter->drawLines(linesY.data(), linesY.size());

  // draw x and y axis
  QPen axisPen{this->axesColor};
  axisPen.setCosmetic(true);
  painter->setPen(axisPen);

  const double marginX = 15 / scaleX;
  const double marginY = 15 / scaleY;

  double xAxis_y = 0;
  double yAxis_x = 0;
  bool xAxisTop = false;
  bool yAxisLeft = false;

  if (sceneViewport.left() > -marginX)
  {
    yAxis_x = sceneViewport.left() + marginX;
    yAxisLeft = true;
  }
  else if (sceneViewport.right() < marginX)
  {
    yAxis_x = sceneViewport.right() - marginX;
  }

  if (sceneViewport.top() > -marginY)
  {
    xAxis_y = sceneViewport.top() + marginY;
    xAxisTop = true;
  }
  else if (sceneViewport.bottom() < marginY)
  {
    xAxis_y = sceneViewport.bottom() - marginY;
  }

  painter->drawLine(
    QLineF{sceneViewport.left(), xAxis_y, sceneViewport.right(), xAxis_y});
  painter->drawLine(
    QLineF{yAxis_x, sceneViewport.bottom(), yAxis_x, sceneViewport.top()});

  for (qreal x = left; x < sceneViewport.right(); x += x_unit)
    painter->drawLine(QLineF(x, xAxis_y - 5 / scaleY, x, xAxis_y + 5 / scaleY));

  for (qreal y = top; y < sceneViewport.bottom(); y += y_unit)
    painter->drawLine(QLineF(yAxis_x - 5 / scaleX, y, yAxis_x + 5 / scaleX, y));

  QPen labelsPen{this->labelsColor};
  axisPen.setCosmetic(true);
  painter->setPen(labelsPen);

  painter->resetTransform();

  auto&& font = painter->font();
  QFontMetrics fm(font);
  qreal txtHeight = fm.height();

  QPointF uy = {uy_line.dx(), uy_line.dy()};
  float sign = 1;
  if (!xAxisTop) sign = -1;
  for (qreal x = left; x < sceneViewport.right(); x += x_unit)
  {
    if (std::abs(x) < x_unit / 2) continue;
    QString txt = QString{"%1"}.arg(x);
    qreal txtWidth = horizontalAdvance(fm, txt) + 4;
    QPointF pos =
      worldTransform.map(QPointF{x, xAxis_y}) +
      sign * QPointF{uy.x() * txtWidth, uy.y() * txtHeight} / scaleY;

    painter->drawText(
      pos.x() - txtWidth / 2, pos.y() - txtHeight / 2, txtWidth, txtHeight,
      Qt::AlignHCenter | Qt::AlignVCenter, txt);
  }

  QPointF ux = {ux_line.dx(), ux_line.dy()};
  sign = 1;
  if (!yAxisLeft) sign = -1;
  for (qreal y = top; y < sceneViewport.bottom(); y += y_unit)
  {
    if (std::abs(y) < y_unit / 2) continue;
    QString txt = QString{"%1"}.arg(y);
    qreal txtWidth = horizontalAdvance(fm, txt) + 4;
    QPointF pos =
      worldTransform.map(QPointF{yAxis_x, y}) +
      sign * QPointF{ux.x() * txtWidth, ux.y() * txtHeight} / scaleX;

    painter->drawText(
      pos.x() - txtWidth / 2, pos.y() - txtHeight / 2, txtWidth, txtHeight,
      Qt::AlignHCenter | Qt::AlignVCenter, txt);
  }

  // revert the painter
  painter->setTransform(worldTransform);
}

QRectF GridGraphicsItem::boundingRect() const { return this->viewportRect(); }

QRectF GridGraphicsItem::viewportRect(QPainter* painter) const
{
  QRectF pixelViewport = painter->viewport();
  return painter->transform().inverted().mapRect(pixelViewport);
}

QRectF GridGraphicsItem::viewportRect() const
{
  QList<QGraphicsView*> views = this->scene()->views();
  if (views.size() == 0) return {};
  // assumes the first view is the right one
  QGraphicsView* view = views.first();
  return view->mapToScene(QRect{0, 0, view->width(), view->height()})
    .boundingRect();
}
