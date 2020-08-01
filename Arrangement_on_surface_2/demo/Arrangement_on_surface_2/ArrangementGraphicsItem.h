// Copyright (c) 2008, 2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>

#ifndef CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
#define CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <QImage>

#include "GraphicsSceneMixin.h"
#include "PointsGraphicsItem.h"

template <typename T>
class ArrangementPainterOstream;

namespace CGAL
{
template <typename T, typename U, typename I, typename S>
class Arr_Bezier_curve_traits_2;
template <typename T, typename U, typename I>
class Arr_conic_traits_2;
template <typename T>
class Arr_linear_traits_2;
template <typename T>
class Arr_algebraic_segment_traits_2;
template <typename T>
class Arr_segment_traits_2;
template <typename T>
class Arr_polyline_traits_2;
template <typename RatK, typename AlgK, typename Nt, typename BoundingTratits>
class Arr_Bezier_curve_traits_2;
} // namespace CGAL

namespace CGAL {
namespace Qt {

class ArrangementGraphicsItemBase :
    public GraphicsItem,
    public QGraphicsSceneMixin
{
public:
  ArrangementGraphicsItemBase();

  const QPen& getVerticesPen() const;
  const QPen& getEdgesPen() const;
  void setVerticesPen(const QPen& pen);
  void setEdgesPen(const QPen& pen);
  void setBackgroundColor(QColor color);

protected:
  CGAL::Bbox_2 bb;
  QPen verticesPen;
  QPen edgesPen;
  QPen facesPen;
  PointsGraphicsItem pointsGraphicsItem;
}; // class ArrangementGraphicsItemBase

template <typename Arr_>
class ArrangementGraphicsItem : public ArrangementGraphicsItemBase
{
  typedef ArrangementGraphicsItemBase                   Superclass;
  typedef Arr_                                          Arrangement;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef typename Arrangement::Edge_iterator           Edge_iterator;
  typedef typename Arrangement::Halfedge                Halfedge;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Face_handle             Face_handle;
  typedef typename Arrangement::Face_iterator           Face_iterator;
  typedef typename Arrangement::Unbounded_face_iterator Unbounded_face_iterator;
  typedef typename Arrangement::Hole_iterator           Holes_iterator;
  typedef typename Arrangement::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Point_2                      Point_2;

public:
  /*! Constructor */
  ArrangementGraphicsItem( Arrangement* t_ );

  /*! Destructor (virtual) */
  ~ArrangementGraphicsItem() {}

public:
  void modelChanged();
  QRectF boundingRect() const override;
  void paint(
    QPainter* painter, const QStyleOptionGraphicsItem* option,
    QWidget* widget) override;

protected:
  void updatePointsItem();

  template <typename TTraits>
  void paint(QPainter*, TTraits);

  template <typename Coefficient_>
  void paint(QPainter*, CGAL::Arr_algebraic_segment_traits_2<Coefficient_>);

  void updateBoundingBox();

  template <typename TTraits>
  void updateBoundingBox(TTraits);

  template <typename RatK, typename AlgK, typename Nt, typename BoundingTratits>
  void updateBoundingBox(
    CGAL::Arr_Bezier_curve_traits_2<RatK, AlgK, Nt, BoundingTratits>);

  template <typename Kernel_>
  void updateBoundingBox(CGAL::Arr_linear_traits_2<Kernel_>);

  void paintFaces(QPainter* painter);

  void paintFaces(QPainter*, QImage&);

  void paintFace(Face_handle f, QPainter* painter);

  void visit_ccb_faces(Face_handle& fh, QPainter* painter);

  /*! antenna - return true if the halfedge and its
   *  twin point to the same face.
   */
  bool antenna(Halfedge_handle h);

  template <typename ArrTraits>
  void paintFace(Face_handle, QPainter*, ArrTraits);

  template <typename Kernel_>
  void paintFace(
    Face_handle f, QPainter* painter, CGAL::Arr_segment_traits_2<Kernel_>);

  template <typename Kernel_>
  void paintFace(
    Face_handle f, QPainter* painter, CGAL::Arr_polyline_traits_2<Kernel_>);

  template <typename RatKernel, typename AlgKernel, typename NtTraits>
  void paintFace(
    Face_handle f, QPainter* painter,
    CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>);

  template <
    typename RatKernel, typename AlgKernel, typename NtTraits,
    typename BoundingTraits>
  void paintFace(
    Face_handle f, QPainter* painter,
    CGAL::Arr_Bezier_curve_traits_2<
      RatKernel, AlgKernel, NtTraits, BoundingTraits>);

  template <typename Coefficient_>
  void paintFace(
    Face_handle f, QPainter* painter,
    CGAL::Arr_algebraic_segment_traits_2<Coefficient_> /* traits */);

  template <typename Kernel_>
  void paintFace(
    Face_handle f, QPainter* painter,
    CGAL::Arr_linear_traits_2<Kernel_> /* traits */);

protected:
  Arrangement* arr;

  // related to painting algebraic faces
  QImage tempImage;
  std::vector<std::pair<uint16_t, uint16_t>> fill_stack;
  static constexpr int margin = 2;
}; // class ArrangementGraphicsItem

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
