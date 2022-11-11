// Copyright (c) 2008, 2012, 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_TYPED_H
#define CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_TYPED_H

#include <CGAL/Object.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <QImage>
#include <QPen>

#include "ArrangementGraphicsItem.h"
#include "FloodFill.h"
#include "ForwardDeclarations.h"
#include "PointsGraphicsItem.h"

namespace CGAL
{
namespace Qt
{

template <typename Arr_>
class ArrangementGraphicsItem : public ArrangementGraphicsItemBase
{
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Arrangement::Halfedge Halfedge;
  typedef typename Arrangement::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement::Face_handle Face_handle;
  typedef typename Arrangement::Hole_iterator Hole_iterator;
  typedef typename Arrangement::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Point_2 Point_2;

public:
  /*! Constructor */
  ArrangementGraphicsItem(Arrangement* t_);

public:
  void modelChanged() override;
  QRectF boundingRect() const override;
  QRectF getInterestingViewport() const override;
  void paint(
    QPainter* painter, const QStyleOptionGraphicsItem* option,
    QWidget* widget) override;

protected:
  void updatePointsItem();

  template <typename TTraits>
  void paint(QPainter*, const TTraits&);

  template <typename Coefficient_>
  void
  paint(QPainter*, const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>&);

  template <typename AlgebraicKernel_d_1>
  void paint(
    QPainter*,
    const CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>&);

  template <typename TTraits>
  void paintWithFloodFill(QPainter*, const TTraits&);

  void updateBoundingBox();

  template <typename TTraits>
  void paintEdges(QPainter*, const TTraits&);

  template <typename Coefficient_>
  void paintEdges(
    QPainter*, const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>&);

  void paintFaces(QPainter* painter);

  void paintFacesFloodFill(QPainter*, QImage&);

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
    Face_handle f, QPainter* painter,
    const CGAL::Arr_segment_traits_2<Kernel_>&);

  template <typename Kernel_>
  void paintFace(
    Face_handle f, QPainter* painter,
    const CGAL::Arr_polyline_traits_2<Kernel_>&);

  template <typename RatKernel, typename AlgKernel, typename NtTraits>
  void paintFace(
    Face_handle f, QPainter* painter,
    const CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>&);

  template <
    typename RatKernel, typename AlgKernel, typename NtTraits,
    typename BoundingTraits>
  void paintFace(
    Face_handle f, QPainter* painter,
    const CGAL::Arr_Bezier_curve_traits_2<
      RatKernel, AlgKernel, NtTraits, BoundingTraits>&);

  template <typename Coefficient_>
  void paintFace(
    Face_handle f, QPainter* painter,
    const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>&);

  template <typename Kernel_>
  void paintFace(
    Face_handle f, QPainter* painter,
    const CGAL::Arr_linear_traits_2<Kernel_>&);

protected:
  Arrangement* arr;
  PointsGraphicsItem* pointsGraphicsItem;
  CGAL::Bbox_2 bb;
  FloodFill flood_fill;
  QImage tempImage;
}; // class ArrangementGraphicsItem

} // namespace Qt
} // namespace CGAL

#endif
