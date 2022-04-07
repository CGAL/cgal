// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>,
//            Saurabh Singh <ssingh@cs.iitr.ac.in>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_H

#include "Callback.h"
#include "ForwardDeclarations.h"

class QEvent;
class PointSnapperBase;

namespace demo_types
{
enum class TraitsType : int;
}

namespace CGAL
{
namespace Qt
{

class CurveInputMethod;
enum class CurveType : int;

class GraphicsViewCurveInputBase : public Callback
{
  Q_OBJECT

public:
  static GraphicsViewCurveInputBase* create(
    demo_types::TraitsType, CGAL::Object arr_obj, QObject* parent,
    QGraphicsScene* scene);

  void setColor(QColor c);
  void reset() override;
  bool eventFilter(QObject* obj, QEvent* event) override;
  virtual void setCurveType(CurveType type) = 0;
  virtual void setPointSnapper(PointSnapperBase*) = 0;
  virtual void generate(CGAL::Object) = 0;

protected:
  GraphicsViewCurveInputBase(QObject* parent, QGraphicsScene* scene);

  void setInputMethod(CurveInputMethod*);

private:
  // active input method
  CurveInputMethod* inputMethod;
}; // class GraphicsViewCurveInputBase

} // namespace Qt
} // namespace CGAL

// TODO: Find a better place for these functions
// (away from ArrangementDemoWindow.cpp, for better compilation speeds)
#ifdef CGAL_USE_CORE
CGAL::Object algebraicCurveFromExpression(
  const CGAL::Object& arr, const std::string&, bool& is_first_curve);
CGAL::Object rationalCurveFromExpression(
  const CGAL::Object& arr, const std::string&, const std::string&,
  bool& is_first_curve);
#endif

#endif // CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
