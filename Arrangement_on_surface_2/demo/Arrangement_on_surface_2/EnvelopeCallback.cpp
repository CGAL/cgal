// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#include "EnvelopeCallback.h"
#include "CurveGraphicsItem.h"
#include "Utils/EnvelopeFunctions.h"
#include "Utils/Utils.h"
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"

#include <CGAL/envelope_2.h>
#include <CGAL/Envelope_diagram_1.h>

#include <vector>

template <typename Arr_>
class EnvelopeCallback : public EnvelopeCallbackBase
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Point_2 Point_2;

  EnvelopeCallback(Arrangement* arr_, QObject* parent);
  void showLowerEnvelope(bool show) override;
  void showUpperEnvelope(bool show) override;
  bool isUpperEnvelopeShown() override;
  bool isLowerEnvelopeShown() override;
  void reset() override;

  /**
     Slot: Update and redraw the envelopes.
  */
  void slotModelChanged() override;

  void setEnvelopeEdgeColor(const QColor& color) override;
  const QColor& getEnvelopeEdgeColor() const override;
  void setEnvelopeEdgeWidth(int width) override;
  int getEnvelopeEdgeWidth() const override;
  void setEnvelopeVertexColor(const QColor& color) override;
  const QColor& getEnvelopeVertexColor() const override;
  void setEnvelopeVertexRadius(int radius) override;
  int getEnvelopeVertexRadius() const override;

protected:
  void setScene(QGraphicsScene* scene_) override;
  void updateEnvelope(bool lower);

  Arrangement* arr;
  CGAL::Qt::CurveGraphicsItem<Traits>* lowerEnvelope;
  CGAL::Qt::CurveGraphicsItem<Traits>* upperEnvelope;
  bool showLower;
  bool showUpper;
}; // class EnvelopeCallback

EnvelopeCallbackBase::EnvelopeCallbackBase( QObject* parent ) :
  CGAL::Qt::Callback( parent )
{ }

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
    res = new EnvelopeCallback<Arrangement>(arr, parent);
  }

  EnvelopeCallbackBase*& res;
  CGAL::Object& arr_obj;
  QObject* parent;
};
} // anonymous namespace

EnvelopeCallbackBase* EnvelopeCallbackBase::create(
  demo_types::TraitsType tt, CGAL::Object arr_obj, QObject* parent)
{
  EnvelopeCallbackBase* res;
  ExplicitLambda explicit_lambda{res, arr_obj, parent};
  demo_types::visitArrangementType(tt, explicit_lambda);
  return res;
}

template <typename Arr_>
EnvelopeCallback<Arr_>::EnvelopeCallback(Arrangement* arr_, QObject* parent) :
    EnvelopeCallbackBase(parent), arr(arr_),
    lowerEnvelope(new CGAL::Qt::CurveGraphicsItem<Traits>()),
    upperEnvelope(new CGAL::Qt::CurveGraphicsItem<Traits>()), showLower(false),
    showUpper(false)
{
  this->lowerEnvelope->hide();
  this->upperEnvelope->hide();
}

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setEnvelopeEdgeColor( const QColor& color )
{
  this->lowerEnvelope->setEdgeColor( color );
  this->upperEnvelope->setEdgeColor( color );
}

template < typename Arr_ >
const QColor& EnvelopeCallback<Arr_>::getEnvelopeEdgeColor( ) const
{
  return this->lowerEnvelope->edgeColor( );
}

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setEnvelopeEdgeWidth( int width )
{
  this->lowerEnvelope->setEdgeWidth( width );
  this->upperEnvelope->setEdgeWidth( width );
}

template < typename Arr_ >
int EnvelopeCallback<Arr_>::getEnvelopeEdgeWidth( ) const
{
  return this->lowerEnvelope->edgeWidth( );
}

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setEnvelopeVertexColor( const QColor& color )
{
  this->lowerEnvelope->setVertexColor( color );
  this->upperEnvelope->setVertexColor( color );
}

template < typename Arr_ >
const QColor& EnvelopeCallback<Arr_>::getEnvelopeVertexColor( ) const
{
  return this->lowerEnvelope->vertexColor( );
}

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setEnvelopeVertexRadius( int radius )
{
  this->lowerEnvelope->setVertexRadius( radius );
  this->upperEnvelope->setVertexRadius( radius );
}

template < typename Arr_ >
int EnvelopeCallback<Arr_>::getEnvelopeVertexRadius( ) const
{
  return this->lowerEnvelope->vertexRadius( );
}

template < typename Arr_ >
void EnvelopeCallback<Arr_>::setScene( QGraphicsScene* scene_ )
{
  EnvelopeCallbackBase::setScene(scene_);
  lowerEnvelope->setScene(scene_);
  upperEnvelope->setScene(scene_);
  this->scene->addItem( this->lowerEnvelope );
  this->scene->addItem( this->upperEnvelope );
}

template < typename Arr_ >
void EnvelopeCallback< Arr_>::slotModelChanged( )
{
  if (showLower) this->updateEnvelope(true);
  if (showUpper) this->updateEnvelope(false);
}

template < typename Arr_ >
void EnvelopeCallback< Arr_>::updateEnvelope( bool lower )
{
  typedef CGAL::Envelope_diagram_1<Traits> Diagram_1;
  typedef typename Diagram_1::Edge_const_handle Edge_const_handle;

  Diagram_1 diagram;
  CGAL::Qt::CurveGraphicsItem<Traits>* envelopeToUpdate;

  EnvelopeFunctions<Arrangement> envelope_functions;
  try
  {
    if (lower)
    {
      envelopeToUpdate = this->lowerEnvelope;
      envelope_functions.lowerEnvelope(this->arr, diagram);
    }
    else
    {
      envelopeToUpdate = this->upperEnvelope;
      envelope_functions.upperEnvelope(this->arr, diagram);
    }
  }
  catch (const std::exception& ex)
  {
    std::cerr << ex.what() << '\n';
    std::cerr << __FILE__ << ':' << __LINE__ << '\n';
    if (lower)
      this->showLowerEnvelope(false);
    else
      this->showUpperEnvelope(false);
    return;
  }

  envelopeToUpdate->clear( );

  auto next_edge = [](const auto& e) -> Edge_const_handle {
    auto&& v = e->right();
    if (v) return v->right();
    else   return nullptr;
  };

  auto construct_x_monotone_subcurve_2 =
    Construct_x_monotone_subcurve_2<Traits>{arr->traits()};
  for (Edge_const_handle e = diagram.leftmost(); e; e = next_edge(e))
  {
    if (!e->is_empty())
    {
      boost::optional<Point_2> leftPoint, rightPoint;
      if (e->left())
        leftPoint = e->left()->point();

      if (e->right())
        rightPoint = e->right()->point();

      X_monotone_curve_2 curve =
        construct_x_monotone_subcurve_2(e->curve(), leftPoint, rightPoint);

      envelopeToUpdate->insert(curve);
      // TODO: visually show leftPoint and rightPoint
      // note: leftPoint and rightPoint are not the actual curve points
    }
  }

  envelopeToUpdate->modelChanged();
}

template < typename Arr_ >
void EnvelopeCallback< Arr_>::showLowerEnvelope( bool show )
{
  this->showLower = show;
  if (this->showLower)
  {
    this->updateEnvelope(true);
    this->lowerEnvelope->show();
  }
  else
  {
    this->lowerEnvelope->hide();
  }
}

template < typename Arr_ >
void EnvelopeCallback< Arr_>::showUpperEnvelope( bool show )
{
  this->showUpper = show;
  if (this->showUpper)
  {
    this->updateEnvelope(false);
    this->upperEnvelope->show();
  }
  else
  {
    this->upperEnvelope->hide();
  }
}

template <typename Arr_>
void EnvelopeCallback<Arr_>::reset()
{
  this->showLowerEnvelope(false);
  this->showUpperEnvelope(false);
}

template <typename Arr_>
bool EnvelopeCallback<Arr_>::isUpperEnvelopeShown()
{
  return this->showUpper;
}

template <typename Arr_>
bool EnvelopeCallback<Arr_>::isLowerEnvelopeShown()
{
  return this->showLower;
}
