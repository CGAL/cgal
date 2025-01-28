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

#ifndef ENVELOPE_CALLBACK_H
#define ENVELOPE_CALLBACK_H

#include "Callback.h"
#include <CGAL/Object.h>

namespace demo_types
{
enum class TraitsType : int;
}

/**
   Updates and draws the lower and upper envelopes of an observed arrangement.
*/
class EnvelopeCallbackBase : public CGAL::Qt::Callback
{
public:
  static EnvelopeCallbackBase*
  create(demo_types::TraitsType, CGAL::Object arr_obj, QObject* parent);

  virtual void setEnvelopeEdgeColor( const QColor& color ) = 0;
  virtual const QColor& getEnvelopeEdgeColor( ) const = 0;
  virtual void setEnvelopeEdgeWidth( int width ) = 0;
  virtual int getEnvelopeEdgeWidth( ) const = 0;
  virtual void setEnvelopeVertexColor( const QColor& color ) = 0;
  virtual const QColor& getEnvelopeVertexColor( ) const = 0;
  virtual void setEnvelopeVertexRadius( int radius ) = 0;
  virtual int getEnvelopeVertexRadius( ) const = 0;
  virtual void showLowerEnvelope( bool b ) = 0;
  virtual void showUpperEnvelope( bool b ) = 0;
  virtual bool isUpperEnvelopeShown() = 0;
  virtual bool isLowerEnvelopeShown() = 0;

protected:
  EnvelopeCallbackBase( QObject* parent );
}; // class EnvelopeCallbackBase

#endif // ENVELOPE_CALLBACK_H
