// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#include "config.h"

#include <QElapsedTimer>
#include <QApplication>

#include "Meshing_thread.h"
#include <CGAL/Three/Three.h>

class Scene_c3t3_item;

Meshing_thread::
Meshing_thread(Mesh_function_interface* f, Scene_c3t3_item* item)
  : f_(f)
  , item_(item)
  , time_(0)
  , timer_(new QTimer(this))
  , timer_period_(1)
{
  connect(timer_, SIGNAL(timeout()),
          this,   SLOT(emit_status()));

  timer_->start(static_cast<int>(timer_period_*1000));
}


Meshing_thread::
~Meshing_thread()
{
  delete f_;
  delete timer_;
  QApplication::restoreOverrideCursor();
}


void
Meshing_thread::
run()
{
  QElapsedTimer timer;
  timer.start();
  CGAL::Three::Three::CursorScopeGuard guard(Qt::BusyCursor);
  f_->launch();
  time_ = double(timer.elapsed()) / 1000;
  Q_EMIT done(this);
}


void
Meshing_thread::
stop()
{
  f_->stop();
  QApplication::setOverrideCursor(Qt::WaitCursor); //restored in mesh_3_plugin lambda
}


void
Meshing_thread::
emit_status()
{
  Q_EMIT (status_report(f_->status(timer_period_)));
}
