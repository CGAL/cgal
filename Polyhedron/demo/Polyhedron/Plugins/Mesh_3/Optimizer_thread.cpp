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

#include "config_mesh_3.h"

#include <QElapsedTimer>
#include <QTimer>
#include "Optimizer_thread.h"
#include "Scene_c3t3_item.h"


Optimizer_thread::
Optimizer_thread(Optimization_function_interface* f, Scene_c3t3_item* item)
  : f_(f)
  , item_(item)
  , rc_()
  , time_(0)
  , timer_(new QTimer(this))
  , timer_period_(1)
{
  connect(timer_, SIGNAL(timeout()),
          this,   SLOT(emit_status()));
  timer_->start(static_cast<int>(timer_period_*1000));
}


Optimizer_thread::~Optimizer_thread()
{
  delete f_;
}


void
Optimizer_thread::
run()
{
  QElapsedTimer timer;
  timer.start();
  //SEGFAULT
  rc_ = f_->launch();
  time_ = double(timer.elapsed()) / 1000;
  Q_EMIT done(this);
}


void
Optimizer_thread::
stop()
{
  f_->stop();
}

void
Optimizer_thread::
emit_status()
{
  Q_EMIT (status_report(f_->status(timer_period_)));
}
