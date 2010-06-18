// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#include <QTime>
#include "Optimizer_thread.h"
#include "Scene_c3t3_item.h"

Optimizer_thread::~Optimizer_thread()
{
  delete f_;
}


void
Optimizer_thread::
run()
{
  QTime timer;
  timer.start();
  
  rc_ = f_->launch();
  time_ = double(timer.elapsed()) / 1000;
  
  emit done(this);
}


void
Optimizer_thread::
stop()
{
  f_->stop();
}


#include "Optimizer_thread.moc"