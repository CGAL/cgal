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

#ifndef DEMO_MESH_3_OPTIMIZER_THREAD_H
#define DEMO_MESH_3_OPTIMIZER_THREAD_H

#include <QThread>
#include <QObject>

#include <CGAL/Mesh_optimization_return_code.h>

class Scene_c3t3_item;

class Optimization_function_interface
{
public:
  virtual ~Optimization_function_interface() {}
  virtual CGAL::Mesh_optimization_return_code launch() const = 0;
  virtual Scene_c3t3_item* item() const = 0;
};

class Optimizer_thread : public QThread
{
  Q_OBJECT
public:
  Optimizer_thread(Optimization_function_interface* f)
    : f_(f), rc_() {}
  
  virtual ~Optimizer_thread();
  
  Scene_c3t3_item* item() const { return f_->item(); }
  CGAL::Mesh_optimization_return_code return_code() const { return rc_; }
  
protected:
  void run();
  
private:
  Optimization_function_interface* f_;
  CGAL::Mesh_optimization_return_code rc_;
};

#endif // DEMO_MESH_3_OPTIMIZER_THREAD_H
