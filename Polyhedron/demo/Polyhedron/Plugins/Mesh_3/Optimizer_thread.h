// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
#include <QString>
#include <QStringList>
#include <QTimer>

#include <CGAL/Mesh_optimization_return_code.h>

class Scene_c3t3_item;

class Optimization_function_interface
{
public:
  virtual ~Optimization_function_interface() {}
  
  // Launch
  virtual CGAL::Mesh_optimization_return_code launch() = 0;
  
  // Stop
  virtual void stop() = 0;
  
  // Logs
  virtual QString name() const = 0;
  virtual QStringList parameters_log() const = 0;
  virtual QString status(double time_period) const = 0;
};


class Optimizer_thread : public QThread
{
  Q_OBJECT
public:
  Optimizer_thread(Optimization_function_interface* f, Scene_c3t3_item* item);
  virtual ~Optimizer_thread();
  
  // Scene item
  Scene_c3t3_item* item() const { return item_; }
  
  // Infos about optimization
  CGAL::Mesh_optimization_return_code return_code() const { return rc_; }
  double time() const { return time_; }
  
  // Logs
  QString optimizer_name() const { return f_->name(); }
  QStringList parameters_log() const { return f_->parameters_log(); }
  
public Q_SLOTS:
  // Stop
  void stop();
  
private Q_SLOTS:
  // emit signal status report
  void emit_status();
  
Q_SIGNALS:
  // Emitted at the end of the process
  void done(Optimizer_thread*);
  // Informs about status of the process
  void status_report(QString);
  
protected:
  // Overload of QThread function
  virtual void run();
  
private:
  Optimization_function_interface* f_;
  Scene_c3t3_item* item_;
  CGAL::Mesh_optimization_return_code rc_;
  double time_; // in seconds
  QTimer* timer_;
  double timer_period_;
};

#endif // DEMO_MESH_3_OPTIMIZER_THREAD_H
