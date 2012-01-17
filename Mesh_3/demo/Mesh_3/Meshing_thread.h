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

#ifndef CGAL_DEMO_MESH_3_MESHING_THREAD_H
#define CGAL_DEMO_MESH_3_MESHING_THREAD_H

#include <QThread>
#include <QObject>
#include <QStringList>
#include <QString>
#include <QTimer>

class Scene_c3t3_item;

class Mesh_function_interface
{
public:
  virtual ~Mesh_function_interface() {}
  
  // Launch
  virtual void launch() = 0;
  
  // Stop
  virtual void stop() = 0;
  
  // Logs
  virtual QStringList parameters_log() const = 0;
  virtual QString status(double time_period) const = 0;
};


class Meshing_thread : public QThread
{
  Q_OBJECT
public:
  // Constructor / Destructor
  Meshing_thread(Mesh_function_interface* f, Scene_c3t3_item* item);
  virtual ~Meshing_thread();
  
  // Scene item
  Scene_c3t3_item* item() const { return item_; }
  
  // Infos about meshing
  double time() const { return time_; }
  
  // Logs
  QStringList parameters_log() const { return f_->parameters_log(); }
  
public slots:
  // Stop
  void stop();
  
private slots:
  // emit signal status report
  void emit_status();
  
signals:
  // Emitted at the end of the process
  void done(Meshing_thread*);
  // Informs about status of meshing
  void status_report(QString);
  
protected:
  // Overload of QThread function
  virtual void run();
  
private:
  Mesh_function_interface* f_;
  Scene_c3t3_item* item_;
  double time_; // in seconds
  QTimer* timer_;
  double timer_period_;
};

#endif // CGAL_DEMO_MESH_3_MESHING_THREAD_H
