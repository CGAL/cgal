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

public Q_SLOTS:
  // Stop
  void stop();

private Q_SLOTS:
  // emit signal status report
  void emit_status();

Q_SIGNALS:
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
