// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://gdamiand@scm.gforge.inria.fr/svn/cgal/branches/features/Linear_cell_complex-gdamiand/Linear_cell_complex/demo/Linear_cell_complex/MainWindow.h $
// $Id: MainWindow.h 65446 2011-09-20 16:55:42Z gdamiand $
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include "typedefs.h"
#include "ui_MainWindow.h"
#include "ui_CreateMesh.h"

#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Random.h>

#include <QDialog>
#include <QSlider>
#include <QLabel>
#include <QFileDialog>
class QWidget;

class DialogMesh : public QDialog, private Ui::createMesh
{
  Q_OBJECT

  public:
  DialogMesh(QWidget* parent)
  { 
    setupUi (this); 
  }

  int getX() { return xvalue->value(); }
  int getY() { return yvalue->value(); }
  int getZ() { return zvalue->value(); }
};


class MainWindow : public CGAL::Qt::DemosMainWindow, private Ui::MainWindow
{
  Q_OBJECT

  public:
  MainWindow(QWidget* parent = 0);

  void connectActions();

  Scene scene;
  Timer timer;

public slots:
  void import_off();
  void add_off();
  void load_off(const QString& fileName, bool clear=true);

  void import_3DTDS();
  void load_3DTDS(const QString& fileName, bool clear=true);
  
  void clear();

  void create_cube();
  void create_3cubes();
  void create_2volumes();
  void create_mesh();

  void subdivide();
  void dual_3();
  void close_volume();
  void remove_current_volume();
  void sew3_same_facets();
  void unsew3_all();
  void triangulate_all_facets();

  void onSceneChanged();   

 signals:
  void sceneChanged();
  
 protected:
  void initVolumeRandomColor(Dart_handle adart);
  void initAllVolumesRandomColor();
	Dart_handle make_iso_cuboid(const Point_3 basepoint, Map::FT lg);
	
 private:
  unsigned int nbcube;
  QLabel* statusMessage;
  Dart_handle tdsdart;
  DialogMesh dialogmesh;
  CGAL::Random random; 
};




#endif
