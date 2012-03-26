// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Kumar Snehasish <kumar.snehasish@gmail.com>
//
#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include "typedefs.h"
#include "ui_MainWindow.h"
#include "ui_CreateMesh.h"
#include "ui_CreateMenger.h"

#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Random.h>

#include <QDialog>
#include <QSlider>
#include <QLabel>
#include <QFileDialog>

#include <QDockWidget>
#include <QTableWidget>
#include <QCheckBox>

class QWidget;

class DialogMesh : public QDialog, public Ui::createMesh
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

class DialogMenger : public QDialog, public Ui::createMenger
{
    Q_OBJECT

    public:

    DialogMenger(QWidget *parent)
    {
        setupUi(this);
    }
};

class MainWindow : public CGAL::Qt::DemosMainWindow, private Ui::MainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget* parent = 0);

public slots:
  // File menu
  void on_actionImportOFF_triggered();
  void on_actionAddOFF_triggered();
  void on_actionImport3DTDS_triggered();
  void on_actionCompute_Voronoi_3D_triggered();
  void on_actionClear_triggered(bool msg=true);
  
  // Creations menu
  Dart_handle on_actionCreate_cube_triggered();
  void on_actionCreate3Cubes_triggered();
  void on_actionCreate2Volumes_triggered();
  void on_actionCreate_mesh_triggered();
  void on_actionCreate_Menger_Sponge_triggered();

  // Operations menu
  void on_actionSubdivide_triggered();
  void on_actionDual_3_triggered();
  void on_actionClose_volume_triggered();
  void on_actionTriangulate_all_facets_triggered();
  void on_actionSew3_same_facets_triggered();
  void on_actionUnsew3_all_triggered();
  void on_actionMerge_all_volumes_triggered();
  void on_actionRemove_selected_volume_triggered();
  void on_actionRemove_filled_volumes_triggered();

  // View menu
  void on_actionExtend_filled_volumes_triggered();
  void on_actionExtend_hidden_volumes_triggered();

  // Other slots
  void load_off(const QString& fileName, bool clear=true);
  void load_3DTDS(const QString& fileName, bool clear=true);
  void extendVolumesSatisfying(char amask, char negatemask);

  void onSceneChanged();

  void connectVolumeListHandlers();
  void onCellChanged(int, int);
  void onItemSelectionChanged();
  void onHeaderClicked(int);

  void onCreateMeshOk();
  
  void onMengerInc();
  void onMengerDec();
  void onMengerChange(int);
  void onMengerOk();
  void onMengerCancel();
    
signals:
  void sceneChanged();
  
protected:
  void onNewVolume(Dart_handle adart);
  void onDeleteVolume(Dart_handle adart);
  void initAllNewVolumes();
  
  Dart_handle make_iso_cuboid(const Point_3 basepoint, LCC::FT lg);

  void connectActions();
  void updateOperationEntry(bool show);

  void check_volume_list();
  bool is_volume_in_list(Dart_handle it);
  void update_volume_list();
  void update_volume_list_add(Dart_handle);
  void update_volume_list_remove(unsigned int);
  void update_volume_list_remove(Dart_handle);

  void split_edge_in_three     (Dart_handle dh);
  void split_face_in_three     (Dart_handle dh);
  void split_face_in_nine      (Dart_handle dh);
  void split_vol_in_three      (Dart_handle dh, bool removecenter);
  void split_vol_in_nine       (Dart_handle dh, bool removecenter);
  void split_vol_in_twentyseven(Dart_handle dh);
  void processFullSlice(Dart_handle init,
                        std::vector<Dart_handle>& faces,
                        int markVols);
  void processInterSlice(Dart_handle init,
                         std::vector<Dart_handle>& faces,
                         int markVols);
  
  Scene scene;
  Timer timer;

  unsigned int nbcube;
  QLabel*      statusMessage;
  DialogMesh   dialogmesh;
  DialogMenger dialogmenger;
  CGAL::Random random;

  unsigned int mengerLevel;
  unsigned int mengerFirstVol;
  Dart_handle mengerDart;
  
  QDockWidget* volumeListDock;
  QTableWidget* volumeList;

  unsigned int volumeUid;
  std::vector< std::pair<unsigned int,Dart_handle> > volumeDartIndex;
  std::vector< char > volumeProperties;
};

#endif
