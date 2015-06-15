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
// Contributor(s): Kumar Snehasish <kumar.snehasish@gmail.com>
//
#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include "typedefs.h"
#include "ui_MainWindow.h"
#include "ui_CreateMesh.h"
#include "ui_CreateMenger.h"
#include "ui_CreateSierpinskiCarpet.h"
#include "ui_CreateSierpinskiTriangle.h"

#include <CGAL/Qt/DemosMainWindow.h>

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
  DialogMesh(QWidget* /*parent*/)
  { setupUi (this); }

  int getX() { return xvalue->value(); }
  int getY() { return yvalue->value(); }
  int getZ() { return zvalue->value(); }
};

class DialogMenger : public QDialog, public Ui::createMenger
{
  Q_OBJECT

public:
  DialogMenger(QWidget* /*parent*/)
  { setupUi(this); }
};

class DialogSierpinskiCarpet : public QDialog,
    public Ui::createSierpinskiCarpet
{
  Q_OBJECT

public:
  DialogSierpinskiCarpet(QWidget* /*parent*/)
  { setupUi(this); }
};

class DialogSierpinskiTriangle : public QDialog,
    public Ui::createSierpinskiTriangle
{
  Q_OBJECT

public:
  DialogSierpinskiTriangle(QWidget* /*parent*/)
  { setupUi(this); }
};


template < class First, class Second, class Third > struct Triplet
{
  First first;
  Second second;
  Third third;

  Triplet(First first, Second second, Third third)
  {
    this->first = first;
    this->second = second;
    this->third = third;
  }

  Triplet()
  {}
};

class MainWindow : public CGAL::Qt::DemosMainWindow, private Ui::MainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget* parent = 0);

public Q_SLOTS:
  // File menu
  void on_actionSave_triggered();
  void on_actionLoad_triggered();  
  void on_actionImportOFF_triggered();
  void on_actionAddOFF_triggered();
  void on_actionImport3DTDS_triggered();
  void on_actionImportMoka_triggered();
  void on_actionCompute_Voronoi_3D_triggered();
  void on_actionClear_triggered();
  
  // Creations menu
  Dart_handle on_actionCreate_cube_triggered();
  void on_actionCreate3Cubes_triggered();
  void on_actionCreate2Volumes_triggered();
  void on_actionCreate_mesh_triggered();
  void on_actionCreate_Menger_Sponge_triggered();
  void on_actionCreate_Sierpinski_Carpet_triggered();
  void on_actionCreate_Sierpinski_Triangle_triggered();

  // Operations menu
  void on_actionSubdivide_triggered();
  void on_actionSubdivide_pqq_triggered();
  void on_actionDual_3_triggered();
  void on_actionClose_volume_triggered();
  void on_actionSew3_same_facets_triggered();
  void on_actionUnsew3_all_triggered();
  void on_actionMerge_coplanar_faces_triggered();
  void on_actionMerge_all_volumes_triggered();
  void on_actionRemove_filled_volumes_triggered();
  void on_actionInsert_center_vertices_triggered();
  void on_actionTriangulate_all_facets_triggered();

  // View menu
  void on_actionExtend_filled_volumes_triggered();
  void on_actionExtend_hidden_volumes_triggered();

  // Other slots
  void load_depend_on_extension(const QString& fileName, bool clear=true);
  void load(const QString& fileName, bool clear=true);
  void save(const QString& fileName);
  void load_off(const QString& fileName, bool clear=true);
  void load_3DTDS(const QString& fileName, bool clear=true);
  void load_moka(const QString& fileName, bool clear=true);

  void onSceneChanged();

  void connectVolumeListHandlers();
  void onCellChanged(int, int);
  void onHeaderClicked(int);

  void onCreateMeshOk();
  
  void onMengerInc();
  void onMengerDec();
  void onMengerChange(int);
  void onMengerOk();
  void onMengerCancel();
  void onMengerUpdateAttributes(bool);

  void onSierpinskiCarpetChangeLevel(int);
/*  void onSierpinskiCarpetNeverUpdateAttributes(bool);
  void onSierpinskiCarpetDuringConstructionUpdateAttributes(bool);
  void onSierpinskiCarpetAfterConstructionUpdateAttributes(bool);
  void onSierpinskiCarpetUpdateAttributesMethodStdMap(bool);
  void onSierpinskiCarpetUpdateAttributesMethodTraversal(bool);
  void onSierpinskiCarpetComputeGeometry(bool);*/
  void onSierpinskiCarpetUpdateAttributes(bool);
  void onSierpinskiCarpetOk();
  void onSierpinskiCarpetCancel();
  void onSierpinskiCarpetInc();
  void onSierpinskiCarpetDec();

  void onSierpinskiTriangleChangeLevel(int);
  void onSierpinskiTriangleUpdateAttributes(bool);
  void onSierpinskiTriangleOk();
  void onSierpinskiTriangleCancel();
  void onSierpinskiTriangleInc();
  void onSierpinskiTriangleDec();

Q_SIGNALS:
  void sceneChanged();
  
protected:
  void clear_all();
  void on_new_volume(Dart_handle adart);
  void on_delete_volume(Dart_handle adart);
  void init_all_new_volumes();
  void mark_all_filled_and_visible_volumes(int amark);

  Dart_handle make_iso_cuboid(const Point_3 basepoint, LCC::FT lg);

  void connect_actions();
  void update_operations_entries(bool show);

  bool is_volume_in_list(LCC::Attribute_handle<3>::type ah);
  void recreate_whole_volume_list();
  void update_volume_list_all_ckeckstates();
  void update_volume_list_add(LCC::Attribute_handle<3>::type ah);
  void update_volume_list_remove(int);
  void update_volume_list_remove(LCC::Attribute_handle<3>::type ah);

  void split_edge_in_three     (Dart_handle dh);
  void split_face_in_three     (Dart_handle dh);
  void split_face_in_nine      (Dart_handle dh);
  void split_vol_in_three      (Dart_handle dh, bool removecenter);
  void split_vol_in_nine       (Dart_handle dh, bool removecenter);
  void split_vol_in_twentyseven(Dart_handle dh);
  void process_full_slice(Dart_handle init,
                          std::vector<Dart_handle>& faces,
                          int markVols);
  void process_inter_slice(Dart_handle init,
                           std::vector<Dart_handle>& faces,
                           int markVols);

  void sierpinski_carpet_copy_attributes_and_embed_vertex(Dart_handle dh,
                                                          LCC::Point& p);
  void sierpinski_carpet_update_geometry();
  void sierpinski_carpet_compute_geometry();
  void sierpinski_carpet_compute_4x4_geometry_matrix(LCC::Point p[4][4],
  LCC::Point& p00, LCC::Point& p03, LCC::Point& p33, LCC::Point& p30);
  void sierpinski_carpet_split_edge_in_three(Dart_handle dh);
  void sierpinski_carpet_split_face_in_three(Dart_handle dh,
                                             bool removecenter);
  void sierpinski_carpet_split_face_in_nine(Dart_handle dh);

  void sierpinski_triangle_split_edge_in_two(Dart_handle dh);
  void sierpinski_triangle_split_face_in_four(Dart_handle dh,
                                              bool removecenter);

  Scene scene;

  unsigned int nbcube;
  QLabel*      statusMessage;
  DialogMesh   dialogmesh;
  DialogMenger dialogmenger;
  DialogSierpinskiCarpet dialogsierpinskicarpet;
  DialogSierpinskiTriangle dialogsierpinskitriangle;

  int mengerLevel;
  bool mengerUpdateAttributes;
  std::vector<Dart_handle> mengerVolumes;

  int sierpinskiCarpetLevel;
  std::size_t nbfacesinit;
  bool sierpinskiCarpetUpdateAttributes;
  bool computeGeometry;
  /*bool neverUpdateAttributes;
  bool duringConstructionUpdateAttributes;
  bool afterConstructionUpdateAttributes;
  bool updateAttributesMethodStdMap;
  bool updateAttributesMethodTraversal;
  bool isComputableGeometry;*/
  std::vector<Dart_handle> sierpinskiCarpetSurfaces;
  // utilisés seulement lorsque pas de mise à jour d'attributs
  std::map<Dart_handle, LCC::Point> dart_map;
  std::vector<Dart_handle> new_darts;

  int sierpinskiTriangleLevel;
  bool sierpinskiTriangleUpdateAttributes;
  std::vector<Dart_handle> sierpinskiTriangleSurfaces;
  std::vector< Triplet<Dart_handle, Dart_handle, Dart_handle> > removedTriangles;

  QDockWidget* volumeListDock;
  QTableWidget* volumeList;
};

#endif
