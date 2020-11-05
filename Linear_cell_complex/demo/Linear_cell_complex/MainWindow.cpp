// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Kumar Snehasish <kumar.snehasish@gmail.com>
//                 Sylvain Brandel <sylvain.brandel@liris.cnrs.fr>
//
#include "MainWindow.h"
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polyhedron_3_to_lcc.h>
#include <CGAL/Triangulation_3_to_lcc.h>
#include <QSettings>
#include <QHeaderView>
#include <CGAL/Timer.h>
#include <CGAL/ipower.h>
#include "import_moka.h"


// Function defined in Linear_cell_complex_3_subivision.cpp
void subdivide_lcc_3 (LCC & m);

// Function defined in Linear_cell_complex_pqq_subivision.cpp
void subdivide_lcc_pqq (LCC & m);

#define DELAY_STATUSMSG 1500

MainWindow::MainWindow (QWidget * parent) : CGAL::Qt::DemosMainWindow (parent),
  nbcube      (0),
  dialogmesh  (this),
  dialogmenger(this),
  dialogsierpinskicarpet(this),
  dialogsierpinskitriangle(this)
{
  setupUi (this);
  scene.lcc = new LCC;

  volumeListDock = new QDockWidget(QString(tr("Volume List")),this);
  volumeListDock->setAllowedAreas(Qt::RightDockWidgetArea |
                                  Qt::LeftDockWidgetArea);
  volumeList = new QTableWidget(0,4,volumeListDock);
  volumeList->verticalHeader()->hide();
  volumeList->setColumnHidden(3,true);
  QObject::connect(this->volumeList, SIGNAL(cellChanged(int,int)),
                   this, SLOT(onCellChanged(int,int)));

  QStringList labels(QString(tr("Color")));
  labels.append(QString(tr("Filled")));
  labels.append(QString(tr("Hidden")));
  volumeList->setHorizontalHeaderLabels(labels);
  //volumeList->resizeColumnsToContents();
  volumeList->setFixedWidth(220);
/*  volumeList->setColumnWidth(0,85);
  volumeList->setColumnWidth(1,35);
  volumeList->setColumnWidth(2,35);*/

  volumeList->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

  volumeList->setSelectionMode(QAbstractItemView::NoSelection);
  //volumeList->setSelectionBehavior(QAbstractItemView::SelectRows);
  volumeListDock->setWidget(volumeList);
  addDockWidget(Qt::RightDockWidgetArea,volumeListDock);
  menuView->addAction(volumeListDock->toggleViewAction());

  QObject::connect(&dialogmesh, SIGNAL(accepted()),
                   this, SLOT(onCreateMeshOk()));
  this->viewer->setScene(&scene, false);

  connect_actions ();
  this->addAboutDemo (":/cgal/help/about_Linear_cell_complex_3.html");
  this->addAboutCGAL ();

  this->addRecentFiles (this->menuFile, this->actionQuit);
  connect (this, SIGNAL (openRecentFile (QString)),
           this, SLOT (load_depend_on_extension(QString)));

  statusMessage = new QLabel
      ("Darts: 0,  Vertices: 0  (Points: 0),  Edges: 0, Facets: 0,"
       " Volume: 0 (Vol color: 0),  Connected components: 0");
  statusBar ()->addWidget (statusMessage);

}

void MainWindow::connect_actions ()
{
  QObject::connect (this->actionQuit, SIGNAL (triggered ()),
                    qApp, SLOT (quit ()));

  QObject::connect (this, SIGNAL (sceneChanged ()),
                    this, SLOT (onSceneChanged ()));

  QObject::connect(this->volumeList->horizontalHeader(),
                   SIGNAL(sectionClicked(int)),
                   this, SLOT(onHeaderClicked(int)));
  QObject::connect(dialogmenger.mengerLevel, SIGNAL(valueChanged(int)),
                   this, SLOT(onMengerChange(int)));
  QObject::connect(dialogmenger.updateAttributes, SIGNAL(clicked(bool)),
                   this, SLOT(onMengerUpdateAttributes(bool)));
  QObject::connect(&dialogmenger, SIGNAL(accepted()),
                   this, SLOT(onMengerOk()));
  QObject::connect(&dialogmenger, SIGNAL(rejected()),
                   this, SLOT(onMengerCancel()));

  QObject::connect(dialogsierpinskicarpet.level, SIGNAL(valueChanged(int)),
                   this, SLOT(onSierpinskiCarpetChangeLevel(int)));
  QObject::connect(dialogsierpinskicarpet.updateAttributes, SIGNAL(clicked(bool)),
                   this, SLOT(onSierpinskiCarpetUpdateAttributes(bool)));
/*  QObject::connect(dialogsierpinskicarpet.never, SIGNAL(clicked(bool)),
                   this, SLOT(onSierpinskiCarpetNeverUpdateAttributes(bool)));
  QObject::connect(dialogsierpinskicarpet.during, SIGNAL(clicked(bool)),
                   this, SLOT(onSierpinskiCarpetDuringConstructionUpdateAttributes(bool)));
  QObject::connect(dialogsierpinskicarpet.after, SIGNAL(clicked(bool)),
                   this, SLOT(onSierpinskiCarpetAfterConstructionUpdateAttributes(bool)));
  QObject::connect(dialogsierpinskicarpet.stdmap, SIGNAL(clicked(bool)),
                   this, SLOT(onSierpinskiCarpetUpdateAttributesMethodStdMap(bool)));
  QObject::connect(dialogsierpinskicarpet.traversal, SIGNAL(clicked(bool)),
                   this, SLOT(onSierpinskiCarpetUpdateAttributesMethodTraversal(bool)));
  QObject::connect(dialogsierpinskicarpet.computeGeometry, SIGNAL(clicked(bool)),
                   this, SLOT(onSierpinskiCarpetComputeGeometry(bool)));*/
  QObject::connect(&dialogsierpinskicarpet, SIGNAL(accepted()),
                   this, SLOT(onSierpinskiCarpetOk()));
  QObject::connect(&dialogsierpinskicarpet, SIGNAL(rejected()),
                   this, SLOT(onSierpinskiCarpetCancel()));

  QObject::connect(dialogsierpinskitriangle.level, SIGNAL(valueChanged(int)),
                   this, SLOT(onSierpinskiTriangleChangeLevel(int)));
  QObject::connect(dialogsierpinskitriangle.updateAttributes, SIGNAL(clicked(bool)),
                   this, SLOT(onSierpinskiTriangleUpdateAttributes(bool)));
  QObject::connect(&dialogsierpinskitriangle, SIGNAL(accepted()),
                   this, SLOT(onSierpinskiTriangleOk()));
  QObject::connect(&dialogsierpinskitriangle, SIGNAL(rejected()),
                   this, SLOT(onSierpinskiTriangleCancel()));
}

void MainWindow::connectVolumeListHandlers()
{
  QObject::connect(this->volumeList, SIGNAL(cellChanged(int,int)),
                   this, SLOT(onCellChanged(int,int)));
}

void MainWindow::update_operations_entries(bool show)
{
  actionImportOFF->setEnabled(show);
  actionImport3DTDS->setEnabled(show);
  actionCompute_Voronoi_3D->setEnabled(show);
  actionClear->setEnabled(show);
  menuCreations->setEnabled(show);
  menuOperations->setEnabled(show);
}

void MainWindow::onSceneChanged ()
{
  QApplication::setOverrideCursor( Qt::WaitCursor );

  LCC::size_type mark = scene.lcc->get_new_mark ();
  scene.lcc->negate_mark (mark);

  std::vector<unsigned int> cells;
  cells.push_back(0);
  cells.push_back(1);
  cells.push_back(2);
  cells.push_back(3);
  cells.push_back(4);

  std::vector<unsigned int> res = scene.lcc->count_cells (cells);

  std::ostringstream os;
  os << "Darts: " << scene.lcc->number_of_darts ()
     << ",  Vertices:" << res[0]
     <<",  (Points:"<<scene.lcc->number_of_attributes<0>()<<")"
    << ",  Edges:" << res[1]
    << ",  Facets:" << res[2]
    << ",  Volumes:" << res[3]
    <<",  (Vol color:"<<scene.lcc->number_of_attributes<3>()<<")"
   << ",  Connected components:" << res[4]
   <<",  Valid:"<<(scene.lcc->is_valid()?"true":"FALSE");

  scene.lcc->negate_mark (mark);
  scene.lcc->free_mark (mark);

  // statusBar()->showMessage (QString ("Update OpenGL lists"), DELAY_STATUSMSG);

  viewer->sceneChanged ();

  statusMessage->setText (os.str().c_str ());
  QApplication::restoreOverrideCursor();
}

void MainWindow::clear_all()
{
  scene.lcc->clear();
  nbcube=0;

  volumeList->clearContents();
  volumeList->setRowCount(0);
}

void MainWindow::on_new_volume(Dart_handle adart)
{
  CGAL_assertion( scene.lcc->attribute<3>(adart)==LCC::null_handle);
  scene.lcc->set_attribute<3>(adart, scene.lcc->create_attribute<3>());
  update_volume_list_add(scene.lcc->attribute<3>(adart));
}

void MainWindow::init_all_new_volumes()
{
  for (LCC::One_dart_per_cell_range<3>::iterator
       it(scene.lcc->one_dart_per_cell<3>().begin());
       it.cont(); ++it)
    if ( scene.lcc->attribute<3>(it)==LCC::null_handle )
    { on_new_volume(it); }
}

void MainWindow::on_actionSave_triggered ()
{
  QString fileName = QFileDialog::getSaveFileName (this,
                                                   tr ("Save"),
                                                   "save.3map",
                                                   tr ("3-map files (*.3map)"));

  if (!fileName.isEmpty ())
  {
     save(fileName);
  }
}

void MainWindow::on_actionLoad_triggered ()
{
  QString fileName = QFileDialog::getOpenFileName (this,
                                                   tr ("Load"),
                                                   "./3map",
                                                   tr ("3-map files (*.3map)"));

  if (!fileName.isEmpty ())
  {
    load(fileName, false);
  }
}

void MainWindow::on_actionImportOFF_triggered ()
{
  QString fileName = QFileDialog::getOpenFileName (this,
                                                   tr ("Import OFF"),
                                                   "./off",
                                                   tr ("off files (*.off)"));

  if (!fileName.isEmpty ())
  {
    load_off (fileName, false);
  }
}

void MainWindow::on_actionImportMoka_triggered()
{
  QString fileName = QFileDialog::getOpenFileName (this,
                                                   tr ("Import Moka"),
                                                   "./moka",
                                                   tr ("Moka files (*.moka)"));

  if (!fileName.isEmpty ())
  {
    load_moka(fileName, false);
  }
}

void MainWindow::on_actionImport3DTDS_triggered ()
{
  QString fileName = QFileDialog::getOpenFileName (this,
                                                   tr ("Import 3DTDS"),
                                                   ".",
                                                   tr ("Data file (*)"));

  if (!fileName.isEmpty ())
  {
    load_3DTDS (fileName, false);
    statusBar ()->showMessage (QString ("Import 3DTDS file") + fileName,
                               DELAY_STATUSMSG);
  }
}

void MainWindow::load_depend_on_extension(const QString & fileName, bool clear)
{
  QString ext = QFileInfo(fileName).suffix();
  if ( ext=="3map")
  {
    load(fileName, clear);
  }
  else if (ext=="off")
  {
    load_off(fileName, clear);
  }
  else if (ext=="moka")
  {
    load_moka(fileName, clear);
  }
  else
  {
    std::cout<<"Extension not considered."<<std::endl;
  }
}

void MainWindow::load(const QString & fileName, bool clear)
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

  if (clear) this->clear_all();

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  bool res = load_combinatorial_map(fileName.toStdString().c_str(), *(scene.lcc));

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to load 3-map "<<qPrintable(fileName)<<": "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  init_all_new_volumes();
  recreate_whole_volume_list();

  this->addToRecentFiles(fileName);
  QApplication::restoreOverrideCursor ();

  if (res)
    statusBar ()->showMessage (QString ("3-map loaded ") + fileName,
                               DELAY_STATUSMSG);
  else
    statusBar ()->showMessage (QString ("Problem: 3-map not loaded ") + fileName,
                               DELAY_STATUSMSG);
  Q_EMIT (sceneChanged ());
}

void MainWindow::save(const QString & fileName)
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  if ( save_combinatorial_map(*(scene.lcc), fileName.toStdString().c_str()) )
    statusBar ()->showMessage (QString ("3-map saved ") + fileName,
                               DELAY_STATUSMSG);
  else
    statusBar ()->showMessage (QString ("Problem: 3-map not saved ") + fileName,
                               DELAY_STATUSMSG);
  QApplication::restoreOverrideCursor ();

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to save 3-map "<<qPrintable(fileName)<<": "
           <<timer.time()<<" seconds."<<std::endl;
#endif
}

void MainWindow::load_off (const QString & fileName, bool clear)
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

  if (clear) this->clear_all();

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  std::ifstream ifs (qPrintable (fileName));

  CGAL::import_from_polyhedron_3_flux < LCC > (*scene.lcc, ifs);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to load off "<<qPrintable(fileName)<<": "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  init_all_new_volumes();

  this->addToRecentFiles (fileName);
  QApplication::restoreOverrideCursor ();

  if (clear)
    statusBar ()->showMessage (QString ("Load off file") + fileName,
                               DELAY_STATUSMSG);
  else
    statusBar ()->showMessage (QString ("Add off file") + fileName,
                               DELAY_STATUSMSG);
  Q_EMIT (sceneChanged ());
}

void MainWindow::load_3DTDS (const QString & fileName, bool clear)
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

  if (clear) this->clear_all();

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  typedef CGAL::Delaunay_triangulation_3 < LCC::Traits > Triangulation;
  Triangulation T;

  std::ifstream ifs (qPrintable (fileName));
  T.insert (std::istream_iterator < Point_3 >(ifs),
            std::istream_iterator < Point_3 >() );

  CGAL::import_from_triangulation_3 < LCC, Triangulation >(*scene.lcc, T);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to import the 3DTDS from "<<qPrintable(fileName)<<": "
           <<timer.time()
           <<" seconds (counting the time to compute denaulay triangulation)."
           <<std::endl;
#endif

  init_all_new_volumes();

  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
}

void MainWindow::load_moka(const QString & fileName, bool clear)
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

  if (clear) this->clear_all();

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  CGAL::import_from_moka < LCC > (*scene.lcc, qPrintable (fileName));

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to load moka "<<qPrintable(fileName)<<": "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  init_all_new_volumes();
  recreate_whole_volume_list();

  this->addToRecentFiles (fileName);
  QApplication::restoreOverrideCursor ();

  if (clear)
    statusBar ()->showMessage (QString ("Load moka file") + fileName,
                               DELAY_STATUSMSG);
  else
    statusBar ()->showMessage (QString ("Add moka file") + fileName,
                               DELAY_STATUSMSG);
  Q_EMIT (sceneChanged ());
}

Dart_handle MainWindow::make_iso_cuboid(const Point_3 basepoint, LCC::FT lg)
{
  Dart_handle dh = scene.lcc->make_hexahedron(basepoint,
                                    LCC::Traits::Construct_translated_point()
                                    (basepoint,LCC::Traits::Vector(lg,0,0)),
                                    LCC::Traits::Construct_translated_point()
                                    (basepoint,LCC::Traits::Vector(lg,lg,0)),
                                    LCC::Traits::Construct_translated_point()
                                    (basepoint,LCC::Traits::Vector(0,lg,0)),
                                    LCC::Traits::Construct_translated_point()
                                    (basepoint,LCC::Traits::Vector(0,lg,lg)),
                                    LCC::Traits::Construct_translated_point()
                                    (basepoint,LCC::Traits::Vector(0,0,lg)),
                                    LCC::Traits::Construct_translated_point()
                                    (basepoint,LCC::Traits::Vector(lg,0,lg)),
                                    LCC::Traits::Construct_translated_point()
                                    (basepoint,LCC::Traits::Vector(lg,lg,lg)));

  scene.lcc->reverse_orientation_connected_component(dh);
  return dh;
}

Dart_handle MainWindow::on_actionCreate_cube_triggered ()
{
  Point_3 basepoint(nbcube%5, (nbcube/5)%5, nbcube/25);

  Dart_handle d = make_iso_cuboid(basepoint, 1);

  on_new_volume(d);

  ++nbcube;

  statusBar ()->showMessage (QString ("Cube created"),DELAY_STATUSMSG);

  Q_EMIT (sceneChanged ());

  return d;
}

void MainWindow::on_actionCreate3Cubes_triggered ()
{
  Dart_handle d1 = make_iso_cuboid (Point_3 (nbcube, nbcube, nbcube),1);
  Dart_handle d2 = make_iso_cuboid (Point_3 (nbcube + 1, nbcube, nbcube),1);
  Dart_handle d3 = make_iso_cuboid (Point_3 (nbcube, nbcube + 1, nbcube), 1);

  on_new_volume(d1);
  on_new_volume(d2);
  on_new_volume(d3);

  scene.lcc->sew<3>(scene.lcc->beta(d1,1,1,2), scene.lcc->beta(d2,2));
  scene.lcc->sew<3>(scene.lcc->beta(d1,2,1,1,2), d3);

  ++nbcube;

  statusBar ()->showMessage (QString ("3 cubes were created"),
                             DELAY_STATUSMSG);

  Q_EMIT (sceneChanged ());
}

void MainWindow::on_actionCreate2Volumes_triggered ()
{
  Dart_handle d1 = make_iso_cuboid(Point_3(nbcube, nbcube, nbcube),1);
  Dart_handle d2 = make_iso_cuboid(Point_3(nbcube + 1, nbcube, nbcube), 1);
  Dart_handle d3 = make_iso_cuboid(Point_3(nbcube, nbcube + 1, nbcube), 1);
  Dart_handle d4 = make_iso_cuboid(Point_3(nbcube + 1, nbcube + 1, nbcube), 1);

  scene.lcc->sew<3>(scene.lcc->beta(d1,1,1,2), scene.lcc->beta(d2,2));
  scene.lcc->sew<3>(scene.lcc->beta(d1,2,1,1,2), d3);

  scene.lcc->sew<3>(scene.lcc->beta(d3,1,1,2), scene.lcc->beta(d4,2));
  scene.lcc->sew<3>(scene.lcc->beta(d2,2,1,1,2), d4);

  scene.lcc->remove_cell<2>(d3);
  scene.lcc->remove_cell<2>(scene.lcc->beta(d2,2));

  on_new_volume(d1);
  on_new_volume(d4);

  ++nbcube;
  statusBar ()->showMessage (QString ("2 volumes were created"),
                             DELAY_STATUSMSG);

  Q_EMIT (sceneChanged());
}

void MainWindow::on_actionCreate_mesh_triggered ()
{
  dialogmesh.show();
}

void MainWindow::onCreateMeshOk()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  for (int x=0; x<dialogmesh.getX(); ++x)
    for (int y=0; y<dialogmesh.getY(); ++y)
      for (int z=0; z<dialogmesh.getZ(); ++z)
      {
        Dart_handle d = make_iso_cuboid
          (Point_3 (x+nbcube, y+nbcube, z+nbcube), 1);
        on_new_volume(d);
      }
  nbcube+=dialogmesh.getX();

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to create mesh ("<<dialogmesh.getX()<<", "
           <<dialogmesh.getY()<<", "<<dialogmesh.getZ()<<"): "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  statusBar ()->showMessage (QString ("Mesh created"),DELAY_STATUSMSG);

  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
}

void MainWindow::on_actionSubdivide_triggered ()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  subdivide_lcc_3 (*(scene.lcc));

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to subdivide the current LCC: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
  statusBar ()->showMessage (QString ("Objects were subdivided"),
                             DELAY_STATUSMSG);
}

void MainWindow::on_actionSubdivide_pqq_triggered ()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  subdivide_lcc_pqq (*(scene.lcc));

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to pqq-subdivide the current LCC: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
  statusBar ()->showMessage (QString ("Objects were subdivided"),
                             DELAY_STATUSMSG);
}


void MainWindow::on_actionClear_triggered()
{
  clear_all();
  statusBar ()->showMessage (QString ("Scene cleared"), DELAY_STATUSMSG);
  Q_EMIT (sceneChanged ());
}

void MainWindow::on_actionCompute_Voronoi_3D_triggered ()
{
  QString fileName = QFileDialog::getOpenFileName (this,
                                                   tr ("Voronoi 3D"),
                                                   ".",
                                                   tr ("Data file (*)"));

  if (fileName.isEmpty ()) return;

  QApplication::setOverrideCursor (Qt::WaitCursor);

  this->clear_all();

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  typedef CGAL::Delaunay_triangulation_3 < LCC::Traits > Triangulation;
  Triangulation T;

  LCC delaunay_lcc;
  Dart_handle dh;

  std::ifstream ifs (qPrintable (fileName));
  T.insert (std::istream_iterator < Point_3 >(ifs),
            std::istream_iterator < Point_3 >() );

  std::map<Triangulation::Cell_handle,
      LCC::Dart_handle > vol_to_dart;

  dh = CGAL::import_from_triangulation_3 < LCC, Triangulation >
      (delaunay_lcc, T, &vol_to_dart);

  Dart_handle ddh=delaunay_lcc.dual(*scene.lcc, dh);

  // We transform all the darts in vol_to_dart into their dual.
  {
    LCC::Dart_range::iterator it1=delaunay_lcc.darts().begin();
    LCC::Dart_range::iterator it2=scene.lcc->darts().begin();

    std::map<LCC::Dart_handle, LCC::Dart_handle> dual;

    for ( ; it1!=delaunay_lcc.darts().end(); ++it1, ++it2 )
    {
      dual[it1]=it2;
    }

    // We update the geometry of dual_lcc by using the std::map face_to_dart.
    for ( std::map<Triangulation::Cell_handle, LCC::Dart_handle>
          ::iterator it=vol_to_dart.begin(), itend=vol_to_dart.end();
          it!=itend; ++it)
    {
      vol_to_dart[it->first]=dual[it->second];
      if ( !T.is_infinite(it->first) )
        scene.lcc->set_vertex_attribute
            (it->second,scene.lcc->create_vertex_attribute(T.dual(it->first)));
      /*       else
                alcc.set_vertex_attribute(it->second,alcc.create_vertex_attribute());*/
    }
  }

  // We remove the infinite volume and all its adjacent volumes.
  {
    std::stack<Dart_handle> toremove;
    LCC::size_type mark_toremove=scene.lcc->get_new_mark();
    toremove.push(ddh);
    CGAL::mark_cell<LCC,3>(*scene.lcc, ddh, mark_toremove);
    for (LCC::Dart_of_cell_range<3>::iterator
         it=scene.lcc->darts_of_cell<3>(ddh).begin(),
         itend=scene.lcc->darts_of_cell<3>(ddh).end(); it!=itend; ++it)
    {
      if ( !scene.lcc->is_marked(scene.lcc->beta(it,3), mark_toremove) )
      {
        CGAL::mark_cell<LCC,3>(*scene.lcc,
                               scene.lcc->beta(it,3), mark_toremove);
        toremove.push(scene.lcc->beta(it,3));
      }
    }
    while( !toremove.empty() )
    {
      scene.lcc->remove_cell<3>(toremove.top());
      toremove.pop();
    }
    CGAL_assertion(scene.lcc->is_whole_map_unmarked(mark_toremove));
    scene.lcc->free_mark(mark_toremove);
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to compute Voronoi 3D from "<<qPrintable(fileName)<<": "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  init_all_new_volumes();
  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
  statusBar ()->showMessage (QString ("Voronoi 3D of points in ") + fileName,
                             DELAY_STATUSMSG);
}

void MainWindow::on_actionDual_3_triggered ()
{
  if ( !scene.lcc->is_without_boundary(3) )
  {
    statusBar()->showMessage
        (QString ("Dual impossible: the lcc has some 3-boundary"),
         DELAY_STATUSMSG);
    return;
  }

  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  LCC* duallcc = new LCC;
  scene.lcc->dual_points_at_barycenter(*duallcc);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to compute the dual: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  this->clear_all();
  delete scene.lcc;
  scene.lcc = duallcc;
  this->viewer->setScene(&scene);
  init_all_new_volumes();

  statusBar ()->showMessage (QString ("Dual_3 computed"), DELAY_STATUSMSG);
  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
}

void MainWindow::on_actionClose_volume_triggered()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  if ( scene.lcc->close<3>() > 0 )
  {
    init_all_new_volumes();
    statusBar ()->showMessage (QString ("All volume(s) closed"),
                               DELAY_STATUSMSG);
    Q_EMIT (sceneChanged ());
  }
  else
    statusBar ()->showMessage
        (QString ("LCC already 3-closed"), DELAY_STATUSMSG);

  QApplication::restoreOverrideCursor ();

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to 3-close the current lcc: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

}

void MainWindow::on_actionSew3_same_facets_triggered()
{
  LCC::size_type mymark = scene.lcc->get_new_mark();
  mark_all_filled_and_visible_volumes(mymark);

  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  if ( scene.lcc->sew3_same_facets(mymark) > 0 )
  {
    statusBar()->showMessage
        (QString ("Same facets of visible and filled volume(s) are 3-sewn"),
         DELAY_STATUSMSG);
    Q_EMIT (sceneChanged ());
  }
  else
    statusBar()->showMessage (QString ("No facets 3-sewn"), DELAY_STATUSMSG);

  scene.lcc->free_mark(mymark);

  QApplication::restoreOverrideCursor ();

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to sew3 all same facets: "
           <<timer.time()<<" seconds."<<std::endl;
#endif
}

void MainWindow::on_actionUnsew3_all_triggered()
{
  unsigned int nb=0;
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  for (LCC::Dart_range::iterator it=scene.lcc->darts().begin();
       it!=scene.lcc->darts().end(); ++it)
  {
    if ( !scene.lcc->is_free(it,3) &&
         scene.lcc->info<3>(it).is_filled_and_visible() &&
         scene.lcc->info<3>(scene.lcc->beta(it,3))
           .is_filled_and_visible())
    { scene.lcc->unsew<3>(it); ++nb; }
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to unsew3 all filled volumes: "
           <<timer.time()<<" seconds."<<std::endl;
#endif
  QApplication::restoreOverrideCursor ();

  if ( nb > 0 )
  {
    statusBar()->showMessage
        (QString ("Darts between visible and filled volume(s) are 3-unsewn"),
         DELAY_STATUSMSG);
    Q_EMIT (sceneChanged ());
  }
  else
    statusBar()->showMessage (QString ("No dart 3-unsewn"), DELAY_STATUSMSG);
}

void MainWindow::on_actionInsideOut_triggered()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  LCC::size_type mymark=scene.lcc->get_new_mark();

  for (LCC::Attribute_range<3>::type::iterator
         it=scene.lcc->attributes<3>().begin(),
         itend=scene.lcc->attributes<3>().end(); it!=itend; )
  {
    LCC::Attribute_handle<3>::type cur = it++;
    if( !scene.lcc->is_marked(scene.lcc->get_attribute<3>(cur).dart(), mymark) &&
        scene.lcc->get_attribute<3>(cur).info().is_filled_and_visible() )
    {
      scene.lcc->reverse_orientation_connected_component
        (scene.lcc->get_attribute<3>(cur).dart(), mymark);
    }
  }

  // unmark all the darts by iterating on all the darts
  // but we cannot do really better
  scene.lcc->free_mark(mymark);

  QApplication::restoreOverrideCursor ();
  Q_EMIT( sceneChanged());

  statusBar()->showMessage
      (QString("Orientation of visible and filled volume(s) reversed"),
       DELAY_STATUSMSG);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to reverse the orientation of all filled volumes: "
           <<timer.time()<<" seconds."<<std::endl;
#endif
}

void MainWindow::on_actionRemove_filled_volumes_triggered()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  unsigned int count = 0;
  for (LCC::Attribute_range<3>::type::iterator
       it=scene.lcc->attributes<3>().begin(),
       itend=scene.lcc->attributes<3>().end(); it!=itend; )
  {
    LCC::Attribute_handle<3>::type cur = it++;
    if( scene.lcc->get_attribute<3>(cur).info().is_filled_and_visible() )
    {
      scene.lcc->remove_cell<3>(scene.lcc->get_attribute<3>(cur).dart());
      ++count;
    }
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to remove all filled volumes: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();
  QApplication::restoreOverrideCursor ();
  Q_EMIT( sceneChanged());

  statusBar()->showMessage
      (QString::number(count)+QString("Visible and filled volume(s) removed"),
       DELAY_STATUSMSG);
}

void MainWindow::on_actionInsert_center_vertices_triggered()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  std::vector<LCC::Dart_handle> v;
  for (LCC::One_dart_per_cell_range<2>::iterator
       it(scene.lcc->one_dart_per_cell<2>().begin()); it.cont(); ++it)
  {
    if ( scene.lcc->info<3>(it).is_filled_and_visible() )
      v.push_back(it);
  }
  for (std::vector<LCC::Dart_handle>::iterator itv(v.begin());
       itv!=v.end(); ++itv)
    scene.lcc->insert_barycenter_in_cell<2>(*itv);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to insert center vertices in all filled faces: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
  statusBar()->showMessage
      (QString ("Vertices are inserted in center of facets of visible and filled volume(s)"),
       DELAY_STATUSMSG);
}

double compute_angle3d(const Vector_3& v1, const Vector_3& v2)
{
  double a = CGAL::to_double( (v1*v2) /
                              ( sqrt(v1.squared_length()) * sqrt(v2.squared_length()) ) ) ;

  if (a < -1.0) return acos(-1.0)/CGAL_PI*180.0;
  else if (a > 1.0) return acos(1.0)/CGAL_PI*180.0;
  else return acos(a)/CGAL_PI*180.0;
}

void MainWindow::on_actionMerge_coplanar_faces_triggered()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  scene.lcc->set_update_attributes(false);

  std::vector<Dart_handle> edges;
  LCC::size_type treated  = scene.lcc->get_new_mark();
  LCC::size_type treated2 = scene.lcc->get_new_mark();

  for ( LCC::Dart_range::iterator it= scene.lcc->darts().begin(),
          itend = scene.lcc->darts().end(); it!=itend; ++it )
  {
    if (!scene.lcc->is_marked(it, treated) )
    {
      if ( scene.lcc->is_removable<1>(it) )
      {
        LCC::Vector normal1 = CGAL::compute_normal_of_cell_2(*scene.lcc,it);
        LCC::Vector normal2 = CGAL::compute_normal_of_cell_2(*scene.lcc, scene.lcc->beta<2>(it) );
        double angle = compute_angle3d(normal1, normal2);

        if ( ((angle<5.0 || angle>355.0) || (angle<185.0 && angle>175.0)) )
        {
          edges.push_back(it);
        }
      }
      CGAL::mark_cell<LCC, 1>(*scene.lcc, it, treated);
    }
  }


  for (std::vector<Dart_handle>::iterator it=edges.begin(),
         itend=edges.end(); it!=itend; ++it)
  {
    CGAL::mark_cell<LCC, 1>(*scene.lcc, *it, treated2);

    if ( scene.lcc->beta<0, 2>(*it)==*it || scene.lcc->beta<1, 2>(*it)==*it)
    { // To process dangling edges

      Dart_handle actu = *it, prev=NULL;
      do
      {
        if ( scene.lcc->beta<0, 2>(actu)==actu ) prev = scene.lcc->beta<1>(actu);
        else prev = scene.lcc->beta<0>(actu);

        if (scene.lcc->is_marked(actu, treated2) &&
            (scene.lcc->beta<0, 2>(actu)!=actu || scene.lcc->beta<1, 2>(actu)!=actu) )
        {
          scene.lcc->remove_cell<1>(actu);
          actu = prev;
        }
        else
          actu = NULL;
      }
      while (actu!=NULL && (scene.lcc->beta<0, 2>(actu)==actu || scene.lcc->beta<1, 2>(actu)==actu));
    }
    else if ( !CGAL::belong_to_same_cell<LCC, 2>(*scene.lcc, *it,
                                                 scene.lcc->beta<2>(*it)) )
      scene.lcc->remove_cell<1>(*it);
  }

  assert(scene.lcc->is_whole_map_marked(treated));
  scene.lcc->free_mark(treated);
  scene.lcc->free_mark(treated2);

  scene.lcc->set_update_attributes(true);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to merge all coplanar faces: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();

  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
  statusBar()->showMessage
      (QString ("Coplanar face(s) merged"), DELAY_STATUSMSG);
}

void MainWindow::on_actionMerge_all_volumes_triggered()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  scene.lcc->set_update_attributes(false);

  Dart_handle prev = scene.lcc->null_handle;
  for (LCC::Dart_range::iterator it(scene.lcc->darts().begin()),
       itend=scene.lcc->darts().end(); it!=itend; )
  {
    if ( !scene.lcc->is_free(it,3) &&
         scene.lcc->info<3>(it).is_filled_and_visible() &&
         scene.lcc->info<3>(scene.lcc->beta(it,3))
          .is_filled_and_visible())
    {
      scene.lcc->remove_cell<2>(it);
      itend=scene.lcc->darts().end();
      if ( prev==scene.lcc->null_handle ) it=scene.lcc->darts().begin();
      else { it=prev; if ( it!=itend ) ++it; }
    }
    else
      ++it;
  }

  scene.lcc->set_update_attributes(true);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to merge all filled volumes: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();

  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
  statusBar()->showMessage
      (QString ("Visible and filled volume(s) merged"), DELAY_STATUSMSG);
}

bool is_external(CDT::Face_handle fh)
{
  return fh->info().is_external;
}

int number_of_existing_edge(CDT::Face_handle fh)
{
  unsigned res=0;
  for (int i=0; i<3; ++i)
    if (fh->info().exist_edge[i]) ++res;
  return res;
}

int get_free_edge(CDT::Face_handle fh)
{
  CGAL_assertion( number_of_existing_edge(fh)==2 );
  for (int i=0; i<3; ++i)
    if (!fh->info().exist_edge[i]) return i;

  CGAL_assertion(false);
  return -1;
}

void constrained_delaunay_triangulation(LCC &lcc, Dart_handle d1)
{
  Vector_3 normal = CGAL::compute_normal_of_cell_2(lcc,d1);
  P_traits cdt_traits(normal);
  CDT cdt(cdt_traits);

  //inserting the constraints edge by edge
  LCC::Dart_of_orbit_range<1>::iterator
    it(lcc.darts_of_orbit<1>(d1).begin());

  CDT::Vertex_handle previous=LCC::null_handle, first=LCC::null_handle,
    vh=LCC::null_handle;

   for (LCC::Dart_of_orbit_range<1>::iterator
          itend(lcc.darts_of_orbit<1>(d1).end()); it!=itend; ++it)
   {
     vh = cdt.insert(lcc.point(it));
     vh->info().dh=it;
     if( first==NULL )
     {
       first=vh;
     }
     if( previous!=NULL)
     {
       CGAL_assertion( previous !=vh );
       cdt.insert_constraint(previous,vh);
     }

     previous=vh;
   }
   cdt.insert_constraint(previous,first);
   CGAL_assertion(cdt.is_valid());

   // sets mark is_external
   for( CDT::All_faces_iterator fit = cdt.all_faces_begin(),
          fitend = cdt.all_faces_end(); fit != fitend; ++fit)
   {
     fit->info().is_external = true;
     fit->info().is_process = false;
     fit->info().exist_edge[0]=false;
     fit->info().exist_edge[1]=false;
     fit->info().exist_edge[2]=false;
   }

   std::queue<CDT::Face_handle> face_queue;
   CDT::Face_handle face_internal = NULL;

   face_queue.push(cdt.infinite_vertex()->face());
   while(! face_queue.empty() )
   {
     CDT::Face_handle fh = face_queue.front();
     face_queue.pop();
     if(!fh->info().is_process)
     {
       fh->info().is_process = true;
       for(int i = 0; i <3; ++i)
       {
         if(!cdt.is_constrained(std::make_pair(fh, i)))
         {
           face_queue.push(fh->neighbor(i));
         }
         else if (face_internal==NULL)
         {
           face_internal = fh->neighbor(i);
         }
       }
     }
   }
   if ( face_internal!=NULL )
     face_queue.push(face_internal);

   while(! face_queue.empty() )
   {
     CDT::Face_handle fh = face_queue.front();
     face_queue.pop();
     if(!fh->info().is_process)
     {
       fh->info().is_process = true;
       fh->info().is_external = false;
       for(int i = 0; i <3; ++i)
       {
         if(!cdt.is_constrained(std::make_pair(fh, i)))
         {
           face_queue.push(fh->neighbor(i));
         }
       }
     }
   }

   for( CDT::Finite_edges_iterator eit = cdt.finite_edges_begin(),
          eitend = cdt.finite_edges_end(); eit != eitend; ++eit)
   {
     CDT::Face_handle fh = eit->first;
     int index = eit->second;
     CDT::Face_handle opposite_fh = fh->neighbor(index);
     if(cdt.is_constrained(std::make_pair(fh, index)))
     {
       fh->info().exist_edge[index]=true;
       opposite_fh->info().exist_edge[cdt.mirror_index(fh,index)]=true;

       if ( !fh->info().is_external && number_of_existing_edge(fh)==2 )
         face_queue.push(fh);
       if ( !opposite_fh->info().is_external &&
            number_of_existing_edge(opposite_fh)==2 )
         face_queue.push(opposite_fh);
     }
   }

   while( !face_queue.empty() )
   {
     CDT::Face_handle fh = face_queue.front();
     face_queue.pop();
     CGAL_assertion( number_of_existing_edge(fh)>=2 ); // i.e. ==2 or ==3
     CGAL_assertion( !fh->info().is_external );

     if (number_of_existing_edge(fh)==2)
     {
       int index = get_free_edge(fh);
       CDT::Face_handle opposite_fh = fh->neighbor(index);

       CGAL_assertion( !fh->info().exist_edge[index] );
       CGAL_assertion( !opposite_fh->info().
                       exist_edge[cdt.mirror_index(fh,index)] );
       // triangle is (vc, vb, va)
       const CDT::Vertex_handle va = fh->vertex(cdt. cw(index));
       const CDT::Vertex_handle vb = fh->vertex(cdt.ccw(index));
       const CDT::Vertex_handle vc = fh->vertex(index);

       Dart_handle dd1 = NULL;
       for (LCC::Dart_of_cell_range<0, 2>::iterator it(lcc.darts_of_cell<0, 2>(va->info().dh).begin());
            dd1==NULL && it.cont(); ++it)
       {
         if (lcc.point(lcc.beta<1>(it))==vc->point())
           dd1=it;
       }

       Dart_handle dd2 = NULL;
       for (LCC::Dart_of_cell_range<0, 2>::iterator it(lcc.darts_of_cell<0, 2>(vb->info().dh).begin());
            dd2==NULL && it.cont(); ++it)
       {
         if (lcc.point(lcc.beta<0>(it))==vc->point())
           dd2=it;
       }

       //       assert(((lcc.beta<0,0>(dd1)==dd2) || lcc.beta<1,1>(dd1)==dd2));

       Dart_handle ndart=lcc.insert_cell_1_in_cell_2(dd1, dd2);
       va->info().dh=lcc.beta<2>(ndart);

       fh->info().exist_edge[index]=true;
       opposite_fh->info().exist_edge[cdt.mirror_index(fh,index)]=true;

       if ( !opposite_fh->info().is_external &&
            number_of_existing_edge(opposite_fh)==2 )
         face_queue.push(opposite_fh);
     }
   }
}

void MainWindow::on_actionTriangulate_all_facets_triggered()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  std::vector<LCC::Dart_handle> v;
  for (LCC::One_dart_per_cell_range<2>::iterator
       it(scene.lcc->one_dart_per_cell<2>().begin()); it.cont(); ++it)
  {
    if ( scene.lcc->info<3>(it).is_filled_and_visible() ||
         (!scene.lcc->is_free<3>(it) &&
          scene.lcc->info<3>(scene.lcc->beta<3>(it)).is_filled_and_visible()) )
      v.push_back(it);
  }

  for (std::vector<LCC::Dart_handle>::iterator itv(v.begin());
       itv!=v.end(); ++itv)
    constrained_delaunay_triangulation(*scene.lcc, *itv);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to triangulate all filled faces: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();

  QApplication::restoreOverrideCursor ();
  Q_EMIT (sceneChanged ());
  statusBar()->showMessage
      (QString ("All visible and filled faces were triangulated"), DELAY_STATUSMSG);
}

bool MainWindow::is_volume_in_list(LCC::Attribute_handle<3>::type ah)
{
  for(int row=0; row < volumeList->rowCount(); ++row)
  {
    LCC::Attribute_type<3>::type* ptr=
        reinterpret_cast<LCC::Attribute_type<3>::type*>
        ( volumeList->item(row,3)->data(Qt::UserRole).value<quintptr>() );

    if(ptr==&(scene.lcc->get_attribute<3>(ah))) return true;
  }

  return false;
}

void MainWindow::update_volume_list_add(LCC::Attribute_handle<3>::type ah)
{
  // CGAL_assertion( !is_volume_in_list(ah) );

  volumeList->disconnect(this);

  int newRow = volumeList->rowCount();
  volumeList->setRowCount(newRow+1);

  QTableWidgetItem* volumeLabel = new QTableWidgetItem
    (QString((scene.lcc->get_attribute<3>(ah).info().color_name().c_str())));
  volumeLabel->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
  volumeLabel->setTextAlignment(Qt::AlignRight|Qt::AlignVCenter);
  volumeList->setItem(newRow,0,volumeLabel);

  QTableWidgetItem* fillCB = new QTableWidgetItem;
  fillCB->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
  if ( scene.lcc->get_attribute<3>(ah).info().is_filled() )
    fillCB->setCheckState(Qt::Checked);
  else
    fillCB->setCheckState(Qt::Unchecked);
  volumeList->setItem(newRow,1, fillCB);

  QTableWidgetItem* hiddenCB = new QTableWidgetItem();
  hiddenCB->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
  if ( scene.lcc->get_attribute<3>(ah).info().is_visible() )
    hiddenCB->setCheckState(Qt::Unchecked);
  else
    hiddenCB->setCheckState(Qt::Checked);
  volumeList->setItem(newRow,2,hiddenCB);

  QTableWidgetItem* attribHandle = new QTableWidgetItem;
  attribHandle->setData
      (Qt::UserRole,
       reinterpret_cast<quintptr>(&scene.lcc->get_attribute<3>(ah)));

  volumeList->setItem(newRow,3,attribHandle);

  connectVolumeListHandlers();
}

void MainWindow::update_volume_list_remove(int i)
{
  CGAL_assertion(i<volumeList->rowCount());
  volumeList->removeRow(i);
}

void MainWindow::update_volume_list_remove(LCC::Attribute_handle<3>::type ah)
{
  for(int row=0; row < volumeList->rowCount(); ++row)
  {
    LCC::Attribute_type<3>::type* ptr=
        reinterpret_cast<LCC::Attribute_type<3>::type*>
        ( volumeList->item(row,3)->data(Qt::UserRole).value<quintptr>() );

    if(ptr==&scene.lcc->get_attribute<3>(ah))
    {
      update_volume_list_remove(row);
      return;
    }
  }
}

void MainWindow::update_volume_list_all_ckeckstates()
{
  volumeList->disconnect(this);

  for(int row=0; row < volumeList->rowCount(); ++row)
  {
    LCC::Attribute_type<3>::type* ptr=
        reinterpret_cast<LCC::Attribute_type<3>::type*>
        ( volumeList->item(row,3)->data(Qt::UserRole).value<quintptr>() );

    if ( ptr->info().is_filled() )
      volumeList->item(row,1)->setCheckState(Qt::Checked);
    else
      volumeList->item(row,1)->setCheckState(Qt::Unchecked);

    if ( !ptr->info().is_visible() )
      volumeList->item(row,2)->setCheckState(Qt::Checked);
    else
      volumeList->item(row,2)->setCheckState(Qt::Unchecked);
  }

  connectVolumeListHandlers();
}

void MainWindow::recreate_whole_volume_list()
{
  volumeList->clearContents();
  volumeList->setRowCount(0);

  for (LCC::Attribute_range<3>::type::iterator
       it=scene.lcc->attributes<3>().begin(),
       itend=scene.lcc->attributes<3>().end(); it!=itend; ++it)
    update_volume_list_add(it);
}

void MainWindow::onCellChanged(int row, int col)
{
  volumeList->disconnect(this);

  LCC::Attribute_type<3>::type* ptr=
      reinterpret_cast<LCC::Attribute_type<3>::type*>
      ( volumeList->item(row,3)->data(Qt::UserRole).value<quintptr>() );

  if ( col==1 )
  {
    ptr->info().negate_filled();
  }
  else if ( col==2 )
  {
    ptr->info().negate_visible();
    if ( !ptr->info().is_visible() )
      volumeList->item(row,1)->setFlags
          (volumeList->item(row,1)->flags()^Qt::ItemIsEnabled);
    else
      volumeList->item(row,1)->setFlags
          (volumeList->item(row,1)->flags()|Qt::ItemIsEnabled);
  }

  connectVolumeListHandlers();
  Q_EMIT( sceneChanged());
}

void MainWindow::onHeaderClicked(int col)
{
  if(col != 0)
  {
    volumeList->disconnect(this);

    for(int i = 0; i < volumeList->rowCount(); ++i)
    {
      LCC::Attribute_type<3>::type* ptr=
          reinterpret_cast<LCC::Attribute_type<3>::type*>
          ( volumeList->item(i,3)->data(Qt::UserRole).value<quintptr>() );

      switch(qApp->keyboardModifiers())
      {
      case(Qt::ShiftModifier):
        if (col==1)
          ptr->info().set_filled(false);
        else if (col==2)
        {
          ptr->info().set_visible(true);
          volumeList->item(i,1)->setFlags
              (volumeList->item(i,1)->flags()|Qt::ItemIsEnabled);
        }
        volumeList->item(i,col)->setCheckState(Qt::Unchecked);
        break;
      case(Qt::ControlModifier):
        if (col==1)
          ptr->info().negate_filled();
        else if (col==2)
        {
          ptr->info().negate_visible();
          if ( !ptr->info().is_visible() )
            volumeList->item(i,1)->setFlags
                (volumeList->item(i,1)->flags()^Qt::ItemIsEnabled);
          else
            volumeList->item(i,1)->setFlags
                (volumeList->item(i,1)->flags()|Qt::ItemIsEnabled);
        }
        volumeList->item(i,col)->
            setCheckState(volumeList->item(i,col)->checkState() ?
                            Qt::Unchecked: Qt::Checked);
        break;
      default:
        if (col==1)
          ptr->info().set_filled(true);
        else if (col==2)
        {
          if ( ptr->info().is_visible() )
          {
            ptr->info().set_visible(false);
            volumeList->item(i,1)->setFlags
              (volumeList->item(i,1)->flags()^Qt::ItemIsEnabled);
          }
        }
        volumeList->item(i,col)->setCheckState(Qt::Checked);
        break;
      }
    }

    connectVolumeListHandlers();
    Q_EMIT( sceneChanged());
  }
}

void MainWindow::mark_all_filled_and_visible_volumes(LCC::size_type amark)
{
  for (LCC::Attribute_range<3>::type::iterator
       it=scene.lcc->attributes<3>().begin(),
       itend=scene.lcc->attributes<3>().end(); it!=itend; ++it)
  {
    if ( scene.lcc->get_attribute<3>(it).info().is_filled_and_visible() &&
         !scene.lcc->is_marked(it->dart(), amark) )
      CGAL::mark_cell<LCC,3>(*scene.lcc,
                             scene.lcc->get_attribute<3>(it).dart(), amark);
  }
}

void MainWindow::on_actionExtend_filled_volumes_triggered()
{
  volumeList->disconnect(this);

  std::vector<LCC::Attribute_handle<3>::type> tofill;

  LCC::size_type mark_volume = scene.lcc->get_new_mark();
  bool already_tofill;

  for (LCC::Attribute_range<3>::type::iterator
       it=scene.lcc->attributes<3>().begin(),
       itend=scene.lcc->attributes<3>().end(); it!=itend; ++it)
  {
    if ( !scene.lcc->is_marked(it->dart(), mark_volume) )
    {
      if ( !scene.lcc->get_attribute<3>(it).info().is_filled() )
      {
        already_tofill = false;
        for (LCC::Dart_of_cell_basic_range<3>::iterator it2=
               scene.lcc->darts_of_cell_basic<3>(it->dart(), mark_volume).begin();
             it2.cont(); ++it2 )
        {
          scene.lcc->mark(it2, mark_volume);
          if ( !scene.lcc->is_free(it2,3) &&
               scene.lcc->info<3>(scene.lcc->beta(it2,3)).
                 is_filled() && !already_tofill)
          {
            tofill.push_back(scene.lcc->attribute<3>(it2));
            already_tofill = true;
          }
        }
      }
      else
        CGAL::mark_cell<LCC,3>(*scene.lcc, it->dart(), mark_volume);
    }
  }

  CGAL_assertion( scene.lcc->is_whole_map_marked(mark_volume) );
  scene.lcc->free_mark(mark_volume);

  if ( tofill.size()>0 )
  {
    for ( std::vector<LCC::Attribute_handle<3>::type>::iterator
            it=tofill.begin(), itend=tofill.end(); it!=itend; ++it )
    {
      scene.lcc->get_attribute<3>(*it).info().set_filled(true);
    }

    update_volume_list_all_ckeckstates();
    Q_EMIT( sceneChanged());
  }

  connectVolumeListHandlers();
}

void MainWindow::on_actionExtend_hidden_volumes_triggered()
{
  volumeList->disconnect(this);

  std::vector<LCC::Attribute_handle<3>::type> tohide;

  LCC::size_type mark_volume = scene.lcc->get_new_mark();
  bool already_tohide;

  for (LCC::Attribute_range<3>::type::iterator
       it=scene.lcc->attributes<3>().begin(),
       itend=scene.lcc->attributes<3>().end(); it!=itend; ++it)
  {
    if ( !scene.lcc->is_marked(it->dart(), mark_volume) )
    {
      if ( scene.lcc->get_attribute<3>(it).info().is_visible() )
      {
        already_tohide = false;
        for (LCC::Dart_of_cell_basic_range<3>::iterator it2=
               scene.lcc->darts_of_cell_basic<3>(it->dart(), mark_volume).begin();
             it2.cont(); ++it2 )
        {
          scene.lcc->mark(it2, mark_volume);
          if ( !scene.lcc->is_free(it2,3) &&
               !scene.lcc->info<3>(scene.lcc->beta(it2,3)).
                  is_visible() && !already_tohide)
          {
            tohide.push_back(scene.lcc->attribute<3>(it2));
            already_tohide = true;
          }
        }
      }
      else
        CGAL::mark_cell<LCC,3>(*scene.lcc, it->dart(), mark_volume);
    }
  }

  CGAL_assertion( scene.lcc->is_whole_map_marked(mark_volume) );
  scene.lcc->free_mark(mark_volume);

  if ( tohide.size()>0 )
  {
    for ( std::vector<LCC::Attribute_handle<3>::type>::iterator
            it=tohide.begin(), itend=tohide.end(); it!=itend; ++it )
    {
      scene.lcc->get_attribute<3>(*it).info().set_visible(false);
    }

    update_volume_list_all_ckeckstates();
    Q_EMIT( sceneChanged());
  }

  connectVolumeListHandlers();
}

void MainWindow::on_actionCreate_Menger_Sponge_triggered ()
{
  dialogmenger.mengerLevel->disconnect(this);

  mengerUpdateAttributes = dialogmenger.updateAttributes->isChecked();

  dialogmenger.mengerLevel->setValue(0);
  mengerLevel=0;
  CGAL_assertion( mengerVolumes.empty() );
  mengerVolumes.push_back(on_actionCreate_cube_triggered());
  update_operations_entries(false);

  QObject::connect(dialogmenger.mengerLevel, SIGNAL(valueChanged(int)),
                   this, SLOT(onMengerChange(int)));

  dialogmenger.show();
}

void MainWindow::onMengerCancel()
{
  for(std::vector<Dart_handle>::iterator it=mengerVolumes.begin();
      it!=mengerVolumes.end(); ++it)
  {
    scene.lcc->remove_cell<3>(*it);
  }

  recreate_whole_volume_list();
  mengerVolumes.clear();
  update_operations_entries(true);
  statusBar()->showMessage (QString ("Menger sponge creation canceled"),
                            DELAY_STATUSMSG);
  Q_EMIT( sceneChanged());
}

void MainWindow::onMengerOk()
{
  update_operations_entries(true);
  mengerVolumes.clear();
}

void MainWindow::onMengerChange(int newLevel)
{
  while ( newLevel > mengerLevel ) onMengerInc();
  while ( newLevel < mengerLevel ) onMengerDec();
}

void MainWindow::onMengerUpdateAttributes(bool newValue)
{
  mengerUpdateAttributes = newValue;
}

void MainWindow::onMengerInc()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  this->mengerLevel++;

  if (!mengerUpdateAttributes)
  {
    scene.lcc->set_update_attributes(false);
  }

  std::vector<Dart_handle> edges;
  std::vector<Dart_handle> faces;
  std::size_t nbvolinit = mengerVolumes.size();

  LCC::size_type markEdges = (scene.lcc)->get_new_mark();
  LCC::size_type markFaces = (scene.lcc)->get_new_mark();
  LCC::size_type markVols  = (scene.lcc)->get_new_mark();

  for(std::vector<Dart_handle>::iterator itvol=mengerVolumes.begin();
        itvol!=mengerVolumes.end(); ++itvol)
  {
    CGAL_assertion( !(scene.lcc)->is_marked(*itvol, markVols) );
    for (LCC::Dart_of_cell_basic_range<3>::iterator
         it=(scene.lcc)->darts_of_cell_basic<3>(*itvol, markVols).begin(),
         itend=(scene.lcc)->darts_of_cell_basic<3>(*itvol, markVols).end();
         it!=itend; ++it)
    {
      if ( !(scene.lcc)->is_marked(it, markEdges) )
      {
        edges.push_back(it);
        CGAL::mark_cell<LCC,1>(*(scene.lcc), it, markEdges);
      }
      if ( !(scene.lcc)->is_marked(it, markFaces) )
      {
        faces.push_back(it);
        CGAL::mark_cell<LCC,2>(*(scene.lcc), it, markFaces);
      }
    }
  }

  (scene.lcc)->negate_mark(markVols);
  for(std::vector<Dart_handle>::iterator itvol=mengerVolumes.begin();
        itvol!=mengerVolumes.end(); ++itvol)
  {
    for (LCC::Dart_of_cell_basic_range<3>::iterator
         it=(scene.lcc)->darts_of_cell_basic<3>(*itvol, markVols).begin(),
         itend=(scene.lcc)->darts_of_cell_basic<3>(*itvol, markVols).end();
         it!=itend; ++it)
    {
      (scene.lcc)->unmark(it, markEdges);
      (scene.lcc)->unmark(it, markFaces);
    }
  }

  (scene.lcc)->negate_mark(markVols);
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markVols) );
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markFaces) );
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markEdges) );

  (scene.lcc)->free_mark(markEdges);
  (scene.lcc)->free_mark(markFaces);
  (scene.lcc)->free_mark(markVols);

  for(std::size_t i = 0; i < edges.size(); i++)
  {
    split_edge_in_three(edges[i]);
  }
  edges.clear();

  for(std::size_t i = 0; i < faces.size(); i++)
  {
    split_face_in_nine(faces[i]);
  }
  faces.clear();

  for(std::size_t i = 0; i < nbvolinit; i++)
  {
    split_vol_in_twentyseven(mengerVolumes[i]);
  }

  if (!mengerUpdateAttributes)
  {
    for(std::size_t i = nbvolinit; i < mengerVolumes.size(); i++)
    {
      LCC::Attribute_handle<3>::type ah = (scene.lcc)->create_attribute<3>();
      scene.lcc->set_attribute<3>(mengerVolumes[i], ah);
      scene.lcc->info<3>(mengerVolumes[i]).color()=
          (CGAL::Color(myrandom.get_int(0,256),
                       myrandom.get_int(0,256),
                       myrandom.get_int(0,256)));

        update_volume_list_add(scene.lcc->attribute<3>(mengerVolumes[i]));
    }

    scene.lcc->set_update_attributes(true);
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Menger sponge "
           <<this->mengerLevel-1<<" -> "<<this->mengerLevel<<", "
          <<"attributes updated "
         <<(mengerUpdateAttributes ? "DURING" : "AFTER")
        << " construction: "
        <<timer.time()<<" seconds."<<std::endl;
#endif

  CGAL_assertion( (scene.lcc)->is_valid() );

  QApplication::restoreOverrideCursor ();
  statusBar()->showMessage(QString ("Menger sponge creation %1 -> %2").
                           arg(this->mengerLevel-1).arg(this->mengerLevel),
                           DELAY_STATUSMSG);

  Q_EMIT( sceneChanged());
}

void MainWindow::split_edge_in_three(Dart_handle dh)
{
  LCC::Point p1 = scene.lcc->point(dh);
  LCC::Point p2 = scene.lcc->point(scene.lcc->other_extremity(dh));

  LCC::Vector v1 = LCC::Traits::Construct_vector() (p1,p2);
  LCC::Vector v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/3);
  LCC::Vector v3 = LCC::Traits::Construct_scaled_vector() (v1,2.0/3);

  LCC::Point p3 = LCC::Traits::Construct_translated_point() (p1,v2);
  LCC::Point p4 = LCC::Traits::Construct_translated_point() (p1,v3);

  (scene.lcc)->insert_point_in_cell<1>(dh,p4);
  (scene.lcc)->insert_point_in_cell<1>(dh,p3);
}

void MainWindow::split_face_in_three(Dart_handle dh)
{
  scene.lcc->insert_cell_1_in_cell_2(scene.lcc->beta(dh,1,1,1),
                                     scene.lcc->beta(dh,0,0));
  scene.lcc->insert_cell_1_in_cell_2(scene.lcc->beta(dh,1,1),
                                     scene.lcc->beta(dh,0));
}

void MainWindow::split_face_in_nine(Dart_handle dh)
{
  Dart_handle d2 = scene.lcc->beta(dh,1,1,1,1,1,1,1);

  Dart_handle e2= scene.lcc->insert_cell_1_in_cell_2
    (scene.lcc->beta(dh,1,1),d2);
  Dart_handle e1= scene.lcc->insert_cell_1_in_cell_2
    (scene.lcc->beta(dh,1),scene.lcc->beta(d2,1));

  split_edge_in_three(e1);
  split_edge_in_three(e2);

  split_face_in_three(dh);
  split_face_in_three(d2);
  split_face_in_three(scene.lcc->beta(e2,0));
}

void MainWindow::split_vol_in_three(Dart_handle dh, bool removecenter)
{
  std::vector<Dart_handle> edges1;
  std::vector<Dart_handle> edges2;

  Dart_handle curd = scene.lcc->beta(dh,2,1,1,2);
  for (unsigned int i=0;i<4;++i)
  {
    edges1.push_back(curd);
    curd=scene.lcc->beta(curd,1,2,1);
  }
  CGAL_assertion( curd==scene.lcc->beta(dh,2,1,1,2) );

  curd = scene.lcc->beta(curd,1,1,2);
  for (unsigned int i=0;i<4;++i)
  {
    edges2.push_back(curd);
    curd=scene.lcc->beta(curd,1,2,1);
  }
  CGAL_assertion( curd==
                  scene.lcc->beta(dh,2,1,1,2,1,1,2) );

  Dart_handle f1=
      scene.lcc->insert_cell_2_in_cell_3(edges1.begin(),edges1.end());

  Dart_handle f2=
      scene.lcc->insert_cell_2_in_cell_3(edges2.begin(),edges2.end());

  if (scene.lcc->are_attributes_automatically_managed())
  {
    scene.lcc->info<3>(f1).color()=
      (CGAL::Color(myrandom.get_int(0,256),
                   myrandom.get_int(0,256),
                   myrandom.get_int(0,256)));
    scene.lcc->info<3>(f2).color()=
      (CGAL::Color(myrandom.get_int(0,256),
                   myrandom.get_int(0,256),
                   myrandom.get_int(0,256)));

    update_volume_list_add(scene.lcc->attribute<3>(dh));
  }

  if ( removecenter )
    scene.lcc->remove_cell<3>(f1);
  else
  {
    mengerVolumes.push_back(f1);

    if (scene.lcc->are_attributes_automatically_managed())
      update_volume_list_add(scene.lcc->attribute<3>(f1));
  }

  mengerVolumes.push_back(f2);
}

void MainWindow::split_vol_in_nine(Dart_handle dh, bool removecenter)
{
  std::vector<Dart_handle> edges1;
  std::vector<Dart_handle> edges2;

  Dart_handle curd = scene.lcc->beta(dh,1,2);
  for (unsigned int i=0;i<8;++i)
  {
    edges1.push_back(curd);
    curd=scene.lcc->beta(curd,1,2,1);
  }
  CGAL_assertion( curd==scene.lcc->beta(dh,1,2) );

  curd = scene.lcc->beta(curd,1,1,2);
  for (unsigned int i=0;i<8;++i)
  {
    edges2.push_back(curd);
    curd=scene.lcc->beta(curd,1,2,1);
  }
  CGAL_assertion( curd==scene.lcc->beta(dh,1,2,1,1,2) );

  Dart_handle f1=
      scene.lcc->insert_cell_2_in_cell_3(edges1.begin(),edges1.end());

  Dart_handle f2=
      scene.lcc->insert_cell_2_in_cell_3(edges2.begin(),edges2.end());

  if (scene.lcc->are_attributes_automatically_managed())
  {
    scene.lcc->info<3>(f1).color()=
      (CGAL::Color(myrandom.get_int(0,256),
                   myrandom.get_int(0,256),
                   myrandom.get_int(0,256)));
    scene.lcc->info<3>(f2).color()=
      (CGAL::Color(myrandom.get_int(0,256),
                   myrandom.get_int(0,256),
                   myrandom.get_int(0,256)));

    update_volume_list_add(scene.lcc->attribute<3>(dh));
    if ( !removecenter)
      update_volume_list_add(scene.lcc->attribute<3>(f1));
  }

  split_face_in_three(f1);
  split_face_in_three(f2);

  split_vol_in_three(dh,removecenter);

  mengerVolumes.push_back(scene.lcc->beta(f2,2,1));
  split_vol_in_three(scene.lcc->beta(f2,2,1),removecenter);

  if ( removecenter )
    scene.lcc->remove_cell<3>(f1);
  else
  {
    mengerVolumes.push_back(scene.lcc->beta(f1,2,1));
    split_vol_in_three(scene.lcc->beta(f1,2,1),true);
  }
}

void MainWindow::split_vol_in_twentyseven(Dart_handle dh)
{
  std::vector<Dart_handle> edges1;
  std::vector<Dart_handle> edges2;

  Dart_handle curd = scene.lcc->beta(dh,1,1,2);
  for (unsigned int i=0;i<12;++i)
  {
    edges1.push_back(curd);
    curd=scene.lcc->beta(curd,1,2,1);
  }
  CGAL_assertion( curd==scene.lcc->beta(dh,1,1,2) );

  curd = scene.lcc->beta(curd,1,1,2);
  for (unsigned int i=0;i<12;++i)
  {
    edges2.push_back(curd);
    curd=scene.lcc->beta(curd,1,2,1);
  }
  CGAL_assertion( curd==scene.lcc->beta(dh,1,1,2,1,1,2) );

  Dart_handle f1=
      scene.lcc->insert_cell_2_in_cell_3(edges1.begin(),edges1.end());

  Dart_handle f2=
      scene.lcc->insert_cell_2_in_cell_3(edges2.begin(),edges2.end());

  if (scene.lcc->are_attributes_automatically_managed())
  {
    scene.lcc->info<3>(f1).color()=
      (CGAL::Color(myrandom.get_int(0,256),
                   myrandom.get_int(0,256),
                   myrandom.get_int(0,256)));
    scene.lcc->info<3>(f2).color()=
      (CGAL::Color(myrandom.get_int(0,256),
                   myrandom.get_int(0,256),
                   myrandom.get_int(0,256)));
    update_volume_list_add(scene.lcc->attribute<3>(dh));
    update_volume_list_add(scene.lcc->attribute<3>(f1));

  }

  mengerVolumes.push_back(scene.lcc->beta(f1,2));
  mengerVolumes.push_back(scene.lcc->beta(f2,2));

  split_face_in_nine(scene.lcc->beta(f1,1));
  split_face_in_nine(scene.lcc->beta(f2,1));

  split_vol_in_nine(dh,false);
  split_vol_in_nine(scene.lcc->beta(f1,2),true);
  split_vol_in_nine(scene.lcc->beta(f2,2),false);
}

void MainWindow::process_full_slice(Dart_handle init,
                                  std::vector<Dart_handle>& faces,
                                  LCC::size_type markVols)
{
  Dart_handle d[12];
  d[0]=scene.lcc->beta(init,1,2);
  d[1]=scene.lcc->beta(d[0],3,1,2,1);
  d[2]=scene.lcc->beta(d[1],1,2,1);
  d[3]=scene.lcc->beta(d[2],3,1,2,1);

  d[4]=scene.lcc->beta(init,1,1,2);
  d[5]=scene.lcc->beta(d[4],3,0,2,0);
  d[6]=scene.lcc->beta(d[5],0,2,0);

  d[7]=scene.lcc->beta(d[6],3,0,2,0);
  d[8]=scene.lcc->beta(d[7],3,0,2,0);
  d[9]=scene.lcc->beta(d[8],0,2,0);

  d[10]=scene.lcc->beta(d[9],3,0,2,0);
  d[11]=scene.lcc->beta(d[10],3,0,2,0);

  for (unsigned int j=0; j<12; ++j)
  {
    if ( !(scene.lcc)->is_marked(d[j], markVols) )
    {
      CGAL::mark_cell<LCC,3>(*(scene.lcc), d[j], markVols);
    }
    faces.push_back(d[j]);
  }
}

void MainWindow::process_inter_slice(Dart_handle init,
                                   std::vector<Dart_handle>& faces,
                                   LCC::size_type markVols)
{
  Dart_handle d[24];
  d[0]=init;
  d[1]=scene.lcc->beta(d[0],0,2,3,2,0);
  d[2]=scene.lcc->beta(d[1],0,2,3,2,0);
  d[3]=scene.lcc->beta(d[2],1,1,2,3,2);
  d[4]=scene.lcc->beta(d[3],1,1,2,3,2);
  d[5]=scene.lcc->beta(d[0],1,1,2,3,2);
  d[6]=scene.lcc->beta(d[5],1,1,2,3,2);
  d[7]=scene.lcc->beta(d[6],0,2,3,2,0);

  init = scene.lcc->beta(init,3,2,1,1,2);
  d[8]=init;
  d[9]=scene.lcc->beta(d[8],3,1,2,3,2,1);
  d[10]=scene.lcc->beta(d[9],1,2,3,2,1,3);
  d[11]=scene.lcc->beta(d[10],3,0,0,2,3,2);
  d[12]=scene.lcc->beta(d[11],0,0,2,3,2,3);
  d[13]=scene.lcc->beta(d[8],3,0,0,2,3,2);
  d[14]=scene.lcc->beta(d[13],0,0,2,3,2,3);
  d[15]=scene.lcc->beta(d[14],3,1,2,3,2,1);

  d[16]=scene.lcc->beta(d[0],3,1,2);
  d[17]=scene.lcc->beta(d[0],3,1,1,2);

  d[18]=scene.lcc->beta(d[4],3,2);
  d[19]=scene.lcc->beta(d[4],3,0,2);

  d[20]=scene.lcc->beta(d[2],3,0,2);
  d[21]=scene.lcc->beta(d[2],3,1,1,2);

  d[22]=scene.lcc->beta(d[6],3,2);
  d[23]=scene.lcc->beta(d[6],3,1,2);

  for (unsigned int j=0; j<24; ++j)
  {
    CGAL_assertion( d[j]!=(scene.lcc)->null_dart_handle );
    if ( !(scene.lcc)->is_marked(d[j], markVols) )
    {
      CGAL::mark_cell<LCC,3>(*(scene.lcc), d[j], markVols);
    }
    faces.push_back(d[j]);
  }
}

void MainWindow::onMengerDec()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  this->mengerLevel--;

  // We know here the number of Menger volume: 20^mengerLevel
  // thus we can directly "cut" the std::vector to the correct size.
  mengerVolumes.resize(CGAL::ipower(20,mengerLevel));

  LCC::size_type markVols     = (scene.lcc)->get_new_mark();
  LCC::size_type markVertices = (scene.lcc)->get_new_mark();

  std::vector<Dart_handle> faces;
  std::vector<Dart_handle> edges;
  std::vector<Dart_handle> vertices;

  // First we remove faces.
  for ( std::vector<Dart_handle>::iterator itvol=mengerVolumes.begin();
        itvol!=mengerVolumes.end(); ++itvol)
  {
    if ( !(scene.lcc)->is_marked(*itvol, markVols) )
    {
      Dart_handle init=*itvol;
      CGAL::mark_cell<LCC,3>(*(scene.lcc), init, markVols);
      process_full_slice(init, faces, markVols);
      init=scene.lcc->beta(init, 2,1,1,2);
      process_inter_slice(init, faces, markVols);
      init=scene.lcc->beta(init, 3,2,1,1,2,3);
      process_full_slice(init, faces, markVols);
    }
  }

  for(std::size_t i = 0; i < faces.size(); i++)
  {
    scene.lcc->remove_cell<2>(faces[i],mengerUpdateAttributes);
  }
  faces.clear();

  // Now we remove edges.
  for ( std::vector<Dart_handle>::iterator itvol=mengerVolumes.begin();
        itvol!=mengerVolumes.end(); ++itvol)
  {
    if ( (scene.lcc)->is_marked(*itvol, markVols) )
      CGAL::unmark_cell<LCC,3>(*(scene.lcc), *itvol, markVols);

    for (LCC::Dart_of_cell_range<3>::iterator
           it=scene.lcc->darts_of_cell<3>(*itvol).begin(),
           itend=scene.lcc->darts_of_cell<3>(*itvol).end();
         it!=itend; ++it)
    {
      if ( scene.lcc->is_free(it,2) &&
           ( scene.lcc->is_free(it,3) || it<scene.lcc->beta(it,3) ) )
        edges.push_back(it);
    }
  }

  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markVols) );

  for(std::size_t i = 0; i < edges.size(); i++)
  {
    scene.lcc->remove_cell<1>(scene.lcc->beta(edges[i],0),mengerUpdateAttributes);
    scene.lcc->remove_cell<1>(scene.lcc->beta(edges[i],1),mengerUpdateAttributes);
    scene.lcc->remove_cell<1>(edges[i],mengerUpdateAttributes);
  }
  edges.clear();

  // Lastly we remove vertices.
  for ( std::vector<Dart_handle>::iterator itvol=mengerVolumes.begin();
        itvol!=mengerVolumes.end(); ++itvol)
  {
    for (LCC::Dart_of_cell_basic_range<3>::iterator
         it=(scene.lcc)->darts_of_cell_basic<3>
           (*itvol, markVols).begin(),
         itend=(scene.lcc)->darts_of_cell_basic<3>
           (*itvol, markVols).end(); it!=itend; ++it)
    {
      if ( !(scene.lcc)->is_marked(it, markVertices) )
      {
        if ( scene.lcc->is_removable<0>(it) )
          vertices.push_back(it);
        CGAL::mark_cell<LCC, 0>(*scene.lcc, it, markVertices);
      }
    }
  }

  (scene.lcc)->negate_mark(markVols);
  for ( std::vector<Dart_handle>::iterator itvol=mengerVolumes.begin();
        itvol!=mengerVolumes.end(); ++itvol)
  {
    for (LCC::Dart_of_cell_basic_range<3>::iterator
         it=(scene.lcc)->darts_of_cell_basic<3>
           (*itvol, markVols).begin(),
         itend=(scene.lcc)->darts_of_cell_basic<3>
           (*itvol, markVols).end(); it!=itend; ++it)
    {
      if ( (scene.lcc)->is_marked(it, markVertices) )
        CGAL::unmark_cell<LCC, 0>(*scene.lcc, it, markVertices);
    }
  }

  (scene.lcc)->negate_mark(markVols);
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markVols) );
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markVertices) );

  for(std::size_t i = 0; i < vertices.size(); i++)
  {
    scene.lcc->remove_cell<0>(vertices[i],mengerUpdateAttributes);
  }
  vertices.clear();

  (scene.lcc)->free_mark(markVols);
  (scene.lcc)->free_mark(markVertices);

  if (!mengerUpdateAttributes)
  {
    scene.lcc->correct_invalid_attributes();
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Menger sponge "
           <<this->mengerLevel+1<<" -> "<<this->mengerLevel<<", "
          <<"attributes updated "
         <<(mengerUpdateAttributes ? "DURING" : "AFTER")
        << " construction: "
        <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();

  statusBar()->showMessage(QString ("Menger sponge creation %1 -> %2").
                           arg(this->mengerLevel+1).arg(this->mengerLevel),
                           DELAY_STATUSMSG);

  QApplication::restoreOverrideCursor ();
  Q_EMIT( sceneChanged());
}

///////////////////////////////////////////////////////////////////////////////////
// SIERPINSKI CARPET
///////////////////////////////////////////////////////////////////////////////////

void MainWindow::on_actionCreate_Sierpinski_Carpet_triggered ()
{
  /*neverUpdateAttributes = dialogsierpinskicarpet.never->isChecked();
  duringConstructionUpdateAttributes = dialogsierpinskicarpet.during->isChecked();
  afterConstructionUpdateAttributes = dialogsierpinskicarpet.after->isChecked();
  updateAttributesMethodStdMap = dialogsierpinskicarpet.stdmap->isChecked();
  updateAttributesMethodTraversal = dialogsierpinskicarpet.traversal->isChecked();
  // By default, the geometry will be computed after the construction
  isComputableGeometry = true;*/

  computeGeometry = false;

  sierpinskiCarpetUpdateAttributes
    = dialogsierpinskicarpet.updateAttributes->isChecked();

  dialogsierpinskicarpet.level->disconnect(this);

  dialogsierpinskicarpet.level->setValue(0);
  sierpinskiCarpetLevel=0;
  CGAL_assertion( sierpinskiCarpetSurfaces.empty() );

  Point_3 basepoint(nbcube%5, (nbcube/5)%5, nbcube/25);

  Dart_handle d = scene.lcc->make_quadrangle(basepoint,
                                             LCC::Traits::Construct_translated_point()
                                             (basepoint,LCC::Traits::Vector(1,0,0)),
                                             LCC::Traits::Construct_translated_point()
                                             (basepoint,LCC::Traits::Vector(1,1,0)),
                                             LCC::Traits::Construct_translated_point()
                                             (basepoint,LCC::Traits::Vector(0,1,0)));

  on_new_volume(d);

  ++nbcube;

  statusBar ()->showMessage (QString ("Square created"),DELAY_STATUSMSG);

  Q_EMIT (sceneChanged ());

  sierpinskiCarpetSurfaces.push_back(d);
  update_operations_entries(false);

  QObject::connect(dialogsierpinskicarpet.level, SIGNAL(valueChanged(int)),
                   this, SLOT(onSierpinskiCarpetChangeLevel(int)));

  dialogsierpinskicarpet.show();
}

void MainWindow::onSierpinskiCarpetCancel()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);
  for(std::vector<Dart_handle>::iterator it=sierpinskiCarpetSurfaces.begin();
      it!=sierpinskiCarpetSurfaces.end(); ++it)
  {
    scene.lcc->remove_cell<2>(*it);
  }

  recreate_whole_volume_list();
  sierpinskiCarpetSurfaces.clear();
  update_operations_entries(true);
  QApplication::restoreOverrideCursor ();
  Q_EMIT( sceneChanged());
}

void MainWindow::onSierpinskiCarpetOk()
{
  update_operations_entries(true);
  sierpinskiCarpetSurfaces.clear();
}

void MainWindow::onSierpinskiCarpetChangeLevel(int newLevel)
{
  while ( newLevel > sierpinskiCarpetLevel ) onSierpinskiCarpetInc();
  while ( newLevel < sierpinskiCarpetLevel ) onSierpinskiCarpetDec();
}

void MainWindow::onSierpinskiCarpetUpdateAttributes(bool newValue)
{
  sierpinskiCarpetUpdateAttributes = newValue;
}

void MainWindow::onSierpinskiCarpetInc()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  this->sierpinskiCarpetLevel++;

/*  if (computeGeometry)
  {
    // Here case where the geometry could be computed after the construction, but it was not updated.
    computeGeometry = false;
    dialogsierpinskicarpet.computeGeometry->setEnabled(false);
    //  => geometry will not be computed later.
    isComputableGeometry = false;
  }*/

  std::vector<Dart_handle> edges;
  nbfacesinit = sierpinskiCarpetSurfaces.size();

  LCC::size_type markEdges = (scene.lcc)->get_new_mark();
  LCC::size_type markFaces = (scene.lcc)->get_new_mark();

  for(std::vector<Dart_handle>::iterator itfaces=sierpinskiCarpetSurfaces.begin();
        itfaces!=sierpinskiCarpetSurfaces.end(); ++itfaces)
  {
    CGAL_assertion( !(scene.lcc)->is_marked(*itfaces, markFaces) );
    for (LCC::Dart_of_cell_basic_range<2>::iterator
         it=(scene.lcc)->darts_of_cell_basic<2>(*itfaces, markFaces).begin(),
         itend=(scene.lcc)->darts_of_cell_basic<2>(*itfaces, markFaces).end();
         it!=itend; ++it)
    {
      if ( !(scene.lcc)->is_marked(it, markEdges) )
      {
        edges.push_back(it);
        CGAL::mark_cell<LCC,1>(*(scene.lcc), it, markEdges);
      }
    }
  }

  (scene.lcc)->negate_mark(markFaces);
  for(std::vector<Dart_handle>::iterator itfaces=sierpinskiCarpetSurfaces.begin();
        itfaces!=sierpinskiCarpetSurfaces.end(); ++itfaces)
  {
    for (LCC::Dart_of_cell_basic_range<2>::iterator
         it=(scene.lcc)->darts_of_cell_basic<2>(*itfaces, markFaces).begin(),
         itend=(scene.lcc)->darts_of_cell_basic<2>(*itfaces, markFaces).end();
         it!=itend; ++it)
    {
      (scene.lcc)->unmark(it, markEdges);
    }
  }

  (scene.lcc)->negate_mark(markFaces);
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markFaces) );
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markEdges) );

  (scene.lcc)->free_mark(markEdges);
  (scene.lcc)->free_mark(markFaces);

/*  if (afterConstructionUpdateAttributes)
  {
    if (updateAttributesMethodStdMap)
    {
      // We create a map to associate embeddings to new darts
      for(std::size_t i = 0; i < edges.size(); i++)
      {
        dart_map.insert(std::pair<Dart_handle, LCC::Point>
                        (edges[i], scene.lcc->point(edges[i])));
        if (!(scene.lcc)->is_free(edges[i],2))
        {
          dart_map.insert(std::pair<Dart_handle, LCC::Point>
                          ((scene.lcc)->beta(edges[i],2),
                           scene.lcc->point((scene.lcc)->beta(edges[i],2))));
        }
      }
    }
  }*/

  for(std::size_t i = 0; i < edges.size(); i++)
  {
    sierpinski_carpet_split_edge_in_three(edges[i]);
  }
  edges.clear();

  for(std::size_t i = 0; i < nbfacesinit; i++)
  {
    sierpinski_carpet_split_face_in_nine(sierpinskiCarpetSurfaces[i]);
  }

  if (!sierpinskiCarpetUpdateAttributes)
  {
    sierpinski_carpet_update_geometry();
  }

/*  if (neverUpdateAttributes)
  {
    scene.lcc->correct_invalid_attributes();

    // Now that the map is valid, we can compute the geometry
    if (isComputableGeometry)
    {
      computeGeometry = true;
      dialogsierpinskicarpet.computeGeometry->setEnabled(true);
    }
  }*/

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Sierpinski carpet "
           <<this->sierpinskiCarpetLevel-1<<" -> "
          <<this->sierpinskiCarpetLevel<<", "
         <<"attributes updated "
        <<(sierpinskiCarpetUpdateAttributes ? "DURING" : "AFTER")
       << " construction: "
       <<timer.time()<<" seconds."<<std::endl;
#endif

  CGAL_assertion( (scene.lcc)->is_valid() );

  statusBar()->showMessage(QString ("Sierpinski carpet creation %1 -> %2").
                           arg(this->sierpinskiCarpetLevel-1).
                           arg(this->sierpinskiCarpetLevel),
                           DELAY_STATUSMSG);

  QApplication::restoreOverrideCursor ();
  Q_EMIT( sceneChanged());
}

void MainWindow::sierpinski_carpet_update_geometry()
{
/*  if (updateAttributesMethodStdMap)
  {
    for(std::size_t i = 0; i < new_darts.size(); i++)
    {
      sierpinski_carpet_copy_attributes_and_embed_vertex(new_darts[i], dart_map[new_darts[i]]);
    }

    dart_map.clear();
    new_darts.clear();
  }

  if (updateAttributesMethodTraversal)*/
  {
    LCC::size_type markVertices = (scene.lcc)->get_new_mark();

    for(std::size_t i = 0; i < nbfacesinit; i++)
    {
      // Geometry of the 4 corners of the current face
      LCC::Point p[4][4];
      Dart_handle d00 = sierpinskiCarpetSurfaces[i];
      Dart_handle d03 = scene.lcc->beta(d00,1,2,1,1,2,1,1);
      Dart_handle d33 = scene.lcc->beta(d03,1,2,1,1,2,1,1);
      Dart_handle d30 = scene.lcc->beta(d33,1,2,1,1,2,1,1);
      sierpinski_carpet_compute_4x4_geometry_matrix(p,
                                                    scene.lcc->point(d00),
                                                    scene.lcc->point(d03),
                                                    scene.lcc->point(d33),
                                                    scene.lcc->point(d30));

      Dart_handle dh = sierpinskiCarpetSurfaces[i];

      // bottom border
      dh = scene.lcc->beta(dh,1,2,1);
      if ( ! (scene.lcc)->is_marked(dh, markVertices) )
      {
        sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[0][1]);
        CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
      }
      dh = scene.lcc->beta(dh,1,2,1);
      if ( ! (scene.lcc)->is_marked(dh, markVertices) )
      {
        sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[0][2]);
        CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
      }

      // right border
      dh = scene.lcc->beta(dh,1,1,2,1);
      if ( ! (scene.lcc)->is_marked(dh, markVertices) )
      {
        sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[1][3]);
        CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
      }
      dh = scene.lcc->beta(dh,1,2,1);
      if ( ! (scene.lcc)->is_marked(dh, markVertices) )
      {
        sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[2][3]);
        CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
      }

      // top border
      dh = scene.lcc->beta(dh,1,1,2,1);
      if ( ! (scene.lcc)->is_marked(dh, markVertices) )
      {
        sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[3][2]);
        CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
      }
      dh = scene.lcc->beta(dh,1,2,1);
      if ( ! (scene.lcc)->is_marked(dh, markVertices) )
      {
        sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[3][1]);
        CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
      }

      // left border
      dh = scene.lcc->beta(dh,1,1,2,1);
      if ( ! (scene.lcc)->is_marked(dh, markVertices) )
      {
        sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[2][0]);
        CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
      }
      dh = scene.lcc->beta(dh,1,2,1);
      if ( ! (scene.lcc)->is_marked(dh, markVertices) )
      {
        sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[1][0]);
        CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
      }

      dh = sierpinskiCarpetSurfaces[i];

      // middle vertex, bottom left
      dh = scene.lcc->beta(dh,1,1);
      sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[1][1]);

      // middle vertex, top left
      dh = scene.lcc->beta(dh,2,1,1);
      sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[2][1]);

      // middle vertex, top right
      dh = scene.lcc->beta(dh,2,1,2,1,1,2,1);
      sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[2][2]);

      // middle vertex, bottom right
      dh = scene.lcc->beta(dh,2,1,1);
      sierpinski_carpet_copy_attributes_and_embed_vertex(dh, p[1][2]);
    }

    scene.lcc->unmark_all(markVertices);
    CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markVertices) );
  }
}

void MainWindow::sierpinski_carpet_compute_geometry()
{
  LCC::size_type markVertices = (scene.lcc)->get_new_mark();

  for(std::size_t i = 0; i < nbfacesinit; i++)
  {
    // on récupère la géométrie des 4 coins de la face courante
    LCC::Point p[4][4];
    Dart_handle d00 = sierpinskiCarpetSurfaces[i];
    Dart_handle d03 = scene.lcc->beta(d00,1,2,1,1,2,1,1);
    Dart_handle d33 = scene.lcc->beta(d03,1,2,1,1,2,1,1);
    Dart_handle d30 = scene.lcc->beta(d33,1,2,1,1,2,1,1);
    sierpinski_carpet_compute_4x4_geometry_matrix(p, scene.lcc->point(d00), scene.lcc->point(d03), scene.lcc->point(d33), scene.lcc->point(d30));

    Dart_handle dh = sierpinskiCarpetSurfaces[i];

    // Geometry of the 4 corners of the current face
    // bottom border
    dh = scene.lcc->beta(dh,1,2,1);
    if ( ! (scene.lcc)->is_marked(dh, markVertices) )
    {
      scene.lcc->point(dh) = p[0][1];
      CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
    }
    dh = scene.lcc->beta(dh,1,2,1);
    if ( ! (scene.lcc)->is_marked(dh, markVertices) )
    {
      scene.lcc->point(dh) = p[0][2];
      CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
    }

    // right border
    dh = scene.lcc->beta(dh,1,1,2,1);
    if ( ! (scene.lcc)->is_marked(dh, markVertices) )
    {
      scene.lcc->point(dh) = p[1][3];
      CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
    }
    dh = scene.lcc->beta(dh,1,2,1);
    if ( ! (scene.lcc)->is_marked(dh, markVertices) )
    {
      scene.lcc->point(dh) = p[2][3];
      CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
    }

    // top border
    dh = scene.lcc->beta(dh,1,1,2,1);
    if ( ! (scene.lcc)->is_marked(dh, markVertices) )
    {
      scene.lcc->point(dh) = p[3][2];
      CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
    }
    dh = scene.lcc->beta(dh,1,2,1);
    if ( ! (scene.lcc)->is_marked(dh, markVertices) )
    {
      scene.lcc->point(dh) = p[3][1];
      CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
    }

    // left border
    dh = scene.lcc->beta(dh,1,1,2,1);
    if ( ! (scene.lcc)->is_marked(dh, markVertices) )
    {
      scene.lcc->point(dh) = p[2][0];
      CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
    }
    dh = scene.lcc->beta(dh,1,2,1);
    if ( ! (scene.lcc)->is_marked(dh, markVertices) )
    {
      scene.lcc->point(dh) = p[1][0];
      CGAL::mark_cell<LCC,0>(*(scene.lcc), dh, markVertices);
    }

    dh = sierpinskiCarpetSurfaces[i];

    // middle vertex, bottom left
    dh = scene.lcc->beta(dh,1,1);
    scene.lcc->point(dh) = p[1][1];

    // middle vertex, top left
    dh = scene.lcc->beta(dh,2,1,1);
    scene.lcc->point(dh) = p[2][1];

    // middle vertex, top right
    dh = scene.lcc->beta(dh,2,1,2,1,1,2,1);
    scene.lcc->point(dh) = p[2][2];

    // middle vertex, bottom right
    dh = scene.lcc->beta(dh,2,1,1);
    scene.lcc->point(dh) = p[1][2];

    }

  scene.lcc->unmark_all(markVertices);
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markVertices) );
}

void MainWindow::sierpinski_carpet_compute_4x4_geometry_matrix
(LCC::Point p[4][4], LCC::Point& p00, LCC::Point& p03,
 LCC::Point& p33, LCC::Point& p30)
{
  p[0][0] = p00;
  p[0][3] = p03;
  p[3][3] = p33;
  p[3][0] = p30;

  LCC::Vector v1, v2, v3;

  // bottom border
  v1 = LCC::Traits::Construct_vector() (p[0][0],p[0][3]);
  v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/3);
  v3 = LCC::Traits::Construct_scaled_vector() (v1,2.0/3);
  p[0][1] = LCC::Traits::Construct_translated_point() (p[0][0],v2);
  p[0][2] = LCC::Traits::Construct_translated_point() (p[0][0],v3);
  // right border
  v1 = LCC::Traits::Construct_vector() (p[0][3],p[3][3]);
  v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/3);
  v3 = LCC::Traits::Construct_scaled_vector() (v1,2.0/3);
  p[1][3] = LCC::Traits::Construct_translated_point() (p[0][3],v2);
  p[2][3] = LCC::Traits::Construct_translated_point() (p[0][3],v3);
  // top border
  v1 = LCC::Traits::Construct_vector() (p[3][3],p[3][0]);
  v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/3);
  v3 = LCC::Traits::Construct_scaled_vector() (v1,2.0/3);
  p[3][2] = LCC::Traits::Construct_translated_point() (p[3][3],v2);
  p[3][1] = LCC::Traits::Construct_translated_point() (p[3][3],v3);
  // left border
  v1 = LCC::Traits::Construct_vector() (p[3][0],p[0][0]);
  v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/3);
  v3 = LCC::Traits::Construct_scaled_vector() (v1,2.0/3);
  p[2][0] = LCC::Traits::Construct_translated_point() (p[3][0],v2);
  p[1][0] = LCC::Traits::Construct_translated_point() (p[3][0],v3);
  // middle, left column
  v1 = LCC::Traits::Construct_vector() (p[0][1],p[3][1]);
  v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/3);
  v3 = LCC::Traits::Construct_scaled_vector() (v1,2.0/3);
  p[1][1] = LCC::Traits::Construct_translated_point() (p[0][1],v2);
  p[2][1] = LCC::Traits::Construct_translated_point() (p[0][1],v3);
  // middle, right column
  v1 = LCC::Traits::Construct_vector() (p[3][2],p[0][2]);
  v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/3);
  v3 = LCC::Traits::Construct_scaled_vector() (v1,2.0/3);
  p[2][2] = LCC::Traits::Construct_translated_point() (p[3][2],v2);
  p[1][2] = LCC::Traits::Construct_translated_point() (p[3][2],v3);
}

void MainWindow::sierpinski_carpet_copy_attributes_and_embed_vertex
(Dart_handle dh, LCC::Point& p)
{
  LCC::Attribute_handle<0>::type ah = (scene.lcc)->create_vertex_attribute(p);
  for ( LCC::Dart_of_cell_range<0>::iterator
        it=(scene.lcc)->darts_of_cell<0>(dh).begin();
        it != (scene.lcc)->darts_of_cell<0>(dh).end(); ++it )
  {
    // We copy all the attributes except for dim=0
    LCC::Helper::Foreach_enabled_attributes_except
      <CGAL::internal::Group_attribute_functor_of_dart<LCC>, 0>::
      run(*(scene.lcc),sierpinskiCarpetSurfaces[0],it);
    // We initialise the 0-atttrib to ah
    scene.lcc->set_dart_attribute<0>(it, ah);
  }
}

void MainWindow::sierpinski_carpet_split_edge_in_three(Dart_handle dh)
{
  if (sierpinskiCarpetUpdateAttributes)
  {
    LCC::Point p1 = scene.lcc->point(dh);
    LCC::Point p2 = scene.lcc->point(scene.lcc->other_extremity(dh));

    LCC::Vector v1 = LCC::Traits::Construct_vector() (p1,p2);
    LCC::Vector v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/3);
    LCC::Vector v3 = LCC::Traits::Construct_scaled_vector() (v1,2.0/3);

    LCC::Point p3 = LCC::Traits::Construct_translated_point() (p1,v2);
    LCC::Point p4 = LCC::Traits::Construct_translated_point() (p1,v3);

    (scene.lcc)->insert_point_in_cell<1>(dh,p4);
    (scene.lcc)->insert_point_in_cell<1>(dh,p3);
  }
  else
  {
    LCC::Point p3, p4;

    /*if (afterConstructionUpdateAttributes && updateAttributesMethodStdMap)
    {
      LCC::Point p1 = dart_map[dh];
      LCC::Point p2 = dart_map[scene.lcc->other_extremity(dh)];
      LCC::Vector v1 = LCC::Traits::Construct_vector() (p1,p2);
      LCC::Vector v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/3);
      LCC::Vector v3 = LCC::Traits::Construct_scaled_vector() (v1,2.0/3);

      p3 = LCC::Traits::Construct_translated_point() (p1,v2);
      p4 = LCC::Traits::Construct_translated_point() (p1,v3);
    }*/

    //    Dart_handle d1=
        scene.lcc->insert_cell_0_in_cell_1
          (dh, LCC::null_handle, false);

    /*if (afterConstructionUpdateAttributes && updateAttributesMethodStdMap)
    {
      dart_map.insert(std::pair<Dart_handle, LCC::Point>(d1, p4));
      if (!(scene.lcc)->is_free(d1,2))
      {
        dart_map.insert(std::pair<Dart_handle, LCC::Point>((scene.lcc)->beta(d1,2,1), p4));
      }
      new_darts.push_back((scene.lcc)->beta(dh,1));
    }*/

    //    Dart_handle d2=
        scene.lcc->insert_cell_0_in_cell_1
          (dh,LCC::null_handle,false);

    /*if (afterConstructionUpdateAttributes && updateAttributesMethodStdMap)
    {
      dart_map.insert(std::pair<Dart_handle, LCC::Point>(d2, p3));
      if (!(scene.lcc)->is_free(d2,2))
      {
        dart_map.insert(std::pair<Dart_handle, LCC::Point>((scene.lcc)->beta(d2,2,1), p3));
      }
      new_darts.push_back((scene.lcc)->beta(dh,1));
    }*/
  }
}

void MainWindow::sierpinski_carpet_split_face_in_three(Dart_handle dh,
                                                       bool removecenter)
{
  Dart_handle d1=
  scene.lcc->insert_cell_1_in_cell_2
    (scene.lcc->beta(dh,1,1,1),scene.lcc->beta(dh,0,0),
     sierpinskiCarpetUpdateAttributes);
  Dart_handle d2=
  scene.lcc->insert_cell_1_in_cell_2
    (scene.lcc->beta(dh,1,1),scene.lcc->beta(dh,0),
     sierpinskiCarpetUpdateAttributes);

  if ( removecenter )
  {
    scene.lcc->remove_cell<2>(d2,sierpinskiCarpetUpdateAttributes);
  }
  else
  {
    sierpinskiCarpetSurfaces.push_back(d2);
  }

  sierpinskiCarpetSurfaces.push_back(d1);
}

void MainWindow::sierpinski_carpet_split_face_in_nine(Dart_handle dh)
{
  Dart_handle d1 = scene.lcc->beta(dh,1,1);
  Dart_handle d2 = scene.lcc->beta(dh,1,1,1,1,1,1,1);
  Dart_handle d3 = scene.lcc->beta(dh,1);
  Dart_handle d4 = scene.lcc->beta(d2,1);

  Dart_handle e2=
  scene.lcc->insert_cell_1_in_cell_2
    (d1,d2,sierpinskiCarpetUpdateAttributes);

  /*if (afterConstructionUpdateAttributes && updateAttributesMethodStdMap)
  {
    dart_map.insert(std::pair<Dart_handle, LCC::Point>(e2, dart_map[d2]));
    dart_map.insert(std::pair<Dart_handle, LCC::Point>((scene.lcc)->beta(e2,2), dart_map[d1]));
  }*/

  Dart_handle e1=
  scene.lcc->insert_cell_1_in_cell_2
    (d3,d4,sierpinskiCarpetUpdateAttributes);

  /*if (afterConstructionUpdateAttributes && updateAttributesMethodStdMap)
  {
    dart_map.insert(std::pair<Dart_handle, LCC::Point>(e1, dart_map[d4]));
    dart_map.insert(std::pair<Dart_handle, LCC::Point>((scene.lcc)->beta(e1,2), dart_map[d3]));
  }*/

  sierpinskiCarpetSurfaces.push_back(e2);
  sierpinskiCarpetSurfaces.push_back(e1);

  // We give the beta2 to not insert in new_darts a dart that will be removed
  // during the removal of the middle face
  sierpinski_carpet_split_edge_in_three(scene.lcc->beta(e1,2));
  sierpinski_carpet_split_edge_in_three(e2);

  sierpinski_carpet_split_face_in_three(dh, false);
  sierpinski_carpet_split_face_in_three(d2, true);
  sierpinski_carpet_split_face_in_three(scene.lcc->beta(e2,0), false);
}

void MainWindow::onSierpinskiCarpetDec()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  this->sierpinskiCarpetLevel--;

  // We know here the number of Sierpinski surfaces: 8^sierpinskiCarpetLevel
  // thus we can directly "cut" the std::vector to the correct size.
  sierpinskiCarpetSurfaces.resize(CGAL::ipower(8,sierpinskiCarpetLevel));

  LCC::size_type markSurfaces = (scene.lcc)->get_new_mark();
  LCC::size_type markVertices = (scene.lcc)->get_new_mark();

  std::vector<Dart_handle> edges;
  std::vector<Dart_handle> vertices;

  // First we remove edges.
  for ( std::vector<Dart_handle>::iterator
        itsurfaces=sierpinskiCarpetSurfaces.begin();
        itsurfaces!=sierpinskiCarpetSurfaces.end(); ++itsurfaces)
  {
    Dart_handle dh = *itsurfaces;
    dh = scene.lcc->beta(dh,1,1,2,1);
    edges.push_back(dh);
    dh = scene.lcc->beta(dh,1,2,1,2,1);
    edges.push_back(dh);
    dh = scene.lcc->beta(dh,1,2,1,2,1);
    edges.push_back(dh);
    dh = scene.lcc->beta(dh,1,2,1,2,1);
    edges.push_back(dh);
  }

  for(std::size_t i = 0; i < edges.size(); i++)
  {
    scene.lcc->remove_cell<1>(scene.lcc->beta(edges[i],0),
                             sierpinskiCarpetUpdateAttributes);
    scene.lcc->remove_cell<1>(scene.lcc->beta(edges[i],1),
                             sierpinskiCarpetUpdateAttributes);
    scene.lcc->remove_cell<1>(edges[i],
                             sierpinskiCarpetUpdateAttributes);
  }
  edges.clear();

  // Lastly we remove vertices.
  for ( std::vector<Dart_handle>::iterator
        itsurfaces=sierpinskiCarpetSurfaces.begin();
        itsurfaces!=sierpinskiCarpetSurfaces.end(); ++itsurfaces)
  {
    Dart_handle dh = scene.lcc->beta(*itsurfaces,1);
    // we proceed side by side
    for (unsigned int i = 0; i < 4; i++)
    {
      if ( !(scene.lcc)->is_marked(dh, markVertices) )
      {
        vertices.push_back(dh);
        CGAL::mark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
      }
      dh = scene.lcc->beta(dh,1);
      if ( !(scene.lcc)->is_marked(dh, markVertices) )
      {
        vertices.push_back(dh);
        CGAL::mark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
      }
      dh = scene.lcc->beta(dh,1,1);
    }
  }

  (scene.lcc)->negate_mark(markSurfaces);
  for ( std::vector<Dart_handle>::iterator
        itsurfaces=sierpinskiCarpetSurfaces.begin();
        itsurfaces!=sierpinskiCarpetSurfaces.end(); ++itsurfaces)
  {
    Dart_handle dh = scene.lcc->beta(*itsurfaces,1);
    for (unsigned int i = 0; i < 4; i++)
    {
      if ( (scene.lcc)->is_marked(dh, markVertices) )
        CGAL::unmark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
      dh = scene.lcc->beta(dh,1);
      if ( (scene.lcc)->is_marked(dh, markVertices) )
        CGAL::unmark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
      dh = scene.lcc->beta(dh,1,1);
    }
  }

  (scene.lcc)->negate_mark(markSurfaces);
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markSurfaces) );
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markVertices) );

  for(std::size_t i = 0; i < vertices.size(); i++)
  {
    scene.lcc->remove_cell<0>(vertices[i],
                              sierpinskiCarpetUpdateAttributes);
  }
  vertices.clear();

  (scene.lcc)->free_mark(markSurfaces);
  (scene.lcc)->free_mark(markVertices);

  if (!sierpinskiCarpetUpdateAttributes)
  {
    scene.lcc->correct_invalid_attributes();
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Sierpinski carpet "
           <<this->sierpinskiCarpetLevel+1<<" -> "
          <<this->sierpinskiCarpetLevel<<", "
         <<"attributes updated "
        <<(sierpinskiCarpetUpdateAttributes ? "DURING" : "AFTER")
       << " construction: "
       <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();
  QApplication::restoreOverrideCursor ();

  statusBar()->showMessage(QString ("Sierpinski carpet creation %1 -> %2").
                           arg(this->sierpinskiCarpetLevel+1).
                           arg(this->sierpinskiCarpetLevel),
                           DELAY_STATUSMSG);
  Q_EMIT( sceneChanged());
}

///////////////////////////////////////////////////////////////////////////////////
// SIERPINSKI TRIANGLE
///////////////////////////////////////////////////////////////////////////////////
void MainWindow::on_actionCreate_Sierpinski_Triangle_triggered ()
{
  sierpinskiTriangleUpdateAttributes
    = dialogsierpinskitriangle.updateAttributes->isChecked();

  dialogsierpinskitriangle.level->disconnect(this);

  dialogsierpinskitriangle.level->setValue(0);
  sierpinskiTriangleLevel=0;
  CGAL_assertion( sierpinskiTriangleSurfaces.empty() );

  CGAL_assertion( removedTriangles.empty() );

  Point_3 basepoint(nbcube%5, (nbcube/5)%5, nbcube/25);

  Dart_handle d = scene.lcc->make_triangle(basepoint,
                                           LCC::Traits::Construct_translated_point()
                                           (basepoint,LCC::Traits::Vector(1,0,0)),
                                           LCC::Traits::Construct_translated_point()
                                           (basepoint,LCC::Traits::Vector(0.5f,CGAL::sqrt(3.f)/2.f,0)));

  on_new_volume(d);

  ++nbcube;

  statusBar ()->showMessage (QString ("Triangle created"),DELAY_STATUSMSG);

  Q_EMIT (sceneChanged ());

  sierpinskiTriangleSurfaces.push_back(d);
  update_operations_entries(false);

  QObject::connect(dialogsierpinskitriangle.level, SIGNAL(valueChanged(int)),
                   this, SLOT(onSierpinskiTriangleChangeLevel(int)));

  dialogsierpinskitriangle.show();
}

void MainWindow::onSierpinskiTriangleCancel()
{
  for(std::vector<Dart_handle>::iterator it=sierpinskiTriangleSurfaces.begin();
      it!=sierpinskiTriangleSurfaces.end(); ++it)
  {
    scene.lcc->remove_cell<2>(*it);
  }

  recreate_whole_volume_list();
  sierpinskiTriangleSurfaces.clear();
  update_operations_entries(true);
  Q_EMIT( sceneChanged());
}

void MainWindow::onSierpinskiTriangleOk()
{
  update_operations_entries(true);
  sierpinskiTriangleSurfaces.clear();
}

void MainWindow::onSierpinskiTriangleChangeLevel(int newLevel)
{
  while ( newLevel > sierpinskiTriangleLevel ) onSierpinskiTriangleInc();
  while ( newLevel < sierpinskiTriangleLevel ) onSierpinskiTriangleDec();
}

void MainWindow::onSierpinskiTriangleUpdateAttributes(bool newValue)
{
  sierpinskiTriangleUpdateAttributes = newValue;
}

void MainWindow::onSierpinskiTriangleInc()
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  this->sierpinskiTriangleLevel++;

  std::vector<Dart_handle> edges;
  nbfacesinit = sierpinskiTriangleSurfaces.size();

  LCC::size_type markEdges = (scene.lcc)->get_new_mark();
  LCC::size_type markFaces = (scene.lcc)->get_new_mark();

  for(std::vector<Dart_handle>::iterator itfaces=sierpinskiTriangleSurfaces.begin();
        itfaces!=sierpinskiTriangleSurfaces.end(); ++itfaces)
  {
    CGAL_assertion( !(scene.lcc)->is_marked(*itfaces, markFaces) );
    for (LCC::Dart_of_cell_basic_range<2>::iterator
         it=(scene.lcc)->darts_of_cell_basic<2>(*itfaces, markFaces).begin(),
         itend=(scene.lcc)->darts_of_cell_basic<2>(*itfaces, markFaces).end();
         it!=itend; ++it)
    {
      if ( !(scene.lcc)->is_marked(it, markEdges) )
      {
        edges.push_back(it);
        CGAL::mark_cell<LCC,1>(*(scene.lcc), it, markEdges);
      }
    }
  }

  (scene.lcc)->negate_mark(markFaces);
  for(std::vector<Dart_handle>::iterator itfaces=sierpinskiTriangleSurfaces.begin();
        itfaces!=sierpinskiTriangleSurfaces.end(); ++itfaces)
  {
    for (LCC::Dart_of_cell_basic_range<2>::iterator
         it=(scene.lcc)->darts_of_cell_basic<2>(*itfaces, markFaces).begin(),
         itend=(scene.lcc)->darts_of_cell_basic<2>(*itfaces, markFaces).end();
         it!=itend; ++it)
    {
      (scene.lcc)->unmark(it, markEdges);
    }
  }

  (scene.lcc)->negate_mark(markFaces);
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markFaces) );
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markEdges) );

  (scene.lcc)->free_mark(markEdges);
  (scene.lcc)->free_mark(markFaces);

  for(std::size_t i = 0; i < edges.size(); i++)
  {
    sierpinski_triangle_split_edge_in_two(edges[i]);
  }
  edges.clear();

  for(std::size_t i = 0; i < nbfacesinit; i++)
  {
    sierpinski_triangle_split_face_in_four(sierpinskiTriangleSurfaces[i],true);
  }

  if (!sierpinskiTriangleUpdateAttributes)
  {
    for(std::size_t i = nbfacesinit; i < sierpinskiTriangleSurfaces.size(); i++)
    {
      LCC::Attribute_handle<3>::type ah = (scene.lcc)->create_attribute<3>();
      scene.lcc->set_attribute<3>(sierpinskiTriangleSurfaces[i], ah);
      scene.lcc->info<3>(sierpinskiTriangleSurfaces[i]).color()=
        (CGAL::Color(myrandom.get_int(0,256),
                     myrandom.get_int(0,256),
                     myrandom.get_int(0,256)));

      update_volume_list_add(scene.lcc->attribute<3>(sierpinskiTriangleSurfaces[i]));
    }

    scene.lcc->correct_invalid_attributes();
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Sierpinski triangle "
           <<this->sierpinskiTriangleLevel-1<<" -> "
          <<this->sierpinskiTriangleLevel<<", "
         <<"attributes updated "
        <<(sierpinskiTriangleUpdateAttributes ? "DURING" : "AFTER")
       << " construction: "
       <<timer.time()<<" seconds."<<std::endl;
#endif

  //CGAL_assertion( (scene.lcc)->is_valid() );
  statusBar()->showMessage(QString ("Sierpinski triangle creation %1 -> %2").
                           arg(this->sierpinskiTriangleLevel-1).
                           arg(this->sierpinskiTriangleLevel),
                           DELAY_STATUSMSG);
  QApplication::restoreOverrideCursor ();

  Q_EMIT( sceneChanged());
}

void MainWindow::sierpinski_triangle_split_edge_in_two(Dart_handle dh)
{
  LCC::Point p1 = scene.lcc->point(dh);
  LCC::Point p2 = scene.lcc->point(scene.lcc->other_extremity(dh));

  LCC::Vector v1 = LCC::Traits::Construct_vector() (p1,p2);
  LCC::Vector v2 = LCC::Traits::Construct_scaled_vector() (v1,1.0/2);

  LCC::Point p3 = LCC::Traits::Construct_translated_point() (p1,v2);

  (scene.lcc)->insert_point_in_cell<1>(dh,p3,sierpinskiTriangleUpdateAttributes);
}

void MainWindow::sierpinski_triangle_split_face_in_four(Dart_handle dh, bool removecenter)
{
  Dart_handle d1=
  scene.lcc->insert_cell_1_in_cell_2
    (scene.lcc->beta(dh,1),scene.lcc->beta(dh,1,1,1),
     sierpinskiTriangleUpdateAttributes);

  Dart_handle d2=
  scene.lcc->insert_cell_1_in_cell_2
    (scene.lcc->beta(d1,2,1),scene.lcc->beta(d1,2,1,1,1),
     sierpinskiTriangleUpdateAttributes);

  Dart_handle d3=
  scene.lcc->insert_cell_1_in_cell_2
    (scene.lcc->beta(d2,2,1),scene.lcc->beta(d2,2,1,1,1),
     sierpinskiTriangleUpdateAttributes);
  if ( removecenter )
  {
    Triplet <Dart_handle, Dart_handle, Dart_handle> triplet(d1,d2,d3);
    removedTriangles.push_back(triplet);

    // at this step, the map is correctly 0-embedded, any other attribute is set
    //  (call of insert_point_in_cell<1> with update_attributes set to true)
    scene.lcc->remove_cell<2>(scene.lcc->beta(d3,2),sierpinskiTriangleUpdateAttributes);

    if (sierpinskiTriangleUpdateAttributes)
    {
      update_volume_list_add(scene.lcc->attribute<3>(scene.lcc->beta(d2,0)));
      update_volume_list_add(scene.lcc->attribute<3>(scene.lcc->beta(d1,0)));
    }
    else
    {
      // we dupplicate all 0-embeddings to set them to the splitted vertices
      (scene.lcc)->set_dart_attribute<0>(scene.lcc->beta(d2,1),(scene.lcc)->create_vertex_attribute(scene.lcc->point(d1)));
      (scene.lcc)->set_dart_attribute<0>(scene.lcc->beta(d3,1),(scene.lcc)->create_vertex_attribute(scene.lcc->point(d2)));
      (scene.lcc)->set_dart_attribute<0>(scene.lcc->beta(d1,1),(scene.lcc)->create_vertex_attribute(scene.lcc->point(d3)));
    }
  }
  else
  {
    sierpinskiTriangleSurfaces.push_back(scene.lcc->beta(d3,2));
  }

  sierpinskiTriangleSurfaces.push_back(scene.lcc->beta(d2,0));
  sierpinskiTriangleSurfaces.push_back(scene.lcc->beta(d1,0));
}

void MainWindow::onSierpinskiTriangleDec()
{
  QApplication::setOverrideCursor( Qt::WaitCursor );

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  this->sierpinskiTriangleLevel--;

  int nbt = CGAL::ipower(3,this->sierpinskiTriangleLevel);

  // First we add triangles removed during construction process
  for ( std::size_t i = removedTriangles.size() - nbt; i < removedTriangles.size(); i++)
  {
    Dart_handle d1 = scene.lcc->create_dart();
    Dart_handle d2 = scene.lcc->create_dart();
    Dart_handle d3 = scene.lcc->create_dart();
    scene.lcc->sew<1>(d1,d2);
    scene.lcc->sew<1>(d2,d3);
    scene.lcc->sew<1>(d3,d1);
    scene.lcc->sew<2>(d1, removedTriangles[i].first);
    scene.lcc->sew<2>(d2, removedTriangles[i].second);
    scene.lcc->sew<2>(d3, removedTriangles[i].third);
  }

  removedTriangles.resize(removedTriangles.size() - nbt);

  // We know here the number of Sierpinski surfaces: 3^sierpinskiTriangleLevel
  // thus we can directly "cut" the std::vector to the correct size.
  sierpinskiTriangleSurfaces.resize(CGAL::ipower(3,sierpinskiTriangleLevel));

  LCC::size_type markSurfaces = (scene.lcc)->get_new_mark();
  LCC::size_type markVertices = (scene.lcc)->get_new_mark();

  std::vector<Dart_handle> edges;
  std::vector<Dart_handle> vertices;

  // Now we remove edges.
  for ( std::vector<Dart_handle>::iterator
        itsurfaces=sierpinskiTriangleSurfaces.begin();
        itsurfaces!=sierpinskiTriangleSurfaces.end(); ++itsurfaces)
  {
    Dart_handle dh = *itsurfaces;
    dh = scene.lcc->beta(dh,1);
    edges.push_back(dh);
    dh = scene.lcc->beta(dh,2,1,2);
    edges.push_back(dh);
    dh = scene.lcc->beta(dh,2,1,2);
    edges.push_back(dh);
  }

  for(std::size_t i = 0; i < edges.size(); i++)
  {
    scene.lcc->remove_cell<1>(edges[i],
                              sierpinskiTriangleUpdateAttributes);
  }
  edges.clear();

  // Lastly we remove vertices.
  for ( std::vector<Dart_handle>::iterator
        itsurfaces=sierpinskiTriangleSurfaces.begin();
        itsurfaces!=sierpinskiTriangleSurfaces.end(); ++itsurfaces)
  {
    Dart_handle dh = scene.lcc->beta(*itsurfaces,1);
    if ( !(scene.lcc)->is_marked(dh, markVertices) )
    {
      vertices.push_back(dh);
      CGAL::mark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
    }
    dh = scene.lcc->beta(dh,1,1);
    if ( !(scene.lcc)->is_marked(dh, markVertices) )
    {
      vertices.push_back(dh);
      CGAL::mark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
    }
    dh = scene.lcc->beta(dh,1,1);
    if ( !(scene.lcc)->is_marked(dh, markVertices) )
    {
      vertices.push_back(dh);
      CGAL::mark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
    }
  }

  (scene.lcc)->negate_mark(markSurfaces);
  for ( std::vector<Dart_handle>::iterator
        itsurfaces=sierpinskiTriangleSurfaces.begin();
        itsurfaces!=sierpinskiTriangleSurfaces.end(); ++itsurfaces)
  {
    Dart_handle dh = scene.lcc->beta(*itsurfaces,1);
    if ( (scene.lcc)->is_marked(dh, markVertices) )
      CGAL::unmark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
    dh = scene.lcc->beta(dh,1,1);
    if ( (scene.lcc)->is_marked(dh, markVertices) )
      CGAL::unmark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
    dh = scene.lcc->beta(dh,1,1);
    if ( (scene.lcc)->is_marked(dh, markVertices) )
      CGAL::unmark_cell<LCC, 0>(*scene.lcc, dh, markVertices);
  }

  (scene.lcc)->negate_mark(markSurfaces);
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markSurfaces) );
  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markVertices) );

  for(std::size_t i = 0; i < vertices.size(); i++)
  {
    scene.lcc->remove_cell<0>(vertices[i],
                              sierpinskiTriangleUpdateAttributes);
  }
  vertices.clear();

  (scene.lcc)->free_mark(markSurfaces);
  (scene.lcc)->free_mark(markVertices);

  if (!sierpinskiTriangleUpdateAttributes)
  {
    scene.lcc->correct_invalid_attributes();
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Sierpinski triangle "
           <<this->sierpinskiTriangleLevel+1<<" -> "
          <<this->sierpinskiTriangleLevel<<", "
         <<"attributes updated "
        <<(sierpinskiTriangleUpdateAttributes ? "DURING" : "AFTER")
       << " construction: "
       <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();
  statusBar()->showMessage(QString ("Sierpinski triangle creation %1 -> %2").
                           arg(this->sierpinskiTriangleLevel+1).
                           arg(this->sierpinskiTriangleLevel),
                           DELAY_STATUSMSG);
  QApplication::restoreOverrideCursor ();

  Q_EMIT( sceneChanged());
}

#undef DELAY_STATUSMSG
