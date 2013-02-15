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
#include "MainWindow.h"
#include <CGAL/Delaunay_triangulation_3.h>
#include <QSettings>
#include "MainWindow.moc"
#include <CGAL/Timer.h>

// Function defined in Linear_cell_complex_3_subivision.cpp
void subdivide_lcc_3 (LCC & m);

// Function defined in Linear_cell_complex_pqq_subivision.cpp
void subdivide_lcc_pqq (LCC & m);

#define DELAY_STATUSMSG 1500

MainWindow::MainWindow (QWidget * parent):CGAL::Qt::DemosMainWindow (parent),
  nbcube      (0),
  dialogmesh  (this),
  dialogmenger(this)
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
  volumeList->setFixedWidth(200);
/*  volumeList->setColumnWidth(0,85);
  volumeList->setColumnWidth(1,35);
  volumeList->setColumnWidth(2,35);*/
  volumeList->horizontalHeader()->setResizeMode(QHeaderView::Stretch);
  volumeList->setSelectionMode(QAbstractItemView::NoSelection);
  //volumeList->setSelectionBehavior(QAbstractItemView::SelectRows);
  volumeListDock->setWidget(volumeList);
  addDockWidget(Qt::RightDockWidgetArea,volumeListDock);
  menuView->addAction(volumeListDock->toggleViewAction());

  QObject::connect(&dialogmesh, SIGNAL(accepted()),
                   this, SLOT(onCreateMeshOk()));

  this->viewer->setScene(&scene);

  connect_actions ();
  this->addAboutDemo (":/cgal/help/about_Linear_cell_complex_3.html");
  this->addAboutCGAL ();

  this->addRecentFiles (this->menuFile, this->actionQuit);
  connect (this, SIGNAL (openRecentFile (QString)),
           this, SLOT (load_off (QString)));

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
  QObject::connect(&dialogmenger, SIGNAL(accepted()),
                   this, SLOT(onMengerOk()));
  QObject::connect(&dialogmenger, SIGNAL(rejected()),
                   this, SLOT(onMengerCancel()));
}

void MainWindow::connectVolumeListHandlers()
{
  QObject::connect(this->volumeList, SIGNAL(cellChanged(int,int)),
                   this, SLOT(onCellChanged(int,int)));
}

void MainWindow::update_operations_entries(bool show)
{
  actionImportOFF->setEnabled(show);
  actionAddOFF->setEnabled(show);
  actionImport3DTDS->setEnabled(show);
  actionCompute_Voronoi_3D->setEnabled(show);
  actionClear->setEnabled(show);
  menuCreations->setEnabled(show);
  menuOperations->setEnabled(show);
}

void MainWindow::onSceneChanged ()
{
  int mark = scene.lcc->get_new_mark ();
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

  viewer->sceneChanged ();

  statusMessage->setText (os.str().c_str ());
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
  CGAL_assertion( adart->attribute<3>()==NULL);
  CGAL::Set_i_attribute_functor<LCC, 3>::
      run(scene.lcc, adart, scene.lcc->create_attribute<3>());
  update_volume_list_add(adart->attribute<3>());
}

void MainWindow::init_all_new_volumes()
{
  for (LCC::One_dart_per_cell_range<3>::iterator
       it(scene.lcc->one_dart_per_cell<3>().begin());
       it.cont(); ++it)
    if ( it->attribute<3>()==NULL )
    { on_new_volume(it); }
}

void MainWindow::on_actionImportOFF_triggered ()
{
  QString fileName = QFileDialog::getOpenFileName (this,
                                                   tr ("Import OFF"),
                                                   "./off",
                                                   tr ("off files (*.off)"));

  if (!fileName.isEmpty ())
  {
    load_off (fileName, true);
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
    load_3DTDS (fileName, true);
    statusBar ()->showMessage (QString ("Import 3DTDS file") + fileName,
                               DELAY_STATUSMSG);
  }
}

void MainWindow::on_actionAddOFF_triggered()
{
  QString fileName = QFileDialog::getOpenFileName (this,
                                                   tr ("Add OFF"),
                                                   "./off",
                                                   tr ("off files (*.off)"));

  if (!fileName.isEmpty ())
  {
    load_off (fileName, false);
  }
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
  emit (sceneChanged ());
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
  emit (sceneChanged ());
}

Dart_handle MainWindow::make_iso_cuboid(const Point_3 basepoint, LCC::FT lg)
{
  return scene.lcc->make_hexahedron(basepoint,
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
}

Dart_handle MainWindow::on_actionCreate_cube_triggered ()
{
  Point_3 basepoint(nbcube%5, (nbcube/5)%5, nbcube/25);

  Dart_handle d = make_iso_cuboid(basepoint, 1);

  on_new_volume(d);

  ++nbcube;

  statusBar ()->showMessage (QString ("Cube created"),DELAY_STATUSMSG);

  emit (sceneChanged ());

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

  scene.lcc->sew<3> (d1->beta(1)->beta(1)->beta(2), d2->beta(2));
  scene.lcc->sew<3> (d1->beta(2)->beta(1)->beta(1)->beta(2), d3);

  ++nbcube;

  statusBar ()->showMessage (QString ("3 cubes were created"),
                             DELAY_STATUSMSG);

  emit (sceneChanged ());
}

void MainWindow::on_actionCreate2Volumes_triggered ()
{
  Dart_handle d1 = make_iso_cuboid(Point_3(nbcube, nbcube, nbcube),1);
  Dart_handle d2 = make_iso_cuboid(Point_3(nbcube + 1, nbcube, nbcube), 1);
  Dart_handle d3 = make_iso_cuboid(Point_3(nbcube, nbcube + 1, nbcube), 1);
  Dart_handle d4 = make_iso_cuboid(Point_3(nbcube + 1, nbcube + 1, nbcube), 1);

  scene.lcc->sew<3>(d1->beta(1)->beta(1)->beta(2), d2->beta (2));
  scene.lcc->sew<3>(d1->beta(2)->beta(1)->beta(1)->beta (2), d3);

  scene.lcc->sew<3>(d3->beta(1)->beta(1)->beta(2), d4->beta (2));
  scene.lcc->sew<3>(d2->beta(2)->beta(1)->beta(1)->beta (2), d4);

  CGAL::remove_cell<LCC,2>(*scene.lcc, d3);
  CGAL::remove_cell<LCC,2>(*scene.lcc, d2->beta (2));

  on_new_volume(d1);
  on_new_volume(d4);

  ++nbcube;
  statusBar ()->showMessage (QString ("2 volumes were created"),
                             DELAY_STATUSMSG);

  emit (sceneChanged());
}

void MainWindow::on_actionCreate_mesh_triggered ()
{
  dialogmesh.show();
}

void MainWindow::onCreateMeshOk()
{
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

  statusBar ()->showMessage (QString ("mesh created"),DELAY_STATUSMSG);

  emit (sceneChanged ());
}

void MainWindow::on_actionSubdivide_triggered ()
{
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

  emit (sceneChanged ());
  statusBar ()->showMessage (QString ("Objects were subdivided"),
                             DELAY_STATUSMSG);
}

void MainWindow::on_actionSubdivide_pqq_triggered ()
{
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

  emit (sceneChanged ());
  statusBar ()->showMessage (QString ("Objects were subdivided"),
                             DELAY_STATUSMSG);
}


void MainWindow::on_actionClear_triggered()
{
  clear_all();
  statusBar ()->showMessage (QString ("Scene cleared"), DELAY_STATUSMSG);
  emit (sceneChanged ());
}

void MainWindow::on_actionCompute_Voronoi_3D_triggered ()
{
  QString fileName = QFileDialog::getOpenFileName (this,
                                                   tr ("Voronoi 3D"),
                                                   ".",
                                                   tr ("Data file (*)"));

  if (fileName.isEmpty ()) return;

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
    int mark_toremove=scene.lcc->get_new_mark();
    toremove.push(ddh);
    CGAL::mark_cell<LCC,3>(*scene.lcc, ddh, mark_toremove);
    for (LCC::Dart_of_cell_range<3>::iterator
         it=scene.lcc->darts_of_cell<3>(ddh).begin(),
         itend=scene.lcc->darts_of_cell<3>(ddh).end(); it!=itend; ++it)
    {
      if ( !scene.lcc->is_marked(it->beta(3), mark_toremove) )
      {
        CGAL::mark_cell<LCC,3>(*scene.lcc, it->beta(3), mark_toremove);
        toremove.push(it->beta(3));
      }
    }
    while( !toremove.empty() )
    {
      CGAL::remove_cell<LCC, 3>(*scene.lcc, toremove.top());
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
  emit (sceneChanged ());
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
  emit (sceneChanged ());
}

void MainWindow::on_actionClose_volume_triggered()
{
#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  if ( scene.lcc->close(3) > 0 )
  {
    init_all_new_volumes();
    statusBar ()->showMessage (QString ("All volume(s) closed"),
                               DELAY_STATUSMSG);
    emit (sceneChanged ());
  }
  else
    statusBar ()->showMessage
        (QString ("LCC already 3-closed"), DELAY_STATUSMSG);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to 3-close the current lcc: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

}

void MainWindow::on_actionSew3_same_facets_triggered()
{
  int mymark = scene.lcc->get_new_mark();
  mark_all_filled_and_visible_volumes(mymark);

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  if ( scene.lcc->sew3_same_facets(mymark) > 0 )
  {
    statusBar()->showMessage
        (QString ("Same facets of visible and filled volume(s) are 3-sewn"),
         DELAY_STATUSMSG);
    emit (sceneChanged ());
  }
  else
    statusBar()->showMessage (QString ("No facets 3-sewn"), DELAY_STATUSMSG);

  scene.lcc->free_mark(mymark);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to sew3 all same facets: "
           <<timer.time()<<" seconds."<<std::endl;
#endif
}

void MainWindow::on_actionUnsew3_all_triggered()
{
  unsigned int nb=0;

#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  for (LCC::Dart_range::iterator it=scene.lcc->darts().begin();
       it!=scene.lcc->darts().end(); ++it)
  {
    if ( !it->is_free(3) &&
         it->attribute<3>()->info().is_filled_and_visible() &&
         it->beta(3)->attribute<3>()->info().is_filled_and_visible())
    { scene.lcc->unsew<3>(it); ++nb; }
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to unsew3 all filled volumes: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  if ( nb > 0 )
  {
    statusBar()->showMessage
        (QString ("Darts between visible and filled volume(s) are 3-unsewn"),
         DELAY_STATUSMSG);
    emit (sceneChanged ());
  }
  else
    statusBar()->showMessage (QString ("No dart 3-unsewn"), DELAY_STATUSMSG);
}

void MainWindow::on_actionRemove_filled_volumes_triggered()
{
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
    if( cur->info().is_filled_and_visible() )
    {
      CGAL::remove_cell<LCC,3>(*scene.lcc,cur->dart());
      ++count;
    }
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to remove all filled volumes: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();
  emit(sceneChanged());

  statusBar()->showMessage
      (QString::number(count)+QString("Visible and filled volume(s) removed"),
       DELAY_STATUSMSG);
}

void MainWindow::on_actionTriangulate_all_facets_triggered()
{
#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  std::vector<LCC::Dart_handle> v;
  for (LCC::One_dart_per_cell_range<2>::iterator
       it(scene.lcc->one_dart_per_cell<2>().begin()); it.cont(); ++it)
  {
    if ( it->attribute<3>()->info().is_filled_and_visible() )
      v.push_back(it);
  }
  for (std::vector<LCC::Dart_handle>::iterator itv(v.begin());
       itv!=v.end(); ++itv)
    scene.lcc->insert_barycenter_in_cell<2>(*itv);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to triangulate all filled faces: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  emit (sceneChanged ());
  statusBar()->showMessage
      (QString ("Facets of visible and filled volume(s) triangulated"),
       DELAY_STATUSMSG);
}

void MainWindow::on_actionMerge_all_volumes_triggered()
{
#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  Dart_handle prev = NULL;
  for (LCC::Dart_range::iterator it(scene.lcc->darts().begin()),
       itend=scene.lcc->darts().end(); it!=itend; )
  {
    if ( !it->is_free(3) &&
         it->attribute<3>()->info().is_filled_and_visible() &&
         it->beta(3)->attribute<3>()->info().is_filled_and_visible() )
    {
      CGAL::remove_cell<LCC,2>(*scene.lcc,it);
      itend=scene.lcc->darts().end();
      if ( prev==NULL ) it=scene.lcc->darts().begin();
      else { it=prev; if ( it!=itend ) ++it; }
    }
    else
      ++it;
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to merge all filled volumes: "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();

  emit (sceneChanged ());
  statusBar()->showMessage
      (QString ("All volume(s) merged"), DELAY_STATUSMSG);
}

bool MainWindow::is_volume_in_list(LCC::Attribute_handle<3>::type ah)
{
  for(int row=0; row < volumeList->rowCount(); ++row)
  {
    LCC::Attribute_type<3>::type* ptr=
        reinterpret_cast<LCC::Attribute_type<3>::type*>
        ( volumeList->item(row,3)->data(Qt::UserRole).value<quintptr>() );

    if(ptr==&*ah) return true;
  }

  return false;
}

void MainWindow::update_volume_list_add(LCC::Attribute_handle<3>::type ah)
{
  CGAL_assertion( !is_volume_in_list(ah) );

  volumeList->disconnect(this);

  int newRow = volumeList->rowCount();
  volumeList->setRowCount(newRow+1);

  QTableWidgetItem* volumeLabel = new QTableWidgetItem
      (QString(ah->info().color_name().c_str()));
  volumeLabel->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
  volumeLabel->setTextAlignment(Qt::AlignRight|Qt::AlignVCenter);
  volumeList->setItem(newRow,0,volumeLabel);

  QTableWidgetItem* fillCB = new QTableWidgetItem;
  fillCB->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
  if ( ah->info().is_filled() ) fillCB->setCheckState(Qt::Checked);
  else                          fillCB->setCheckState(Qt::Unchecked);
  volumeList->setItem(newRow,1, fillCB);

  QTableWidgetItem* hiddenCB = new QTableWidgetItem();
  hiddenCB->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
  if ( ah->info().is_visible() ) hiddenCB->setCheckState(Qt::Unchecked);
  else                           hiddenCB->setCheckState(Qt::Checked);
  volumeList->setItem(newRow,2,hiddenCB);

  QTableWidgetItem* attribHandle = new QTableWidgetItem;
  attribHandle->setData
      (Qt::UserRole, reinterpret_cast<quintptr>(&*(ah)));

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

    if(ptr==&*ah)
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

  emit(sceneChanged());
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
        if (col==1) ptr->info().set_filled(false);
        else if (col==2)
        {
          ptr->info().set_visible(true);
          volumeList->item(i,1)->setFlags
              (volumeList->item(i,1)->flags()|Qt::ItemIsEnabled);
        }
        volumeList->item(i,col)->setCheckState(Qt::Unchecked);
        break;
      case(Qt::ControlModifier):
        if (col==1) ptr->info().negate_filled();
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
        if (col==1) ptr->info().set_filled(true);
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
    emit(sceneChanged());
  }
}

void MainWindow::mark_all_filled_and_visible_volumes(int amark)
{
  for (LCC::Attribute_range<3>::type::iterator
       it=scene.lcc->attributes<3>().begin(),
       itend=scene.lcc->attributes<3>().end(); it!=itend; ++it)
  {
    if ( it->info().is_filled_and_visible() &&
         !scene.lcc->is_marked(it->dart(), amark) )
      CGAL::mark_cell<LCC,3>(*scene.lcc, it->dart(), amark);
  }
}

void MainWindow::on_actionExtend_filled_volumes_triggered()
{
  volumeList->disconnect(this);

  std::vector<LCC::Attribute_handle<3>::type> tofill;

  int mark_volume = scene.lcc->get_new_mark();
  bool already_tofill;

  for (LCC::Attribute_range<3>::type::iterator
       it=scene.lcc->attributes<3>().begin(),
       itend=scene.lcc->attributes<3>().end(); it!=itend; ++it)
  {
    if ( !scene.lcc->is_marked(it->dart(), mark_volume) )
    {
      if ( !it->info().is_filled() )
      {
        already_tofill = false;
        for (LCC::Dart_of_cell_basic_range<3>::iterator it2=
               scene.lcc->darts_of_cell_basic<3>(it->dart(), mark_volume).begin();
             it2.cont(); ++it2 )
        {
          scene.lcc->mark(it2, mark_volume);
          if ( !it2->is_free(3) &&
               it2->beta(3)->attribute<3>()->info().is_filled() &&
               !already_tofill)
          {
            tofill.push_back(it2->attribute<3>());
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
      (*it)->info().set_filled(true);
    }

    update_volume_list_all_ckeckstates();
    emit(sceneChanged());
  }

  connectVolumeListHandlers();
}

void MainWindow::on_actionExtend_hidden_volumes_triggered()
{
  volumeList->disconnect(this);

  std::vector<LCC::Attribute_handle<3>::type> tohide;

  int mark_volume = scene.lcc->get_new_mark();
  bool already_tohide;

  for (LCC::Attribute_range<3>::type::iterator
       it=scene.lcc->attributes<3>().begin(),
       itend=scene.lcc->attributes<3>().end(); it!=itend; ++it)
  {
    if ( !scene.lcc->is_marked(it->dart(), mark_volume) )
    {
      if ( it->info().is_visible() )
      {
        already_tohide = false;
        for (LCC::Dart_of_cell_basic_range<3>::iterator it2=
               scene.lcc->darts_of_cell_basic<3>(it->dart(), mark_volume).begin();
             it2.cont(); ++it2 )
        {
          scene.lcc->mark(it2, mark_volume);
          if ( !it2->is_free(3) &&
               !it2->beta(3)->attribute<3>()->info().is_visible() &&
               !already_tohide)
          {
            tohide.push_back(it2->attribute<3>());
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
      (*it)->info().set_visible(false);
    }

    update_volume_list_all_ckeckstates();
    emit(sceneChanged());
  }

  connectVolumeListHandlers();
}

void MainWindow::on_actionCreate_Menger_Sponge_triggered ()
{
  dialogmenger.mengerLevel->disconnect(this);

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
    CGAL::remove_cell<LCC,3>(*scene.lcc, *it);
  }

  recreate_whole_volume_list();
  mengerVolumes.clear();
  emit(sceneChanged());
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

void MainWindow::onMengerInc()
{
#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  this->mengerLevel++;

  std::vector<Dart_handle> edges;
  std::vector<Dart_handle> faces;
  unsigned int nbvolinit = (unsigned int)mengerVolumes.size();

  int markEdges = (scene.lcc)->get_new_mark();
  int markFaces = (scene.lcc)->get_new_mark();
  int markVols  = (scene.lcc)->get_new_mark();

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

  for(unsigned int i = 0; i < (unsigned int)edges.size(); i++)
  {
    split_edge_in_three(edges[i]);
  }
  edges.clear();

  for(unsigned int i = 0; i < (unsigned int)faces.size(); i++)
  {
    split_face_in_nine(faces[i]);
  }
  faces.clear();

  for(unsigned int i = 0; i < nbvolinit; i++)
  {
    split_vol_in_twentyseven(mengerVolumes[i]);
  }

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to increase the level of menger sponge ("
           <<this->mengerLevel-1<<" -> "<<this->mengerLevel<<"): "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  CGAL_assertion( (scene.lcc)->is_valid() );

  emit(sceneChanged());
}

void MainWindow::split_edge_in_three(Dart_handle dh)
{
  LCC::Point p1 = LCC::point(dh);
  LCC::Point p2 = LCC::point(dh->other_extremity());

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
  CGAL::insert_cell_1_in_cell_2(*(scene.lcc),dh->beta(1)->beta(1)->beta(1),
                                dh->beta(0)->beta(0));
  CGAL::insert_cell_1_in_cell_2(*(scene.lcc),dh->beta(1)->beta(1),dh->beta(0));
}

void MainWindow::split_face_in_nine(Dart_handle dh)
{
  Dart_handle d2 = dh->beta(1)->beta(1)->beta(1)->beta(1)
      ->beta(1)->beta(1)->beta(1);

  Dart_handle e2= CGAL::insert_cell_1_in_cell_2(*(scene.lcc),
                                                dh->beta(1)->beta(1),d2);
  Dart_handle e1= CGAL::insert_cell_1_in_cell_2(*(scene.lcc),
                                                dh->beta(1),d2->beta(1));

  split_edge_in_three(e1);
  split_edge_in_three(e2);

  split_face_in_three(dh);
  split_face_in_three(d2);
  split_face_in_three(e2->beta(0));
}

void MainWindow::split_vol_in_three(Dart_handle dh, bool removecenter)
{
  std::vector<Dart_handle> edges1;
  std::vector<Dart_handle> edges2;

  Dart_handle curd = dh->beta(2)->beta(1)->beta(1)->beta(2);
  for (unsigned int i=0;i<4;++i)
  {
    edges1.push_back(curd);
    curd=curd->beta(1)->beta(2)->beta(1);
  }
  CGAL_assertion( curd==dh->beta(2)->beta(1)->beta(1)->beta(2) );

  curd = curd->beta(1)->beta(1)->beta(2);
  for (unsigned int i=0;i<4;++i)
  {
    edges2.push_back(curd);
    curd=curd->beta(1)->beta(2)->beta(1);
  }
  CGAL_assertion( curd==
          dh->beta(2)->beta(1)->beta(1)->beta(2)->beta(1)->beta(1)->beta(2) );

  Dart_handle f1=
      insert_cell_2_in_cell_3(*(scene.lcc),edges1.begin(),edges1.end());

  Dart_handle f2=
      insert_cell_2_in_cell_3(*(scene.lcc),edges2.begin(),edges2.end());

  f1->attribute<3>()->info().color()=
    (CGAL::Color(myrandom.get_int(0,256),
                 myrandom.get_int(0,256),
                 myrandom.get_int(0,256)));
  f2->attribute<3>()->info().color()=
      (CGAL::Color(myrandom.get_int(0,256),
                   myrandom.get_int(0,256),
                   myrandom.get_int(0,256)));

  update_volume_list_add(dh->attribute<3>());

  if ( removecenter )
    CGAL::remove_cell<LCC,3>(*scene.lcc,f1);
  else
  {
    mengerVolumes.push_back(f1);
    update_volume_list_add(f1->attribute<3>());
  }

  mengerVolumes.push_back(f2);
}

void MainWindow::split_vol_in_nine(Dart_handle dh, bool removecenter)
{
  std::vector<Dart_handle> edges1;
  std::vector<Dart_handle> edges2;

  Dart_handle curd = dh->beta(1)->beta(2);
  for (unsigned int i=0;i<8;++i)
  {
    edges1.push_back(curd);
    curd=curd->beta(1)->beta(2)->beta(1);
  }
  CGAL_assertion( curd==dh->beta(1)->beta(2) );

  curd = curd->beta(1)->beta(1)->beta(2);
  for (unsigned int i=0;i<8;++i)
  {
    edges2.push_back(curd);
    curd=curd->beta(1)->beta(2)->beta(1);
  }
  CGAL_assertion( curd==dh->beta(1)->beta(2)->beta(1)->beta(1)->beta(2) );

  Dart_handle f1=
      insert_cell_2_in_cell_3(*(scene.lcc),edges1.begin(),edges1.end());

  Dart_handle f2=
      insert_cell_2_in_cell_3(*(scene.lcc),edges2.begin(),edges2.end());

  f1->attribute<3>()->info().color()=
    (CGAL::Color(myrandom.get_int(0,256),
                 myrandom.get_int(0,256),
                 myrandom.get_int(0,256)));
  f2->attribute<3>()->info().color()=
    (CGAL::Color(myrandom.get_int(0,256),
                 myrandom.get_int(0,256),
                 myrandom.get_int(0,256)));

  update_volume_list_add(dh->attribute<3>());
  if ( !removecenter)
    update_volume_list_add(f1->attribute<3>());

  split_face_in_three(f1);
  split_face_in_three(f2);

  split_vol_in_three(dh,removecenter);

  mengerVolumes.push_back(f2->beta(2)->beta(1));
  split_vol_in_three(f2->beta(2)->beta(1),removecenter);

  if ( removecenter )
    CGAL::remove_cell<LCC,3>(*scene.lcc,f1);
  else
  {
    mengerVolumes.push_back(f1->beta(2)->beta(1));
    split_vol_in_three(f1->beta(2)->beta(1),true);
  }
}

void MainWindow::split_vol_in_twentyseven(Dart_handle dh)
{
  std::vector<Dart_handle> edges1;
  std::vector<Dart_handle> edges2;

  Dart_handle curd = dh->beta(1)->beta(1)->beta(2);
  for (unsigned int i=0;i<12;++i)
  {
    edges1.push_back(curd);
    curd=curd->beta(1)->beta(2)->beta(1);
  }
  CGAL_assertion( curd==dh->beta(1)->beta(1)->beta(2) );

  curd = curd->beta(1)->beta(1)->beta(2);
  for (unsigned int i=0;i<12;++i)
  {
    edges2.push_back(curd);
    curd=curd->beta(1)->beta(2)->beta(1);
  }
  CGAL_assertion( curd==dh->beta(1)->beta(1)->beta(2)->beta(1)->beta(1)->beta(2) );

  Dart_handle f1=
      insert_cell_2_in_cell_3(*(scene.lcc),edges1.begin(),edges1.end());

  Dart_handle f2=
      insert_cell_2_in_cell_3(*(scene.lcc),edges2.begin(),edges2.end());

  f1->attribute<3>()->info().color()=
    (CGAL::Color(myrandom.get_int(0,256),
                 myrandom.get_int(0,256),
                 myrandom.get_int(0,256)));
  f2->attribute<3>()->info().color()=
    (CGAL::Color(myrandom.get_int(0,256),
                 myrandom.get_int(0,256),
                 myrandom.get_int(0,256)));

  update_volume_list_add(dh->attribute<3>());
  update_volume_list_add(f1->attribute<3>());

  mengerVolumes.push_back(f1->beta(2));
  mengerVolumes.push_back(f2->beta(2));

  split_face_in_nine(f1->beta(1));
  split_face_in_nine(f2->beta(1));

  split_vol_in_nine(dh,false);
  split_vol_in_nine(f1->beta(2),true);
  split_vol_in_nine(f2->beta(2),false);
}

void MainWindow::process_full_slice(Dart_handle init,
                                  std::vector<Dart_handle>& faces,
                                  int markVols)
{
  Dart_handle d[12];
  d[0]=init->beta(1)->beta(2);
  d[1]=d[0]->beta(3)->beta(1)->beta(2)->beta(1);
  d[2]=d[1]->beta(1)->beta(2)->beta(1);
  d[3]=d[2]->beta(3)->beta(1)->beta(2)->beta(1);

  d[4]=init->beta(1)->beta(1)->beta(2);
  d[5]=d[4]->beta(3)->beta(0)->beta(2)->beta(0);
  d[6]=d[5]->beta(0)->beta(2)->beta(0);

  d[7]=d[6]->beta(3)->beta(0)->beta(2)->beta(0);
  d[8]=d[7]->beta(3)->beta(0)->beta(2)->beta(0);
  d[9]=d[8]->beta(0)->beta(2)->beta(0);

  d[10]=d[9]->beta(3)->beta(0)->beta(2)->beta(0);
  d[11]=d[10]->beta(3)->beta(0)->beta(2)->beta(0);

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
                                   int markVols)
{
  Dart_handle d[24];
  d[0]=init;
  d[1]=d[0]->beta(0)->beta(2)->beta(3)->beta(2)->beta(0);
  d[2]=d[1]->beta(0)->beta(2)->beta(3)->beta(2)->beta(0);
  d[3]=d[2]->beta(1)->beta(1)->beta(2)->beta(3)->beta(2);
  d[4]=d[3]->beta(1)->beta(1)->beta(2)->beta(3)->beta(2);
  d[5]=d[0]->beta(1)->beta(1)->beta(2)->beta(3)->beta(2);
  d[6]=d[5]->beta(1)->beta(1)->beta(2)->beta(3)->beta(2);
  d[7]=d[6]->beta(0)->beta(2)->beta(3)->beta(2)->beta(0);

  init = init->beta(3)->beta(2)->beta(1)->beta(1)->beta(2);
  d[8]=init;
  d[9]=d[8]->beta(3)->beta(1)->beta(2)->beta(3)->beta(2)->beta(1);
  d[10]=d[9]->beta(1)->beta(2)->beta(3)->beta(2)->beta(1)->beta(3);
  d[11]=d[10]->beta(3)->beta(0)->beta(0)->beta(2)->beta(3)->beta(2);
  d[12]=d[11]->beta(0)->beta(0)->beta(2)->beta(3)->beta(2)->beta(3);
  d[13]=d[8]->beta(3)->beta(0)->beta(0)->beta(2)->beta(3)->beta(2);
  d[14]=d[13]->beta(0)->beta(0)->beta(2)->beta(3)->beta(2)->beta(3);
  d[15]=d[14]->beta(3)->beta(1)->beta(2)->beta(3)->beta(2)->beta(1);

  d[16]=d[0]->beta(3)->beta(1)->beta(2);
  d[17]=d[0]->beta(3)->beta(1)->beta(1)->beta(2);

  d[18]=d[4]->beta(3)->beta(2);
  d[19]=d[4]->beta(3)->beta(0)->beta(2);

  d[20]=d[2]->beta(3)->beta(0)->beta(2);
  d[21]=d[2]->beta(3)->beta(1)->beta(1)->beta(2);

  d[22]=d[6]->beta(3)->beta(2);
  d[23]=d[6]->beta(3)->beta(1)->beta(2);

  for (unsigned int j=0; j<24; ++j)
  {
    CGAL_assertion( d[j]!=LCC::null_dart_handle );
    if ( !(scene.lcc)->is_marked(d[j], markVols) )
    {
      CGAL::mark_cell<LCC,3>(*(scene.lcc), d[j], markVols);
    }
    faces.push_back(d[j]);
  }
}

void MainWindow::onMengerDec()
{
#ifdef CGAL_PROFILE_LCC_DEMO
  CGAL::Timer timer;
  timer.start();
#endif

  this->mengerLevel--;

  // We know here the number of Menger volume: 20^mengerLevel
  // thus we can directly "cut" the std::vector to the correct size.
  mengerVolumes.resize(CGAL::ipower(20,mengerLevel));

  int markVols     = (scene.lcc)->get_new_mark();
  int markVertices = (scene.lcc)->get_new_mark();

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
      init=init->beta(2)->beta(1)->beta(1)->beta(2);
      process_inter_slice(init, faces, markVols);
      init=init->beta(3)->beta(2)->beta(1)->beta(1)->beta(2)->beta(3);
      process_full_slice(init, faces, markVols);
    }
  }

  for(unsigned int i = 0; i < faces.size(); i++)
  {
    CGAL::remove_cell<LCC,2>(*scene.lcc, faces[i]);
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
      if ( it->is_free(2) && ( it->is_free(3) || &*it<&*it->beta(3) ) )
        edges.push_back(it);
    }
  }

  CGAL_assertion( (scene.lcc)->is_whole_map_unmarked(markVols) );

  for(unsigned int i = 0; i < edges.size(); i++)
  {
    CGAL::remove_cell<LCC,1>(*scene.lcc, edges[i]->beta(0));
    CGAL::remove_cell<LCC,1>(*scene.lcc, edges[i]->beta(1));
    CGAL::remove_cell<LCC,1>(*scene.lcc, edges[i]);
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
        if ( CGAL::is_removable<LCC, 0>(*scene.lcc, it) )
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

  for(unsigned int i = 0; i < vertices.size(); i++)
  {
    CGAL::remove_cell<LCC,0>(*scene.lcc, vertices[i]);
  }
  vertices.clear();

  (scene.lcc)->free_mark(markVols);
  (scene.lcc)->free_mark(markVertices);

#ifdef CGAL_PROFILE_LCC_DEMO
  timer.stop();
  std::cout<<"Time to decrease the level of menger sponge ("
           <<this->mengerLevel+1<<" -> "<<this->mengerLevel<<"): "
           <<timer.time()<<" seconds."<<std::endl;
#endif

  recreate_whole_volume_list();

  statusBar ()->showMessage (QString ("Menger Dec"),DELAY_STATUSMSG);
  emit(sceneChanged());
}


#undef DELAY_STATUSMSG
