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

// Function defined in Linear_cell_complex_3_subivision.cpp
void subdivide_lcc_3 (LCC & m);

#define DELAY_STATUSMSG 1500

bool isVisibleAndFilled(char property)
{ return (property & LCC_DEMO_FILLED) && (property & LCC_DEMO_VISIBLE); }
  
bool isVisible(char property)
{ return (property & LCC_DEMO_VISIBLE); }
  
bool isFilled(char property)
{ return (property & LCC_DEMO_FILLED); }

char setVisible(char property)
{ return property | LCC_DEMO_VISIBLE; }

char setHidden(char property)
{
  if ( !isVisible(property) ) return property;
  return property ^ LCC_DEMO_VISIBLE;
}

char setFilled(char property)
{ return property | LCC_DEMO_FILLED; }
  
char setWireframe(char property)
{
  if ( !isFilled(property) ) return property;
  return property ^ LCC_DEMO_FILLED;
}
  
char setVisibleAndFilled(char property)
{ return property | LCC_DEMO_FILLED | LCC_DEMO_VISIBLE; }

char negateVisible(char property)
{ return property ^ LCC_DEMO_VISIBLE; }

char negateFilled(char property)
{ return property ^ LCC_DEMO_FILLED; }


MainWindow::MainWindow (QWidget * parent):CGAL::Qt::DemosMainWindow (parent),
					  nbcube (0),
                                          dialogmesh(this),
                                          volumeUid(1)
{
  setupUi (this);
  scene.lcc = new LCC;
  
  volumeListDock = new QDockWidget(QString(tr("Volume List")),this);
  volumeListDock->setAllowedAreas(Qt::RightDockWidgetArea |
                                  Qt::LeftDockWidgetArea);
  volumeList = new QTableWidget(0,3,volumeListDock);
  volumeList->verticalHeader()->hide();
  QObject::connect(this->volumeList, SIGNAL(cellChanged(int,int)),
                   this, SLOT(onCellChanged(int,int)));

  QStringList labels(QString(tr("Volume")));
  labels.append(QString(tr("Filled")));
  labels.append(QString(tr("Hidden")));
  volumeList->setHorizontalHeaderLabels(labels);
  //  volumeList->resizeColumnsToContents();
  volumeList->setFixedWidth(170);  
  volumeList->horizontalHeader()->setResizeMode(QHeaderView::Stretch);
  volumeList->setSelectionMode(QAbstractItemView::SingleSelection);
  volumeList->setSelectionBehavior(QAbstractItemView::SelectRows);
  volumeListDock->setWidget(volumeList);
  addDockWidget(Qt::RightDockWidgetArea,volumeListDock);
  menuView->addAction(volumeListDock->toggleViewAction());

  this->viewer->setScene (&scene);
  this->viewer->setVectorPointers(&volumeDartIndex,&volumeProperties);

  connectActions ();
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

void MainWindow::connectActions ()
{
  QObject::connect (this->actionImportOFF, SIGNAL (triggered ()),
                    this, SLOT (import_off ()));

  QObject::connect (this->actionAddOFF, SIGNAL (triggered ()),
                    this, SLOT (add_off ()));

  QObject::connect (this->actionCompute_Voronoi_3D, SIGNAL (triggered ()),
                    this, SLOT (voronoi_3 ()));

  QObject::connect (this->actionImport3DTDS, SIGNAL (triggered ()),
                    this, SLOT (import_3DTDS ()));

  QObject::connect (this->actionQuit, SIGNAL (triggered ()),
                    qApp, SLOT (quit ()));

  QObject::connect (this->actionSubdivide, SIGNAL (triggered ()),
                    this, SLOT (subdivide ()));

  QObject::connect (this->actionCreate_cube, SIGNAL (triggered ()),
                    this, SLOT (create_cube ()));

  QObject::connect (this->actionCreate_mesh, SIGNAL (triggered ()),
                    this, SLOT (create_mesh ()));

  QObject::connect (this->actionCreate3Cubes, SIGNAL (triggered ()),
                    this, SLOT (create_3cubes ()));

  QObject::connect (this->actionCreate2Volumes, SIGNAL (triggered ()),
                    this, SLOT (create_2volumes ()));

  QObject::connect (this, SIGNAL (sceneChanged ()),
                    this, SLOT (onSceneChanged ()));

  QObject::connect (this->actionClear, SIGNAL (triggered ()),
                    this, SLOT (clear ()));

  QObject::connect (this->actionDual_3, SIGNAL (triggered ()),
                    this, SLOT (dual_3 ()));

  QObject::connect (this->actionClose_volume, SIGNAL (triggered ()),
                    this, SLOT (close_volume ()));

  QObject::connect (this->actionRemove_filled_volumes, SIGNAL (triggered ()),
                    this, SLOT (remove_filled_volumes ()));

  QObject::connect (this->actionRemove_selected_volume, SIGNAL (triggered ()),
                    this, SLOT (remove_selected_volume ()));

  QObject::connect (this->actionSew3_same_facets, SIGNAL (triggered ()),
                    this, SLOT (sew3_same_facets ()));

  QObject::connect (this->actionUnsew3_all, SIGNAL (triggered ()),
                    this, SLOT (unsew3_all ()));

  QObject::connect (this->actionTriangulate_all_facets, SIGNAL (triggered ()),
                    this, SLOT (triangulate_all_facets ()));

  QObject::connect (this->actionExtend_filled_volumes, SIGNAL (triggered ()),
                    this, SLOT (extendFilledVolumes()));
  QObject::connect (this->actionExtend_hidden_volumes, SIGNAL (triggered ()),
                    this, SLOT (extendHiddenVolumes()));

  QObject::connect(this->volumeList->horizontalHeader(),
                   SIGNAL(sectionClicked(int)),
                   this, SLOT(onHeaderClicked(int)));
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

void MainWindow::import_off ()
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

void MainWindow::import_3DTDS ()
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

void MainWindow::add_off ()
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

  if (clear)
    this->clear(false);

  std::ifstream ifs (qPrintable (fileName));

  CGAL::import_from_polyhedron_3_flux < LCC > (*scene.lcc, ifs);
  initAllNewVolumes();

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

void MainWindow::onNewVolume(Dart_handle adart)
{
  scene.lcc->set_attribute<3>(adart,scene.lcc->create_attribute<3>
                              (CGAL::Color(random.get_int(0,256),
                                           random.get_int(0,256),
                                           random.get_int(0,256))));
  update_volume_list_add(adart);
}

void MainWindow::initAllNewVolumes()
{
  for (LCC::One_dart_per_cell_range<3>::iterator 
	 it(scene.lcc->one_dart_per_cell<3>().begin());
       it.cont(); ++it)
    if ( it->attribute<3>()==NULL )
    { onNewVolume(it); }
}

void MainWindow::load_3DTDS (const QString & fileName, bool clear)
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

  if (clear)
    this->clear(false);

  typedef CGAL::Delaunay_triangulation_3 < LCC::Traits > Triangulation;
  Triangulation T;

  std::ifstream ifs (qPrintable (fileName));
  std::istream_iterator < Point_3 > begin (ifs), end;
  T.insert (begin, end);

  CGAL::import_from_triangulation_3 < LCC, Triangulation >(*scene.lcc, T);
  initAllNewVolumes();

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

void MainWindow::create_cube ()
{
  Point_3 basepoint(nbcube%5, (nbcube/5)%5, nbcube/25);
	
  Dart_handle d = make_iso_cuboid(basepoint, 1);
	
  onNewVolume(d);

  ++nbcube;

  statusBar ()->showMessage (QString ("Cube created"),DELAY_STATUSMSG);

  emit (sceneChanged ());
}

void MainWindow::create_3cubes ()
{
  Dart_handle d1 = make_iso_cuboid (Point_3 (nbcube, nbcube, nbcube),1);
  Dart_handle d2 = make_iso_cuboid (Point_3 (nbcube + 1, nbcube, nbcube),1);
  Dart_handle d3 = make_iso_cuboid (Point_3 (nbcube, nbcube + 1, nbcube), 1);

  onNewVolume(d1);
  onNewVolume(d2);
  onNewVolume(d3);

  scene.lcc->sew<3> (d1->beta(1)->beta(1)->beta(2), d2->beta(2));
  scene.lcc->sew<3> (d1->beta(2)->beta(1)->beta(1)->beta(2), d3);

  ++nbcube;

  statusBar ()->showMessage (QString ("3 cubes were created"),
                             DELAY_STATUSMSG);

  emit (sceneChanged ());
}

void MainWindow::create_2volumes ()
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

  onNewVolume(d1);
  onNewVolume(d4);

  ++nbcube;
  statusBar ()->showMessage (QString ("2 volumes were created"),
                             DELAY_STATUSMSG);

  emit (sceneChanged ());
}

void MainWindow::create_mesh ()
{
  if ( dialogmesh.exec()==QDialog::Accepted )
  {
    for (int x=0; x<dialogmesh.getX(); ++x)
      for (int y=0; y<dialogmesh.getY(); ++y)
        for (int z=0; z<dialogmesh.getZ(); ++z)
        {
          Dart_handle d = make_iso_cuboid
            (Point_3 (x+nbcube, y+nbcube, z+nbcube), 1);
          onNewVolume(d);
        }
    ++nbcube;
      
    statusBar ()->showMessage (QString ("mesh created"),DELAY_STATUSMSG);
      
    emit (sceneChanged ());
  }
}

void MainWindow::subdivide ()
{
  subdivide_lcc_3 (*(scene.lcc));
  emit (sceneChanged ());
  statusBar ()->showMessage (QString ("Objects were subdivided"),
			     DELAY_STATUSMSG);
}

void MainWindow::clear (bool msg)
{
  scene.lcc->clear ();
  volumeUid = 1;
  nbcube=0;
  if (msg)
  {
    statusBar ()->showMessage (QString ("Scene was cleared"), DELAY_STATUSMSG);
    emit (sceneChanged ());
  }
  
  volumeDartIndex.clear();
  volumeProperties.clear();
  volumeList->clearContents();
  volumeList->setRowCount(0);
}

void MainWindow::voronoi_3 ()
{
  QString fileName = QFileDialog::getOpenFileName (this,
                                                   tr ("Voronoi 3D"),
                                                   ".",
                                                   tr ("Data file (*)"));

  if (fileName.isEmpty ()) return;
  
  this->clear(false);
  typedef CGAL::Delaunay_triangulation_3 < LCC::Traits > Triangulation;
  Triangulation T;

  LCC delaunay_lcc;
  Dart_handle dh;
  
  std::ifstream ifs (qPrintable (fileName));
  std::istream_iterator < Point_3 > begin (ifs), end;
  T.insert (begin, end);

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
  
  initAllNewVolumes();
  emit (sceneChanged ());
  statusBar ()->showMessage (QString ("Voronoi 3D of points in ") + fileName,
                             DELAY_STATUSMSG);
}

void MainWindow::dual_3 ()
{
  if ( !scene.lcc->is_without_boundary(3) )
  {
    statusBar()->showMessage
      (QString ("Dual impossible: the lcc has some 3-boundary"), 
       DELAY_STATUSMSG);
    return;
  }

  LCC* duallcc = new LCC;
  scene.lcc->dual_points_at_barycenter(*duallcc);

  this->clear(false);
  delete scene.lcc;
  scene.lcc = duallcc;
  this->viewer->setScene(&scene);
  initAllNewVolumes();
  
  statusBar ()->showMessage (QString ("Dual_3 computed"), DELAY_STATUSMSG);
  emit (sceneChanged ());
}

void MainWindow::close_volume()
{
  if ( scene.lcc->close(3) > 0 )
  {  
    initAllNewVolumes();
    statusBar ()->showMessage (QString ("Volume are closed"), DELAY_STATUSMSG);
    emit (sceneChanged ());
  }
  else
    statusBar ()->showMessage
      (QString ("LCC already 3-closed"), DELAY_STATUSMSG);
}

void MainWindow::sew3_same_facets()
{
  //  timer.reset();
  //  timer.start();
  if ( scene.lcc->sew3_same_facets() > 0 )
  {
    statusBar()->showMessage
      (QString ("Same facets are 3-sewn"), DELAY_STATUSMSG);
    emit (sceneChanged ());
  }
  else
    statusBar()->showMessage (QString ("No facets 3-sewn"), DELAY_STATUSMSG);
  //  timer.stop();
  //  std::cout<<"sew3_same_facets in "<<timer.time()<<" seconds."<<std::endl;
}

void MainWindow::unsew3_all()
{
  unsigned int nb=0;

  for (LCC::Dart_range::iterator it=scene.lcc->darts().begin();
       it!=scene.lcc->darts().end(); ++it)
  {
    if ( !it->is_free(3) )
    { scene.lcc->unsew<3>(it); ++nb; }
  }

  if ( nb > 0 )
  {
    statusBar()->showMessage
      (QString ("All darts are 3-unsewn"), DELAY_STATUSMSG);
    emit (sceneChanged ());
  }
  else
    statusBar()->showMessage (QString ("No dart 3-unsewn"), DELAY_STATUSMSG);
}

void MainWindow::remove_filled_volumes()
{
  unsigned int count = 0;
  if(volumeDartIndex.size() > 0)
  {
    for(unsigned int i = 0; i < volumeDartIndex.size();)
    {
      if( isVisibleAndFilled(volumeProperties[i]) )
      {
        CGAL::remove_cell<LCC,3>(*scene.lcc,volumeDartIndex[i].second);
        update_volume_list_remove(i);
        ++count;
      }
      else
      {
        i++;
      }
    }
    emit(sceneChanged());
  }
  statusBar()->showMessage
    (QString::number(count)+QString(" visible volume(s) removed"),
     DELAY_STATUSMSG);
}

void MainWindow::remove_selected_volume()
{
  bool nothingSelected = true;
  int row = 0;
  for(; row < volumeList->rowCount(); row++)
  {
    if(volumeList->item(row,0)->isSelected())
    {
      nothingSelected = false;
      break;
    }
  }

  if(nothingSelected)
    statusBar()->showMessage (QString("Nothing Selected"), DELAY_STATUSMSG);
  else
  {
    if(::isVisible(volumeProperties[row]))
    {
      CGAL::remove_cell<LCC,3>(*scene.lcc,volumeDartIndex[row].second);
      update_volume_list_remove(row);
      emit(sceneChanged());
    }
    else
    {
      statusBar()->showMessage (QString("Volume is hidden"), DELAY_STATUSMSG);
    }
  }
}

void MainWindow::triangulate_all_facets()
{
  std::vector<LCC::Dart_handle> v;
  for (LCC::One_dart_per_cell_range<2>::iterator 
	 it(scene.lcc->one_dart_per_cell<2>().begin()); it.cont(); ++it)
  {
    v.push_back(it);
  }
  for (std::vector<LCC::Dart_handle>::iterator itv(v.begin());
       itv!=v.end(); ++itv)
    scene.lcc->insert_barycenter_in_cell<2>(*itv);

  emit (sceneChanged ());
  statusBar()->showMessage
    (QString ("All facets were triangulated"), DELAY_STATUSMSG);    
}

void MainWindow::update_volume_list_add(Dart_handle it)
{
  volumeDartIndex.push_back(std::pair<int, Dart_handle>(volumeUid,it));
  volumeProperties.push_back(setVisibleAndFilled(0));
  int newRow = volumeList->rowCount();
  volumeList->setRowCount(newRow+1);

  volumeList->disconnect(this);
    
  QTableWidgetItem* volumeLabel =
    new QTableWidgetItem(QString::number(volumeUid));
  volumeLabel->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
  volumeLabel->setTextAlignment(Qt::AlignRight|Qt::AlignVCenter);
  volumeList->setItem(newRow,0,volumeLabel);
    
  QTableWidgetItem* fillCB = new QTableWidgetItem;
  fillCB->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
  fillCB->setCheckState(Qt::Checked);
  volumeList->setItem(newRow,1, fillCB);
    
  QTableWidgetItem* hiddenCB = new QTableWidgetItem();
  hiddenCB->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
  hiddenCB->setCheckState(Qt::Unchecked);    
  volumeList->setItem(newRow,2,hiddenCB);
    
  volumeUid++;

  connectVolumeListHandlers();
                     
}

void MainWindow::connectVolumeListHandlers()
{
  QObject::connect(this->volumeList, SIGNAL(cellChanged(int,int)),
                   this, SLOT(onCellChanged(int,int)));
  QObject::connect(this->volumeList, SIGNAL(itemSelectionChanged()),
                   this, SLOT(onItemSelectionChanged()));


}

void MainWindow::update_volume_list_remove(int i)
{
  volumeDartIndex.erase(volumeDartIndex.begin()+i);
  volumeProperties.erase(volumeProperties.begin()+i);
  volumeList->removeRow(i);

  this->viewer->setSelectedVolumeIndex(-1);
  if(volumeList->rowCount() > i)
    volumeList->item(i,0)->setSelected(false);

}

void MainWindow::onCellChanged(int row, int col)
{

  volumeProperties[row] = 
    (volumeList->item(row,2)->checkState() == Qt::Unchecked ? ::setVisible(0) : 0) |
    (volumeList->item(row,1)->checkState() == Qt::Checked ? ::setFilled(0) : 0 );

  // Change selection when toggling any checkbox?
  // volumeList->item(row,0)->setSelected(true);

  emit(sceneChanged());
}

void MainWindow::onItemSelectionChanged()
{
  // This gives a QList index out of range assert fail if volumeUid is selected
  // and some other volume checkbox is clicked 
  // QItemSelectionModel* selectionModel = volumeList->selectionModel();
  // QModelIndexList selectedRow = selectionModel->selectedIndexes();

  int row = 0;
  for(; row < volumeList->rowCount(); row++)
  {
    if(volumeList->item(row,0)->isSelected())
      break;
  }

  this->viewer->setSelectedVolumeIndex(row);
  emit(sceneChanged());
  // statusBar()->showMessage (QString::number(row)+QString(" is selected"), DELAY_STATUSMSG);
}

void MainWindow::onHeaderClicked(int col)
{
  if(col != 0)
  {
    volumeList->disconnect(this);

    for(unsigned int i = 0; i < volumeProperties.size(); i++)
    {
      switch(qApp->keyboardModifiers())
      {
      case(Qt::ShiftModifier):
        volumeProperties[i] =
          (col == 1 ? ::setWireframe(volumeProperties[i]) :
           ::setVisible(volumeProperties[i]));          
        volumeList->item(i,col)->setCheckState( Qt::Unchecked);
        break;
      case(Qt::ControlModifier):
        volumeProperties[i] =
          ( col == 1 ? ::negateFilled(volumeProperties[i]) :
            ::negateVisible(volumeProperties[i]));
        volumeList->item(i,col)->
          setCheckState(volumeList->item(i,col)->checkState() ?
                        Qt::Unchecked: Qt::Checked);
        break;
      default:
        volumeProperties[i] =
          (col == 1 ? ::setFilled(volumeProperties[i]) :
           ::setHidden(volumeProperties[i]));          
        volumeList->item(i,col)->setCheckState(Qt::Checked);
        break;
      }
    }

    connectVolumeListHandlers();
    emit(sceneChanged());
  }
}

void MainWindow::extendVolumesSatisfying(char amask, char negatemask)
{
  bool changed = false;
  
  volumeList->disconnect(this);

  int mark_volume = scene.lcc->get_new_mark();
  
  for(unsigned int i = 0; i < volumeProperties.size(); i++)
  {
    if ( ((volumeProperties[i] & amask) == amask) &&
         ((volumeProperties[i] & negatemask) ==0 ) )
    {
      for (LCC::Dart_of_cell_range<3>::iterator
             it=scene.lcc->darts_of_cell<3>(volumeDartIndex[i].second).begin(),
             itend=scene.lcc->darts_of_cell<3>(volumeDartIndex[i].second).end();
           it!=itend; ++it )
      {
        scene.lcc->mark(it, mark_volume);
        if ( !it->is_free(3) &&
             !scene.lcc->is_marked( it->beta(3), mark_volume) )
        {
          CGAL::mark_cell<LCC,3>(*scene.lcc, it->beta(3), mark_volume);
          changed = true;
        }
      }
    }
  }

  for(unsigned int i = 0; i < volumeProperties.size(); i++)
  {
    if ( scene.lcc->is_marked(volumeDartIndex[i].second, mark_volume) )
    {
      volumeProperties[i] |= amask;
      volumeProperties[i] ^= (volumeProperties[i] & negatemask);
      CGAL::unmark_cell<LCC,3>(*scene.lcc, volumeDartIndex[i].second,
                               mark_volume);
      
      volumeList->item(i,1)->setCheckState
        ( (::isFilled(volumeProperties[i])? Qt::Checked : Qt::Unchecked) );
      volumeList->item(i,2)->setCheckState
        ( (::isVisible(volumeProperties[i])? Qt::Unchecked: Qt::Checked) );      
    }
  }

  CGAL_assertion( scene.lcc->is_whole_map_unmarked(mark_volume) );
  scene.lcc->free_mark(mark_volume);
  
  connectVolumeListHandlers();
  if ( changed ) emit(sceneChanged());
}

void MainWindow::extendFilledVolumes()
{ extendVolumesSatisfying( LCC_DEMO_VISIBLE | LCC_DEMO_FILLED, 0 ); }
void MainWindow::extendHiddenVolumes()
{ extendVolumesSatisfying( 0, LCC_DEMO_VISIBLE ); }

#undef DELAY_STATUSMSG

#include "MainWindow.moc"
