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
// $URL: svn+ssh://gdamiand@scm.gforge.inria.fr/svn/cgal/branches/features/Linear_cell_complex-gdamiand/Linear_cell_complex/demo/Linear_cell_complex/MainWindow.cpp $
// $Id: MainWindow.cpp 65446 2011-09-20 16:55:42Z gdamiand $
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#include "MainWindow.h"
#include <CGAL/Delaunay_triangulation_3.h>

// Function defined in map_3_subivision.cpp
void subdivide_map_3 (Map & m);

#define DELAY_STATUSMSG 1500

MainWindow::MainWindow (QWidget * parent):CGAL::Qt::DemosMainWindow (parent),
					  nbcube (0),
					  tdsdart(NULL),
					  dialogmesh(this)
{
  setupUi (this);
  scene.map = new Map;
  
  this->viewer->setScene (&scene);
  connectActions ();
  this->addAboutDemo (":/cgal/help/about_Combinatorial_map_3.html");
  this->addAboutCGAL ();

  this->addRecentFiles (this->menuFile, this->actionQuit);
  connect (this, SIGNAL (openRecentFile (QString)),
           this, SLOT (load_off (QString)));

  statusMessage = new QLabel ("Darts: 0,  Vertices: 0  (Points: 0),  Edges: 0, Facets: 0,"
               " Volume: 0 (Vol color: 0),  Connected components: 0");
  statusBar ()->addWidget (statusMessage);
}


void MainWindow::connectActions ()
{
  QObject::connect (this->actionImportOFF, SIGNAL (triggered ()),
          this, SLOT (import_off ()));

  QObject::connect (this->actionAddOFF, SIGNAL (triggered ()),
          this, SLOT (add_off ()));

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

  QObject::connect (this->actionRemove_current_volume, SIGNAL (triggered ()),
          this, SLOT (remove_current_volume ()));

  QObject::connect (this->actionSew3_same_facets, SIGNAL (triggered ()),
          this, SLOT (sew3_same_facets ()));

  QObject::connect (this->actionUnsew3_all, SIGNAL (triggered ()),
          this, SLOT (unsew3_all ()));

  QObject::connect (this->actionTriangulate_all_facets, SIGNAL (triggered ()),
          this, SLOT (triangulate_all_facets ()));
}

void MainWindow::onSceneChanged ()
{
  int mark = scene.map->get_new_mark ();
  scene.map->negate_mark (mark);

  std::vector<unsigned int> cells;
  cells.push_back(0);
  cells.push_back(1);
  cells.push_back(2);
  cells.push_back(3);
  cells.push_back(4);
  
  std::vector<unsigned int> res = scene.map->count_cells (cells);

  std::ostringstream os;
  os << "Darts: " << scene.map->number_of_darts ()
    << ",  Vertices:" << res[0]
     <<",  (Points:"<<scene.map->number_of_attributes<0>()<<")"
    << ",  Edges:" << res[1]
    << ",  Facets:" << res[2]
    << ",  Volumes:" << res[3]
#ifdef COLOR_VOLUME
     <<",  (Vol color:"<<scene.map->number_of_attributes<3>()<<")"
#endif
    << ",  Connected components:" << res[4]
     <<",  Valid:"<<(scene.map->is_valid()?"true":"FALSE");

  scene.map->negate_mark (mark);
  scene.map->free_mark (mark);

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
    scene.map->clear ();

  std::ifstream ifs (qPrintable (fileName));

  CGAL::import_from_polyhedron_flux < Map > (*scene.map, ifs);
  initAllVolumesRandomColor();

  this->addToRecentFiles (fileName);
  QApplication::restoreOverrideCursor ();

  if (clear)
    statusBar ()->showMessage (QString ("Load off file") + fileName,
                DELAY_STATUSMSG);
  else
    statusBar ()->showMessage (QString ("Add off file") + fileName,
                DELAY_STATUSMSG);
  tdsdart = NULL;

  emit (sceneChanged ());
}

void MainWindow::initVolumeRandomColor(Dart_handle adart)
{
#ifdef COLOR_VOLUME
  scene.map->set_attribute<3>(adart,scene.map->create_attribute<3>(CGAL::Color(random.get_int(0,256),
									       random.get_int(0,256),
									       random.get_int(0,256))));
#endif
}

void MainWindow::initAllVolumesRandomColor()
{
#ifdef COLOR_VOLUME
  for (Map::One_dart_per_cell_range<3>::iterator 
	 it(scene.map->one_dart_per_cell<3>().begin());
       it.cont(); ++it)
    if ( it->attribute<3>()==NULL ) initVolumeRandomColor(it);
#endif
}

void MainWindow::load_3DTDS (const QString & fileName, bool clear)
{
  QApplication::setOverrideCursor (Qt::WaitCursor);

  if (clear)
    scene.map->clear ();

  typedef CGAL::Delaunay_triangulation_3 < Map::Traits > Triangulation;
  Triangulation T;

  std::ifstream ifs (qPrintable (fileName));
  std::istream_iterator < Point_3 > begin (ifs), end;
  T.insert (begin, end);

  tdsdart = CGAL::import_from_triangulation_3 < Map, Triangulation > (*scene.map, T);
  initAllVolumesRandomColor();

  QApplication::restoreOverrideCursor ();
  emit (sceneChanged ());
}

Dart_handle MainWindow::make_iso_cuboid(const Point_3 basepoint, Map::FT lg)
{
  return make_hexahedron(*scene.map,
												 basepoint,
												 Map::Construct_translated_point()(basepoint,Map::Vector(lg,0,0)),
												 Map::Construct_translated_point()(basepoint,Map::Vector(lg,lg,0)),
												 Map::Construct_translated_point()(basepoint,Map::Vector(0,lg,0)),
												 Map::Construct_translated_point()(basepoint,Map::Vector(0,lg,lg)),
												 Map::Construct_translated_point()(basepoint,Map::Vector(0,0,lg)),
												 Map::Construct_translated_point()(basepoint,Map::Vector(lg,0,lg)),
												 Map::Construct_translated_point()(basepoint,Map::Vector(lg,lg,lg)));
}

void MainWindow::create_cube ()
{
	Point_3 basepoint(nbcube%5, (nbcube/5)%5, nbcube/25);
	
  Dart_handle d = make_iso_cuboid(basepoint, 1);
	
  //  scene.map->display_characteristics(std::cout)<<std::endl;

  initVolumeRandomColor(d);

  ++nbcube;

  tdsdart = NULL;
  statusBar ()->showMessage (QString ("Cube created"),DELAY_STATUSMSG);

  emit (sceneChanged ());
}

void MainWindow::create_3cubes ()
{
  Dart_handle d1 = make_iso_cuboid (Point_3 (nbcube, nbcube, nbcube),1);
  Dart_handle d2 = make_iso_cuboid (Point_3 (nbcube + 1, nbcube, nbcube),1);
  Dart_handle d3 = make_iso_cuboid (Point_3 (nbcube, nbcube + 1, nbcube), 1);

  initVolumeRandomColor(d1);
  initVolumeRandomColor(d2);
  initVolumeRandomColor(d3);

  scene.map->sew<3> (d1->beta(1)->beta(1)->beta(2), d2->beta(2));
  scene.map->sew<3> (d1->beta(2)->beta(1)->beta(1)->beta(2), d3);

  ++nbcube;

  tdsdart = NULL;
  statusBar ()->showMessage (QString ("3 cubes were created"),
              DELAY_STATUSMSG);

  emit (sceneChanged ());
}

void MainWindow::create_2volumes ()
{
  Dart_handle d1 = make_iso_cuboid (Point_3 (nbcube, nbcube, nbcube),1);
  Dart_handle d2 = make_iso_cuboid (Point_3 (nbcube + 1, nbcube, nbcube), 1);
  Dart_handle d3 = make_iso_cuboid (Point_3 (nbcube, nbcube + 1, nbcube), 1);
  Dart_handle d4 = make_iso_cuboid (Point_3 (nbcube + 1, nbcube + 1, nbcube), 1);
  
  initVolumeRandomColor(d1);
  initVolumeRandomColor(d2);
  initVolumeRandomColor(d3);
  initVolumeRandomColor(d4);

  scene.map->sew<3>(d1->beta(1)->beta(1)->beta(2), d2->beta (2));
  scene.map->sew<3>(d1->beta(2)->beta(1)->beta(1)->beta (2), d3);

  scene.map->sew<3>(d3->beta(1)->beta(1)->beta(2), d4->beta (2));
  scene.map->sew<3>(d2->beta(2)->beta(1)->beta(1)->beta (2), d4);

  /*  scene.map->display_characteristics(std::cout)
    <<" is_valid="<<scene.map->is_valid()<<std::endl;

  std::cout<<"AVANT"<<std::endl;
  scene.map->display_darts(std::cout)<<std::endl;
  std::cout<<" is_valid="<<scene.map->is_valid()<<std::endl;*/

  CGAL::remove_cell<Map,2>(*scene.map, d3);
  CGAL::remove_cell<Map,2>(*scene.map, d2->beta (2));

  /*  std::cout<<"APRES"<<std::endl;
  scene.map->display_darts(std::cout)<<std::endl;
  std::cout<<" is_valid="<<scene.map->is_valid()<<std::endl;
  scene.map->display_characteristics(std::cout);*/

  tdsdart = NULL;
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
	      Dart_handle d = make_iso_cuboid (Point_3 (x+nbcube, y+nbcube, z+nbcube), 1);
	      initVolumeRandomColor(d);
	    }
      ++nbcube;
      
      tdsdart = NULL;
      statusBar ()->showMessage (QString ("mesh created"),DELAY_STATUSMSG);
      
      emit (sceneChanged ());
    }
}

void MainWindow::subdivide ()
{
  subdivide_map_3 (*(scene.map));
  tdsdart = NULL;
  emit (sceneChanged ());
  statusBar ()->showMessage (QString ("Objects were subdivided"),
			     DELAY_STATUSMSG);
}

void MainWindow::clear ()
{
  scene.map->clear ();
  tdsdart = NULL;
  statusBar ()->showMessage (QString ("Scene was cleared"), DELAY_STATUSMSG);
  emit (sceneChanged ());
}

void MainWindow::dual_3 ()
{
  if ( !scene.map->is_without_boundary(3) )
    {
      statusBar()->showMessage (QString ("Dual impossible: the map has some 3-boundary"), 
				DELAY_STATUSMSG);
      return;
    }

  Map* dualmap = new Map;
  Dart_handle infinitevolume = CGAL::dual<Map>(*scene.map,*dualmap,tdsdart);

  if ( tdsdart!=NULL )
    CGAL::remove_cell<Map,3>(*dualmap,infinitevolume);

  delete scene.map;
  scene.map = dualmap;
  this->viewer->setScene (&scene);
  initAllVolumesRandomColor();
  
  statusBar ()->showMessage (QString ("Dual_3 computed"), DELAY_STATUSMSG);
  emit (sceneChanged ());
}

void MainWindow::close_volume()
{
  tdsdart = NULL;
  if ( scene.map->close(3) > 0 )
    {  
      initAllVolumesRandomColor();
      statusBar ()->showMessage (QString ("Volume are closed"), DELAY_STATUSMSG);
      emit (sceneChanged ());
    }
  else
    statusBar ()->showMessage (QString ("Map already 3-closed"), DELAY_STATUSMSG);
}

void MainWindow::sew3_same_facets()
{
  tdsdart = NULL;
  //  timer.reset();
  //  timer.start();
  if ( scene.map->sew3_same_facets() > 0 )
    {
      statusBar()->showMessage (QString ("Same facets are 3-sewn"), DELAY_STATUSMSG);
      emit (sceneChanged ());
    }
  else
    statusBar()->showMessage (QString ("No facets 3-sewn"), DELAY_STATUSMSG);
  //  timer.stop();
  //  std::cout<<"sew3_same_facets in "<<timer.time()<<" seconds."<<std::endl;
}

void MainWindow::unsew3_all()
{
  tdsdart = NULL;
  unsigned int nb=0;

  for (Map::Dart_range::iterator it=scene.map->darts().begin();
       it!=scene.map->darts().end(); ++it)
    {
      if ( !it->is_free(3) )
	{ scene.map->unsew<3>(it); ++nb; }
    }

  if ( nb > 0 )
    {
      statusBar()->showMessage (QString ("All darts are 3-unsewn"), DELAY_STATUSMSG);
      emit (sceneChanged ());
    }
  else
    statusBar()->showMessage (QString ("No dart 3-unsewn"), DELAY_STATUSMSG);
}

void MainWindow::remove_current_volume()
{
  if ( this->viewer->getCurrentDart()!=scene.map->darts().end() )
    {
      CGAL::remove_cell<Map,3>(*scene.map,this->viewer->getCurrentDart());
      emit (sceneChanged ());
      statusBar()->showMessage (QString ("Current volume removed"), DELAY_STATUSMSG);    
    }
  else
    statusBar()->showMessage (QString ("No volume removed"), DELAY_STATUSMSG);
}

void MainWindow::triangulate_all_facets()
{
  std::vector<Map::Dart_handle> v;
  for (Map::One_dart_per_cell_range<2>::iterator 
	 it(scene.map->one_dart_per_cell<2>().begin()); it.cont(); ++it)
    {
      v.push_back(it);
    }
  for (std::vector<Map::Dart_handle>::iterator itv(v.begin());
       itv!=v.end(); ++itv)
    CGAL::insert_center_cell_0_in_cell_2(*scene.map,*itv);

  emit (sceneChanged ());
  statusBar()->showMessage (QString ("All facets were triangulated"), DELAY_STATUSMSG);    
}

#undef DELAY_STATUSMSG

#include "MainWindow.moc"
