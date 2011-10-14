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
// $URL$
// $Id$
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_MAP_VIEWER_VTK_H
#define CGAL_MAP_VIEWER_VTK_H

#include <QApplication>
#include <QVBoxLayout>
#include <QMenuBar>
#include <QMainWindow>

#include <QVTKWidget.h>

#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <CGAL/Combinatorial_map_with_points_operations.h>

class SimpleViewVtk : public QMainWindow
{
public:

  // Constructor/Destructor
  SimpleViewVtk(QWidget* p = 0) : QMainWindow(p)
  {
    //setupUi(this);
    centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    vtkWidget = new QVTKWidget(centralWidget);

    vboxLayout = new QVBoxLayout(centralWidget);
    vboxLayout->addWidget(vtkWidget);

    QAction* a_fileExit = new QAction(tr("&Exit"), this);
    a_fileExit->setShortcut(tr("Ctrl+Q"));
    a_fileExit->setStatusTip(tr("Exit"));
    connect(a_fileExit, SIGNAL(triggered()), this, SLOT(close()));
    
    QMenu* file_menu = this->menuBar()->addMenu(tr("&File"));
    file_menu->addAction(a_fileExit);
    
    // QT/VTK interact
    ren = vtkRenderer::New();
    vtkWidget->GetRenderWindow()->AddRenderer(ren);

    resize(500, 450);
  }
  
  ~SimpleViewVtk() {};
		  
protected:
  QWidget     *centralWidget;
  QVBoxLayout *vboxLayout;
  QVTKWidget  *vtkWidget;
  
  vtkPolyDataMapper *mapper;
  vtkActor          *actor;
  vtkRenderer       *ren;   
};

template<class Map>
class SimpleMapViewerVtk : public SimpleViewVtk
{
public:
  SimpleMapViewerVtk(Map& m) : SimpleViewVtk()
  {
    setWindowTitle("3D map viewer");	

    typedef typename Map::Dart_handle Dart_handle;

    vtkPolyData *polydata = vtkPolyData::New();
  
    unsigned int facettreated = m.get_new_mark();
    unsigned int vertextreated = m.get_new_mark();

    vtkCellArray* polygons = vtkCellArray::New();
    vtkCellArray* vertices = vtkCellArray::New();  
    vtkPoints* points = vtkPoints::New();
    unsigned nbpoints=0;

    for(typename Map::Dart_iterator_of_all it(m); it.cont(); ++it)
      {
	if (!m.is_marked(*it,facettreated))
	  {
	    unsigned int nb=0;

	    for (typename Map::Dart_iterator_of_beta1 it2(m,*it); it2.cont(); ++it2)
	      {
		++nb;
		m.set_mark(*it2,facettreated);
		if ( !it2->is_free(3) ) m.set_mark(it2->beta(3),facettreated);
	      }

	    polygons->InsertNextCell(nb);
	    for (typename Map::Dart_iterator_of_beta1 it2(m,*it); it2.cont(); ++it2)
	      {
		typename Map::Point p = it2->vertex()->point();
		vtkIdType id=points->InsertNextPoint(p.x(),p.y(),p.z());
		++nbpoints;

		if ( !m.is_marked(*it2,vertextreated) )
		  {
		    vertices->InsertNextCell(1);
		    vertices->InsertCellPoint(id);

		    mark_orbit(m,*it2,Map::VERTEX_ORBIT,vertextreated);
		  }
	      
		polygons->InsertCellPoint(id);
	      }
	  }
	polydata->SetPoints(points);
	polydata->SetVerts(vertices);
	polydata->SetPolys(polygons);
      }

    assert(m.is_whole_map_marked(vertextreated));
    assert(m.is_whole_map_marked(facettreated));
  
    m.free_mark(vertextreated);
    m.free_mark(facettreated);  
  
  
    // Mapper
    mapper = vtkPolyDataMapper::New();
    mapper->ImmediateModeRenderingOn();
    mapper->SetInput(polydata);

    // Actor in scene
    actor = vtkActor::New();
    actor->SetMapper(mapper);

    // Add Actor to renderer
    ren->AddActor(actor);

    // Reset camera
    ren->ResetCamera();

    ren->GetRenderWindow()->Render();
  }
};

template<class Map>
void display_map(Map& amap)
{
  int argc=1;
  typedef char* s;
  
  const char* argv[2]={"mapviewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));
  
  SimpleMapViewerVtk<Map> mainwindow(amap);
  mainwindow.show();

  app.exec();
};

#endif // CGAL_MAP_VIEWER_VTK_H

