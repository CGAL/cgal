// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
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
#ifndef CGAL_LCC_3_VIEWER_VTK_H
#define CGAL_LCC_3_VIEWER_VTK_H

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
#include <CGAL/Linear_cell_complex.h>

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
  
protected:
  QWidget     *centralWidget;
  QVBoxLayout *vboxLayout;
  QVTKWidget  *vtkWidget;
  
  vtkPolyDataMapper *mapper;
  vtkActor          *actor;
  vtkRenderer       *ren;   
};

template<class LCC>
class SimpleLCCViewerVtk : public SimpleViewVtk
{
public:
  SimpleLCCViewerVtk(LCC& lcc) : SimpleViewVtk()
  {
    setWindowTitle("3D lcc viewer");	

    typedef typename LCC::Dart_handle Dart_handle;

    vtkPolyData *polydata = vtkPolyData::New();
  
    int facettreated = lcc.get_new_mark();
    int vertextreated = lcc.get_new_mark();

    vtkCellArray* polygons = vtkCellArray::New();
    vtkCellArray* vertices = vtkCellArray::New();  
    vtkPoints* points = vtkPoints::New();
    unsigned int nbpoints=0;

    for(typename LCC::Dart_range::iterator it=lcc.darts().begin(),
        itend=lcc.darts().end(); it!=itend; ++it)
    {
      if (!lcc.is_marked(it,facettreated))
      {
        unsigned int nb=0;

        for (typename LCC::template Dart_of_orbit_range<1>::iterator
               it2=lcc.template darts_of_orbit<1>(it).begin();
             it2.cont(); ++it2)
        {
          ++nb;
          lcc.mark(it2, facettreated);
          if ( !it2->is_free(3) ) lcc.mark(it2->beta(3), facettreated);
        }

        polygons->InsertNextCell(nb);
        for (typename LCC::template Dart_of_orbit_range<1>::iterator
               it2=lcc.template darts_of_orbit<1>(it).begin();
             it2.cont(); ++it2)
        {
          typename LCC::Point p =  lcc.vertex_attribute(it2)->point();
          vtkIdType id=points->InsertNextPoint(p.x(),p.y(),p.z());
          ++nbpoints;

          if ( !lcc.is_marked(it2,vertextreated) )
          {
            vertices->InsertNextCell(1);
            vertices->InsertCellPoint(id);
            
            CGAL::mark_cell<LCC, 0>(lcc, it2, vertextreated);
          }
	      
          polygons->InsertCellPoint(id);
        }
      }
      polydata->SetPoints(points);
      polydata->SetVerts(vertices);
      polydata->SetPolys(polygons);
    }

    assert(lcc.is_whole_map_marked(vertextreated));
    assert(lcc.is_whole_map_marked(facettreated));
  
    lcc.free_mark(vertextreated);
    lcc.free_mark(facettreated);  
  
  
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

template<class LCC>
void display_lcc(LCC& alcc)
{
  int argc=1;
  typedef char* s;
  
  const char* argv[2]={"lccviewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));
  
  SimpleLCCViewerVtk<LCC> mainwindow(alcc);
  mainwindow.show();

  app.exec();
};

#endif // CGAL_LCC_3_VIEWER_VTK_H
