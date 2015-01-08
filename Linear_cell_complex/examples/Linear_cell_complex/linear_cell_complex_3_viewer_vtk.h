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
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_converter.h>

typedef CGAL::Cartesian<double> Local_kernel;
typedef typename Local_kernel::Point_3  Local_point;
typedef typename Local_kernel::Vector_3 Local_vector;

template<class LCC, int dim=LCC::ambient_dimension>
struct Geom_utils;

template<class LCC>
struct Geom_utils<LCC,3>
{
  Local_point get_point(LCC& lcc, typename LCC::Vertex_attribute_const_handle vh)
  { return converter(lcc.point_of_vertex_attribute(vh)); }

  Local_point get_point(LCC& lcc, typename LCC::Dart_const_handle dh)
  { return converter(lcc.point(dh)); }
  
  Local_vector get_facet_normal(LCC& lcc, typename LCC::Dart_const_handle dh)
  {
    Local_vector n = converter(CGAL::compute_normal_of_cell_2<LCC>(lcc,dh));
    n = n/(CGAL::sqrt(n*n));
    return n;
  }

  Local_vector get_vertex_normal(LCC& lcc, typename LCC::Dart_const_handle dh)
  {
    Local_vector n = converter(CGAL::compute_normal_of_cell_0<LCC>(lcc,dh));
    n = n/(CGAL::sqrt(n*n));
    return n;
  }
protected:
  CGAL::Cartesian_converter<typename LCC::Traits, Local_kernel> converter;
};

template<class LCC>
struct Geom_utils<LCC,2>
{
  Local_point get_point(LCC& lcc, typename LCC::Vertex_attribute_const_handle vh)
  {
    Local_point p(converter(lcc.point_of_vertex_attribute(vh).x()),0,
                  converter(lcc.point_of_vertex_attribute(vh).y()));
    return p;
  }

  Local_point get_point(LCC& lcc, typename LCC::Dart_const_handle dh)
  { return get_point(lcc, lcc.vertex_attribute(dh)); }

  Local_vector get_facet_normal(LCC&, typename LCC::Dart_const_handle)
  {
    Local_vector n(0,1,0);
    return n;
  }

  Local_vector get_vertex_normal(LCC&, typename LCC::Dart_const_handle)
  {
    Local_vector n(0,1,0);
    return n;    
  }
protected:
  CGAL::Cartesian_converter<typename LCC::Traits, Local_kernel> converter;  
};

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
  typedef typename LCC::Dart_handle Dart_handle;
  Geom_utils<LCC> geomutils;

public:
  SimpleLCCViewerVtk(LCC& lcc) : SimpleViewVtk()
  {
    setWindowTitle("3D lcc viewer");	

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
          if ( lcc.dimension>=3 && !lcc.is_free(it2, 3) )
            lcc.mark(lcc.beta(it2, 3), facettreated);
        }

        polygons->InsertNextCell(nb);
        for (typename LCC::template Dart_of_orbit_range<1>::iterator
               it2=lcc.template darts_of_orbit<1>(it).begin();
             it2.cont(); ++it2)
        {
          Local_point p =  geomutils.get_point(lcc, it2);
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

    CGAL_assertion(lcc.is_whole_map_marked(vertextreated));
    CGAL_assertion(lcc.is_whole_map_marked(facettreated));
  
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
