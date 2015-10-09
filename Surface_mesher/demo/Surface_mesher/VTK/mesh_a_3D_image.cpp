// This demo has been inspired from the file
// Examples/Medical/Cxx/Medical1.cxx in VTK-5.0.

// Adaptation to CGAL_ImageIO made by: Laurent Rineau

/*
  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
*/

#ifdef CGAL_USE_VTK

#include <QApplication>
#include <iostream>
#include <cstdlib>
#include <sstream>

#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkOutlineFilter.h>
#include <vtkCamera.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkContourFilter.h>

#include <QVTKWidget.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_to_vtk.h>
#include <fstream>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

int main(int argc, char** argv) {
  QApplication app(argc, argv);

  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // the 'function' is a 3D gray level image
  Gray_level_image image("../../../examples/Surface_mesher/data/skull_2.9.inr", 2.9);

  // Carefully choosen bounding sphere: the center must be inside the
  // surface defined by 'image' and the radius must be high enough so that
  // the sphere actually bounds the whole image.
  GT::Point_3 bounding_sphere_center(122., 102., 117.);
  GT::FT bounding_sphere_squared_radius = 200.*200.*2.;
  GT::Sphere_3 bounding_sphere(bounding_sphere_center,
                                   bounding_sphere_squared_radius);

  // definition of the surface, with 10^-2 as relative precision
  Surface_3 surface(image, bounding_sphere, 1e-5);

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,
                                                     5.,
                                                     1.);

  // meshing surface, with the "manifold without boundary" algorithm
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());

  QVTKWidget widget;
  widget.resize(256,256);

//   vtkImageData* vtk_image = CGAL::vtk_image_sharing_same_data_pointer(image);

  vtkRenderer *aRenderer = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(aRenderer);

  widget.SetRenderWindow(renWin);

//   vtkContourFilter *skinExtractor = vtkContourFilter::New();
//     skinExtractor->SetInput(vtk_image);
//     skinExtractor->SetValue(0, isovalue);
//     skinExtractor->SetComputeNormals(0);
  vtkPolyDataNormals *skinNormals = vtkPolyDataNormals::New();
//     skinNormals->SetInputConnection(skinExtractor->GetOutputPort());
  vtkPolyData* polydata = CGAL::output_c2t3_to_vtk_polydata(c2t3);
  skinNormals->SetInputData(polydata);
    skinNormals->SetFeatureAngle(60.0);
  vtkPolyDataMapper *skinMapper = vtkPolyDataMapper::New();
//     skinMapper->SetInputConnection(skinExtractor->GetOutputPort());
  skinMapper->SetInputData(polydata);
    skinMapper->ScalarVisibilityOff();
  vtkActor *skin = vtkActor::New();
    skin->SetMapper(skinMapper);

  // An outline provides context around the data.
  //
//   vtkOutlineFilter *outlineData = vtkOutlineFilter::New();
//     outlineData->SetInput(vtk_image);
//   vtkPolyDataMapper *mapOutline = vtkPolyDataMapper::New();
//     mapOutline->SetInputConnection(outlineData->GetOutputPort());
//   vtkActor *outline = vtkActor::New();
//     outline->SetMapper(mapOutline);
//     outline->GetProperty()->SetColor(0,0,0);

  // It is convenient to create an initial view of the data. The FocalPoint
  // and Position form a vector direction. Later on (ResetCamera() method)
  // this vector is used to position the camera to look at the data in
  // this direction.
  vtkCamera *aCamera = vtkCamera::New();
    aCamera->SetViewUp (0, 0, -1);
    aCamera->SetPosition (0, 1, 0);
    aCamera->SetFocalPoint (0, 0, 0);
    aCamera->ComputeViewPlaneNormal();

  // Actors are added to the renderer. An initial camera view is created.
  // The Dolly() method moves the camera towards the FocalPoint,
  // thereby enlarging the image.
//   aRenderer->AddActor(outline);
  aRenderer->AddActor(skin);
  aRenderer->SetActiveCamera(aCamera);
  aRenderer->ResetCamera ();
  aCamera->Dolly(1.5);

  // Set a background color for the renderer and set the size of the
  // render window (expressed in pixels).
  aRenderer->SetBackground(1,1,1);
  renWin->SetSize(640, 480);

  // Note that when camera movement occurs (as it does in the Dolly()
  // method), the clipping planes often need adjusting. Clipping planes
  // consist of two planes: near and far along the view direction. The 
  // near plane clips out objects in front of the plane; the far plane
  // clips out objects behind the plane. This way only what is drawn
  // between the planes is actually rendered.
  aRenderer->ResetCameraClippingRange ();

  // Initialize the event loop and then start it.
//   iren->Initialize();
//   iren->Start(); 

  // It is important to delete all objects created previously to prevent
  // memory leaks. In this case, since the program is on its way to
  // exiting, it is not so important. But in applications it is
  // essential.
//   vtk_image->Delete();
//   skinExtractor->Delete();
  skinNormals->Delete();
  skinMapper->Delete();
  skin->Delete();
//   outlineData->Delete();
//   mapOutline->Delete();
//   outline->Delete();
  aCamera->Delete();
//   iren->Delete();
  renWin->Delete();
  aRenderer->Delete();
  polydata->Delete();

  widget.show();

  app.exec();
  
  return 0;
}

#else // #ifdef CGAL_USE_VTK

#include <iostream>

int main()
{
  std::cerr << "That demo needs VTK support.\n";
  return 0;
}

#endif // CGAL_USE_VTK
