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

#include <CGAL/Image_3.h>
#include <CGAL/Image_3_vtk_interface.h>

void usage_and_exit(char *argv0)
{
  std::cerr << "Usage:\n  "
            << argv0 << " <filename> <iso-value>\n";
  std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
  QApplication app(argc, argv);

  if(argc != 3)
    usage_and_exit(argv[0]);

  QVTKWidget widget;
  widget.resize(256,256);
 
  CGAL::Image_3 image;
  if(!image.read(argv[1]))
  {
    std::cerr << "Cannot read image file \"" << argv[1] << "\"!\n";
    usage_and_exit(argv[0]);
  }

  std::stringstream argv2;
  argv2 << argv[2];
  double isovalue;
  if(!(argv2 >> isovalue))
  {
    std::cerr << "Invalid iso-value \"" << argv[2] << "\"!\n";
    usage_and_exit(argv[0]);
  }

  vtkImageData* vtk_image = CGAL::vtk_image_sharing_same_data_pointer(image);

  vtkRenderer *aRenderer = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(aRenderer);

  widget.SetRenderWindow(renWin);

  vtkContourFilter *skinExtractor = vtkContourFilter::New();
    skinExtractor->SetInputData(vtk_image);
    skinExtractor->SetValue(0, isovalue);
//     skinExtractor->SetComputeNormals(0);
  vtkPolyDataNormals *skinNormals = vtkPolyDataNormals::New();
    skinNormals->SetInputConnection(skinExtractor->GetOutputPort());
    skinNormals->SetFeatureAngle(60.0);
  vtkPolyDataMapper *skinMapper = vtkPolyDataMapper::New();
    skinMapper->SetInputConnection(skinExtractor->GetOutputPort());
    skinMapper->ScalarVisibilityOff();
  vtkActor *skin = vtkActor::New();
    skin->SetMapper(skinMapper);

  // An outline provides context around the data.
  //
  vtkOutlineFilter *outlineData = vtkOutlineFilter::New();
    outlineData->SetInputData(vtk_image);
  vtkPolyDataMapper *mapOutline = vtkPolyDataMapper::New();
    mapOutline->SetInputConnection(outlineData->GetOutputPort());
  vtkActor *outline = vtkActor::New();
    outline->SetMapper(mapOutline);
    outline->GetProperty()->SetColor(0,0,0);

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
  aRenderer->AddActor(outline);
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
  vtk_image->Delete();
  skinExtractor->Delete();
  skinNormals->Delete();
  skinMapper->Delete();
  skin->Delete();
  outlineData->Delete();
  mapOutline->Delete();
  outline->Delete();
  aCamera->Delete();
//   iren->Delete();
  renWin->Delete();
  aRenderer->Delete();

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
