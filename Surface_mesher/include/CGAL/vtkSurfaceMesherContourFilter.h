// Copyright (c) 2008  GeometryFactory, Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_VTK_SURFACE_MESHER_CONTOUR_FILTER_H
#define CGAL_VTK_SURFACE_MESHER_CONTOUR_FILTER_H

#include <CGAL/license/Surface_mesher.h>


#ifdef CGAL_USE_VTK

#include <CGAL/config.h>
#include <vtkPolyDataAlgorithm.h>

class vtkCGALSurfaceMesherContourFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkCGALSurfaceMesherContourFilter *New();
  vtkTypeMacro(vtkCGALSurfaceMesherContourFilter,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Methods to set contour values
  vtkSetMacro(Value,double);
  vtkGetMacro(Value,double);

protected:
  vtkCGALSurfaceMesherContourFilter();
  ~vtkCGALSurfaceMesherContourFilter();

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  double Value;

private:
  vtkCGALSurfaceMesherContourFilter(const vtkCGALSurfaceMesherContourFilter&);  // Not implemented.
  void operator=(const vtkCGALSurfaceMesherContourFilter&);  // Not implemented.
};
  
// IMPLEMENTATION

#include "vtkCellArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"

#ifndef CGAL_USE_VTK
#  define CGAL_USE_VTK
#endif

#include <CGAL/Image_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_to_vtk.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>

vtkStandardNewMacro(vtkCGALSurfaceMesherContourFilter);

vtkCGALSurfaceMesherContourFilter::vtkCGALSurfaceMesherContourFilter()
{
  Value = 0.;
}
  
vtkCGALSurfaceMesherContourFilter::~vtkCGALSurfaceMesherContourFilter()
{
}

int 
vtkCGALSurfaceMesherContourFilter::
FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

void vtkCGALSurfaceMesherContourFilter::PrintSelf(ostream& os,
                                                  vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Value: " << this->Value << "\n";
}

int vtkCGALSurfaceMesherContourFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and ouptut
  vtkImageData *inData = vtkImageData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

  // c2t3
  typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

  typedef Tr::Geom_traits GT;
  typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
  typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  CGAL::Image_3 image;
  if(!image.read_vtk_image_data(inData))
    return 0;
  Gray_level_image gray_level_image(image, Value);

  GT::FT radius = std::max(image.xdim() * image.vx(),
                           std::max(image.ydim() * image.vy(),
                                    image.zdim() * image.vz())
                           );
  GT::Sphere_3 bounding_sphere(GT::Point_3(image.xdim() * image.vx()/2.,
                                           image.ydim() * image.vy()/2.,
                                           image.zdim() * image.vz()/2.),
                               radius*radius);
  // definition of the surface, with 10^-2 as relative precision
  Surface_3 surface(gray_level_image, bounding_sphere, 1e-5);
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,
                                                     radius/50.,
                                                     radius/500.);
  // meshing surface, with the "manifold without boundary" algorithm
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());

  CGAL::output_c2t3_to_vtk_polydata(c2t3, output);
  output->Squeeze();

  return 1;
}
#endif

#endif // CGAL_VTK_SURFACE_MESHER_CONTOUR_FILTER_H
