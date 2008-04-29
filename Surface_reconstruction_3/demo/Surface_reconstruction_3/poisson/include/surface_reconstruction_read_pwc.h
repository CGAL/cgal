#ifndef SURFACE_RECONSTRUCTION_READ_PWC_H
#define SURFACE_RECONSTRUCTION_READ_PWC_H

#include <CGAL/basic.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/surface_reconstruction_assertions.h>

#include <stdio.h>
#include <vector>


/// Read 3D points + cameras from a Gyroviz .pwc file.
/// @return true on success.
template <typename OutputIterator> ///< OutputIterator value_type must be Gyroviz_point_3
bool surface_reconstruction_read_pwc(const char* pFilename, OutputIterator output)
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  typedef typename CGAL::value_type_traits<OutputIterator>::type Gyroviz_point;
  
  typedef typename Gyroviz_point::Geom_traits Geom_traits;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Vector_3 Vector;

  CGAL_precondition(pFilename != NULL);

  FILE *pFile = fopen(pFilename,"rt");
  if(pFile == NULL)
  {
    std::cerr << "Error: cannot open " << pFilename;
    return false;
  }

  long pointsCount = 0, camCount = 0; // number of points and number of cameras in the file
  int lineNumber = 0; // current line number
  char pLine[4096]; // current line buffer
  std::vector<Point> list_of_camera_coordinates; // container of cameras read
  while (fgets(pLine,4096,pFile))
  {
    lineNumber++;

    // Read file signature on first line
    if (lineNumber == 1)
    {
      char signature[4096];
      if ( (sscanf(pLine,"%s",signature) != 1) || (strcmp(signature, "CamOFF") != 0) )
      {
        // if unsupported file format
        std::cerr << "Incorrect file format line " << lineNumber << " of " << pFilename;
        return false;
      }
    }    

    // Read number of cameras and points on 2nd line    
    else if (lineNumber == 2)
    {
      if (sscanf(pLine,"%ld %ld",&camCount,&pointsCount) != 2)
      {
        std::cerr << "Error line " << lineNumber << " of " << pFilename;
        return false;
      }
    }       

    // Read cameras on next lines
    else if (list_of_camera_coordinates.size() < camCount)
    {
      // Read position...
      double Cx,Cy,Cz;   
      if(sscanf(pLine,"%lg\t%lg\t%lg",&Cx,&Cy,&Cz) == 3)
      {
        Point cam_coord(Cx,Cy,Cz);
        list_of_camera_coordinates.push_back(cam_coord);
      }
      // ...or skip comment line
    }

    // Read 3D points + camera indices on next lines
    else 
    {
      // Read position + camera indices...
      std::istringstream iss(pLine);
      Point position;
      if (iss >> position)
      {
        int value;
        std::vector<Point> list_of_cameras;
        while (iss >> value)
        {
          Point cam_coord = list_of_camera_coordinates[value-1];
          list_of_cameras.push_back(cam_coord);
        } 

        // TEMPORARY? Skip 3D points with no cameras.
        if (list_of_cameras.begin() != list_of_cameras.end())
        {
          *output = Gyroviz_point(position, 
                                  CGAL::NULL_VECTOR,  
                                  list_of_cameras.begin(), list_of_cameras.end());
          output++;
        }
      }
      // ...or skip comment line
    }
  }

  fclose(pFile);
  return true;
}


#endif // SURFACE_RECONSTRUCTION_READ_PWC_H
