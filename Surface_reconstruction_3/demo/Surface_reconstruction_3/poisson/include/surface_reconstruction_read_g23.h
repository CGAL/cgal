#ifndef SURFACE_RECONSTRUCTION_READ_G23_H
#define SURFACE_RECONSTRUCTION_READ_G23_H

#include <CGAL/basic.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/surface_reconstruction_assertions.h>

#include <stdio.h>
#include <vector>
#include <string>


/// .g23 format:
///
/// G23
/// movie_file_name
/// number_of_cameras number_of_3D_points number_of_2D_points
/// camera_index Cx Cy Cz                 for each camera
/// X Y Z                                 for each 3D point
/// camera_index point_3D_index x y       for each 2D point
/// with:
/// - camera_index is the camera's index (aka frame's index).
///   Camera indices may not be consecutive.
/// - (Cx, Cy, Cz) is the camera's position.
/// - point_3D_index is the 3D point's index.
///   3D point indices must start at 1 and be consecutive.
/// - (X,Y,Z) is 3D the point's position.
/// - (x,y) is 2D the point's position.


/// Read 3D points + cameras from a Gyroviz .g23 file.
///
/// @heading Parameters:
/// @param GyrovizPointOutputIterator value_type must be Gyroviz_point_3.
/// @param PointOutputIterator value_type must be Point_3.
///
/// @return true on success.
template <typename GyrovizPointOutputIterator,
          typename PointOutputIterator>
bool surface_reconstruction_read_g23(const char* pFilename, 
                                     GyrovizPointOutputIterator gyroviz_point_output,
                                     PointOutputIterator camera_output,
                                     std::string* movie_file_name)
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  typedef typename CGAL::value_type_traits<GyrovizPointOutputIterator>::type 
                                              Gyroviz_point;

  typedef typename Gyroviz_point::Geom_traits Geom_traits;
  typedef typename Geom_traits::Point_3       Point_3;
  typedef typename Geom_traits::Point_2       Point_2;

  CGAL_precondition(pFilename != NULL);
  CGAL_precondition(movie_file_name != NULL);

  FILE *pFile = fopen(pFilename,"rt");
  if(pFile == NULL)
  {
    std::cerr << "Error: cannot open " << pFilename << std::endl;
    return false;
  }

  long positions_2D_count = -1, positions_3D_count = -1, cameras_count = -1; // number of points and cameras in the file
  int lineNumber = 0; // current line number
  char pLine[4096]; // current line buffer
  std::map<int, Point_3> cameras; // container of (indexed) cameras
  std::vector<Gyroviz_point> gyroviz_points; // container of 3D points + (camera, 2D point) pairs
  
  *movie_file_name = "";
  
  while(fgets(pLine,sizeof(pLine),pFile) != NULL)
  {
    lineNumber++;
    
    // Skip blank line or comment
    if (strlen(pLine) == 0)
      continue;
    if (pLine[0] == '\r' || pLine[0] == '\n')
      continue;
    if (pLine[0] == '#')
      continue;

    // Read file signature on first line
    if (lineNumber == 1)
    {
      char signature[512];
      if ( (sscanf(pLine,"%s",signature) != 1) || (strcmp(signature, "G23") != 0) )
      {
        // if unsupported file format
        std::cerr << "Incorrect file format line " << lineNumber << " of " << pFilename << std::endl;
        return false;
      }
    }    

    // Read movie file name on 2nd (significant) line    
    else if (movie_file_name->size() == 0)
    {
      char file_name[512];
      if (sscanf(pLine,"%s",file_name) != 1)
      {
        std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
        return false;
      }
      *movie_file_name = file_name;
    }       

    // Read number of cameras and points on 3rd (significant) line    
    else if (cameras_count == -1)
    {
      if (sscanf(pLine,"%ld %ld %ld",&cameras_count,&positions_3D_count,&positions_2D_count) != 3)
      {
        std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
        return false;
      }
    }       

    // Read cameras on next lines
    else if (cameras.size() < cameras_count)
    {
      // Read index + position...
      int camera_index;
      double Cx,Cy,Cz;   
      if(sscanf(pLine, "%d %lg %lg %lg", &camera_index, &Cx,&Cy,&Cz) == 4)
      {
        Point_3 camera(Cx,Cy,Cz);
        cameras[camera_index] = camera;
      }
    }

    // Read 3D points on next lines
    else if (gyroviz_points.size() < positions_3D_count)
    {
      // Read position...
      double X,Y,Z;   
      if(sscanf(pLine, "%lg %lg %lg", &X,&Y,&Z) == 3)
      {
        Point_3 position_3D(X,Y,Z);
        gyroviz_points.push_back(position_3D);
      }
    }

    // Read 2D points + camera indices + 3D points indices on remaining lines
    else 
    {
      // Read camera index, 3D point index and 2D position...
      int camera_index, point_3D_index;
      double x,y;   
      if(sscanf(pLine, "%d %d %lg %lg", &camera_index, &point_3D_index, &x,&y) == 4)
      {
        if (cameras.find(camera_index) == cameras.end() || point_3D_index-1 >= gyroviz_points.size())
        {
          std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
          return false;
        }
        Point_3 camera = cameras[camera_index];
        Point_2 position_2D(x, y);
        gyroviz_points[point_3D_index-1].add_camera_point2_pair( std::make_pair(camera, position_2D) );
      }
    }
  }

  fclose(pFile);

  // Copy gyroviz_points[] to gyroviz_point_output
  //std::copy(gyroviz_points.begin(), gyroviz_points.end(), gyroviz_point_output);
  for (std::vector<Gyroviz_point>::iterator it=gyroviz_points.begin(); it != gyroviz_points.end(); it++)
  {
    // TEMPORARY? Skip 3D points with no cameras.
    if (it->cameras_begin() != it->cameras_end())
      *gyroviz_point_output++ = *it;
    //else
    //  std::cerr << "Skip (" << (Point_3)*it << ") line " << lineNumber << " of " << pFilename << std::endl;
  }
  
  // Copy cameras[] to camera_output
  for (std::map<int,Point_3>::iterator it=cameras.begin(); it != cameras.end(); it++)
    *camera_output++ = it->second;
  
  return true;
}


#endif // SURFACE_RECONSTRUCTION_READ_G23_H
