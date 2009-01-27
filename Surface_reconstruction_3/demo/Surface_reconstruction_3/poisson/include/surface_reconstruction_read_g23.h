#ifndef SURFACE_RECONSTRUCTION_READ_G23_H
#define SURFACE_RECONSTRUCTION_READ_G23_H

#include <CGAL/value_type_traits.h>

#include <stdio.h>
#include <vector>
#include <string>


// .g23 version 1 format:
// ----------------------
// G23
// movie_file_name
// number_of_cameras number_of_3D_points number_of_2D_points
// frame_index Cx Cy Cz                  for each camera
// point_3D_label X Y Z                  for each 3D point
// frame_index point_3D_label x y        for each 2D point
// with:
// - movie_file_name is the tracked movie's name, with path
//   (using % notation for image sequences, e.g. foo/bar%3d.tga).
// - frame_index is the frame's index (== number in image's name).
//   Camera indices may not be consecutive.
// - (Cx, Cy, Cz) is the camera's position.
// - point_3D_label is the 3D point's label (between double quotes).
// - (X,Y,Z) is the 3D point's position.
// - (x,y) is the 2D point's position.
//
// .g23 version 2 format:
// ----------------------
// G23 v2
// movie_file_name
// number_of_cameras number_of_3D_points number_of_2D_points
// frame_index M11 M12...M34             for each camera
// point_3D_label X Y Z                  for each 3D point
// frame_index point_3D_label x y        for each 2D point
// with:
// - movie_file_name is the tracked movie's name, with path
//   (using % notation for image sequences, e.g. foo/bar%3d.png).
// - frame_index is the frame's index (== number in image's name).
//   Camera indices may not be consecutive.
// - (M11 M12...M34) is the camera's 3x4 matrix.
// - point_3D_label is the 3D point's label (between double quotes).
// - (X,Y,Z) is the 3D point's position.
// - (x,y) is the 2D point's position.



/// Read 3D points + cameras from a Gyroviz .g23 file.
///
/// @heading Parameters:
/// @param GyrovizPointOutputIterator value_type must be Gyroviz_point_3.
///
/// @return true on success.
template <typename Point_3,
          typename GyrovizPointOutputIterator>
bool surface_reconstruction_read_g23(
  const char* pFilename, 
  GyrovizPointOutputIterator gyroviz_point_output,
  std::map<int, Point_3>* cameras, // container of (indexed) cameras
  std::string* movie_file_name)
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  typedef typename CGAL::value_type_traits<GyrovizPointOutputIterator>::type 
                                              Gyroviz_point;

  typedef typename Gyroviz_point::Geom_traits Geom_traits;
  //typedef typename Geom_traits::Point_3               Point_3;
  typedef typename Geom_traits::Vector_3              Vector_3;
  typedef typename Geom_traits::Point_2               Point_2;
  typedef typename Geom_traits::Aff_transformation_3  Aff_transformation_3;

  CGAL_precondition(pFilename != NULL);
  CGAL_precondition(movie_file_name != NULL);

  FILE *pFile = fopen(pFilename,"rt");
  if(pFile == NULL)
  {
    std::cerr << "Error: cannot open " << pFilename << std::endl;
    return false;
  }

  int version;
  long positions_2D_count = -1, positions_3D_count = -1, cameras_count = -1; // number of points and cameras in the file
  int lineNumber = 0; // current line number
  char pLine[4096]; // current line buffer
  std::map<std::string, Gyroviz_point> gyroviz_points; // container of (labelled) 3D points + camera/2D point pairs
  
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
      version = 1; // version == 1 if lacking in file
      if ( (sscanf(pLine,"%s v%d", signature, &version) == 0) 
        || (strcmp(signature, "G23") != 0) 
        || (version != 1 && version != 2) )
      {
        // if incorrect file format
        std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
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
    else if (cameras->size() < cameras_count)
    {
      // If version 1, read frame index + camera's 3D position
      if (version == 1)
      {
        int camera_index;
        double Cx,Cy,Cz;   
        if(sscanf(pLine, "%d %lg %lg %lg", &camera_index, &Cx,&Cy,&Cz) != 4)
        {
          std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
          return false;
        }
        Point_3 camera(Cx,Cy,Cz);
        (*cameras)[camera_index] = camera;
      }
      // If version 2, read frame index + 3x4 matrix and compute camera's 3D position
      else if (version == 2)
      {
        // read frame index + 3x4 matrix
        int camera_index;
        double m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23;   
        if(sscanf(pLine, 
                  "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", 
                  &camera_index, &m00,&m01,&m02,&m03,&m10,&m11,&m12,&m13,&m20,&m21,&m22,&m23) != 13)
        {
          std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
          return false;
        }
        // compute camera's 3D position = - transpose(R) * T if R/T are the matrix rotation/translation
        Aff_transformation_3 R(m00, m01, m02, m10, m11, m12, m20, m21, m22);
        Vector_3 T(m03, m13, m23);
        Point_3 camera = CGAL::ORIGIN - (R.inverse())(T);
        (*cameras)[camera_index] = camera;
      }
    }

    // Read 3D points on next lines
    else if (gyroviz_points.size() < positions_3D_count)
    {
      // Read label (with double quotes) + position...
      // WARNING: this code does not support spaces in labels.
      char point_3D_label[512];
      double X,Y,Z;   
      if (sscanf(pLine, "\"%[^\"]\" %lg %lg %lg", point_3D_label, &X,&Y,&Z) != 4)
      {
          std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
          return false;
      }
      Point_3 position_3D(X,Y,Z);
      gyroviz_points[point_3D_label] = position_3D;
    }

    // Read 2D points + camera indices + 3D points indices on remaining lines
    else 
    {
      // Read camera index, 3D point label (with double quotes) and 2D position...
      // WARNING: this code does not support spaces in labels.
      int camera_index;
      char point_3D_label[512];
      double x,y;   
      if (sscanf(pLine, "%d \"%[^\"]\" %lg %lg", &camera_index, point_3D_label, &x,&y) != 4 ||
          cameras->find(camera_index) == cameras->end())
      {
          std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
          return false;
      }
      // TEMPORARY? Skip 2D points not reconstructed in 3D.
      if (gyroviz_points.find(point_3D_label) == gyroviz_points.end())
      {
          //std::cerr << "Skip incorrect 2D point on line " << lineNumber << " of " << pFilename << std::endl;
      }
      else
      {
          Point_3 camera = (*cameras)[camera_index];
          Point_2 position_2D(x, y);
          gyroviz_points[point_3D_label].add_camera_point2_pair( std::make_pair(camera, position_2D) );
      }
    }
  }

  fclose(pFile);

  // Copy gyroviz_points[] to gyroviz_point_output
  for (std::map<std::string,Gyroviz_point>::iterator it=gyroviz_points.begin(); it != gyroviz_points.end(); it++)
  {
    Gyroviz_point gpt = it->second;
    
    // TEMPORARY? Skip 3D points with no cameras.
    if (gpt.cameras_begin() != gpt.cameras_end())
      *gyroviz_point_output++ = gpt;
    //else
    //  std::cerr << "Skip (" << (Point_3)gpt << ") line " << lineNumber << " of " << pFilename << std::endl;
  }
  
  return true;
}


#endif // SURFACE_RECONSTRUCTION_READ_G23_H
