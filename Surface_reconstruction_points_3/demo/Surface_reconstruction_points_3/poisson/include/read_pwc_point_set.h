#ifndef READ_PWC_POINT_SET_H
#define READ_PWC_POINT_SET_H

#include <CGAL/value_type_traits.h>

#include <stdio.h>
#include <vector>


///
/// WARNING: .pwc format is deprecated.
///


/// Read 3D points + cameras from a Gyroviz .pwc file.
/// @return true on success.
template <typename GyrovizPointOutputIterator, ///< GyrovizPointOutputIterator value_type must be Gyroviz_point_3
          typename PointOutputIterator> ///< PointOutputIterator value_type must be Point_3
bool read_pwc_point_set(const char* pFilename,
                        GyrovizPointOutputIterator gyroviz_point_output,
                        PointOutputIterator camera_output)
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  typedef typename CGAL::value_type_traits<GyrovizPointOutputIterator>::type
                                              Gyroviz_point;

  typedef typename Gyroviz_point::Geom_traits Geom_traits;
  typedef typename Geom_traits::Point_3       Point_3;
  typedef typename Geom_traits::Point_2       Point_2;

  CGAL_precondition(pFilename != NULL);

  FILE *pFile = fopen(pFilename,"rt");
  if(pFile == NULL)
  {
    std::cerr << "Error: cannot open " << pFilename << std::endl;
    return false;
  }

  long positions_3D_count = 0, cameras_count = 0; // number of points and number of cameras in the file
  int lineNumber = 0; // current line number
  char pLine[4096]; // current line buffer
  std::vector<Point_3> cameras; // container of cameras read

  // TEST //
  // std::ofstream fout("../scores.txt");
  // TEST//

  while(fgets(pLine,sizeof(pLine),pFile))
  {
    lineNumber++;

    // Read file signature on first line
    if (lineNumber == 1)
    {
      char signature[512];
      if ( (sscanf(pLine,"%s",signature) != 1) || (strcmp(signature, "CamOFF") != 0) )
      {
        // if unsupported file format
        std::cerr << "Incorrect file format line " << lineNumber << " of " << pFilename << std::endl;
        return false;
      }
    }

    // Read number of cameras and points on 2nd line
    else if (lineNumber == 2)
    {
      if (sscanf(pLine,"%ld %ld",&cameras_count,&positions_3D_count) != 2)
      {
        std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
        return false;
      }
      cameras.reserve(cameras_count);
    }

    // Read cameras on next lines
    else if (cameras.size() < cameras_count)
    {
      // Read position...
      double Cx,Cy,Cz;
      if(sscanf(pLine,"%lg\t%lg\t%lg",&Cx,&Cy,&Cz) == 3)
      {
        Point_3 camera(Cx,Cy,Cz);
        cameras.push_back(camera);
      }
      // ...or skip comment line
    }

    // Read 3D points + camera indices on next lines
    else
    {
      // Read position + camera indices...
      std::istringstream iss(pLine);
      Point_3 position_3D;
      if (iss >> position_3D)
      {
        Gyroviz_point gyroviz_point = position_3D;

        int camera_index;
        while (iss >> camera_index)
        {
          if (camera_index-1 >= cameras.size())
          {
            std::cerr << "Error line " << lineNumber << " of " << pFilename << std::endl;
            return false;
          }
          Point_3 camera = cameras[camera_index-1];
          Point_2 fake(-1, -1); // TEMPORARY: fake 2D position
          gyroviz_point.add_camera_point2_pair( std::make_pair(camera, fake) );
        }

        // TEMPORARY? Skip 3D points with no cameras.
        if (gyroviz_point.cameras_begin() != gyroviz_point.cameras_end())
        {
          // TEST //
          //fout << position_3D << '\t' << list_of_cameras.size() << '\t' << vertex_score << std::endl;
          // TEST //

          *gyroviz_point_output++ = gyroviz_point;
        }
        //else
        //{
        //  std::cerr << "Skip (" << position_3D << ") line " << lineNumber << " of " << pFilename << std::endl;
        //}
      }
      // ...or skip comment line
    }
  }

  fclose(pFile);

  // TEST //
  // fout.close();
  // TEST //

  // Copy cameras[] to camera_output
  std::copy(cameras.begin(), cameras.end(), camera_output);

  return true;
}


#endif // READ_PWC_POINT_SET_H
