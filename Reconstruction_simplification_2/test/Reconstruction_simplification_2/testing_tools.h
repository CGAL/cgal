#include <fstream>


template <class PointMassList, class Point>
void load_xy_file(const std::string& filename, PointMassList& points)
{
   std::ifstream ifs(filename);
   Point point;
   while (ifs >> point)
      points.push_back(std::make_pair(point, 1));

   ifs.close();
}

template <class Point>
void load_xy_file_points(const std::string& fileName, std::list<Point>& points)
{
   std::ifstream ifs(fileName);
   Point point;
   while (ifs >> point)
   {
       points.push_back(point);
   }
   ifs.close();
}
