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
