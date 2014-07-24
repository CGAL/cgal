#include <fstream>


template <class PointMassList, class Point>
void load_xy_file(const std::string& fileName, PointMassList& points)
{
   std::ifstream ifs(fileName);
   std::cerr << "read xy...";
   Point point;
   unsigned int nb = 0;
   while (ifs >> point)
   {
	   points.push_back(std::make_pair(point, 1));
   }
   std::cerr << "done (" << nb << " points)" << std::endl;
   ifs.close();
}
