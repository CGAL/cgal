#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 400)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.0 or higher installed!\n";
 std::cout << "A LEDA version >= 4.0 is required to run GeoWin!\n";
 return 0;
}
#else 

#include <CGAL/Cartesian.h>
#include <CGAL/geowin_support.h>

int main()
{

  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointList"));
  geowin_init_default_type((CGALSegmentlist*)0, leda_string("CGALSegmentList"));
  geowin_init_default_type((CGALCirclelist*)0, leda_string("CGALCircleList"));
  geowin_init_default_type((CGALLinelist*)0, leda_string("CGALLineList"));
  geowin_init_default_type((CGALRaylist*)0, leda_string("CGALRayList"));
  geowin_init_default_type((CGALTrianglelist*)0, leda_string("CGALTriangleList"));
  geowin_init_default_type((CGALRectanglelist*)0, leda_string("CGALRectangleList"));

  geowin_init_default_type((CGALPolygonlist*)0, leda_string("CGALPolygonList"));

  geowin_init_default_type((CGALPoint_3_list*)0, leda_string("CGALPoint_3_List")); 

  GeoWin GW("GeoWin Demo using CGAL objects");
  GW.edit();
  return 0;  
}

#endif
