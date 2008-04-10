// Author     : Nader Salman

 /// PAS BESOIN DE CA ///

// an object used as info in the CGAL::Triangulation_vertex_base_with_info_3
// to be revised !!

#ifndef _Gyroviz_info_for_dt3_
#define _Gyroviz_info_for_dt3_

#include <vector>

class	Gyroviz_info_for_dt3
{

protected:

  // List of images
  std::vector<Point_3> list_of_cameras;

public:
  Gyroviz_info_for_dt3(){}
  Gyroviz_info_for_dt3(std::vector<Point_3> v):list_of_cameras(v){}

  // accessors
  const	std::vector<Point_3>	  get_list_of_cameras() const { return	list_of_cameras; }
};

#endif // _Gyroviz_info_for_dt3_