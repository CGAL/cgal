// Author     : Nader Salman


// an object used as info in the CGAL::Triangulation_vertex_base_with_info_2

#ifndef _Gyroviz_info_for_cdt2_
#define _Gyroviz_info_for_cdt2_

class	Gyroviz_info_for_cdt2
{

protected:
	
	// 3D Point
	Point_3 point;

	// Camera number
	int camera_number;

	// Flag := on border or not 
	bool flag;

public:
  Gyroviz_info_for_cdt2(){}
  Gyroviz_info_for_cdt2(Point_3 p, int cam_num, bool f):point(p),camera_number(cam_num),flag(f){}

  // accessors
  const	Point_3	    get_point3() const { return	point; }
  const int         get_camera_number() const { return camera_number; }
  const	bool	    get_flag() const { return	flag;	}
  
  // modificators
  void	set_point3(Point_3 point3) { point = point3; }
  void set_camera_number(int cam_num) { camera_number = cam_num; }
  void	set_flag(bool f) { flag = f; }

};
#endif // _Gyroviz_info_for_cdt2_