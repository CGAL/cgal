#include <CGAL/basic.h>
#if !defined(CGAL_USE_FLTK) || !defined(CGAL_USE_OPENGL)

#include <iostream>

int main()
{
 std::cout << "No FLTK or OpenGL or Mesa version installed!\n";
 return 0;
}

#else 
#include <CGAL/Cartesian.h>
#include <CGAL/Viewer_stream.h>

typedef CGAL::Cartesian<double> rep_t;
typedef CGAL::Point_3<rep_t> point_t;

int main()
{
 CGAL::Viewer_3 W(500);
 point_t p(100,100,100);
 CGAL::Drawable_point_3<point_t> dp(p,CGAL::RED,CGAL::FILL,25,50);
 W.add_drawable(&dp);
 return 0;
}
#endif
