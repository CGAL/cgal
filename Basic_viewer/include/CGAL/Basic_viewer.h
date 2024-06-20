#include <CGAL/Graphics_scene.h>

#ifdef CGAL_USE_BASIC_VIEWER || CGAL_USE_BASIC_VIEWER_QT 
#include <CGAL/Qt/Basic_viewer.h>
#elif CGAL_USE_BASIC_VIEWER_GLFW
#include <CGAL/GLFW/Basic_viewer.h>
#else 
namespace CGAL
{
  inline
  void draw_graphics_scene(const Graphics_scene&,
                           const char* ="CGAL Basic Viewer")
  {
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }
} // End namespace CGAL
#endif
