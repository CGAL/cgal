#include <CGAL/Graphics_scene.h>

#if defined(CGAL_USE_BASIC_VIEWER_QT) || defined(CGAL_USE_BASIC_VIEWER)
#include <CGAL/Qt/Basic_viewer.h>
#elif defined(CGAL_USE_BASIC_VIEWER_GLFW)
#include <CGAL/GLFW/Basic_viewer.h>
#else 
namespace CGAL
{
  inline
  void draw_graphics_scene(const Graphics_scene&,
                           const char* ="CGAL Basic Viewer (GLFW)")
  {
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }
} // End namespace CGAL
#endif
