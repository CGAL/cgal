#ifndef SCENE_IMPLICIT_FUNCTION_ITEM_H
#define SCENE_IMPLICIT_FUNCTION_ITEM_H

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include "Scene_implicit_function_item_config.h"
#include "implicit_functions/Implicit_function_interface.h"
#include "Color_ramp.h"
#include <QGLViewer/manipulatedFrame.h>
#include <Viewer.h>

#define SCENE_IMPLICIT_GRID_SIZE 120

class Viewer_interface;

class Texture{
private:
     int Width;
     int Height;
     int size;
    GLubyte *data;
public:
    Texture(int w, int h)
    {
        Width = w;
        Height = h;
        size = 3*Height*Width;
        data = new GLubyte[size];
    }
    int getWidth() const {return Width;}
    int getHeight() const {return Height;}
    int getSize() const {return size;}
    void setData(int i, int j, int r, int g, int b){
        data[j*Width*3 +i*3] = r;
        data[j*Width*3 +i*3+1] = g;
        data[j*Width*3 +i*3+2] = b;}
    GLubyte* getData(){return data; }

};
class SCENE_IMPLICIT_FUNCTION_ITEM_EXPORT Scene_implicit_function_item 
  : public CGAL::Three::Scene_item
{
  Q_OBJECT
  
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;
  
public:
  Scene_implicit_function_item(Implicit_function_interface*);
  virtual ~Scene_implicit_function_item();
  
  Implicit_function_interface* function() const { return function_; }

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  void compute_bbox() const;

  Scene_implicit_function_item* clone() const { return NULL; }

  // rendering mode
  virtual bool supportsRenderingMode(RenderingMode m) const;
  virtual bool manipulatable() const { return true; }
  virtual ManipulatedFrame* manipulatedFrame() { return frame_; }
  


  // actually draw() is also overloaded to detect when the cut plane is moved
  virtual void draw()const {}
  virtual void draw(CGAL::Three::Viewer_interface*) const;
  virtual void draw_edges(CGAL::Three::Viewer_interface*) const;

  virtual QString toolTip() const;
  virtual void invalidateOpenGLBuffers();
public Q_SLOTS:
  void plane_was_moved() { need_update_ = true; }
  void compute_function_grid() const;
  void timerEvent(QTimerEvent*);

private:
  typedef qglviewer::Vec                  Point;
  typedef std::pair <Point,double>        Point_value;
  void compute_min_max();
  
private:
  Implicit_function_interface* function_;
  ManipulatedFrame* frame_;
  
  mutable bool need_update_;
  int grid_size_;
  double max_value_;
  double min_value_;
  mutable Point_value implicit_grid_[SCENE_IMPLICIT_GRID_SIZE][SCENE_IMPLICIT_GRID_SIZE];
  
  Color_ramp blue_color_ramp_;
  Color_ramp red_color_ramp_;

  enum VAOs {
      Plane = 0,
      BBox,
      Grid,
      NbOfVaos = Grid +1
  };
  enum VBOs {
      Quad_vertices = 0,
      TexMap,
      Cube_vertices,
      Grid_vertices,
      NbOfVbos = Grid_vertices +1
  };

  std::vector<float> positions_cube;
  std::vector<float> positions_grid;
  std::vector<float> positions_tex_quad;
  std::vector<float> texture_map;
  Texture *texture;


  mutable QOpenGLShaderProgram *program;
  mutable GLuint textureId;



  GLuint vao;
  GLuint buffer[4];
  using CGAL::Three::Scene_item::initialize_buffers;
  void initialize_buffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_vertices_and_texmap(void);
  void compute_texture(int, int);
};

#endif // SCENE_IMPLICIT_FUNCTION_ITEM
