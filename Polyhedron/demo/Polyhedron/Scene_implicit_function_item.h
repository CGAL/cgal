#ifndef SCENE_IMPLICIT_FUNCTION_ITEM_H
#define SCENE_IMPLICIT_FUNCTION_ITEM_H

#include <Scene_item.h>
#include <CGAL_demo/Scene_interface.h>
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
  : public Scene_item
{
  Q_OBJECT
  
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;
  
public:
  Scene_implicit_function_item(Implicit_function_interface*);
  virtual ~Scene_implicit_function_item();
  
  Implicit_function_interface* function() const { return function_; }

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  Bbox bbox() const;

  Scene_implicit_function_item* clone() const { return NULL; }

  // rendering mode
  virtual bool supportsRenderingMode(RenderingMode m) const;
  virtual bool manipulatable() const { return true; }
  virtual ManipulatedFrame* manipulatedFrame() { return frame_; }
  


  // actually draw() is also overloaded to detect when the cut plane is moved
  virtual void draw()const {}
  virtual void draw(Viewer_interface*) const;
  virtual void draw_edges(Viewer_interface*) const;

  virtual QString toolTip() const;
  virtual void contextual_changed();
  virtual void invalidate_buffers();
public Q_SLOTS:
  void plane_was_moved() { need_update_ = true; }
  void compute_function_grid() const;

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

  std::vector<float> positions_cube;
  std::vector<float> positions_grid;
  std::vector<float> positions_tex_quad;
  std::vector<float> texture_map;
  Texture *texture;


  mutable QOpenGLShaderProgram *program;
  mutable GLuint textureId;



  GLuint vao;
  GLuint buffer[4];
  using Scene_item::initialize_buffers;
  void initialize_buffers(Viewer_interface *viewer) const;
  void compute_vertices_and_texmap(void);
  void compute_texture(int, int);
};

#endif // SCENE_IMPLICIT_FUNCTION_ITEM
