#ifndef SCENE_IMPLICIT_FUNCTION_ITEM_H
#define SCENE_IMPLICIT_FUNCTION_ITEM_H

#include <CGAL_demo/Scene_item.h>
#include <CGAL_demo/Scene_interface.h>
#include "Scene_implicit_function_item_config.h"
#include "implicit_functions/Implicit_function_interface.h"
#include "Color_ramp.h"

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

#define SCENE_IMPLICIT_GRID_SIZE 120

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
  
  // draw (overload only direct_draw() to use display list of base class)
  virtual void draw(Viewer* viewer) const;
  
  virtual QString toolTip() const;

public Q_SLOTS:
  void compute_function_grid();

private:
  typedef qglviewer::Vec                  Point;
  typedef std::pair <Point,double>        Point_value;
  
  void compute_min_max();
  
private:
  Implicit_function_interface* function_;
  ManipulatedFrame* frame_;
  
  int grid_size_;
  double max_value_;
  double min_value_;
  Point_value implicit_grid_[SCENE_IMPLICIT_GRID_SIZE][SCENE_IMPLICIT_GRID_SIZE];
  
  Color_ramp blue_color_ramp_;
  Color_ramp red_color_ramp_;

  mutable bool are_buffers_initialized;

  static const int vaoSize = 2;
  static const int vboSize = 3;
  mutable int vertexLocation[2];
  mutable int mvpLocation[2];
  mutable int colorLocation[2];
  mutable int tex_Location;
  mutable int f_Location;


   std::vector<float> v_cube;
   std::vector<float> v_plan;

   std::vector<float> texture_map;
   Texture *texture;
   mutable GLuint textureId;
   mutable bool texture_initialized;
   GLint sampler_location;

  mutable QOpenGLBuffer buffers[vboSize];
  mutable QOpenGLVertexArrayObject vao[vaoSize];
  mutable QOpenGLShaderProgram rendering_program;
  mutable QOpenGLShaderProgram tex_rendering_program;
  void initialize_buffers(Viewer*) const;
  void compute_elements();
  void attrib_buffers(Viewer*) const;
  void compile_shaders();
  void compute_texture(int, int);

public Q_SLOTS:
    void changed();
};

#endif // SCENE_IMPLICIT_FUNCTION_ITEM
