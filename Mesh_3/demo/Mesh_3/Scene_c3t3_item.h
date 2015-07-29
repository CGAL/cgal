#ifndef SCENE_C3T3_ITEM_H
#define SCENE_C3T3_ITEM_H

#include "Scene_c3t3_item_config.h"
#include "C3t3_type.h"
#include <CGAL_demo/Scene_item.h>

#include <QVector>
#include <QColor>
#include <set>

#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

struct Scene_c3t3_item_priv;

class SCENE_C3T3_ITEM_EXPORT Scene_c3t3_item
  : public Scene_item
{
  Q_OBJECT
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_c3t3_item();
  Scene_c3t3_item(const C3t3& c3t3);
  ~Scene_c3t3_item();

  const C3t3& c3t3() const;
  C3t3& c3t3();
  
  bool manipulatable() const {
    return true;
  }

  ManipulatedFrame* manipulatedFrame() {
    return frame;
  }

  void setPosition(float x, float y, float z) {
    frame->setPosition(x, y, z);
  }

  void setNormal(float x, float y, float z) {
    frame->setOrientation(x, y, z, 0.f);
  }

  Kernel::Plane_3 plane() const;

  bool isFinite() const { return true; }
  bool isEmpty() const {
    return c3t3().triangulation().number_of_vertices() == 0;
  }

  Bbox bbox() const;

  Scene_c3t3_item* clone() const {
    return 0;
  }

  QString toolTip() const;
  virtual QPixmap graphicalToolTip() const;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m != Gouraud); // CHECK THIS!
  }

  void draw(Viewer* viewer) const;
  void draw_edges(Viewer* viewer) const;
  
  // data item
  inline const Scene_item* data_item() const;
  inline void set_data_item(const Scene_item* data_item);
  
  // rebuild histogram
  inline void update_histogram();
  
  // Call this if c3t3 has been modified
  void c3t3_changed();
  void contextual_changed();

public Q_SLOTS:
  inline void data_item_destroyed();
  virtual void setColor(QColor c);
  
private:
  void build_histogram();
  void compute_color_map(const QColor& c);
  QColor get_histogram_color(const double v) const;
  
protected:
  Scene_c3t3_item_priv* d;

  qglviewer::ManipulatedFrame* frame;
  
private:
  QPixmap histogram_;
  const Scene_item* data_item_;
  
  typedef std::set<int> Indices;
  Indices indices_;

  static const int vaoSize = 2;
  static const int vboSize = 4;
  mutable bool are_buffers_initialized;
  mutable int poly_vertexLocation[2];
  mutable int normalsLocation[2];
  mutable int mvpLocation[2];
  mutable int mvLocation[2];
  mutable int fmatLocation;
  mutable int colorLocation[2];
  mutable int lightLocation[5*2];
  mutable int twosideLocation;

  std::vector<float> v_poly;
  std::vector<float> normal;
  std::vector<float> color_triangles;
  std::vector<float> v_grid;



  mutable QOpenGLBuffer buffers[vboSize];
  mutable QOpenGLVertexArrayObject vao[vaoSize];
  mutable QOpenGLShaderProgram rendering_program;
  mutable QOpenGLShaderProgram rendering_program_grid;
  void initialize_buffers() const;
  void compute_elements();
  void attrib_buffers(Viewer*) const;
  void compile_shaders();
  void draw_grid(float diag, std::vector<float> *positions_grid);
};

inline
const Scene_item*
Scene_c3t3_item::data_item() const
{
  return data_item_;
}

inline
void
Scene_c3t3_item::set_data_item(const Scene_item* data_item)
{
  data_item_ = data_item;
  
  if ( NULL != data_item )
  {
    connect(data_item, SIGNAL(aboutToBeDestroyed()),
            this, SLOT(data_item_destroyed()));
  }
}

inline
void
Scene_c3t3_item::update_histogram()
{
  build_histogram();
}

inline
void
Scene_c3t3_item::data_item_destroyed()
{
  set_data_item(NULL);
}



#endif // SCENE_C3T3_ITEM_H
