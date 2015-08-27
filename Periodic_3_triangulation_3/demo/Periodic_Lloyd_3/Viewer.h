#ifndef VIEWER_H
#define VIEWER_H

#include "Scene.h"
#include <QGLViewer/qglviewer.h>
#include <QOpenGLFunctions_2_1>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <CGAL/Qt/CreateOpenGLContext.h>


class Viewer : public QGLViewer, QOpenGLFunctions_2_1 {

  typedef qglviewer::Vec Vec;

  Q_OBJECT

  CGAL::Timer timer;
  Scene* scene;

  int nr_of_facets;
public:
  Viewer(QWidget* parent)
    : QGLViewer(CGAL::Qt::createOpenGLContext(), parent)
  {}
  ~Viewer()
  {
   for(int i=0; i<4; i++)
   {
    buffers[i].destroy();
    vao[i].destroy();
   }
  }

  void setScene(Scene* scene_)
  {
    scene = scene_;
  }

  void init();
  void clear();

public:
  void draw();


public slots :
  void changed();
  void sceneChanged();
  void render_video();
  
signals:
  void valueChanged(int i);

private:
  Vec next_around_circle(const float& phi, const Vec& pos, const Vec& ori);

      int vertexLocation[3];
      int mvpLocation[3];
      int colorLocation[3];

      bool are_buffers_initialized;
      std::vector<float> pos_points;
      std::vector<float> pos_lines;
      std::vector<float> pos_8lines2D;
      std::vector<float> pos_8lines;


      QOpenGLBuffer buffers[4];
      QOpenGLVertexArrayObject vao[4];
      QOpenGLShaderProgram rendering_program;
      void initialize_buffers();
      void compute_elements();
      void attrib_buffers(QGLViewer*);
      void compile_shaders();
};

#endif
