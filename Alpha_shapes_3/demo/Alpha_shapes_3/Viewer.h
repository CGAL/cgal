#ifndef VIEWER_H
#define VIEWER_H

#include "typedefs.h"
#include <QGLViewer/qglviewer.h>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>



class Viewer : public QGLViewer, protected QOpenGLFunctions_3_3_Core{
  Q_OBJECT

  CGAL::Timer timer;
  Scene* scene;

  int nr_of_facets;
public:
  Viewer(QWidget* parent)
    : QGLViewer(createContext(),parent)
  {
    are_buffers_initialized = false;
  }
  ~Viewer()
  {
    buffers[0].destroy();
    buffers[1].destroy();
    buffers[2].destroy();
    vao[0].destroy();
    vao[1].destroy();
  }
  void setScene(Scene* scene_)
  {
    scene = scene_;
  }

  void clear();

public:
  void draw();

private:

  static QGLContext* createContext()
  {
      QOpenGLContext *context = new QOpenGLContext();
      QSurfaceFormat format;
      format.setVersion(3,3);
      format.setProfile(QSurfaceFormat::CompatibilityProfile);
      context->setFormat(format);
      return QGLContext::fromOpenGLContext(context);
  }
  bool are_buffers_initialized;
  //Shaders elements
    int poly_vertexLocation;
    int points_vertexLocation;
    int normalsLocation;
    int mvpLocation;
    int mvpLocation_points;
    int mvLocation;
    int colorLocation;
    int colorLocation_points;
    int lightLocation[5];


    std::vector<float> pos_points;
    std::vector<float> pos_poly;
    std::vector<float> normals;

    QOpenGLBuffer buffers[3];
    QOpenGLVertexArrayObject vao[2];
    QOpenGLShaderProgram rendering_program;
    QOpenGLShaderProgram rendering_program_points;
    void initialize_buffers();
    void compute_elements();
    void attrib_buffers(QGLViewer*);
    void compile_shaders();
 public Q_SLOTS:
    void initializeGL();
    void sceneChanged();
    void changed(){
        compute_elements();
        are_buffers_initialized = false;
    }
    void alphaChanged();

};

#endif
