#ifndef TRIANGULATION_ON_SPHERE_2_VIEWER_H
#define TRIANGULATION_ON_SPHERE_2_VIEWER_H

#include <QMap>
#include <CGAL/Qt/qglviewer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <QOpenGLFunctions_2_1>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Cartesian_converter.h>

#include "Circular_arc_3_subsampling.h"



class Viewer : public CGAL::QGLViewer
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel          Kernel;
  typedef CGAL::Exact_spherical_kernel_3 SK;
  typedef Kernel::Point_3 Point_3;
  typedef std::list< std::vector<Point_3> > Subsampled_arcs;

  template <class Triangulation_on_sphere>
  void build_the_boundary(const Triangulation_on_sphere& T)
  {
    for (typename Triangulation_on_sphere::All_edges_iterator
         it=T.all_edges_begin();it!=T.all_edges_end();++it)
    {
      if ( it->first->is_ghost() &&
           it->first->neighbor(it->second)->is_ghost() )
        continue;

      Point_3 source=it->first->vertex( (it->second+1)%3 )->point();
      Point_3 target=it->first->vertex( (it->second+2)%3 )->point();
      subsampled_arcs.push_front(std::vector<Kernel::Point_3>());

      Kernel::Plane_3  plane(source,target,center_);
      Kernel::Circle_3 circle(center_,radius_,plane);
      subsample_circular_arc_3<Kernel>(circle,plane,source,target,std::back_inserter(*subsampled_arcs.begin()),min_edge_size);
    }
  }

public:
  Viewer(QWidget* parent = 0);
  template <class Triangulation_on_sphere,class Iterator>
  void open(Iterator begin, Iterator end,
            Triangulation_on_sphere& T,
            Point_3 center,
            double scale)
  {
    center_ = center;
    radius_ = scale;
    draw_balls = true;
    draw_inputs = false;
    min_edge_size = scale/100.;
    subsampled_arcs.clear();
    inputs.clear();
    pos_points.clear();
    pos_lines.clear();
    pos_sphere_inter.clear();
    normals_inter.clear();
    normals_lines.clear();
    std::copy(begin,end,std::back_inserter(inputs));
    build_the_boundary(T);
    compute_elements();
    initialize_buffers();
    qDebug()<<"done loading.";
  }

  GLuint dl_nb;

protected :
  virtual void draw();
  virtual void init();
  virtual QString helpString() const;
  virtual void keyPressEvent(QKeyEvent *e);
private:
  //Shaders elements

  int vertexLocation[3];
  int normalsLocation[3];
  int centerLocation;
  int trivialCenterLocation;
  int mvpLocation;
  int mvLocation;
  int colorLocation;
  int lightLocation[5];


  std::vector<float> pos_points;
  std::vector<float> pos_lines;
  std::vector<float> pos_sphere;
  std::vector<float> pos_sphere_inter;
  std::vector<float> normals;
  std::vector<float> normals_inter;
  std::vector<float> trivial_center;
  std::vector<float> normals_lines;
  Subsampled_arcs subsampled_arcs;
  std::vector<Point_3> inputs;
  Point_3 center_;
  double radius_;
  bool draw_balls;
  bool draw_inputs;
  double min_edge_size;

  enum VBO {
    SPHERE_POINTS =0,
    SPHERE_NORMALS,
    SPHERE_CENTER,
    POINTS_POINTS,
    POINTS_NORMALS,
    POINTS_CENTERS,
    EDGES_POINTS,
    EDGES_NORMALS,
    EDGES_CENTER,
    SIZE_OF_VBO
  };
  enum VAO {
    SPHERE_VAO = 0,
    POINTS_VAO,
    EDGES_VAO,
    SIZE_OF_VAO
  };
  QOpenGLBuffer buffers[SIZE_OF_VBO];
  QOpenGLVertexArrayObject vao[SIZE_OF_VAO];
  QOpenGLShaderProgram rendering_program;
  QOpenGLShaderProgram rendering_program_no_ext;
  typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
  typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
  PFNGLDRAWARRAYSINSTANCEDARBPROC glDrawArraysInstanced;
  PFNGLVERTEXATTRIBDIVISORARBPROC glVertexAttribDivisor;


  void initialize_buffers();
  void compute_elements();
  void attrib_buffers(CGAL::QGLViewer*);
  void compile_shaders();
};
#endif //TRIANGULATION_ON_SPHERE_2_VIEWER_H
