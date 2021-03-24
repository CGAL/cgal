#ifndef TRIANGULATION_ON_SPHERE_2_VIEWER_H
#define TRIANGULATION_ON_SPHERE_2_VIEWER_H

#include <QMap>
#include <CGAL/Qt/qglviewer.h>

#include <QOpenGLFunctions_2_1>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>

#include <list>
#include <vector>

class Viewer
  : public CGAL::QGLViewer
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel          Kernel;
  typedef Kernel::FT                                                   FT;
  typedef Kernel::Point_3                                              Point_3;
  typedef Kernel::Segment_3                                            Segment_3;

  typedef CGAL::Projection_on_sphere_traits_3<Kernel>                  Projection_traits;
  typedef CGAL::Delaunay_triangulation_on_sphere_2<Projection_traits>  Projected_DToS2;

  typedef std::list<std::vector<Point_3> >                             Subsampled_arcs;

  template <class TOS>
  void build_the_boundary(const TOS& tos)
  {
    typedef typename TOS::Geom_traits::SK                            Spherical_kernel;

    if(tos.dimension() < 2)
      return;

    typename TOS::Solid_edges_iterator it = tos.solid_edges_begin();
    for (; it != tos.solid_edges_end(); ++it)
    {
      typename TOS::Edge e = *it;

      const bool diametral_edge = CGAL::collinear(tos.construct_point(tos.point(e.first, (e.second+1)%3)),
                                                  tos.construct_point(tos.point(e.first, (e.second+2)%3)),
                                                  tos.geom_traits().center());

      subsampled_arcs.emplace_front();

      // primal
      if(!diametral_edge)
      {
        typename TOS::Arc_on_sphere_2 as = tos.segment_on_sphere(e);
        std::vector<typename Spherical_kernel::Point_3> discretization_points;
        CGAL::Triangulations_on_sphere_2::internal::subsample_arc_on_sphere_2<Spherical_kernel>(
              as, std::back_inserter(discretization_points), min_edge_size);

        subsampled_arcs.begin()->reserve(discretization_points.size());
        for(const typename Spherical_kernel::Point_3& spt : discretization_points)
          subsampled_arcs.begin()->push_back(Point_3(spt.x(), spt.y(), spt.z()));
      }
      else
      {
        const Segment_3 s = tos.segment(e);
        subsampled_arcs.begin()->push_back(s.source());
        subsampled_arcs.begin()->push_back(s.target());
      }
    }
  }

public:
  Viewer(QWidget* parent = 0);

  template <typename TOS, typename Iterator>
  void open(Iterator begin, Iterator end,
            const TOS& tos)
  {
    draw_balls = true;
    draw_inputs = false;

    subsampled_arcs.clear();
    inputs.clear();

    pos_points.clear();
    pos_lines.clear();
    pos_sphere_inter.clear();
    normals_inter.clear();
    normals_lines.clear();

    radius_ = tos.geom_traits().radius();
    center_ = tos.geom_traits().center();
    min_edge_size = 0.01 * radius_;

    std::copy(begin, end, std::back_inserter(inputs));
    build_the_boundary(tos);
    compute_elements();
    initialize_buffers();

    qDebug() << "Finished loading";
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

  double min_edge_size;
  Subsampled_arcs subsampled_arcs;
  std::vector<Point_3> inputs;

  bool draw_inputs;
  bool draw_balls;

  enum VBO {
    SPHERE_POINTS =0,
    SPHERE_NORMALS,
    SPHERE_CENTER,
    POINTS_POINTS,
    POINTS_NORMALS,
    POINTS_CENTERS,
    EDGES_POINTS,
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
  QOpenGLShaderProgram rendering_program_edges;
  typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
  typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
  PFNGLDRAWARRAYSINSTANCEDARBPROC glDrawArraysInstanced;
  PFNGLVERTEXATTRIBDIVISORARBPROC glVertexAttribDivisor;
  double radius_;
  Point_3 center_;
  void initialize_buffers();
  void compute_elements();
  void attrib_buffers(CGAL::QGLViewer*);
  void compile_shaders();
};

#endif // TRIANGULATION_ON_SPHERE_2_VIEWER_H
