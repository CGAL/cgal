#include "Scene.h"
#include "Viewer.h"

#include <iostream>
#include <fstream>

#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QInputDialog>

#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/subdivision_method_3.h>

#include <CGAL/centroid.h>
#include <CGAL/linear_least_squares_fitting_3.h>

Scene::Scene()
{
  is_gl_init = false;
    m_pPolyhedron = NULL;

    // view options
    m_view_polyhedron = true;

}

Scene::~Scene()
{
    delete m_pPolyhedron;
}

int Scene::open(QString filename)
{
    QTextStream cerr(stderr);
    cerr << QString("Opening file \"%1\"\n").arg(filename);
    QApplication::setOverrideCursor(QCursor(::Qt::WaitCursor));

    QFileInfo fileinfo(filename);
    std::ifstream in(filename.toUtf8());

    if(!in || !fileinfo.isFile() || ! fileinfo.isReadable())
    {
        std::cerr << "unable to open file" << std::endl;
        QApplication::restoreOverrideCursor();
        return -1;
    }

    if(m_pPolyhedron != NULL)
        delete m_pPolyhedron;

    // allocate new polyhedron
    m_pPolyhedron = new Polyhedron;
    in >> *m_pPolyhedron;
    if(!in)
    {
        std::cerr << "invalid OFF file" << std::endl;
        QApplication::restoreOverrideCursor();

        delete m_pPolyhedron;
        m_pPolyhedron = NULL;

        return -1;
    }

    QApplication::restoreOverrideCursor();
    return 0;
}

void Scene::update_bbox()
{
    std::cout << "Compute bbox...";
    m_bbox = Bbox();

    if(m_pPolyhedron == NULL)
    {
        std::cout << "failed (no polyhedron)." << std::endl;
        return;
    }

    if(m_pPolyhedron->empty())
    {
        std::cout << "failed (empty polyhedron)." << std::endl;
        return;
    }

    Polyhedron::Point_iterator it = m_pPolyhedron->points_begin();
    m_bbox = (*it).bbox();
    for(; it != m_pPolyhedron->points_end();it++)
        m_bbox = m_bbox + (*it).bbox();
    std::cout << "done (" << m_pPolyhedron->size_of_facets()
        << " facets)" << std::endl;
}

void Scene::draw(Viewer* viewer)
{
  if(!is_gl_init)
    gl_init();
  rendering_program.bind();
  //setup mvp matrix
  QMatrix4x4 mvp_matrix;
  GLdouble mat[16];
  viewer->camera()->getModelViewProjectionMatrix(mat);
  for(int i=0; i<16; ++i)
    mvp_matrix.data()[i] = (GLfloat)mat[i];
  rendering_program.setUniformValue("mvp_matrix", mvp_matrix);

  if(m_view_polyhedron){
    rendering_program.setUniformValue("color", QColor(Qt::black));
    render_polyhedron(viewer);
  }
  //draw points and lines
  rendering_program.setUniformValue("color", QColor(Qt::darkBlue));
  render_line(viewer);
  rendering_program.setUniformValue("color", QColor(Qt::red));
  render_plane(viewer);
  rendering_program.setUniformValue("color", QColor(Qt::darkGreen));
  render_centroid(viewer);
  rendering_program.release();
}

void Scene::render_plane(Viewer* viewer)
{
    Point o = m_plane.projection(m_centroid);
    Point a = o + normalize(m_plane.base1()) + normalize(m_plane.base2());
    Point b = o + normalize(m_plane.base1()) - normalize(m_plane.base2());
    Point c = o - normalize(m_plane.base1()) - normalize(m_plane.base2());
    Point d = o - normalize(m_plane.base1()) + normalize(m_plane.base2());
    GLfloat verts[12];
    verts[0]= a.x();verts[1]=a.y();verts[2]=a.z();
    verts[3]= b.x();verts[4]=b.y();verts[5]=b.z();
    verts[6]= c.x();verts[7]=c.y();verts[8]=c.z();
    verts[9]= d.x();verts[10]=d.y();verts[11]=d.z();
    vao[0].bind();
    buffers[0].bind();
    buffers[0].allocate(verts, 12*sizeof(GLfloat));
    rendering_program.setAttributeBuffer("vertex", GL_FLOAT, 0,3);
    rendering_program.enableAttributeArray("vertex");
    buffers[0].release();

    viewer->glDrawArrays(GL_LINE_LOOP, 0, 4);
    vao[0].release();


}

void Scene::render_line(Viewer* viewer)
{
    Point o = m_line.projection(m_centroid);
    Point a = o + normalize(m_line.to_vector());
    Point b = o - normalize(m_line.to_vector());
    GLfloat verts[6];
    verts[0]=a.x();verts[1]=a.y();verts[2]=a.z();
    verts[3]=b.x();verts[4]=b.y();verts[5]=b.z();


    vao[1].bind();
    buffers[1].bind();
    buffers[1].allocate(verts, 6*sizeof(GLfloat));
    rendering_program.setAttributeBuffer("vertex", GL_FLOAT, 0,3);
    rendering_program.enableAttributeArray("vertex");
    buffers[1].release();
    viewer->glDrawArrays(GL_LINES, 0, 2);
    vao[1].release();
}

void Scene::render_centroid(Viewer* viewer)
{
    GLfloat verts[3];
    verts[0]=m_centroid.x();verts[1]=m_centroid.y();verts[2]=m_centroid.z();

    vao[2].bind();
    buffers[2].bind();
    buffers[2].allocate(verts, 3*sizeof(GLfloat));
    rendering_program.setAttributeBuffer("vertex", GL_FLOAT, 0,3);
    rendering_program.enableAttributeArray("vertex");
    buffers[2].release();
    viewer->glDrawArrays(GL_POINTS, 0, 1);
    vao[2].release();
}


Vector Scene::normalize(const Vector& v)
{
    return v / std::sqrt(v*v);
}

void Scene::render_polyhedron(Viewer *viewer)
{
    // draw black edges
    if(m_pPolyhedron != NULL)
    {
      typedef Kernel::Point_3 Point;

      std::vector<GLfloat> verts;

      Polyhedron::Edge_iterator he;
      for(he = m_pPolyhedron->edges_begin();
          he != m_pPolyhedron->edges_end();
          he++)
      {
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        verts.push_back(a.x());verts.push_back(a.y());verts.push_back(a.z());
        verts.push_back(b.x());verts.push_back(b.y());verts.push_back(b.z());
      }

      vao[3].bind();
      buffers[3].bind();
      buffers[3].allocate(verts.data(), static_cast<int>(verts.size()*sizeof(GLfloat)));
      rendering_program.setAttributeBuffer("vertex", GL_FLOAT, 0,3);
      rendering_program.enableAttributeArray("vertex");
      buffers[3].release();
      viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(verts.size()/3));
      vao[3].release();
    }
}

void Scene::refine_loop()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }
    std::cout << "Loop subdivision...";
    CGAL::Subdivision_method_3::Loop_subdivision(*m_pPolyhedron);
    std::cout << "done (" << m_pPolyhedron->size_of_facets() << " facets)" << std::endl;
}

void Scene::toggle_view_poyhedron()
{
    m_view_polyhedron = !m_view_polyhedron;
}

void Scene::fit_triangles()
{
    std::cout << "Fit triangles...";

    std::list<Triangle> triangles;
    Polyhedron::Facet_iterator it;
    for(it = m_pPolyhedron->facets_begin();
        it != m_pPolyhedron->facets_end();
        it++)
    {
        Polyhedron::Halfedge_handle he = it->halfedge();
        const Point& a = he->vertex()->point();
        const Point& b = he->next()->vertex()->point();
        const Point& c = he->next()->next()->vertex()->point();
        Triangle triangle(a,b,c);
        triangles.push_back(triangle);
    }

    m_centroid = CGAL::centroid(triangles.begin(),triangles.end());
    CGAL::linear_least_squares_fitting_3(triangles.begin(),
        triangles.end(), m_line, CGAL::Dimension_tag<2>());
    CGAL::linear_least_squares_fitting_3(triangles.begin(),
        triangles.end(), m_plane, CGAL::Dimension_tag<2>());

    std::cout << "done" << std::endl;
}

void Scene::fit_edges()
{
    std::cout << "Fit edges...";

    std::list<Segment> segments;
    Polyhedron::Edge_iterator he;
    for(he = m_pPolyhedron->edges_begin();
        he != m_pPolyhedron->edges_end();
        he++)
    {
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        Segment segment(a,b);
        segments.push_back(segment);
    }

    m_centroid = CGAL::centroid(segments.begin(),segments.end());
    CGAL::linear_least_squares_fitting_3(segments.begin(),
        segments.end(), m_line, CGAL::Dimension_tag<1>());
    CGAL::linear_least_squares_fitting_3(segments.begin(),
        segments.end(), m_plane, CGAL::Dimension_tag<1>());

    std::cout << "done" << std::endl;
}

void Scene::fit_vertices()
{
    std::cout << "Fit vertices...";

    std::list<Point> points;
    Polyhedron::Vertex_iterator v;
    for(v = m_pPolyhedron->vertices_begin();
        v != m_pPolyhedron->vertices_end();
        v++)
    {
        const Point& p = v->point();
        points.push_back(p);
    }

    m_centroid = CGAL::centroid(points.begin(),points.end());
    CGAL::linear_least_squares_fitting_3(points.begin(),
        points.end(), m_line, CGAL::Dimension_tag<0>());
    CGAL::linear_least_squares_fitting_3(points.begin(),
        points.end(), m_plane, CGAL::Dimension_tag<0>());

    std::cout << "done" << std::endl;
}

void Scene::gl_init()
{
  //Vertex source code
  const char vertex_source[] =
  {
      "#version 120 \n"
      "attribute highp vec4 vertex;\n"
      "uniform highp mat4 mvp_matrix;\n"
      "void main(void)\n"
      "{\n"
      "   gl_PointSize = 10.0;\n"
      "   gl_Position = mvp_matrix * vertex;\n"
      "}"
  };
  //Vertex source code
  const char fragment_source[] =
  {
      "#version 120 \n"
      "uniform highp vec4 color; \n"
      "void main(void) { \n"
      "  gl_FragColor = color; \n"
      "} \n"
      "\n"
  };
  QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
  if(!vertex_shader->compileSourceCode(vertex_source))
  {
      std::cerr<<"Compiling vertex source FAILED"<<std::endl;
  }

  QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
  if(!fragment_shader->compileSourceCode(fragment_source))
  {
      std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
  }

  if(!rendering_program.addShader(vertex_shader))
  {
      std::cerr<<"adding vertex shader FAILED"<<std::endl;
  }
  if(!rendering_program.addShader(fragment_shader))
  {
      std::cerr<<"adding fragment shader FAILED"<<std::endl;
  }
  if(!rendering_program.link())
  {
      std::cerr<<"linking Program FAILED"<<std::endl;
  }

  for(int i=0; i< 4; ++i)
  {
    buffers[i].create();
    vao[i].create();
  }
  is_gl_init = true;
}

