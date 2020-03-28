#include "Viewer.h"
#include <vector>
#include <CGAL/bounding_box.h>
#include <CGAL/Qt/vec.h>


void
Viewer::init()
{
  initializeOpenGLFunctions();
  setBackgroundColor(::Qt::white);
  this->camera()->setSceneBoundingBox(
      CGAL::qglviewer::Vec(-1.,-1.,-1.),
      CGAL::qglviewer::Vec( 1., 1., 1.));
  glEnable(GL_LINE_SMOOTH);
  compile_shaders();
  are_buffers_initialized = false;
}

void
Viewer::sceneChanged()
{
  changed();
  this->showEntireScene();
}

void
Viewer::draw()
{
    glEnable(GL_DEPTH_TEST);
  if(!are_buffers_initialized)
      initialize_buffers();
  QColor color;
 //the points
  glEnable(GL_BLEND);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(3.0f,-3.0f);


  attrib_buffers(this);
  vao[0].bind();
  color.setRgbF(1.0f, 0.72f, 0.06f);
  rendering_program.bind();
  rendering_program.setUniformValue(colorLocation[0], color);
  glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(pos_points.size()/3));
  rendering_program.release();
  vao[0].release();

  //The Lines
  glDisable(GL_POLYGON_OFFSET_FILL);

  vao[1].bind();
  color.setRgbF(0.27f, 0.51f, 0.7f);
  rendering_program.bind();
  rendering_program.setUniformValue(colorLocation[0], color);
  glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_lines.size()/3));
  rendering_program.release();
  vao[1].release();

  if (scene->eight_copies) {
      vao[2].bind();
      color.setRgbF(0.69f, 0.77f, 0.87f);
      rendering_program.bind();
      rendering_program.setUniformValue(colorLocation[0], color);
      glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_8lines2D.size()/3));
      rendering_program.release();
      vao[2].release();


      if (!scene->two_dimensional) {
          vao[3].bind();
          color.setRgbF(0.69f, 0.77f, 0.87f);
          rendering_program.bind();
          rendering_program.setUniformValue(colorLocation[0], color);
          glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_8lines.size()/3));
          rendering_program.release();
          vao[3].release();
      }
  }
}


Viewer::Vec
Viewer::next_around_circle(const float& phi, const Vec& pos, const Vec& ori) {
  Vec cam = pos-ori;
  Vec cam_norm = cam/cam.norm();

  Vec y(cam_norm.z, 0, -cam_norm.x);
  Vec y_norm = y/y.norm();

  Vec new_cam = ori + (cam_norm*cos(phi) + y_norm*sin(phi)) * cam.norm();
  return new_cam;
}

void Viewer::initialize_buffers()
{
    rendering_program.bind();

    vao[0].bind();
    buffers[0].bind();
    buffers[0].allocate(pos_points.data(), static_cast<int>(pos_points.size()*sizeof(float)));
    vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(vertexLocation[0]);
    rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
    buffers[0].release();

    vao[0].release();

    vao[1].bind();
    buffers[1].bind();
    buffers[1].allocate(pos_lines.data(), static_cast<int>(pos_lines.size()*sizeof(float)));
    vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(vertexLocation[0]);
    rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
    buffers[1].release();
    vao[1].release();

    vao[2].bind();
    buffers[2].bind();
    buffers[2].allocate(pos_8lines2D.data(), static_cast<int>(pos_8lines2D.size()*sizeof(float)));
    vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(vertexLocation[0]);
    rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
    buffers[2].release();
    vao[2].release();

    vao[3].bind();
    buffers[3].bind();
    buffers[3].allocate(pos_8lines.data(), static_cast<int>(pos_8lines.size()*sizeof(float)));
    vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(vertexLocation[0]);
    rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
    buffers[3].release();
    vao[3].release();
    rendering_program.release();
    are_buffers_initialized = true;
}

void Viewer::compute_elements()
{
    //glColor3f(1.0f, 0.72f, 0.06f);
    pos_points.resize(0);
    pos_lines.resize(0);
    pos_8lines.resize(0);
    pos_8lines2D.resize(0);

    for(Periodic_point_iterator ppit
      = scene->periodic_triangulation.periodic_points_begin(
          P3DT3::UNIQUE) ;
        ppit != scene->periodic_triangulation.periodic_points_end(P3DT3::UNIQUE);
        ++ppit){
      Point_3 p(scene->periodic_triangulation.point(*ppit));
      pos_points.push_back(p.x()); pos_points.push_back(p.y()); pos_points.push_back(p.z());
    }

    //glColor3f(0.27f, 0.51f, 0.7f);
    for (Periodic_triangle_iterator ptit
       = scene->periodic_triangulation.periodic_triangles_begin(
           P3DT3::UNIQUE);
         ptit != scene->periodic_triangulation.periodic_triangles_end(
         P3DT3::UNIQUE);
        ++ptit) {
      for (int i=0 ; i<4 ; i++) {
        Segment_3 dual = scene->periodic_triangulation.construct_segment(
        scene->periodic_triangulation.dual(*(ptit.get_facet())));

        FT sz = dual.source().z();
        FT tz = dual.target().z();

        if (scene->two_dimensional && ((sz-tz > 0.5) || (sz-tz < -0.5))) continue;

        if (scene->two_dimensional) { sz = 0.; tz = 0.; }
        FT sx = dual.source().x();
        FT tx = dual.target().x();
        FT sy = dual.source().y();
        FT ty = dual.target().y();
        pos_lines.push_back(sx); pos_lines.push_back(sy); pos_lines.push_back(sz);
        pos_lines.push_back(tx); pos_lines.push_back(ty); pos_lines.push_back(tz);
        Iso_cuboid_3 d = scene->periodic_triangulation.domain();

        pos_8lines2D.push_back(sx+0.); pos_8lines2D.push_back(sy+d.ymax()-d.ymin()); pos_8lines2D.push_back(sz+0.);
        pos_8lines2D.push_back(tx+0.); pos_8lines2D.push_back(ty+d.ymax()-d.ymin()); pos_8lines2D.push_back(tz+0.);
        pos_8lines2D.push_back(sx+d.xmax()-d.xmin()); pos_8lines2D.push_back(sy+0.); pos_8lines2D.push_back(sz+0.);
        pos_8lines2D.push_back(tx+d.xmax()-d.xmin()); pos_8lines2D.push_back(ty+0.); pos_8lines2D.push_back(tz+0.);
        pos_8lines2D.push_back(sx+d.xmax()-d.xmin()); pos_8lines2D.push_back(sy+d.ymax()-d.ymin()); pos_8lines2D.push_back(sz+0.);
        pos_8lines2D.push_back(tx+d.xmax()-d.xmin()); pos_8lines2D.push_back(ty+d.ymax()-d.ymin()); pos_8lines2D.push_back(tz+0.);

        pos_8lines.push_back(sx+0.); pos_8lines.push_back(sy+0.); pos_8lines.push_back(sz+d.zmax()-d.zmin());
        pos_8lines.push_back(tx+0.); pos_8lines.push_back(ty+0.); pos_8lines.push_back(tz+d.zmax()-d.zmin());
        pos_8lines.push_back(sx+0.); pos_8lines.push_back(sy+d.ymax()-d.ymin()); pos_8lines.push_back(sz+d.zmax()-d.zmin());
        pos_8lines.push_back(tx+0.); pos_8lines.push_back(ty+d.ymax()-d.ymin()); pos_8lines.push_back(tz+d.zmax()-d.zmin());
        pos_8lines.push_back(sx+d.xmax()-d.xmin()); pos_8lines.push_back(sy+0.); pos_8lines.push_back(sz+d.zmax()-d.zmin());
        pos_8lines.push_back(tx+d.xmax()-d.xmin()); pos_8lines.push_back(ty+0.); pos_8lines.push_back(tz+d.zmax()-d.zmin());
        pos_8lines.push_back(sx+d.xmax()-d.xmin()); pos_8lines.push_back(sy+d.ymax()-d.ymin()); pos_8lines.push_back(sz+d.zmax()-d.zmin());
        pos_8lines.push_back(tx+d.xmax()-d.xmin()); pos_8lines.push_back(ty+d.ymax()-d.ymin()); pos_8lines.push_back(tz+d.zmax()-d.zmin());
      }

    }
}

void Viewer::attrib_buffers(CGAL::QGLViewer* viewer)
{
    QMatrix4x4 mvpMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvpMatrix.data()[i] = (float)mat[i];
    }

    rendering_program.bind();
    mvpLocation[0] = rendering_program.uniformLocation("mvp_matrix");
    colorLocation[0] = rendering_program.uniformLocation("color");
    rendering_program.setUniformValue(mvpLocation[0], mvpMatrix);

    rendering_program.release();
}

void Viewer::compile_shaders()
{
    if(! buffers[0].create() || !buffers[1].create() || !buffers[2].create() || !buffers[3].create())
    {
        std::cerr<<"VBO Creation FAILED"<<std::endl;
    }

    if(!vao[0].create() || !vao[1].create() || !vao[2].create() || !vao[3].create())
    {
        std::cerr<<"VAO Creation FAILED"<<std::endl;
    }

    //Vertex source code
    const char vertex_source[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "uniform highp mat4 mvp_matrix;\n"
        "void main(void)\n"
        "{\n"
        "   gl_PointSize = 5.0; \n"
        "   gl_Position = mvp_matrix * vertex; \n"
        "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
        "#version 120 \n"
        "uniform highp vec4 color; \n"
        "void main(void) { \n"
        "gl_FragColor = color; \n"
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
    rendering_program.bind();

}

void Viewer::changed()
{
    compute_elements();
    are_buffers_initialized = false;
}
