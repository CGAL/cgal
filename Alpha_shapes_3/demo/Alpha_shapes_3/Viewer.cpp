#include "Viewer.h"
#include <vector>
#include <CGAL/bounding_box.h>
#include <QGLViewer/vec.h>
#include "CGAL/Qt/CreateOpenGLContext.h"

Viewer::Viewer(QWidget* parent)
  : QGLViewer(CGAL::Qt::createOpenGLContext(),parent)
{
  are_buffers_initialized = false;
}


void Viewer::compile_shaders()
{
    initializeOpenGLFunctions();
    if(! buffers[0].create() || !buffers[1].create() || !buffers[2].create()  )
    {
        std::cerr<<"VBO Creation FAILED"<<std::endl;
    }

    if(!vao[0].create() || !vao[1].create() )
    {
        std::cerr<<"VAO Creation FAILED"<<std::endl;
    }

    //Facets

    //Vertex source code
    const char vertex_source[] =
    {
         "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec3 normal;\n"

        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 mv_matrix; \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "void main(void)\n"
        "{\n"
        "   fP = mv_matrix * vertex; \n"
        "   fN = mat3(mv_matrix)* normal; \n"
        "   gl_Position = mvp_matrix * vertex;\n"
        "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
        "#version 120 \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "uniform highp vec4 color; \n"
        "uniform vec4 light_pos;  \n"
        "uniform vec4 light_diff; \n"
        "uniform vec4 light_spec; \n"
        "uniform vec4 light_amb;  \n"
        "uniform float spec_power ; \n"

        "void main(void) { \n"

        "   vec3 L = light_pos.xyz - fP.xyz; \n"
        "   vec3 V = -fP.xyz; \n"

        "   vec3 N = normalize(fN); \n"
        "   L = normalize(L); \n"
        "   V = normalize(V); \n"

        "   vec3 R = reflect(-L, N); \n"
        "   vec4 diffuse = abs(dot(N,L)) * light_diff * color; \n"
        "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

        "gl_FragColor = light_amb*color + diffuse + specular ; \n"
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

//Points

//Vertex source code
const char vertex_source_points[] =
{
    "#version 120 \n"
    "attribute highp vec4 vertex;\n"

    "uniform highp mat4 mvp_matrix;\n"
    "void main(void)\n"
    "{\n"
    "   gl_Position = mvp_matrix * vertex;\n"
    "}"
};
//Vertex source code
const char fragment_source_points[] =
{
    "#version 120 \n"
    "uniform highp vec4 color; \n"

    "void main(void) { \n"
    "gl_FragColor = color; \n"
    "} \n"
    "\n"
};
vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
if(!vertex_shader->compileSourceCode(vertex_source_points))
{
std::cerr<<"Compiling vertex source FAILED"<<std::endl;
}

fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
if(!fragment_shader->compileSourceCode(fragment_source_points))
{
std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
}

if(!rendering_program_points.addShader(vertex_shader))
{
std::cerr<<"adding vertex shader FAILED"<<std::endl;
}
if(!rendering_program_points.addShader(fragment_shader))
{
std::cerr<<"adding fragment shader FAILED"<<std::endl;
}
if(!rendering_program_points.link())
{
std::cerr<<"linking Program FAILED"<<std::endl;
}
rendering_program_points.bind();


}

void Viewer::initialize_buffers()
{

    vao[0].bind();
    buffers[0].bind();
    buffers[0].allocate(pos_poly.data(),
                        static_cast<int>(pos_poly.size()*sizeof(float)));
    poly_vertexLocation = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.enableAttributeArray(poly_vertexLocation);
    rendering_program.setAttributeBuffer(poly_vertexLocation,GL_FLOAT,0,3);
    buffers[0].release();

    buffers[1].bind();
    buffers[1].allocate(normals.data(), 
                        static_cast<int>(normals.size()*sizeof(float)));
    normalsLocation = rendering_program.attributeLocation("normal");
    rendering_program.bind();
    rendering_program.enableAttributeArray(normalsLocation);
    rendering_program.setAttributeBuffer(normalsLocation,GL_FLOAT,0,3);
    buffers[1].release();


    rendering_program.release();
    vao[0].release();

    vao[1].bind();
    buffers[2].bind();
    buffers[2].allocate(pos_points.data(),
                        static_cast<int>(pos_points.size()*sizeof(float)));
    points_vertexLocation = rendering_program_points.attributeLocation("vertex");
    rendering_program_points.bind();
    rendering_program_points.enableAttributeArray(points_vertexLocation);
    rendering_program_points.setAttributeBuffer(points_vertexLocation,GL_FLOAT,0,3);
    buffers[2].release();


    rendering_program_points.release();
    vao[0].release();

    are_buffers_initialized = true;
}


void Viewer::attrib_buffers(QGLViewer* viewer)
{
    QMatrix4x4 mvpMatrix;
    QMatrix4x4 mvMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvpMatrix.data()[i] = (float)mat[i];
    }
    viewer->camera()->getModelViewMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvMatrix.data()[i] = (float)mat[i];
    }
    // define material
     QVector4D	ambient(0.25f, 0.20725f, 0.20725f, 0.922f);
     QVector4D	diffuse( 1.0f,
                            0.829f,
                            0.829f,
                            0.922f );

    QVector4D	specular(  0.6f,
                            0.6f,
                            0.6f,
                            1.0f );

    QVector4D	position(0.0f,0.0f,1.0f,1.0f );
     GLfloat shininess =  11.264f;


    rendering_program.bind();
    mvpLocation = rendering_program.uniformLocation("mvp_matrix");
    mvLocation = rendering_program.uniformLocation("mv_matrix");
    colorLocation = rendering_program.uniformLocation("color");
    lightLocation[0] = rendering_program.uniformLocation("light_pos");
    lightLocation[1] = rendering_program.uniformLocation("light_diff");
    lightLocation[2] = rendering_program.uniformLocation("light_spec");
    lightLocation[3] = rendering_program.uniformLocation("light_amb");
    lightLocation[4] = rendering_program.uniformLocation("spec_power");

    rendering_program.setUniformValue(lightLocation[0], position);
    rendering_program.setUniformValue(lightLocation[1], diffuse);
    rendering_program.setUniformValue(lightLocation[2], specular);
    rendering_program.setUniformValue(lightLocation[3], ambient);
    rendering_program.setUniformValue(lightLocation[4], shininess);
    rendering_program.setUniformValue(mvpLocation, mvpMatrix);
    rendering_program.setUniformValue(mvLocation, mvMatrix);


    rendering_program.release();

    rendering_program_points.bind();

    mvpLocation_points = rendering_program_points.uniformLocation("mvp_matrix");
    colorLocation_points = rendering_program_points.uniformLocation("color");
    rendering_program_points.setUniformValue(mvpLocation_points, mvpMatrix);

    rendering_program_points.release();
}

void Viewer::initializeGL()
{
  QGLViewer::initializeGL();
  compile_shaders();
}


void
Viewer::sceneChanged()
{

  Iso_cuboid_3 bb = CGAL::bounding_box(scene->points.begin(), scene->points.end());
   
  this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(), bb.ymin(), bb.zmin()),
				      qglviewer::Vec(bb.xmax(),
						     bb.ymax(),
						     bb.zmax()));

  this->showEntireScene();

}

void Viewer::compute_elements()
{

    normals.resize(0);
    pos_points.resize(0);
    pos_poly.resize(0);
    for(std::list<Point_3>::iterator it = scene->points.begin();
        it != scene->points.end();
        ++it){
        pos_points.push_back(it->x()); pos_points.push_back(it->y()); pos_points.push_back(it->z());
    }

    std::list<Facet> facets;
    scene->alpha_shape.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);

    for(std::list<Facet>::iterator fit = facets.begin();
        fit != facets.end();
        ++fit) {
      const Cell_handle& ch = fit->first;
      const int index = fit->second;

      //const Vector_3& n = ch->normal(index); // must be unit vector

      const Point_3& a = ch->vertex((index+1)&3)->point();
      const Point_3& b = ch->vertex((index+2)&3)->point();
      const Point_3& c = ch->vertex((index+3)&3)->point();

      Vector_3 v = CGAL::unit_normal(a,b,c);

      normals.push_back(v.x()); normals.push_back(v.y()); normals.push_back(v.z());
      normals.push_back(v.x()); normals.push_back(v.y()); normals.push_back(v.z());
      normals.push_back(v.x()); normals.push_back(v.y()); normals.push_back(v.z());
      pos_poly.push_back(a.x()); pos_poly.push_back(a.y()); pos_poly.push_back(a.z());
      pos_poly.push_back(b.x()); pos_poly.push_back(b.y()); pos_poly.push_back(b.z());
      pos_poly.push_back(c.x()); pos_poly.push_back(c.y()); pos_poly.push_back(c.z());

    }


}

void
Viewer::draw()
{
    glEnable(GL_DEPTH_TEST);
    if(!are_buffers_initialized)
        initialize_buffers();
    QColor color;
    //points
    vao[1].bind();
    attrib_buffers(this);
    rendering_program_points.bind();
    color.setRgbF(1.0f, 0.0f, 0.0f);
    glPointSize(5);
    ::glEnable(GL_POINT_SMOOTH);
    rendering_program_points.setUniformValue(colorLocation_points, color);
    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(pos_points.size()/3));
    rendering_program_points.release();
    vao[1].release();
    //facets
    vao[0].bind();
    attrib_buffers(this);
    rendering_program.bind();
    color.setRgbF(0.5f, 1.0f, 0.5f);
    rendering_program.setUniformValue(colorLocation, color);
    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(pos_poly.size()/3));
    rendering_program.release();
    vao[0].release();

}

void Viewer::alphaChanged()
{

    normals.resize(0);
    pos_poly.resize(0);

    std::list<Facet> facets;
    scene->alpha_shape.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);

    for(std::list<Facet>::iterator fit = facets.begin();
        fit != facets.end();
        ++fit) {
      const Cell_handle& ch = fit->first;
      const int index = fit->second;

      //const Vector_3& n = ch->normal(index); // must be unit vector

      const Point_3& a = ch->vertex((index+1)&3)->point();
      const Point_3& b = ch->vertex((index+2)&3)->point();
      const Point_3& c = ch->vertex((index+3)&3)->point();

      Vector_3 v = CGAL::unit_normal(a,b,c);

      normals.push_back(v.x()); normals.push_back(v.y()); normals.push_back(v.z());
      normals.push_back(v.x()); normals.push_back(v.y()); normals.push_back(v.z());
      normals.push_back(v.x()); normals.push_back(v.y()); normals.push_back(v.z());
      pos_poly.push_back(a.x()); pos_poly.push_back(a.y()); pos_poly.push_back(a.z());
      pos_poly.push_back(b.x()); pos_poly.push_back(b.y()); pos_poly.push_back(b.z());
      pos_poly.push_back(c.x()); pos_poly.push_back(c.y()); pos_poly.push_back(c.z());

    }

    initialize_buffers();

}
