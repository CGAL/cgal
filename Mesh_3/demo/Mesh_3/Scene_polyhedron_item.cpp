#include "config.h"

#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <QObject>


Scene_polyhedron_item::Scene_polyhedron_item()
  : Scene_item(),
    poly(new Polyhedron)
{
    are_buffers_initialized = false;
    compile_shaders();
}

Scene_polyhedron_item::Scene_polyhedron_item(Polyhedron* const p)
  : Scene_item(),
    poly(p)
{
    are_buffers_initialized = false;
    compile_shaders();
}

Scene_polyhedron_item::Scene_polyhedron_item(const Polyhedron& p)
  : Scene_item(),
    poly(new Polyhedron(p))
{
    are_buffers_initialized = false;
    compile_shaders();
}

Scene_polyhedron_item::~Scene_polyhedron_item()
{
  for(int i=0; i<vboSize; i++)
      buffers[i].destroy();
  for(int i=0; i<vaoSize; i++)
      vao[i].destroy();
  delete poly;
}

/**************************************************
****************SHADER FUNCTIONS******************/

void Scene_polyhedron_item::compile_shaders()
{

    for(int i=0; i< vboSize; i++)
        buffers[i].create();
    for(int i=0; i< vaoSize; i++)
        vao[i].create();

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
        "   gl_Position = mvp_matrix * vertex; \n"
        "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
        "#version 120 \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "uniform vec4 color; \n"
        "uniform bool is_two_side; \n"
        "uniform highp vec4 light_pos;  \n"
        "uniform highp vec4 light_diff; \n"
        "uniform highp vec4 light_spec; \n"
        "uniform highp vec4 light_amb;  \n"
        "uniform float spec_power ; \n"

        "void main(void) { \n"

        "   vec3 L = light_pos.xyz - fP.xyz; \n"
        "   vec3 V = -fP.xyz; \n"

        "   vec3 N; \n"
        "   if(fN == vec3(0.0,0.0,0.0)) \n"
        "       N = vec3(0.0,0.0,0.0); \n"
        "   else \n"
        "       N = normalize(fN); \n"
        "   L = normalize(L); \n"
        "   V = normalize(V); \n"

        "   vec3 R = reflect(-L, N); \n"
        "   vec4 diffuse; \n"
        "   if(!is_two_side) \n"
        "       diffuse = max(dot(N,L),0) * light_diff*color; \n"
        "   else \n"
        "       diffuse = max(abs(dot(N,L)),0) * light_diff*color; \n"
        "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

         "gl_FragColor = color*light_amb + diffuse + specular; \n"
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
}

void Scene_polyhedron_item::compute_elements()
{

    v_poly.resize(0);
    v_edge.resize(0);
    normal_flat.resize(0);
    normal_smooth.resize(0);
    Polyhedron& polyhedron =*poly;
    //FACETS

    typedef Polyhedron::Traits	    Kernel;
    typedef Kernel::Point_3	    Point;
    typedef Kernel::Vector_3	    Vector;
    typedef Polyhedron::Facet_iterator Facet_iterator;
    typedef Polyhedron::Halfedge_around_facet_circulator HF_circulator;



    Facet_iterator f;
    for(f = polyhedron.facets_begin();
      f != polyhedron.facets_end();
      f++)
    {

      // If Flat shading: 1 normal per polygon

        Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(f, polyhedron);

        normal_flat.push_back(n.x()); normal_flat.push_back(n.y()); normal_flat.push_back(n.z());
        normal_flat.push_back(n.x()); normal_flat.push_back(n.y()); normal_flat.push_back(n.z());
        normal_flat.push_back(n.x()); normal_flat.push_back(n.y()); normal_flat.push_back(n.z());


      // revolve around current face to get vertices
      HF_circulator he = f->facet_begin();
      HF_circulator end = he;
      CGAL_For_all(he,end)
      {

          Vector n = CGAL::Polygon_mesh_processing::compute_vertex_normal(he->vertex(), polyhedron);
          normal_smooth.push_back(n.x()); normal_smooth.push_back(n.y()); normal_smooth.push_back(n.z());

        const Point& p = he->vertex()->point();
        v_poly.push_back(p.x()); v_poly.push_back(p.y()); v_poly.push_back(p.z());
      }
    }

    //EDGES

    typedef Polyhedron::Edge_iterator	Edge_iterator;


    Edge_iterator he;
    for(he = polyhedron.edges_begin();
      he != polyhedron.edges_end();
      he++)
    {
      const Point& a = he->vertex()->point();
      const Point& b = he->opposite()->vertex()->point();
      v_edge.push_back(a.x()); v_edge.push_back(a.y()); v_edge.push_back(a.z());
      v_edge.push_back(b.x()); v_edge.push_back(b.y()); v_edge.push_back(b.z());
    }

}

void Scene_polyhedron_item::initialize_buffers() const
{
    rendering_program.bind();

        vao[0].bind();
        buffers[0].bind();
        buffers[0].allocate(v_poly.data(), static_cast<int>(v_poly.size()*sizeof(float)));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[0].release();

        buffers[1].bind();
        buffers[1].allocate(normal_flat.data(), static_cast<int>(normal_flat.size()*sizeof(float)));
        normalsLocation[0] = rendering_program.attributeLocation("normal");
        rendering_program.enableAttributeArray(normalsLocation[0]);
        rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[1].release();

        vao[0].release();

        vao[1].bind();
        buffers[2].bind();
        buffers[2].allocate(v_poly.data(), static_cast<int>(v_poly.size()*sizeof(float)));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[2].release();

        buffers[3].bind();
        buffers[3].allocate(normal_smooth.data(), static_cast<int>(normal_smooth.size()*sizeof(float)));
        normalsLocation[0] = rendering_program.attributeLocation("normal");
        rendering_program.enableAttributeArray(normalsLocation[0]);
        rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[3].release();

        vao[1].release();

        vao[2].bind();
        buffers[4].bind();
        buffers[4].allocate(v_edge.data(), static_cast<int>(v_edge.size()*sizeof(float)));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[4].release();
        std::vector<float> empty_array;
        empty_array.resize(0);
        for(std::size_t i=0; i<v_edge.size(); i++)
            empty_array.push_back(0.0);
        buffers[5].bind();
        buffers[5].allocate(empty_array.data(), static_cast<int>(empty_array.size()*sizeof(float)));
        normalsLocation[0] = rendering_program.attributeLocation("normal");
        rendering_program.enableAttributeArray(normalsLocation[0]);
        rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[5].release();

        vao[2].release();
        rendering_program.release();
        are_buffers_initialized = true;

}

void Scene_polyhedron_item::attrib_buffers(Viewer* viewer) const
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
    QVector4D	position(0.0f,0.0f,1.0f,1.0f );
    GLboolean isTwoSide;
    viewer->glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE,&isTwoSide);
    // define material
     QVector4D	ambient;
     QVector4D	diffuse;
     QVector4D	specular;
     GLfloat      shininess ;
    // Ambient
    ambient[0] = 0.29225f;
    ambient[1] = 0.29225f;
    ambient[2] = 0.29225f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.50754f;
    diffuse[1] = 0.50754f;
    diffuse[2] = 0.50754f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.0f;
    specular[1] = 0.0f;
    specular[2] = 0.0f;
    specular[3] = 0.0f;
    // Shininess
    shininess = 51.2f;


    rendering_program.bind();
    colorLocation[0] = rendering_program.uniformLocation("color");
    twosideLocation = rendering_program.uniformLocation("is_two_side");
    mvpLocation[0] = rendering_program.uniformLocation("mvp_matrix");
    mvLocation[0] = rendering_program.uniformLocation("mv_matrix");
    lightLocation[0] = rendering_program.uniformLocation("light_pos");
    lightLocation[1] = rendering_program.uniformLocation("light_diff");
    lightLocation[2] = rendering_program.uniformLocation("light_spec");
    lightLocation[3] = rendering_program.uniformLocation("light_amb");
    lightLocation[4] = rendering_program.uniformLocation("spec_power");

    rendering_program.setUniformValue(lightLocation[0], position);
    rendering_program.setUniformValue(twosideLocation, isTwoSide);
    rendering_program.setUniformValue(mvpLocation[0], mvpMatrix);
    rendering_program.setUniformValue(mvLocation[0], mvMatrix);
    rendering_program.setUniformValue(lightLocation[1], diffuse);
    rendering_program.setUniformValue(lightLocation[2], specular);
    rendering_program.setUniformValue(lightLocation[3], ambient);
    rendering_program.setUniformValue(lightLocation[4], shininess);

    rendering_program.release();
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_polyhedron_item::draw(Viewer* viewer) const {
    if(!are_buffers_initialized)
        initialize_buffers();
    QColor color;
    GLint shading;
    viewer->glGetIntegerv(GL_SHADE_MODEL, &shading);
    if(shading == GL_FLAT)
        vao[0].bind();
    else if(shading == GL_SMOOTH)
        vao[1].bind();
    else
        return;
    attrib_buffers(viewer);

    rendering_program.bind();
    float current_color[4];
    viewer->glGetFloatv(GL_CURRENT_COLOR, current_color);
    color.setRgbF(current_color[0],current_color[1],current_color[2],current_color[3]);
    rendering_program.setUniformValue(colorLocation[0], color);
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(v_poly.size()/3));
    rendering_program.release();
    if(shading == GL_FLAT)
        vao[0].release();
    else if(shading == GL_SMOOTH)
        vao[1].release();


}
void Scene_polyhedron_item::draw_edges(Viewer* viewer) const {
    if(!are_buffers_initialized)
        initialize_buffers();
    QColor color;

    vao[2].bind();
    attrib_buffers(viewer);
    rendering_program.bind();
    float current_color[4];
    viewer->glGetFloatv(GL_CURRENT_COLOR, current_color);
    color.setRgbF(current_color[0],current_color[1],current_color[2],current_color[3]);
    rendering_program.setUniformValue(colorLocation[0], color);
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(v_edge.size()/3));
    rendering_program.release();
    vao[2].release();

}
void Scene_polyhedron_item::draw_points(Viewer* viewer) const {
    if(!are_buffers_initialized)
        initialize_buffers();
    QColor color;

    vao[0].bind();
    attrib_buffers(viewer);
    rendering_program.bind();
    float current_color[4];
    viewer->glGetFloatv(GL_CURRENT_COLOR, current_color);
    color.setRgbF(current_color[0],current_color[1],current_color[2],current_color[3]);
    rendering_program.setUniformValue(colorLocation[0], color);
    viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(v_poly.size()/3));
    rendering_program.release();
    vao[0].release();
}

Scene_polyhedron_item* 
Scene_polyhedron_item::clone() const {
  return new Scene_polyhedron_item(*poly);
}

// Load polyhedron from .OFF file
bool
Scene_polyhedron_item::load(std::istream& in)
{
  in >> *poly;
  changed();
  return in && !isEmpty();
}

// Write polyhedron to .OFF file
bool 
Scene_polyhedron_item::save(std::ostream& out) const
{
  out << *poly;
  return (bool) out;
}

QString 
Scene_polyhedron_item::toolTip() const
{
  if(!poly)
    return QString();

  return QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4</p>")
    .arg(this->name())
    .arg(poly->size_of_vertices())
    .arg(poly->size_of_halfedges()/2)
    .arg(poly->size_of_facets())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}



Polyhedron* 
Scene_polyhedron_item::polyhedron()       { return poly; }
const Polyhedron* 
Scene_polyhedron_item::polyhedron() const { return poly; }

bool
Scene_polyhedron_item::isEmpty() const {
  return (poly == 0) || poly->empty();
}

Scene_polyhedron_item::Bbox
Scene_polyhedron_item::bbox() const {
  const Point& p = *(poly->points_begin());
  CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
  for(Polyhedron::Point_iterator it = poly->points_begin();
      it != poly->points_end();
      ++it) {
    bbox = bbox + it->bbox();
  }
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

void Scene_polyhedron_item::changed()
{
    compute_elements();
    are_buffers_initialized = false;
}

