#include <CGAL/Three/Triangle_container.h>
#include <QOpenGLFramebufferObject>

typedef Viewer_interface VI;
using namespace CGAL::Three;

struct D{

  D():
    shrink_factor(1.0f),
    comparing(false),
    plane(QVector4D()),
    width(0.0f),
    height(0.0f),
    near(0.0f),
    far(0.0f),
    writing(false),
    alpha(1.0f)
  {}

  Triangle_container* container;
  float shrink_factor;

  bool comparing;

  QVector4D plane;

  float width;

  float height;

  float near;

  float far;

  bool writing;
  float alpha;
};

Triangle_container::Triangle_container(int program, bool indexed)
  : Primitive_container(program, indexed),
    d(new D())
{
  std::vector<Vbo*> vbos(NbOfVbos, NULL);
  setVbos(vbos);
}

void Triangle_container::initGL( Viewer_interface* viewer)const
{
  viewer->makeCurrent();
  if(isDataIndexed())
  {
    switch(getProgram())
    {
    case VI::PROGRAM_WITHOUT_LIGHT:
    case VI::PROGRAM_WITH_LIGHT:
    {
      if(!getVbo(Smooth_vertices))
        setVbo(Smooth_vertices,
            new Vbo("vertex", Vbo::GEOMETRY));
      if(!getVbo(Vertex_indices))
        setVbo(Vertex_indices,
            new Vbo("indices",
                    Vbo::GEOMETRY,
                    QOpenGLBuffer::IndexBuffer));
      if(!getVao(viewer))
        setVao(viewer, new Vao(viewer->getShaderProgram(getProgram())));
      getVao(viewer)->addVbo(getVbo(Smooth_vertices));
      getVao(viewer)->addVbo(getVbo(Vertex_indices));
    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }
    switch(getProgram())
    {
    case VI::PROGRAM_WITH_LIGHT:
    {
      if(!getVbo(Smooth_normals))
        setVbo(Smooth_normals,
            new Vbo("normals",
                    Vbo::NORMALS));
      if(!getVbo(VColors))
        setVbo(VColors,
            new Vbo("colors",
                    Vbo::COLORS,
                    QOpenGLBuffer::VertexBuffer, GL_FLOAT, 0, 4));
      getVao(viewer)->addVbo(getVbo(Smooth_normals));
      getVao(viewer)->addVbo(getVbo(VColors));
    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }
  }
  else
  {
    switch(getProgram())
    {

    case VI::PROGRAM_WITH_LIGHT:
    case VI::PROGRAM_C3T3:
    //case VI::PROGRAM_C3T3_TETS:
    case VI::PROGRAM_SPHERES:
    case VI::PROGRAM_CUTPLANE_SPHERES:
    {
      if(!getVbo(Flat_vertices))
      {
        setVbo(Flat_vertices,
            new Vbo("vertex",
                    Vbo::GEOMETRY));
      }
      if(!getVbo(Flat_normals))
        setVbo(Flat_normals,
            new Vbo("normals",
                    Vbo::NORMALS));
      if(!getVbo(FColors))
        setVbo(FColors,
            new Vbo("colors",
                    Vbo::COLORS,
                    QOpenGLBuffer::VertexBuffer, GL_FLOAT, 0, 4));
      setVao(viewer, new Vao(viewer->getShaderProgram(getProgram())));
      getVao(viewer)->addVbo(getVbo(Flat_vertices));
      getVao(viewer)->addVbo(getVbo(Flat_normals));
      getVao(viewer)->addVbo(getVbo(FColors));

    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }
    switch(getProgram())
    {
    case VI::PROGRAM_C3T3:
      //case VI::PROGRAM_C3T3_TETS:
    case VI::PROGRAM_SPHERES:
    case VI::PROGRAM_CUTPLANE_SPHERES:
    {
      if(!getVbo(Facet_barycenters))
        setVbo(Facet_barycenters,
            new Vbo("barycenter",
                    Vbo::GEOMETRY));
      getVao(viewer)->addVbo(getVbo(Facet_barycenters));
    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }
    if(getProgram() == VI::PROGRAM_SPHERES
       || getProgram() == VI::PROGRAM_CUTPLANE_SPHERES)
    {
      if(!getVbo(Radius))
        setVbo(Radius,
            new Vbo("radius",
                    Vbo::GEOMETRY,
                    QOpenGLBuffer::VertexBuffer, GL_FLOAT, 0, 1));
      getVao(viewer)->addVbo(getVbo(Radius));
      viewer->glVertexAttribDivisor(getVao(viewer)->program->attributeLocation("barycenter"), 1);
      viewer->glVertexAttribDivisor(getVao(viewer)->program->attributeLocation("radius"), 1);
      viewer->glVertexAttribDivisor(getVao(viewer)->program->attributeLocation("colors"), 1);
    }
  }
  setGLInit(viewer, true);
}

void Triangle_container::draw(Viewer_interface* viewer,
                              bool is_color_uniform,
                              QOpenGLFramebufferObject* fbo) const
{

  bindUniformValues(viewer);

  if(isDataIndexed())
  {
    getVao(viewer)->bind();
    if(is_color_uniform)
      getVao(viewer)->program->setAttributeValue("colors", getColor());
    getVbo(Vertex_indices)->bind();
    if(getProgram() == VI::PROGRAM_WITH_LIGHT)
    {
      getVao(viewer)->program->setUniformValue("comparing", d->comparing);
      getVao(viewer)->program->setUniformValue("width", d->width);
      getVao(viewer)->program->setUniformValue("height", d->height);
      getVao(viewer)->program->setUniformValue("near", d->near);
      getVao(viewer)->program->setUniformValue("far", d->far);
      getVao(viewer)->program->setUniformValue("writing", d->writing);
      getVao(viewer)->program->setUniformValue("alpha", d->alpha);
      if( fbo)
        viewer->glBindTexture(GL_TEXTURE_2D, fbo->texture());
    }
    viewer->glDrawElements(GL_TRIANGLES, static_cast<unsigned int>(getIdxSize()),
                           GL_UNSIGNED_INT, 0 );
    getVbo(Vertex_indices)->release();
    getVao(viewer)->release();
  }
  else
  {
    getVao(viewer)->bind();
    if(getProgram() == VI::PROGRAM_C3T3)
      getVao(viewer)->program->setUniformValue("shrink_factor", d->shrink_factor);
    if(getProgram() == VI::PROGRAM_C3T3
       || getProgram() == VI::PROGRAM_CUTPLANE_SPHERES)
      getVao(viewer)->program->setUniformValue("cutplane", d->plane);
    if(is_color_uniform)
      getVao(viewer)->program->setAttributeValue("colors", getColor());
    getVao(viewer)->program->setUniformValue("alpha", d->alpha);
    if(getProgram() == VI::PROGRAM_SPHERES
       || getProgram() == VI::PROGRAM_CUTPLANE_SPHERES)
    {
      viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                    static_cast<GLsizei>(getFlatDataSize()/3),
                                    static_cast<GLsizei>(getCenterSize()/3));
    }
    else
    {
      if(getProgram() == VI::PROGRAM_WITH_LIGHT)
      {
        getVao(viewer)->program->setUniformValue("comparing", d->comparing);
        getVao(viewer)->program->setUniformValue("width", d->width);
        getVao(viewer)->program->setUniformValue("height", d->height);
        getVao(viewer)->program->setUniformValue("near", d->near);
        getVao(viewer)->program->setUniformValue("far", d->far);
        getVao(viewer)->program->setUniformValue("writing", d->writing);
        if( fbo)
          viewer->glBindTexture(GL_TEXTURE_2D, fbo->texture());
      }
      viewer->glDrawArrays(GL_TRIANGLES,0,static_cast<GLsizei>(getFlatDataSize()/3));
    }

    getVao(viewer)->release();
  }
}



float     Triangle_container::getShrinkFactor()const { return d->shrink_factor ; }
bool      Triangle_container::isComparing()const     { return d->comparing; }
QVector4D Triangle_container::getPlane()const        { return d->plane; }
float     Triangle_container::getWidth()const        { return d->width; }
float     Triangle_container::getHeight()const       { return d->height; }
float     Triangle_container::getNear()const         { return d->near; }
float     Triangle_container::getFar()const          { return d->far; }
bool      Triangle_container::isDepthWriting()const  { return d->writing; }
float     Triangle_container::getAlpha()const        { return d->alpha; }

void Triangle_container::setShrinkFactor(const float& f)const { d->shrink_factor = f; }
void Triangle_container::setComparing   (const bool& b)const     { d->comparing = b; }
void Triangle_container::setPlane       (const QVector4D& p)const    { d->plane = p; }
void Triangle_container::setWidth       (const float& f)const        { d->width = f; }
void Triangle_container::setHeight      (const float& f)const       { d->height = f; }
void Triangle_container::setNear        (const float& f)const         { d->near = f; }
void Triangle_container::setFar         (const float& f)const          { d->far = f; }
void Triangle_container::setDepthWriting(const bool& b)const  { d->writing = b; }
void Triangle_container::setAlpha       (const float& f)const        { d->alpha = f ; }
