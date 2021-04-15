#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Three.h>


typedef Viewer_interface VI;
using namespace CGAL::Three;

struct Edge_d{

  Edge_d()
    :
      plane(QVector4D())
  {
    f_matrix.setToIdentity();
  }

  QVector4D plane;
  QMatrix4x4 f_matrix;
  GLfloat width;
  QVector2D viewport;
  bool is_surface;
};

Edge_container::Edge_container(int program, bool indexed)
  :Primitive_container(program, indexed),
    d(new Edge_d())
{
  std::vector<Vbo*> vbos(NbOfVbos, NULL);
  setVbos(vbos);

}

void Edge_container::initGL(Viewer_interface *viewer)
{
  viewer->makeCurrent();
  if(viewer->isSharing())
  {
    if(!getVao(viewer))
      setVao(viewer, new Vao(getVao(Three::mainViewer()),
                             viewer->getShaderProgram(getProgram())));
  }
  else
  {
    if(!getVao(viewer))
      setVao(viewer, new Vao(viewer->getShaderProgram(getProgram())));
    if(isDataIndexed())
    {
      if(!getVbo(Vertices))
        setVbo(Vertices,
               new Vbo("vertex",
                       Vbo::GEOMETRY));
      if(!getVbo(Indices))
        setVbo(Indices,
               new Vbo("indices",
                       Vbo::GEOMETRY,
                       QOpenGLBuffer::IndexBuffer));
      getVao(viewer)->addVbo(getVbo(Vertices));
      getVao(viewer)->addVbo(getVbo(Indices));
    }
    else
    {
      if(!getVbo(Vertices))
        setVbo(Vertices,
               new Vbo("vertex",
                       Vbo::GEOMETRY));
      if(!getVbo(Colors))
        setVbo(Colors,
               new Vbo("colors",
                       Vbo::COLORS));
      setVao(viewer, new Vao(viewer->getShaderProgram(getProgram())));
      getVao(viewer)->addVbo(getVbo(Vertices));
      getVao(viewer)->addVbo(getVbo(Colors));

      if(viewer->getShaderProgram(getProgram())->property("hasNormals").toBool())
      {
        if(!getVbo(Normals))
          setVbo(Normals,
                 new Vbo("normals",
                         Vbo::NORMALS));
        getVao(viewer)->addVbo(getVbo(Normals));
      }
      if(viewer->getShaderProgram(getProgram())->property("hasRadius").toBool())
      {
        if(!getVbo(Radius))
          setVbo(Radius,
                 new Vbo("radius",
                         Vbo::NOT_INSTANCED,
                         QOpenGLBuffer::VertexBuffer, GL_FLOAT, 0, 1));
        getVao(viewer)->addVbo(getVbo(Radius));
      }
      if(viewer->getShaderProgram(getProgram())->property("hasCenter").toBool())
      {
        if(!getVbo(Centers))
          setVbo(Centers,
                 new Vbo("center",
                  viewer->getShaderProgram(getProgram())->property("isInstanced").toBool()
                  ? Vbo::NOT_INSTANCED
                  : Vbo::GEOMETRY));
        getVao(viewer)->addVbo(getVbo(Centers));
      }
      if(viewer->getShaderProgram(getProgram())->property("hasTexture").toBool())
      {
        if(!getVbo(Texture_map))
          setVbo(Texture_map,
                 new Vbo("v_texCoord",
                         Vbo::GEOMETRY,
                         QOpenGLBuffer::VertexBuffer,
                         GL_FLOAT, 0, 2));
        getVao(viewer)->addVbo(getVbo(Texture_map));
        if(!getTexture())
          setTexture(new Texture());
      }
    }
  }
  setGLInit(viewer, true);
}

void Edge_container::initializeBuffers(Viewer_interface *viewer)
{
  Primitive_container::initializeBuffers(viewer);
  if(viewer->getShaderProgram(getProgram())->property("isInstanced").toBool())
  {
    getVao(viewer)->bind();
    if(viewer->getShaderProgram(getProgram())->property("hasCenter").toBool())
      viewer->glVertexAttribDivisor(getVao(viewer)->program->attributeLocation("center"), 1);
    if(viewer->getShaderProgram(getProgram())->property("hasRadius").toBool())
      viewer->glVertexAttribDivisor(getVao(viewer)->program->attributeLocation("radius"), 1);
    viewer->glVertexAttribDivisor(getVao(viewer)->program->attributeLocation("colors"), 1);
    getVao(viewer)->release();
  }
}

void Edge_container::draw(Viewer_interface *viewer,
                          bool is_color_uniform)

{
  bindUniformValues(viewer);
  if(isDataIndexed())
  {
    getVao(viewer)->bind();
    if(is_color_uniform)
      getVao(viewer)->program->setAttributeValue("colors", getColor());
    if(getVao(viewer)->program->property("hasFMatrix").toBool())
      getVao(viewer)->program->setUniformValue("f_matrix", getFrameMatrix());
    getVbo(Indices)->bind();
    viewer->glDrawElements(GL_LINES, static_cast<GLuint>(getIdxSize()),
                           GL_UNSIGNED_INT, 0);
    getVbo(Indices)->release();
    getVao(viewer)->release();
  }
  else
  {
    getVao(viewer)->bind();
    if(getVao(viewer)->program->property("hasTexture").toBool()){
      viewer->glActiveTexture(GL_TEXTURE0);
      viewer->glBindTexture(GL_TEXTURE_2D, getTextureId());
    }
    if(is_color_uniform)
      getVao(viewer)->program->setAttributeValue("colors", getColor());
    if(viewer->getShaderProgram(getProgram())->property("hasCutPlane").toBool())
      getVao(viewer)->program->setUniformValue("cutplane", d->plane);
    if(getVao(viewer)->program->property("hasSurfaceMode").toBool())
      getVao(viewer)->program->setUniformValue("is_surface", d->is_surface);
    if(getVao(viewer)->program->property("hasFMatrix").toBool())
      getVao(viewer)->program->setUniformValue("f_matrix", getFrameMatrix());

    if(getVao(viewer)->program->property("hasViewport").toBool())
    {
      getVao(viewer)->program->setUniformValue("viewport", getViewport());
      getVao(viewer)->program->setUniformValue("near",(GLfloat)viewer->camera()->zNear());
      getVao(viewer)->program->setUniformValue("far",(GLfloat)viewer->camera()->zFar());
    }
    if(getVao(viewer)->program->property("hasWidth").toBool())
      getVao(viewer)->program->setUniformValue("width", getWidth());

    if(viewer->getShaderProgram(getProgram())->property("isInstanced").toBool())
    {
      viewer->glDrawArraysInstanced(GL_LINES, 0,
                                    static_cast<GLsizei>(getFlatDataSize()/getTupleSize()),
                                    static_cast<GLsizei>(getCenterSize()/3));
    }
    else
    {
      viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(getFlatDataSize()/getTupleSize()));
    }
    getVao(viewer)->release();

  }
}

QVector4D Edge_container::getPlane() const { return d->plane; }
QMatrix4x4 Edge_container::getFrameMatrix() const { return d->f_matrix; }
GLfloat Edge_container::getWidth() const { return d->width; }
QVector2D Edge_container::getViewport() const { return d->viewport; }

void Edge_container::setPlane(const QVector4D& p) { d->plane = p; }
void Edge_container::setFrameMatrix(const QMatrix4x4& m) { d->f_matrix = m; }
void Edge_container::setWidth(const GLfloat & w) { d->width = w; }
void Edge_container::setViewport(const QVector2D & v) { d-> viewport = v; }
void Edge_container::setIsSurface  (const bool b) { d->is_surface = b; }
