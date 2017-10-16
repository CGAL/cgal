#include <CGAL/Three/Edge_container.h>


typedef Viewer_interface VI;
using namespace CGAL::Three;

struct D{

  D()
    :
      plane(QVector4D())
  {
    f_matrix.setToIdentity();
  }

  QVector4D plane;
  QMatrix4x4 f_matrix;
};

Edge_container::Edge_container(int program, bool indexed)
  :Primitive_container(program, indexed),
    d(new D)
{
  std::vector<Vbo*> vbos(NbOfVbos, NULL);
  setVbos(vbos);

}

void Edge_container::initGL(Viewer_interface *viewer)
{
  viewer->makeCurrent();
  if(isDataIndexed())
  {
    switch(getProgram())
    {
    case VI::PROGRAM_WITHOUT_LIGHT:
    case VI::PROGRAM_NO_SELECTION:
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
      if(getVao(viewer))
        delete getVao(viewer);
      setVao(viewer, new Vao(viewer->getShaderProgram(getProgram())));
      getVao(viewer)->addVbo(getVbo(Vertices));
      getVao(viewer)->addVbo(getVbo(Indices));
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
    case VI::PROGRAM_SPHERES:
    case VI::PROGRAM_CUTPLANE_SPHERES:
    case VI::PROGRAM_NO_SELECTION:
    case VI::PROGRAM_WITHOUT_LIGHT:
    case VI::PROGRAM_C3T3_EDGES:
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
    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }

    switch(getProgram())
    {
    case VI::PROGRAM_SPHERES:
    case VI::PROGRAM_CUTPLANE_SPHERES:
    {
      if(!getVbo(Normals))
        setVbo(Normals,
            new Vbo("normals",
                    Vbo::NORMALS));
      if(!getVbo(Radius))
        setVbo(Radius,
            new Vbo("radius",
                    Vbo::GEOMETRY,
                    QOpenGLBuffer::VertexBuffer, GL_FLOAT, 0, 1));
      if(!getVbo(Barycenters))
        setVbo(Barycenters,
            new Vbo("barycenter",
                    Vbo::GEOMETRY));
      getVao(viewer)->addVbo(getVbo(Normals));
      getVao(viewer)->addVbo(getVbo(Radius));
      getVao(viewer)->addVbo(getVbo(Barycenters));
    }
      break;
    default:
      break;
    }
  }
  setGLInit(viewer, true);
}

void Edge_container::initializeBuffers(Viewer_interface *viewer)
{
  Primitive_container::initializeBuffers(viewer);
  if(getProgram() == VI::PROGRAM_SPHERES
     || getProgram() == VI::PROGRAM_CUTPLANE_SPHERES)
  {
    getVao(viewer)->bind();
    viewer->glVertexAttribDivisor(getVao(viewer)->program->attributeLocation("barycenter"), 1);
    viewer->glVertexAttribDivisor(getVao(viewer)->program->attributeLocation("radius"), 1);
    viewer->glVertexAttribDivisor(getVao(viewer)->program->attributeLocation("colors"), 1);
    getVao(viewer)->release();
  }
}

void Edge_container::draw(Viewer_interface *viewer,
                          bool is_color_uniform, QOpenGLFramebufferObject *)

{
  bindUniformValues(viewer);

  if(isDataIndexed())
  {
    getVao(viewer)->bind();
    if(is_color_uniform)
      getVao(viewer)->program->setAttributeValue("colors", getColor());
    getVbo(Indices)->bind();
    viewer->glDrawElements(GL_LINES, static_cast<GLuint>(getIdxSize()),
                           GL_UNSIGNED_INT, 0);
    getVbo(Indices)->release();
    getVao(viewer)->release();
  }
  else
  {
    getVao(viewer)->bind();
    if(is_color_uniform)
      getVao(viewer)->program->setAttributeValue("colors", getColor());
    if(getProgram() == VI::PROGRAM_WITHOUT_LIGHT)
      getVao(viewer)->program->setUniformValue("f_matrix", d->f_matrix);
    if(getProgram() == VI::PROGRAM_C3T3_EDGES
       || getProgram() == VI::PROGRAM_CUTPLANE_SPHERES)
      getVao(viewer)->program->setUniformValue("cutplane", d->plane);
    if(getProgram() == VI::PROGRAM_SPHERES
       || getProgram() == VI::PROGRAM_CUTPLANE_SPHERES)
    {
      viewer->glDrawArraysInstanced(GL_LINES, 0,
                                    static_cast<GLsizei>(getFlatDataSize()/3),
                                    static_cast<GLsizei>(getCenterSize()/3));
    }
    else
    {
      viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(getFlatDataSize()/3));
    }
    getVao(viewer)->release();

  }
}

QVector4D Edge_container::getPlane() const { return d->plane; }
QMatrix4x4 Edge_container::getFrameMatrix() const { return d->f_matrix; }

void Edge_container::setPlane(const QVector4D& p) { d->plane = p; }
void Edge_container::setFrameMatrix(const QMatrix4x4& m) { d->f_matrix = m; }
