#include <CGAL/Three/Point_container.h>


typedef Viewer_interface VI;
using namespace CGAL::Three;

struct Point_d{

  Point_d()
    :
      plane(QVector4D())
  {
    f_matrix=QMatrix4x4();
  }

  QVector4D plane;
  QMatrix4x4 f_matrix;
};

Point_container::Point_container(int program, bool indexed)
  :Primitive_container(program, indexed),
    d(new Point_d())
{
  std::vector<Vbo*> vbos(NbOfVbos, NULL);
  setVbos(vbos);

}

void Point_container::initGL(Viewer_interface *viewer)
{
  viewer->makeCurrent();
  if(getVao(viewer))
    delete getVao(viewer);
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

  }
  setGLInit(viewer, true);
}

void Point_container::draw(Viewer_interface *viewer,
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
    viewer->glDrawElements(GL_POINTS, static_cast<GLuint>(getIdxSize()),
                           GL_UNSIGNED_INT, 0);
    getVbo(Indices)->release();
    getVao(viewer)->release();
  }
  else
  {
    getVao(viewer)->bind();
    if(is_color_uniform)
      getVao(viewer)->program->setAttributeValue("colors", getColor());
    if(getVao(viewer)->program->property("hasFMatrix").toBool())
      getVao(viewer)->program->setUniformValue("f_matrix", getFrameMatrix());
    viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(getFlatDataSize()/3));
    getVao(viewer)->release();

  }
}

QVector4D Point_container::getPlane() const { return d->plane; }
QMatrix4x4 Point_container::getFrameMatrix() const { return d->f_matrix; }

void Point_container::setFrameMatrix(const QMatrix4x4& m) { d->f_matrix = m; }
