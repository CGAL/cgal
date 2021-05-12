#include <CGAL/Three/Primitive_container.h>
using namespace CGAL::Three;
struct D
{
  D(int program, bool indexed)
    :
      program_id(program),
      indexed(indexed),
      is_selected(false),
      clipping(true),
      flat_size(0),
      idx_size(0),
      tuple_size(3),
      color(QColor()),
      textureId(-1),
      texture(nullptr)
  {}
  ~D()
  {
    if(texture)
      delete texture;
  }
  QMap<CGAL::Three::Viewer_interface*, Vao*> VAOs;
  std::vector<Vbo*> VBOs;
  int program_id;
  bool indexed;
  QMap<CGAL::Three::Viewer_interface*, bool> is_init;
  QMap<CGAL::Three::Viewer_interface*, bool> is_gl_init;
  bool is_selected;
  bool clipping;
  std::size_t flat_size;
  std::size_t center_size;
  std::size_t idx_size;
  int tuple_size;
  QColor color;
  GLuint textureId;
  Texture *texture;
};

Primitive_container::Primitive_container(int program, bool indexed)
  :d(new D(program, indexed))
{}

Primitive_container::~Primitive_container()
{
  Q_FOREACH(Vbo* vbo, d->VBOs)
    if(vbo)
      delete vbo;
  Q_FOREACH(CGAL::Three::Viewer_interface*viewer, d->VAOs.keys())
  {
    removeViewer(viewer);
  }
}

void Primitive_container::bindUniformValues(CGAL::Three::Viewer_interface* viewer)
{
  viewer->attribBuffers(d->program_id);
  viewer->getShaderProgram(d->program_id)->bind();
  if(d->is_selected)
    viewer->getShaderProgram(d->program_id)->setUniformValue("is_selected", true);
  else
    viewer->getShaderProgram(d->program_id)->setUniformValue("is_selected", false);
  //override clipping if necessary
  if(!d->clipping)
    viewer->getShaderProgram(d->program_id)->setUniformValue("is_clipbox_on", false);
  viewer->getShaderProgram(d->program_id)->release();
}

void Primitive_container::initializeBuffers(CGAL::Three::Viewer_interface* viewer)
{
  if(!d->VAOs[viewer])
    return;
  viewer->makeCurrent();
  d->VAOs[viewer]->bind();
  Q_FOREACH(CGAL::Three::Vbo* vbo, d->VAOs[viewer]->vbos)
  {
    vbo->bind();
    if(vbo->dataSize !=0)
    {
      if(!vbo->allocated)
      {
        if(vbo->vbo_type == QOpenGLBuffer::IndexBuffer)
          vbo->vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);

        vbo->vbo.allocate(vbo->data, vbo->dataSize);
        vbo->allocated = true;
      }
      if(vbo->vbo_type == QOpenGLBuffer::VertexBuffer)
      {
        d->VAOs[viewer]->program->enableAttributeArray(vbo->attribute);
        d->VAOs[viewer]->program->setAttributeBuffer(
              vbo->attribute, vbo->data_type, vbo->offset, vbo->tupleSize, vbo->stride);
      }
    }
    else if(vbo->vbo_type == QOpenGLBuffer::VertexBuffer)
    {
      d->VAOs[viewer]->program->disableAttributeArray(vbo->attribute);
    }
    vbo->release();
  }
  if(getVao(viewer)->program->property("hasTexture").toBool())
  {
    getVao(viewer)->bind();
    if(GLuint(-1) == d->textureId) {
      viewer->glGenTextures(1, &d->textureId);
    }
    viewer->glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    viewer->glBindTexture(GL_TEXTURE_2D, d->textureId);
    viewer->glTexImage2D(GL_TEXTURE_2D,
                         0,
                         GL_RGB,
                         d->texture->Width,
                         d->texture->Height,
                         0,
                         GL_RGB,
                         GL_UNSIGNED_BYTE,
                         d->texture->data);
    viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    getVao(viewer)->release();
  }
  d->VAOs[viewer]->release();

  d->is_init[viewer] = true;
}

void Primitive_container::removeViewer(CGAL::Three::Viewer_interface* viewer)
{
  setGLInit(viewer, false);
  delete d->VAOs[viewer];
  d->VAOs.remove(viewer);
}

void Primitive_container::reset_vbos(Scene_item_rendering_helper::Gl_data_names name)
{
  Q_FOREACH(CGAL::Three::Vbo* vbo, d->VBOs)
  {
    if(!vbo)
      continue;
    if(
       (name.testFlag(Scene_item_rendering_helper::GEOMETRY) && vbo->flag == Vbo::GEOMETRY)
       || (name.testFlag(Scene_item_rendering_helper::NORMALS) && vbo->flag == Vbo::NORMALS)
       || (name.testFlag(Scene_item_rendering_helper::COLORS) && vbo->flag == Vbo::COLORS)
       || (name.testFlag(Scene_item_rendering_helper::NOT_INSTANCED) && vbo->flag == Vbo::NOT_INSTANCED)){
       vbo->allocated = false;
       vbo->dataSize=0;
    }
  }
}

void Primitive_container::allocate(std::size_t vbo_id, void* data, int datasize)
{
  d->VBOs[vbo_id]->allocate(data, datasize);
}

void Primitive_container::setVao(Viewer_interface* viewer, Vao* vao)
{
  d->VAOs[viewer]=vao;
}

void Primitive_container::setVbos(std::vector<Vbo*> vbos)
{
  d->VBOs = vbos;
}

void Primitive_container::setInit(Viewer_interface* viewer, bool b)
{
  d->is_init[viewer] = b;
}

void Primitive_container::setGLInit(Viewer_interface* viewer, bool b)
{
  d->is_gl_init[viewer] = b;
}

void Primitive_container::setSelected(bool b)
{
  d->is_selected = b;
}

void Primitive_container::setFlatDataSize(std::size_t s)
{
  d->flat_size = s;
}

void Primitive_container::setCenterSize(std::size_t s)
{
  d->center_size = s;
}

void Primitive_container::setIdxSize(std::size_t s)
{
   d->idx_size = s;
}

void Primitive_container::setColor(QColor c)
{
  d->color = c;
}

Vao* Primitive_container::getVao(Viewer_interface* viewer) const  {

  return d->VAOs[viewer];
}

Vbo* Primitive_container::getVbo(std::size_t id) const
{

  return d->VBOs[id];
}

void Primitive_container::setVbo(std::size_t vbo_id, Vbo* vbo) { d->VBOs[vbo_id] = vbo; }

void Primitive_container::setStride(std::size_t id, int stride) {
  d->VBOs[id]->stride = stride;
}

void Primitive_container::setOffset(std::size_t id, int offset) {
  d->VBOs[id]->offset = offset;
}

void Primitive_container::setTupleSize(int ts) { d->tuple_size = ts; }
int Primitive_container::getProgram() const  { return d->program_id; }

bool Primitive_container::isDataIndexed()  { return d->indexed; }

bool Primitive_container::isInit(Viewer_interface* viewer)const  {
  return d->is_init[viewer];
}

bool Primitive_container::isGLInit(Viewer_interface* viewer)const  {
  return d->is_gl_init[viewer];
}

bool Primitive_container::isSelected() const  { return d->is_selected; }

std::size_t Primitive_container::getFlatDataSize() const  {return d->flat_size; }

std::size_t Primitive_container::getCenterSize() const  { return d->center_size; }

std::size_t Primitive_container::getIdxSize() const  { return d->idx_size; }

QColor Primitive_container::getColor()const  { return d->color; }

int Primitive_container::getTupleSize() const { return d->tuple_size; }

Texture* Primitive_container::getTexture() const { return d->texture; }

void Primitive_container::setTexture(Texture * t) { d->texture = t; }

GLuint Primitive_container::getTextureId() const { return d->textureId; }

void Primitive_container::setTextureSize(const QSize &size) {
  if(!d->texture)
    d->texture = new Texture();
  d->texture->Width = size.width();
  d->texture->Height = size.height();
  d->texture->size = 3 * size.width()*size.height();
  d->texture->data = new GLubyte[d->texture->size];
}

void Primitive_container::setTextureData(int i, int j, int r, int g, int b)
{
  d->texture->setData(i,j,r,g,b);
}

QSize Primitive_container::getTextureSize() const
{
  return QSize(d->texture->Width, d->texture->Height);
}

bool Primitive_container::getClipping() const
{
  return d->clipping;
}

void Primitive_container::setClipping(bool b)
{
  d->clipping = b;
}
