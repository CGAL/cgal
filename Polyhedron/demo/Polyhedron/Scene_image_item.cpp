#include "config.h"

#include "Scene_image_item.h"
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include "Image_type.h"
#include <QColor>
#include <map>
#include <CGAL/ImageIO.h>
#include <CGAL/use.h>
#include <QApplication>

// -----------------------------------
// Internal classes
// -----------------------------------
using namespace  CGAL::Three;
typedef Triangle_container Tc;
typedef Edge_container Ec;
typedef Viewer_interface Vi;
namespace internal {

class Image_accessor
{
public:
  Image_accessor(const Image& im, int dx=1, int dy=1, int dz=1);

  bool is_vertex_active(std::size_t i, std::size_t j, std::size_t k) const;
  const QColor& vertex_color(std::size_t i, std::size_t j, std::size_t k) const;
  void normal(std::size_t i, std::size_t j, std::size_t k,
              float& x, float& y, float& z) const;

  int dx() const { return dx_; }
  int dy() const { return dy_; }
  int dz() const { return dz_; }
  std::size_t xdim() const { return im_.xdim(); }
  std::size_t ydim() const { return im_.ydim(); }
  std::size_t zdim() const { return im_.zdim(); }
  double vx() const { return im_.vx(); }
  double vy() const { return im_.vy(); }
  double vz() const { return im_.vz(); }

  double tx() const { return im_.image()->tx; }
  double ty() const { return im_.image()->ty; }
  double tz() const { return im_.image()->tz; }

private:
  unsigned char non_null_neighbor_data(std::size_t i,
                                       std::size_t j,
                                       std::size_t k) const;

  unsigned char image_data(std::size_t i, std::size_t j, std::size_t k) const;

  void add_to_normal(unsigned char v,
                     float& x, float& y, float& z,
                     int dx, int dy, int dz) const;

private:
  const Image& im_;
  int dx_, dy_, dz_;
  const QColor default_color_;
  std::map<unsigned char, QColor> colors_;
};


Image_accessor::Image_accessor(const Image& im, int dx, int dy, int dz)
: im_(im)
, dx_(dx)
, dy_(dy)
, dz_(dz)
, default_color_()
, colors_()
{
  const std::size_t xdim = im_.xdim();
  const std::size_t ydim = im_.ydim();
  const std::size_t zdim = im_.zdim();

  for(std::size_t i=0 ; i<xdim ; i+=dx_)
  {
    for(std::size_t j=0 ; j<ydim ; j+=dy_)
    {
      for(std::size_t k=0 ; k<zdim ; k+=dz_)
      {
        unsigned char c = image_data(i,j,k);
        if ( 0 != c ) { colors_.insert(std::make_pair(c,QColor())); }
      }
    }
  }

  double i=0;
  const double starting_hue = 45./360.; // magenta
  for ( std::map<unsigned char, QColor>::iterator it = colors_.begin(),
       end = colors_.end() ; it != end ; ++it, i += 1.)
  {
    double hue =  starting_hue + 1./double(colors_.size()) * i;
    if ( hue > 1. ) { hue -= 1.; }
    it->second = QColor::fromHsvF(hue, .75, .75);
  }
}

bool
Image_accessor::
is_vertex_active(std::size_t i, std::size_t j, std::size_t k) const
{
  unsigned char v1 = image_data(i-dx_, j-dy_, k-dz_);
  unsigned char v2 = image_data(i-dx_, j-dy_, k  );
  unsigned char v3 = image_data(i-dx_, j    , k-dz_);
  unsigned char v4 = image_data(i-dx_, j    , k  );
  unsigned char v5 = image_data(i    , j-dy_, k-dz_);
  unsigned char v6 = image_data(i    , j-dy_, k  );
  unsigned char v7 = image_data(i    , j    , k-dz_);
  unsigned char v8 = image_data(i    , j    , k  );

  // don't draw interior vertices
  if ( v1 != 0 && v2 != 0 && v3 != 0 && v4 != 0 &&
       v5 != 0 && v6 != 0 && v7 != 0 && v8 != 0 )
  {
    return false;
  }

  return ( v1 != 0 || v2 != 0 || v3 != 0 || v4 != 0 ||
           v5 != 0 || v6 != 0 || v7 != 0 || v8 != 0 );
}

const QColor&
Image_accessor::vertex_color(std::size_t i, std::size_t j, std::size_t k) const
{
  unsigned char c = non_null_neighbor_data(i,j,k);
  if ( 0 == c ) { return default_color_; }

  std::map<unsigned char, QColor>::const_iterator color = colors_.find(c);
  if ( colors_.end() == color ) { return default_color_; }

  return color->second;
}

unsigned char
Image_accessor::image_data(std::size_t i, std::size_t j, std::size_t k) const
{
  if ( i<im_.xdim() && j<im_.ydim() && k<im_.zdim() )
    return CGAL::IMAGEIO::static_evaluate<unsigned char>(im_.image(),i,j,k);
  else
    return 0;
}

unsigned char
Image_accessor::
non_null_neighbor_data(std::size_t i, std::size_t j, std::size_t k) const
{
  unsigned char v = image_data(i-dx_, j-dy_, k-dz_);
  if ( v != 0 ) { return v; }

  v = image_data(i-dx_, j-dy_, k  );
  if ( v != 0 ) { return v; }

  v = image_data(i-dx_, j    , k-dz_);
  if ( v != 0 ) { return v; }

  v = image_data(i-dx_, j    , k  );
  if ( v != 0 ) { return v; }

  v = image_data(i    , j-dy_, k-dz_);
  if ( v != 0 ) { return v; }

  v = image_data(i    , j-dy_, k  );
  if ( v != 0 ) { return v; }

  v = image_data(i    , j    , k-dz_);
  if ( v != 0 ) { return v; }

  v = image_data(i    , j    , k  );
  if ( v != 0 ) { return v; }

  return 0;
}

void
Image_accessor::
normal(std::size_t i, std::size_t j, std::size_t k,
       float& x, float& y, float& z) const
{
  unsigned char v = image_data(i-dx_, j-dy_, k-dz_);
  add_to_normal(v,x,y,z,       1    , 1    , 1);

  v = image_data(        i-dx_, j-dy_, k);
  add_to_normal(v,x,y,z, 1    , 1    , -1);

  v = image_data(        i-dx_, j    , k-dz_);
  add_to_normal(v,x,y,z, 1    , -1   , 1);

  v = image_data(        i-dx_, j    , k  );
  add_to_normal(v,x,y,z, 1    , -1   , -1);

  v = image_data(        i    , j-dy_, k-dz_);
  add_to_normal(v,x,y,z, -1   , 1    , 1);

  v = image_data(        i    , j-dy_, k  );
  add_to_normal(v,x,y,z, -1   , 1    , -1);

  v = image_data(        i    , j    , k-dz_);
  add_to_normal(v,x,y,z, -1   , -1   , 1);

  v = image_data(        i    , j    , k);
  add_to_normal(v,x,y,z, -1   , -1   , -1);
}

void
Image_accessor::
add_to_normal(unsigned char v,
              float& x, float& y, float& z,
              int dx, int dy, int dz) const
{
  if ( 0 != v )
  {
    x += float(dx);
    y += float(dy);
    z += float(dz);
  }
}



class Vertex_buffer_helper
{
public:
  Vertex_buffer_helper(const Image_accessor& data, bool is_ogl_4_3);

  void fill_buffer_data();

  float* colors() { return colors_.data(); }
  float* normals() { return normals_.data(); }
  float* vertices() { return vertices_.data();}
  unsigned int* quads() { return quads_.data(); }

  std::size_t color_size() const { return colors_.size(); }
  std::size_t normal_size() const { return normals_.size(); }
  std::size_t vertex_size() const { return vertices_.size(); }
  std::size_t quad_size() const { return quads_.size(); }

private:
  void treat_vertex(std::size_t i, std::size_t j, std::size_t k);

  void push_color(std::size_t i, std::size_t j, std::size_t k);
  void push_normal(std::size_t i, std::size_t j, std::size_t k);
  void push_vertex(std::size_t i, std::size_t j, std::size_t k);
  void push_quads(std::size_t i, std::size_t j, std::size_t k);
  void push_quad(int pos1, int pos2, int pos3, int pos4);

  int compute_position(std::size_t i, std::size_t j, std::size_t k) const;
  int vertex_index(std::size_t i, std::size_t j, std::size_t k) const;

  int dx() const { return data_.dx(); }
  int dy() const { return data_.dy(); }
  int dz() const { return data_.dz(); }

private:
  static int vertex_not_found_;

  const Image_accessor& data_;
  typedef std::map<int, std::size_t> Indices;
  Indices indices_;
  std::vector<float> colors_, normals_, vertices_;
  std::vector<unsigned int> quads_;
  bool is_ogl_4_3;
};

int Vertex_buffer_helper::vertex_not_found_ = -1;

Vertex_buffer_helper::
Vertex_buffer_helper(const Image_accessor& data, bool b)
  : data_(data), is_ogl_4_3(b)
{}


void
Vertex_buffer_helper::
fill_buffer_data()
{
  std::size_t i,j,k;
  for ( i = 0 ; i <= data_.xdim() ; i+=dx() )
  {
    for ( j = 0 ; j <= data_.ydim() ; j+=dy() )
    {
      for ( k = 0 ; k <= data_.zdim() ; k+=dz() )
      {
        treat_vertex(i,j,k);
      }
    }
  }
}

void
Vertex_buffer_helper::treat_vertex(std::size_t i, std::size_t j, std::size_t k)
{
  if ( data_.is_vertex_active(i,j,k) )
  {
    push_vertex(i,j,k);
    push_color(i,j,k);
    push_normal(i,j,k);
    push_quads(i,j,k);
  }
}

void
Vertex_buffer_helper::push_color(std::size_t i, std::size_t j, std::size_t k)
{
  const QColor& color = data_.vertex_color(i,j,k);
  if ( ! color.isValid() ) { return; }
  colors_.push_back(float(color.red())/255.f);
  colors_.push_back(float(color.green())/255.f);
  colors_.push_back(float(color.blue())/255.f);
}

void
Vertex_buffer_helper::push_normal(std::size_t i, std::size_t j, std::size_t k)
{
  float x=0.f, y=0.f, z=0.f;
  data_.normal(i,j,k,x,y,z);

  float norm = std::sqrt(x*x+y*y+z*z);
  x = x / norm;
  y = y / norm;
  z = z / norm;

  normals_.push_back(x);
  normals_.push_back(y);
  normals_.push_back(z);
}

void
Vertex_buffer_helper::push_vertex(std::size_t i, std::size_t j, std::size_t k)
{
  indices_.insert(std::make_pair(compute_position(i,j,k),
                                 vertices_.size()/3));
  //resize the "border vertices"
  double di = double(i),dj = double(j),dk = double(k);
  if (di == 0)
    di = 0.5;
  if (dj == 0)
    dj = 0.5;
  if (dk == 0)
    dk = 0.5;
  if (di == double(data_.xdim()))
    di = double(data_.xdim())-0.5;
  if (dj == double(data_.ydim()))
    dj = double(data_.ydim())-0.5;
  if (dk == double(data_.zdim()))
    dk = double(data_.zdim())-0.5;

  vertices_.push_back( (di - 0.5) * data_.vx() + data_.tx());
  vertices_.push_back( (dj - 0.5) * data_.vy() + data_.ty());
  vertices_.push_back( (dk - 0.5) * data_.vz() + data_.tz());
}

void
Vertex_buffer_helper::push_quads(std::size_t i, std::size_t j, std::size_t k)
{
  int pos1 = vertex_index(i-dx(), j     , k);
  int pos2 = vertex_index(i-dx(), j-dy(), k);
  int pos3 = vertex_index(i     , j-dy(), k) ;
  int pos4 = vertex_index(i     ,j      , k);
  push_quad(pos1, pos2, pos3, pos4);

  pos1 = vertex_index(i-dx(), j, k);
  pos2 = vertex_index(i-dx(), j, k-dz());
  pos3 = vertex_index(i     , j, k-dz());
  push_quad(pos1, pos2, pos3, pos4);

  pos1 = vertex_index(i, j-dy(), k);
  pos2 = vertex_index(i, j-dy(), k-dz());
  pos3 = vertex_index(i, j     , k-dz());
  push_quad(pos1, pos2, pos3, pos4);
}

void
Vertex_buffer_helper::push_quad(int pos1, int pos2, int pos3, int pos4)
{
  if (   pos1 != vertex_not_found_
      && pos2 != vertex_not_found_
      && pos3 != vertex_not_found_ )
  {
    quads_.push_back(pos1);
    quads_.push_back(pos2);
    quads_.push_back(pos3);
    if(!is_ogl_4_3)
    {
      quads_.push_back(pos1);
      quads_.push_back(pos3);
    }
    quads_.push_back(pos4);
  }
}

int
Vertex_buffer_helper::
compute_position(std::size_t i, std::size_t j, std::size_t k) const
{
  return  static_cast<int>(
    i/dx() * (data_.ydim()/dy()+1) * (data_.zdim()/dz()+1)
         + j/dy() * (data_.zdim()/dz()+1)
         + k/dz());
}

int
Vertex_buffer_helper::
vertex_index(std::size_t i, std::size_t j, std::size_t k) const
{
  if ( i > data_.xdim() || j > data_.ydim() || k > data_.zdim() )
  {
    return vertex_not_found_;
  }

  int vertex_key = compute_position(i,j,k);
  Indices::const_iterator it = indices_.find(vertex_key);
  if ( it != indices_.end() )
  {
    return static_cast<int>(it->second);
  }

  return vertex_not_found_;
}


} // namespace internal

struct Scene_image_item_priv
{

  Scene_image_item_priv(int display_scale, bool hidden, Scene_image_item* parent)
    : m_initialized(false)
    , m_voxel_scale(display_scale)
  {
    item = parent;
    is_ogl_4_3 = Three::mainViewer()->isOpenGL_4_3();
    is_hidden = hidden;
    box_size = 0;
    idx_size = 0;
    helper = nullptr;
  }

  void draw_bbox();
  void draw_Bbox(Scene_item::Bbox bbox, std::vector<float> *vertices);

  bool m_initialized;
//#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  int m_voxel_scale;
  std::vector<float> v_box;
  std::size_t idx_size, box_size;
  bool is_hidden;
  bool is_ogl_4_3;
  Scene_image_item* item;
  internal::Vertex_buffer_helper* helper;

//#endif // SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
};
// -----------------------------------
// Scene_image_item
// -----------------------------------
Scene_image_item::Scene_image_item(Image* im, int display_scale, bool hidden = false)
  : m_image(im)
{
  CGAL_USE(display_scale);
  d = new Scene_image_item_priv(display_scale, hidden, this);
  setRenderingMode(Flat);
  setTriangleContainer(0,new Tc(
                         Three::mainViewer()->isOpenGL_4_3()
                         ? Vi::PROGRAM_NO_INTERPOLATION
                         : Vi::PROGRAM_WITH_LIGHT,
                         true));
  setEdgeContainer(0, new Ec(
                     Three::mainViewer()->isOpenGL_4_3()
                     ? Vi::PROGRAM_SOLID_WIREFRAME
                     : Vi::PROGRAM_NO_SELECTION,
                     false));

}


Scene_image_item::~Scene_image_item()
{
   delete d;
  delete m_image;
}

void
Scene_image_item::compute_bbox() const
{
  if(!m_image)
    setBbox(Bbox());
  else
   setBbox(Bbox(m_image->image()->tx,
                m_image->image()->ty,
                m_image->image()->tz,
                m_image->image()->tx+(m_image->xdim()-1) * m_image->vx(),
                m_image->image()->ty+(m_image->ydim()-1) * m_image->vy(),
                m_image->image()->tz+(m_image->zdim()-1) * m_image->vz()));
}

void
Scene_image_item::draw(Viewer_interface* viewer) const
{
  if(m_image)
  {
    if(!isInit(viewer))
      initGL(viewer);
    if ( getBuffersFilled() &&
         ! getBuffersInit(viewer))
    {
      initializeBuffers(viewer);
      setBuffersInit(viewer, true);
    }
    if(!getBuffersFilled())
    {
      computeElements();
      initializeBuffers(viewer);
    }
    if(!d->is_hidden)
    {
      getTriangleContainer(0)->draw(viewer, false);
    }
  }
  drawEdges(viewer);
}

template <typename T> const char* whatType(T) { return "unknown"; }    // default
template <> const char* whatType(float) { return "float"; }
template <> const char* whatType(double) { return "double"; }
template <> const char* whatType(char) { return "int8_t (char)"; }
template <> const char* whatType(boost::uint8_t) { return "uint8_t (unsigned char)"; }
template <> const char* whatType(boost::int16_t) { return "int16_t (short)"; }
template <> const char* whatType(boost::uint16_t) { return "uint16_t (unsigned short)"; }
template <> const char* whatType(boost::int32_t) { return "int32_t (int)"; }
template <> const char* whatType(boost::uint32_t) { return "uint32_t (unsigned int)"; }

template<typename Word>
QString explicitWordType()
{
  return QString(whatType(Word(0)));
}

QString
Scene_image_item::toolTip() const
{
  QString w_type;
  CGAL_IMAGE_IO_CASE(image()->image(), w_type = explicitWordType<Word>())
  return tr("<p>Image <b>%1</b></p>"
            "<p>Word type: %2</p>"
            "<p>Dimensions: %3 x %4 x %5</p>"
            "<p>Spacings: ( %6 , %7 , %8 )</p>")
    .arg(this->name())
    .arg(w_type)
    .arg(m_image->xdim())
    .arg(m_image->ydim())
    .arg(m_image->zdim())
    .arg(m_image->vx())
    .arg(m_image->vy())
    .arg(m_image->vz());
}

bool
Scene_image_item::supportsRenderingMode(RenderingMode m) const
{
  switch ( m )
  {
  case Wireframe:
  case Flat:
  case FlatPlusEdges:
    return true;

  default:
    return false;
  }

  return false;
}


void Scene_image_item_priv::draw_Bbox(Scene_item::Bbox bbox, std::vector<float> *vertices)
{
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmin()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmin()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymax()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

    vertices->push_back(bbox.xmax()+offset.x);
    vertices->push_back(bbox.ymin()+offset.y);
    vertices->push_back(bbox.zmax()+offset.z);

}

void Scene_image_item::drawEdges(Viewer_interface *viewer) const
{
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }
  if(viewer->isOpenGL_4_3())
  {
    QVector2D vp(viewer->width(), viewer->height());
    getEdgeContainer(0)->setViewport(vp);
    getEdgeContainer(0)->setWidth(3.0f);
  }
  getEdgeContainer(0)->setColor(QColor(Qt::black));
  viewer->glDepthRangef(0.00001f, 0.99999f);
  getEdgeContainer(0)->draw(viewer, true);
  viewer->glDepthRangef(0.0f, 1.0f);
}

bool Scene_image_item::isGray() { return d->is_hidden;}

void Scene_image_item::invalidateOpenGLBuffers()
{
  setBuffersFilled(false);
  getTriangleContainer(0)->reset_vbos(ALL);
  getEdgeContainer(0)->reset_vbos(ALL);
}

void Scene_image_item::initializeBuffers(Viewer_interface *v) const
{
  getTriangleContainer(0)->initializeBuffers(v);
  getTriangleContainer(0)->setIdxSize(d->idx_size);

  getEdgeContainer(0)->initializeBuffers(v);
  getEdgeContainer(0)->setFlatDataSize(d->box_size);

  d->v_box.clear();
  d->v_box.shrink_to_fit();
  if(d->helper)
  {
    delete d->helper;
    d->helper = nullptr;
  }
}

void Scene_image_item::computeElements() const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  d->draw_Bbox(bbox(), &d->v_box);
  if(!d->is_hidden)
  {
    internal::Image_accessor image_data_accessor (*m_image,
                                                  d->m_voxel_scale,
                                                  d->m_voxel_scale,
                                                  d->m_voxel_scale);
    d->helper = new internal::Vertex_buffer_helper(image_data_accessor, d->is_ogl_4_3);
    d->helper->fill_buffer_data();
    getTriangleContainer(0)->allocate(
          Tc::Smooth_vertices,
          d->helper->vertices(), static_cast<int>(d->helper->vertex_size()*sizeof(float)));
    getTriangleContainer(0)->allocate(
          Tc::Smooth_normals,
          d->helper->normals(), static_cast<int>(d->helper->normal_size()*sizeof(float)));

    getTriangleContainer(0)->allocate(
          Tc::VColors,
          d->helper->colors(), static_cast<int>(d->helper->color_size()*sizeof(float)));
    getTriangleContainer(0)->allocate(
          Tc::Vertex_indices,
          d->helper->quads(), static_cast<int>(d->helper->quad_size()*sizeof(unsigned int)));
    d->idx_size = d->helper->quad_size();
  }
  getEdgeContainer(0)->allocate(
        Ec::Vertices,
        d->v_box.data(), static_cast<int>(d->v_box.size()*sizeof(float)));
  d->box_size = d->v_box.size();
  setBuffersFilled(true);
  QApplication::restoreOverrideCursor();
}
