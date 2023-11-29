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

template <typename Word_type>
class Image_accessor
{
public:
  Image_accessor(const Image& im, int dx=1, int dy=1, int dz=1, const QColor& default_color = QColor());

  bool is_vertex_active(std::size_t i, std::size_t j, std::size_t k) const;
  const QColor& vertex_color(std::size_t i, std::size_t j, std::size_t k) const;
  void normal(std::size_t i, std::size_t j, std::size_t k,
              float& x, float& y, float& z) const;

  int dx() const { return dx_; }
  int dy() const { return dy_; }
  int dz() const { return dz_; }
  std::size_t xdim() const { return im_->xdim(); }
  std::size_t ydim() const { return im_->ydim(); }
  std::size_t zdim() const { return im_->zdim(); }
  double vx() const { return im_->vx(); }
  double vy() const { return im_->vy(); }
  double vz() const { return im_->vz(); }

  double tx() const { return im_->image()->tx; }
  double ty() const { return im_->image()->ty; }
  double tz() const { return im_->image()->tz; }

private:
  Word_type non_null_neighbor_data(std::size_t i,
                                       std::size_t j,
                                       std::size_t k) const;

  Word_type image_data(std::size_t i, std::size_t j, std::size_t k) const;

  void add_to_normal(Word_type v,
                     float& x, float& y, float& z,
                     int dx, int dy, int dz) const;

private:
  const Image* im_;
  int dx_, dy_, dz_;
  QColor default_color_;
  std::map<Word_type, QColor> colors_;
};

template <typename Word_type>
Image_accessor<Word_type>::Image_accessor(const Image& im, int dx, int dy, int dz, const QColor& default_color)
: im_(&im)
, dx_(dx)
, dy_(dy)
, dz_(dz)
, default_color_(default_color)
, colors_()
{
  const std::size_t xdim = im_->xdim();
  const std::size_t ydim = im_->ydim();
  const std::size_t zdim = im_->zdim();

  for(std::size_t i=0 ; i<xdim ; i+=dx_)
  {
    for(std::size_t j=0 ; j<ydim ; j+=dy_)
    {
      for(std::size_t k=0 ; k<zdim ; k+=dz_)
      {
        Word_type c = image_data(i,j,k);
        if ( 0 != c ) { colors_.insert(std::make_pair(c,QColor())); }
      }
    }
  }

  const double nb_Colors = colors_.size()+1;
  double i=0;
  const double starting_hue = default_color.hueF();
  for ( auto it = colors_.begin(),
       end = colors_.end() ; it != end ; ++it, i += 1.)
  {
    double hue =  starting_hue + 1./nb_Colors * i;
    if ( hue > 1. ) { hue -= 1.; }
    it->second = QColor::fromHsvF(hue, default_color.saturationF(), default_color.valueF());
  }
}

template <typename Word_type>
bool
Image_accessor<Word_type>::
is_vertex_active(std::size_t i, std::size_t j, std::size_t k) const
{
  Word_type v1 = image_data(i-dx_, j-dy_, k-dz_);
  Word_type v2 = image_data(i-dx_, j-dy_, k  );
  Word_type v3 = image_data(i-dx_, j    , k-dz_);
  Word_type v4 = image_data(i-dx_, j    , k  );
  Word_type v5 = image_data(i    , j-dy_, k-dz_);
  Word_type v6 = image_data(i    , j-dy_, k  );
  Word_type v7 = image_data(i    , j    , k-dz_);
  Word_type v8 = image_data(i    , j    , k  );

  return v1 != v2 || v1 != v3 || v1 != v4 || v1 != v5 || v1 != v6 || v1 != v7 || v1 != v8;
}
template <typename Word_type>
const QColor&
Image_accessor<Word_type>::vertex_color(std::size_t i, std::size_t j, std::size_t k) const
{
  Word_type c = non_null_neighbor_data(i,j,k);
  if ( 0 == c ) { return default_color_; }

  auto color = colors_.find(c);
  if ( colors_.end() == color ) { return default_color_; }

  return color->second;
}

template <typename Word_type>
Word_type
Image_accessor<Word_type>::image_data(std::size_t i, std::size_t j, std::size_t k) const
{
  if ( i<im_->xdim() && j<im_->ydim() && k<im_->zdim() )
    return CGAL::IMAGEIO::static_evaluate<Word_type>(im_->image(),i,j,k);
  else
    return 0;
}

template <typename Word_type>
Word_type
Image_accessor<Word_type>::
non_null_neighbor_data(std::size_t i, std::size_t j, std::size_t k) const
{
  Word_type v = image_data(i-dx_, j-dy_, k-dz_);
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

template <typename Word_type>
void
Image_accessor<Word_type>::
normal(std::size_t i, std::size_t j, std::size_t k,
       float& x, float& y, float& z) const
{
  Word_type v = image_data(i-dx_, j-dy_, k-dz_);
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

template <typename Word_type>
void
Image_accessor<Word_type>::
add_to_normal(Word_type v,
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
  virtual void fill_buffer_data() = 0;
  virtual ~Vertex_buffer_helper() {}

  float* colors() { return colors_.data(); }
  float* normals() { return normals_.data(); }
  float* vertices() { return vertices_.data();}
  unsigned int* quads() { return quads_.data(); }

  std::size_t color_size() const { return colors_.size(); }
  std::size_t normal_size() const { return normals_.size(); }
  std::size_t vertex_size() const { return vertices_.size(); }
  std::size_t quad_size() const { return quads_.size(); }
protected:
  std::vector<float> colors_, normals_, vertices_;
  std::vector<unsigned int> quads_;
};

template <typename Word_type>
class Vertex_buffer_helper_impl : public Vertex_buffer_helper
{
public:
  Vertex_buffer_helper_impl(Image_accessor<Word_type>&& data, bool is_ogl_4_3);

  ~Vertex_buffer_helper_impl() override {}
  void fill_buffer_data() override;

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
  static constexpr int vertex_not_found_ = -1;

  const Image_accessor<Word_type> data_;
  typedef std::map<int, std::size_t> Indices;
  Indices indices_;
  bool is_ogl_4_3;
};

// template <typename Word_type>
// int Vertex_buffer_helper_impl<Word_type>::vertex_not_found_ = -1;


template <typename Word_type>
Vertex_buffer_helper_impl<Word_type>::
Vertex_buffer_helper_impl(Image_accessor<Word_type>&& data, bool b)
  : data_(std::move(data)), is_ogl_4_3(b)
{}


template <typename Word_type>
void
Vertex_buffer_helper_impl<Word_type>::
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

template <typename Word_type>
void
Vertex_buffer_helper_impl<Word_type>::treat_vertex(std::size_t i, std::size_t j, std::size_t k)
{
  if ( data_.is_vertex_active(i,j,k) )
  {
    push_vertex(i,j,k);
    push_color(i,j,k);
    push_normal(i,j,k);
    push_quads(i,j,k);
  }
}

template <typename Word_type>
void
Vertex_buffer_helper_impl<Word_type>::push_color(std::size_t i, std::size_t j, std::size_t k)
{
  const QColor& color = data_.vertex_color(i,j,k);
  if ( ! color.isValid() ) { return; }
  colors_.push_back(float(color.red())/255.f);
  colors_.push_back(float(color.green())/255.f);
  colors_.push_back(float(color.blue())/255.f);
}

template <typename Word_type>
void
Vertex_buffer_helper_impl<Word_type>::push_normal(std::size_t i, std::size_t j, std::size_t k)
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

template <typename Word_type>
void
Vertex_buffer_helper_impl<Word_type>::push_vertex(std::size_t i, std::size_t j, std::size_t k)
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

template <typename Word_type>
void
Vertex_buffer_helper_impl<Word_type>::push_quads(std::size_t i, std::size_t j, std::size_t k)
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

template <typename Word_type>
void
Vertex_buffer_helper_impl<Word_type>::push_quad(int pos1, int pos2, int pos3, int pos4)
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

template <typename Word_type>
int
Vertex_buffer_helper_impl<Word_type>::
compute_position(std::size_t i, std::size_t j, std::size_t k) const
{
  return  static_cast<int>(
    i/dx() * (data_.ydim()/dy()+1) * (data_.zdim()/dz()+1)
         + j/dy() * (data_.zdim()/dz()+1)
         + k/dz());
}

template <typename Word_type>
int
Vertex_buffer_helper_impl<Word_type>::
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
    , m_image_weights()
    , m_sigma_weights(0.f)
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
  Image m_image_weights;
  float m_sigma_weights;

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

Image*
Scene_image_item::image_weights() const
{
  return &d->m_image_weights;
}
void
Scene_image_item::set_image_weights(const Image& img, const float sigma)
{
  d->m_image_weights = img;
  d->m_sigma_weights = sigma;
}
float
Scene_image_item::sigma_weights() const
{
  return d->m_sigma_weights;
}
float
Scene_image_item::default_sigma_weights() const
{
  if(!m_image)
    return 0.f;
  else
    return (std::max)(m_image->image()->vx,
           (std::max)(m_image->image()->vy,
                      m_image->image()->vz));
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
template <> const char* whatType(std::uint8_t) { return "uint8_t (unsigned char)"; }
template <> const char* whatType(std::int16_t) { return "int16_t (short)"; }
template <> const char* whatType(std::uint16_t) { return "uint16_t (unsigned short)"; }
template <> const char* whatType(std::int32_t) { return "int32_t (int)"; }
template <> const char* whatType(std::uint32_t) { return "uint32_t (unsigned int)"; }

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
  getEdgeContainer(0)->setClipping(false);
  viewer->glDepthRangef(0.00001f, 0.99999f);
  getEdgeContainer(0)->draw(viewer, true);
  viewer->glDepthRangef(0.0f, 1.0f);
}

bool Scene_image_item::isGray() { return d->is_hidden;}

void Scene_image_item::setColor(QColor c) {
  color_ = c;
  invalidateOpenGLBuffers();
}

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

template <typename Word_type>
internal::Vertex_buffer_helper*
init_helper(const Image &im, int dx, int dy, int dz, bool is_ogl_4_3, const QColor& default_color = QColor())
{
  internal::Image_accessor<Word_type> image_data_accessor(im, dx, dy, dz, default_color);
  return new internal::Vertex_buffer_helper_impl<Word_type>(std::move(image_data_accessor),
                                                            is_ogl_4_3);
}

void Scene_image_item::computeElements() const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  d->draw_Bbox(bbox(), &d->v_box);
  if(!d->is_hidden)
  {
    // internal::Image_accessor image_data_accessor (*m_image,
    //                                               d->m_voxel_scale,
    //                                               d->m_voxel_scale,
    //                                               d->m_voxel_scale);
    // d->helper = new internal::Vertex_buffer_helper(image_data_accessor, d->is_ogl_4_3);
    CGAL_IMAGE_IO_CASE(m_image->image(),
                       d->helper = init_helper<Word>(*m_image,
                                                     d->m_voxel_scale, d->m_voxel_scale, d->m_voxel_scale,
                                                     d->is_ogl_4_3,
                                                     color_));
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
