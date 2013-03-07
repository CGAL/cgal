#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
#  include <GL/glew.h>
#endif
#include "Scene_segmented_image_item.h"
#include "Image_type.h"
#include <QColor>
#include <map>
#include <CGAL/gl.h>
#include <CGAL/ImageIO.h>

//#define SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE

#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
bool gl_vbo_available() {
  return  glewIsSupported("GL_VERSION_1_4");
}
#else
bool gl_vbo_available() {
  return false;
}
#endif



// -----------------------------------
// Internal classes
// -----------------------------------
namespace internal {

class Image_accessor
{
public:
  Image_accessor(const Image& im, int dx=1, int dy=1, int dz=1);
  
  bool is_vertex_active(unsigned int i, unsigned int j, unsigned int k) const;
  const QColor& vertex_color(unsigned int i, unsigned int j, unsigned int k) const;
  void normal(unsigned int i, unsigned int j, unsigned int k,
              float& x, float& y, float& z) const;
  
  int dx() const { return dx_; }
  int dy() const { return dy_; }
  int dz() const { return dz_; }
  unsigned int xdim() const { return im_.xdim(); }
  unsigned int ydim() const { return im_.ydim(); }
  unsigned int zdim() const { return im_.zdim(); }
  double vx() const { return im_.vx(); }
  double vy() const { return im_.vy(); }
  double vz() const { return im_.vz(); }
  
private:
  unsigned char non_null_neighbor_data(unsigned int i,
                                       unsigned int j,
                                       unsigned int k) const;
  
  unsigned char image_data(unsigned int i, unsigned int j, unsigned int k) const;
  
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
  const unsigned int xdim = im_.xdim();
  const unsigned int ydim = im_.ydim();
  const unsigned int zdim = im_.zdim();  
  
  for(unsigned int i=0 ; i<xdim ; i+=dx_)
  { 
    for(unsigned int j=0 ; j<ydim ; j+=dy_)
    { 
      for(unsigned int k=0 ; k<zdim ; k+=dz_)
      {
        unsigned char c = image_data(i,j,k);
        if ( 0 != c ) { colors_.insert(std::make_pair(c,QColor())); }
      }
    }
  }
  
  int i=0;
  const double starting_hue = 45./360.; // magenta
  for ( std::map<unsigned char, QColor>::iterator it = colors_.begin(),
       end = colors_.end() ; it != end ; ++it, ++i )
  {
    double hue =  starting_hue + 1./colors_.size() * i;
    if ( hue > 1. ) { hue -= 1.; }
    it->second = QColor::fromHsvF(hue, .75, .75);
  }
}

bool 
Image_accessor::
is_vertex_active(unsigned int i, unsigned int j, unsigned int k) const
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
Image_accessor::vertex_color(unsigned int i, unsigned int j, unsigned int k) const
{
  unsigned char c = non_null_neighbor_data(i,j,k);
  if ( 0 == c ) { return default_color_; }
  
  std::map<unsigned char, QColor>::const_iterator color = colors_.find(c);
  if ( colors_.end() == color ) { return default_color_; }
  
  return color->second;
}

unsigned char
Image_accessor::image_data(unsigned int i, unsigned int j, unsigned int k) const
{
  if ( i<im_.xdim() && j<im_.ydim() && k<im_.zdim() )
    return CGAL::IMAGEIO::static_evaluate<unsigned char>(im_.image(),i,j,k);
  else
    return 0;
}

unsigned char
Image_accessor::
non_null_neighbor_data(unsigned int i, unsigned int j, unsigned int k) const
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
normal(unsigned int i, unsigned int j, unsigned int k,
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
    x += dx;
    y += dy;
    z += dz;    
  }
}



class Vertex_buffer_helper
{
public:
  Vertex_buffer_helper(const Image_accessor& data);
  ~Vertex_buffer_helper();
  
  void fill_buffer_data();

  const GLfloat* colors() const { return color_array_; }
  const GLfloat* normals() const { return normal_array_; }
  const GLfloat* vertices() const { return vertex_array_; }
  const GLuint* quads() const { return quad_array_; }
  
  std::size_t color_size() const { return color_size_*sizeof(GLfloat); }
  std::size_t normal_size() const { return normal_size_*sizeof(GLfloat); }
  std::size_t vertex_size() const { return vertex_size_*sizeof(GLfloat); }
  std::size_t quad_size() const { return quad_size_*sizeof(GLuint); }
  
private:
  void treat_vertex(unsigned int i, unsigned int j, unsigned int k);
  
  void push_color(unsigned int i, unsigned int j, unsigned int k);
  void push_normal(unsigned int i, unsigned int j, unsigned int k);
  void push_vertex(unsigned int i, unsigned int j, unsigned int k);
  void push_quads(unsigned int i, unsigned int j, unsigned int k);
  void push_quad(int pos1, int pos2, int pos3, int pos4);
  
  int compute_position(unsigned int i, unsigned int j, unsigned int k) const;
  int vertex_index(unsigned int i, unsigned int j, unsigned int k) const;
  
  void create_arrays();
  
  template <typename T>
  void create_array(T*& destination, std::size_t& size, const std::vector<T>& source);
  
  int dx() const { return data_.dx(); }
  int dy() const { return data_.dy(); }
  int dz() const { return data_.dz(); }
  
private:
  static int vertex_not_found_;
  
  const Image_accessor& data_;
  typedef std::map<int, std::size_t> Indices;
  Indices indices_;
  std::vector<GLfloat> colors_, normals_, vertices_;
  std::vector<GLuint> quads_;
  
  GLfloat *color_array_, *normal_array_, *vertex_array_;
  GLuint* quad_array_;
  std::size_t color_size_, normal_size_, vertex_size_, quad_size_;
};

int Vertex_buffer_helper::vertex_not_found_ = -1;

Vertex_buffer_helper::
Vertex_buffer_helper(const Image_accessor& data)
: data_(data)
, color_array_(NULL)
, normal_array_(NULL)
, vertex_array_(NULL)
, quad_array_(NULL)
, color_size_(0)
, normal_size_(0)
, vertex_size_(0)
, quad_size_(0)
{
}

Vertex_buffer_helper::
~Vertex_buffer_helper()
{
  delete[] color_array_;
  delete[] normal_array_;
  delete[] vertex_array_;
  delete[] quad_array_;
}

void
Vertex_buffer_helper::
fill_buffer_data()
{
  unsigned int i,j,k;
  
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
  
  create_arrays();
}

void
Vertex_buffer_helper::treat_vertex(unsigned int i, unsigned int j, unsigned int k)
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
Vertex_buffer_helper::push_color(unsigned int i, unsigned int j, unsigned int k)
{
  const QColor& color = data_.vertex_color(i,j,k);
  if ( ! color.isValid() ) { return; }
  
  colors_.push_back(color.red()/255.f);
  colors_.push_back(color.green()/255.f);
  colors_.push_back(color.blue()/255.f);
}

void
Vertex_buffer_helper::push_normal(unsigned int i, unsigned int j, unsigned int k)
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
Vertex_buffer_helper::push_vertex(unsigned int i, unsigned int j, unsigned int k)
{
  indices_.insert(std::make_pair(compute_position(i,j,k),
                                 vertices_.size()/3)); 
  
  vertices_.push_back(i*data_.vx());
  vertices_.push_back(j*data_.vy());
  vertices_.push_back(k*data_.vz());
}

void
Vertex_buffer_helper::push_quads(unsigned int i, unsigned int j, unsigned int k)
{
  int pos1 = vertex_index(i-dx(), j     , k);
  int pos2 = vertex_index(i-dx(), j-dy(), k);
  int pos3 = vertex_index(i     , j-dy(), k);
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
    quads_.push_back(pos4);    
  }
}

int
Vertex_buffer_helper::
compute_position(unsigned int i, unsigned int j, unsigned int k) const
{
  return   i/dx() * (data_.ydim()/dy()+1) * (data_.zdim()/dz()+1)
         + j/dy() * (data_.zdim()/dz()+1)
         + k/dz();
}

int
Vertex_buffer_helper::
vertex_index(unsigned int i, unsigned int j, unsigned int k) const
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


void
Vertex_buffer_helper::
create_arrays()
{
  create_array(color_array_, color_size_, colors_);
  create_array(normal_array_, normal_size_, normals_);
  create_array(vertex_array_, vertex_size_, vertices_);
  create_array(quad_array_, quad_size_, quads_);
}

template <typename T>
void
Vertex_buffer_helper::
create_array(T*& destination, std::size_t& size, const std::vector<T>& source)
{
  size = source.size();
  destination = new T[size];
  
  int i=0;
  for ( typename std::vector<T>::const_iterator it = source.begin(),
       end = source.end() ; it != end ; ++it )
  {
    destination[i++] = *it;
  }
}
  
} // namespace internal



// -----------------------------------
// Scene_segmented_image_item
// -----------------------------------
Scene_segmented_image_item::Scene_segmented_image_item(Image* im,
                                                       int display_scale)
  : m_image(im)
  , m_initialized(false)
  , m_voxel_scale(display_scale)
{
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  if(gl_vbo_available()) {
    ::glGenBuffers(3,m_vbo);
    ::glGenBuffers(1,&m_ibo);
  }
#endif // SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  
  initialize_buffers();
  setRenderingMode(Flat);
}


Scene_segmented_image_item::~Scene_segmented_image_item()
{
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  if(gl_vbo_available()) {
    ::glDeleteBuffers(3,m_vbo);
    ::glDeleteBuffers(1,&m_ibo);
  }
#endif // SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
}


Scene_segmented_image_item::Bbox
Scene_segmented_image_item::bbox() const
{
  if(!m_image) return Bbox();
  return Bbox(0, 0, 0,
              m_image->xdim() * m_image->vx(),
              m_image->ydim() * m_image->vy(),
              m_image->zdim() * m_image->vz());
}


void
Scene_segmented_image_item::draw() const
{
  if(m_image)
  {
    m_image->gl_draw_bbox(3.0f,0,0,0);
    draw_gl();
  }
}


QString
Scene_segmented_image_item::toolTip() const
{
  return tr("<p>Image <b>%1</b></p>"
            "<p>Word type: %2</p>"
            "<p>Dimensions: %3 x %4 x %5</p>"
            "<p>Spacings: ( %6 , %7 , %8 )</p>")
    .arg(this->name())
    .arg("...")
    .arg(m_image->xdim()) 
    .arg(m_image->ydim())
    .arg(m_image->zdim())
    .arg(m_image->vx())
    .arg(m_image->vy())
    .arg(m_image->vz());
}

bool
Scene_segmented_image_item::supportsRenderingMode(RenderingMode m) const
{ 
  switch ( m )
  {
    case Gouraud:
      return false;
      
    case Points:
    case Wireframe:
    case Flat:
    case FlatPlusEdges:
      return true;
      
    default:
      return false;
  }
  
  return false;
}

void
Scene_segmented_image_item::initialize_buffers() 
{
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  if(!gl_vbo_available()) {
    m_initialized = true;
    return;
  }

  internal::Image_accessor image_data_accessor (*m_image,
                                                m_voxel_scale,
                                                m_voxel_scale,
                                                m_voxel_scale);
  
  internal::Vertex_buffer_helper helper (image_data_accessor);
  helper.fill_buffer_data();
  
  // Vertex buffer
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[0]);
  ::glBufferData(GL_ARRAY_BUFFER, helper.vertex_size(), helper.vertices(), GL_STATIC_DRAW);
  
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[1]);
  ::glBufferData(GL_ARRAY_BUFFER, helper.normal_size(), helper.normals(), GL_STATIC_DRAW);
  
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[2]);
  ::glBufferData(GL_ARRAY_BUFFER, helper.color_size(), helper.colors(), GL_STATIC_DRAW);
  
  // Indices buffer
  ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
  ::glBufferData(GL_ELEMENT_ARRAY_BUFFER, helper.quad_size(), helper.quads(), GL_STATIC_DRAW);
  
  // Close buffers
  ::glBindBuffer(GL_ARRAY_BUFFER, 0);
  ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  
#endif // SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE

  m_initialized = true;
}


void
Scene_segmented_image_item::draw_gl() const
{
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  if(!gl_vbo_available()) return;

  ::glShadeModel(GL_SMOOTH);
  
  // Draw faces
  ::glEnableClientState( GL_VERTEX_ARRAY );
  ::glEnableClientState( GL_NORMAL_ARRAY );
  ::glEnableClientState( GL_COLOR_ARRAY );
  
  // Get buffers
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[0]);
  ::glVertexPointer(3, GL_FLOAT, 0, 0);
  
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[1]);
  ::glNormalPointer(GL_FLOAT, 0, 0);
  
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[2]);
  ::glColorPointer(3, GL_FLOAT, 0, 0);
  
  // Render
  ::glDrawElements(GL_QUADS, ibo_size(), GL_UNSIGNED_INT, 0);
  
  // Cleanup
  ::glBindBuffer(GL_ARRAY_BUFFER, 0);
  ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  
  ::glDisableClientState( GL_COLOR_ARRAY );
  ::glDisableClientState( GL_NORMAL_ARRAY );
  ::glDisableClientState( GL_VERTEX_ARRAY );
#endif // SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
}


void
Scene_segmented_image_item::draw_gl_edges() const
{
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  if(!gl_vbo_available()) return;

  // Ensure edges are drawn in black
  ::glColor3f( 0.f, 0.f, 0.f );
  
  // Draw edges
  ::glEnableClientState( GL_VERTEX_ARRAY );

  // Get buffers
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[0]);
  ::glVertexPointer(3, GL_FLOAT, 0, 0);

  ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);

  // Render
  ::glDrawElements(GL_QUADS, ibo_size(), GL_UNSIGNED_INT, 0);

  // Cleanup
  ::glBindBuffer(GL_ARRAY_BUFFER, 0);
  ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  ::glDisableClientState( GL_VERTEX_ARRAY ); 
#endif // SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
}


GLint
Scene_segmented_image_item::ibo_size() const
{
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  if(gl_vbo_available()) {
    GLint nb_elts = 0;
    ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
    ::glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &nb_elts);

    return nb_elts/sizeof(GLuint);
  }
#endif // SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE

  return 0;
}


#include "Scene_segmented_image_item.moc"

