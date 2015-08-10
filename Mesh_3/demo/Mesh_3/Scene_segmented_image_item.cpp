#include "config.h"

#include "Scene_segmented_image_item.h"
#include "Image_type.h"
#include <QColor>
#include <map>
#include <CGAL/gl.h>
#include <CGAL/ImageIO.h>
#include <CGAL/use.h>



// -----------------------------------
// Internal classes
// -----------------------------------
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
  void treat_vertex(std::size_t i, std::size_t j, std::size_t k);
  
  void push_color(std::size_t i, std::size_t j, std::size_t k);
  void push_normal(std::size_t i, std::size_t j, std::size_t k);
  void push_vertex(std::size_t i, std::size_t j, std::size_t k);
  void push_quads(std::size_t i, std::size_t j, std::size_t k);
  void push_quad(int pos1, int pos2, int pos3, int pos4);
  
  int compute_position(std::size_t i, std::size_t j, std::size_t k) const;
  int vertex_index(std::size_t i, std::size_t j, std::size_t k) const;
  
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
  
  create_arrays();
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
  
  colors_.push_back(color.red()/255.f);
  colors_.push_back(color.green()/255.f);
  colors_.push_back(color.blue()/255.f);
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
  
  vertices_.push_back(i*data_.vx());
  vertices_.push_back(j*data_.vy());
  vertices_.push_back(k*data_.vz());

}

void
Vertex_buffer_helper::push_quads(std::size_t i, std::size_t j, std::size_t k)
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
    quads_.push_back(pos1);
    quads_.push_back(pos3);
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
  CGAL_USE(display_scale);

  v_box = new std::vector<float>();
  compile_shaders();
  initialize_buffers();
  setRenderingMode(Flat);
}


Scene_segmented_image_item::~Scene_segmented_image_item()
{
    for(int i=0; i<vboSize; i++)
        m_vbo[i].destroy();
    for(int i=0; i<vaoSize; i++)
        vao[i].destroy();
}

/**************************************************
****************SHADER FUNCTIONS******************/

void Scene_segmented_image_item::compile_shaders()
{

    for(int i=0; i< vboSize; i++)
        m_vbo[i].create();
    for(int i=0; i< vaoSize; i++)
        vao[i].create();
    m_ibo = new QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
    m_ibo->create();
    //Vertex source code
    const char vertex_source[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec3 normal;\n"
        "attribute highp vec4 inColor;\n"

        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 mv_matrix; \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "varying highp vec4 color; \n"
        "void main(void)\n"
        "{\n"
        "   color=inColor; \n"
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
        "varying highp vec4 color; \n"
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

void Scene_segmented_image_item::attrib_buffers(Viewer* viewer) const
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
Scene_segmented_image_item::draw(Viewer* viewer) const
{
  if(m_image)
  {

    //m_image->gl_draw_bbox(3.0f,0,0,0);
    draw_gl(viewer);
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
  internal::Image_accessor image_data_accessor (*m_image,
                                                m_voxel_scale,
                                                m_voxel_scale,
                                                m_voxel_scale);
  
  internal::Vertex_buffer_helper helper (image_data_accessor);
  helper.fill_buffer_data();

  draw_Bbox(bbox(), v_box);
  std::vector<float> nul_vec(0);
  for(std::size_t i=0; i<v_box->size(); i++)
      nul_vec.push_back(0.0);

  rendering_program.bind();
  vao[0].bind();
  m_vbo[0].bind();
  m_vbo[0].allocate(helper.vertices(), static_cast<int>(helper.vertex_size()));
  poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
  rendering_program.enableAttributeArray(poly_vertexLocation[0]);
  rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
  m_vbo[0].release();

  m_vbo[1].bind();
  m_vbo[1].allocate(helper.normals(), static_cast<int>(helper.normal_size()));
  normalsLocation[0] = rendering_program.attributeLocation("normal");
  rendering_program.enableAttributeArray(normalsLocation[0]);
  rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
  m_vbo[1].release();

  m_vbo[2].bind();
  m_vbo[2].allocate(helper.colors(), static_cast<int>(helper.color_size()));
  colorLocation[0] = rendering_program.attributeLocation("inColor");
  rendering_program.enableAttributeArray(colorLocation[0]);
  rendering_program.setAttributeBuffer(colorLocation[0],GL_FLOAT,0,3);
  m_vbo[2].release();

  m_ibo->bind();
  m_ibo->allocate(helper.quads(), static_cast<int>(helper.quad_size()));
  vao[0].release();

  color.resize(0);
  for(std::size_t i=0; i<helper.color_size()/sizeof(GLuint); i++)
      color.push_back(0.0);

  vao[1].bind();
  m_vbo[3].bind();
  m_vbo[3].allocate(v_box->data(), static_cast<int>(v_box->size()*sizeof(float)));
  poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
  rendering_program.enableAttributeArray(poly_vertexLocation[0]);
  rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
  m_vbo[3].release();

  m_vbo[4].bind();
  m_vbo[3].allocate(nul_vec.data(), static_cast<int>(nul_vec.size()*sizeof(float)));
  normalsLocation[0] = rendering_program.attributeLocation("normal");
  rendering_program.enableAttributeArray(normalsLocation[0]);
  rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
  m_vbo[4].release();

  m_vbo[5].bind();
  m_vbo[5].allocate(nul_vec.data(), static_cast<int>(nul_vec.size()*sizeof(float)));
  colorLocation[0] = rendering_program.attributeLocation("inColor");
  rendering_program.enableAttributeArray(colorLocation[0]);
  rendering_program.setAttributeBuffer(colorLocation[0],GL_FLOAT,0,3);
  m_vbo[5].release();

  m_ibo->bind();
  vao[1].release();
  rendering_program.release();

  m_initialized = true;
}


void
Scene_segmented_image_item::draw_gl(Viewer* viewer) const
{
  attrib_buffers(viewer);
  rendering_program.bind();
  vao[0].bind();
  viewer->glDrawElements(GL_TRIANGLES, m_ibo->size()/sizeof(GLuint), GL_UNSIGNED_INT, 0);
  vao[0].release();

  vao[1].bind();
  viewer->glLineWidth(3);
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(v_box->size()/3));
  vao[1].release();
  rendering_program.release();
}

GLint
Scene_segmented_image_item::ibo_size() const
{
      m_ibo->bind();
    GLint nb_elts = m_ibo->size();
    m_ibo->release();

    return nb_elts/sizeof(GLuint);

  return 0;
}

void Scene_segmented_image_item::changed()
{
    initialize_buffers();
}

void Scene_segmented_image_item::draw_Bbox(Bbox bbox, std::vector<float> *vertices)
{
    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmin);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmin);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymax);
    vertices->push_back(bbox.zmax);

    vertices->push_back(bbox.xmax);
    vertices->push_back(bbox.ymin);
    vertices->push_back(bbox.zmax);

}


