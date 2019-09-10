#include "config.h"

#include "Scene_image_item.h"
#include "Image_type.h"
#include <QColor>
#include <map>
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

  const GLfloat* colors() const { return colors_.data(); }
  const GLfloat* normals() const { return normals_.data(); }
  const GLfloat* vertices() const { return vertices_.data(); }
  const GLuint* quads() const { return quads_.data(); }
  
  std::size_t color_size() const { return colors_.size()*sizeof(GLfloat); }
  std::size_t normal_size() const { return normals_.size()*sizeof(GLfloat); }
  std::size_t vertex_size() const { return vertices_.size()*sizeof(GLfloat); }
  std::size_t quad_size() const { return quads_.size()*sizeof(GLuint); }
  
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
  std::vector<GLfloat> colors_, normals_, vertices_;
  std::vector<GLuint> quads_;
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
    v_box = new std::vector<float>();
    is_ogl_4_3 = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->isOpenGL_4_3();
    is_hidden = hidden;
    compile_shaders();
    initializeBuffers();
  }

  ~Scene_image_item_priv()
  {
    for(int i=0; i<vboSize; i++)
        m_vbo[i].destroy();
    for(int i=0; i<vaoSize; i++)
        vao[i].destroy();
  }
  void draw_gl(CGAL::Three::Viewer_interface* viewer) const;

  void initializeBuffers();
  GLint ibo_size() const;
  void draw_bbox();
  void attribBuffers(CGAL::Three::Viewer_interface*) const;
  void compile_shaders();
  void draw_Bbox(Scene_item::Bbox bbox, std::vector<float> *vertices);

  bool m_initialized;
//#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  int m_voxel_scale;
  std::vector<float> *v_box;
  std::vector<float> color;
  static const int vaoSize = 2;
  static const int vboSize = 6;
  mutable int poly_vertexLocation[1];
  mutable int normalsLocation[1];
  mutable int mvpLocation[1];
  mutable int mvLocation[1];
  mutable int colorLocation[1];
  mutable int lightLocation[5];
  mutable int twosideLocation;
  mutable QOpenGLBuffer m_vbo[vboSize];
  mutable QOpenGLBuffer *m_ibo;
  mutable QOpenGLVertexArrayObject vao[vaoSize];
  mutable QOpenGLShaderProgram rendering_program;
  bool is_hidden;
  bool is_ogl_4_3;
  Scene_image_item* item;

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
}


Scene_image_item::~Scene_image_item()
{
   delete d;
}

/**************************************************
****************SHADER FUNCTIONS******************/

void Scene_image_item_priv::compile_shaders()
{
  for(int i=0; i< vboSize; i++)
    m_vbo[i].create();
  for(int i=0; i< vaoSize; i++)
    vao[i].create();
  m_ibo = new QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
  m_ibo->create();
  if(!is_ogl_4_3)
  {
    //Vertex source code
    const char vertex_source[] =
    {
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
      "   mat3 mv_matrix_3;                    "
      "   mv_matrix_3[0] = mv_matrix[0].xyz;   "
      "   mv_matrix_3[1] = mv_matrix[1].xyz;   "
      "   mv_matrix_3[2] = mv_matrix[2].xyz;   "
      "   fN = mv_matrix_3* normal;           "
      "   gl_Position = mvp_matrix * vertex; \n"
      "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
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
      "       diffuse = max(dot(N,L),0.0) * light_diff*color; \n"
      "   else \n"
      "       diffuse = max(abs(dot(N,L)),0.0) * light_diff*color; \n"
      "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

      "gl_FragColor = color*light_amb + diffuse + specular; \n"
      "} \n"
      "\n"
    };
    QOpenGLShader vertex_shader(QOpenGLShader::Vertex);
    if(!vertex_shader.compileSourceCode(vertex_source))
    {
      std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader fragment_shader(QOpenGLShader::Fragment);
    if(!fragment_shader.compileSourceCode(fragment_source))
    {
      std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!rendering_program.addShader(&vertex_shader))
    {
      std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }

    if(!rendering_program.addShader(&fragment_shader))
    {
      std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program.link())
    {
      std::cerr<<"linking Program FAILED"<<std::endl;
    }
  }
  else
  {
    if(!rendering_program.addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/no_interpolation_shader.v"))
    {
      std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }

    if(!rendering_program.addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/no_interpolation_shader.f"))
    {
      std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program.addShaderFromSourceFile(QOpenGLShader::Geometry,":/cgal/Polyhedron_3/resources/no_interpolation_shader.g"))
    {
      std::cerr<<"adding geometry shader FAILED"<<std::endl;
    }
    rendering_program.link();
  }
  rendering_program.bindAttributeLocation("colors", 1);
}

void Scene_image_item_priv::attribBuffers(Viewer_interface* viewer) const
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

void
Scene_image_item::compute_bbox() const
{
  if(!m_image)
    _bbox = Bbox();
  else
   _bbox = Bbox(m_image->image()->tx,
                m_image->image()->ty,
                m_image->image()->tz,
                m_image->image()->tx+(m_image->xdim()-1) * m_image->vx(),
                m_image->image()->ty+(m_image->ydim()-1) * m_image->vy(),
                m_image->image()->tz+(m_image->zdim()-1) * m_image->vz());
}

void
Scene_image_item::draw(Viewer_interface* viewer) const
{
  if(m_image)
  {
    d->draw_gl(viewer);
  }
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

void
Scene_image_item_priv::initializeBuffers()
{
  draw_Bbox(item->bbox(), v_box);
  if(!is_hidden)
  {
    internal::Image_accessor image_data_accessor (*item->m_image,
                                                  m_voxel_scale,
                                                  m_voxel_scale,
                                                  m_voxel_scale);
    internal::Vertex_buffer_helper helper (image_data_accessor, is_ogl_4_3);
    helper.fill_buffer_data();
    rendering_program.bind();
    vao[0].bind();
    m_vbo[0].bind();
    m_vbo[0].allocate(helper.vertices(), static_cast<int>(helper.vertex_size()));
    rendering_program.enableAttributeArray("vertex");
    rendering_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
    m_vbo[0].release();
    m_vbo[1].bind();
    m_vbo[1].allocate(helper.normals(), static_cast<int>(helper.normal_size()));
    rendering_program.enableAttributeArray("normal");
    rendering_program.setAttributeBuffer("normal",GL_FLOAT,0,3);
    m_vbo[1].release();

    m_vbo[2].bind();
    m_vbo[2].allocate(helper.colors(), static_cast<int>(helper.color_size()));
    rendering_program.enableAttributeArray("inColor");
    rendering_program.setAttributeBuffer("inColor",GL_FLOAT,0,3);
    m_vbo[2].release();

    m_ibo->bind();
    m_ibo->allocate(helper.quads(), static_cast<int>(helper.quad_size()));
    vao[0].release();

    color.resize(0);
    for(std::size_t i=0; i<helper.color_size()/sizeof(GLuint); i++)
      color.push_back(0.0);
    rendering_program.release();
  }
  QOpenGLShaderProgram* line_program = item->getShaderProgram(Scene_item::PROGRAM_NO_SELECTION);
  line_program->bind();
  vao[1].bind();
  m_vbo[3].bind();
  m_vbo[3].allocate(v_box->data(), static_cast<int>(v_box->size()*sizeof(float)));
  line_program->enableAttributeArray("vertex");
  line_program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  m_vbo[3].release();

  line_program->disableAttributeArray("colors");

  vao[1].release();
  line_program->release();
  m_initialized = true;
}

void
Scene_image_item_priv::draw_gl(Viewer_interface* viewer) const
{
  attribBuffers(viewer);
  rendering_program.bind();
  if(!is_hidden)
  {
    vao[0].bind();
    if(!is_ogl_4_3)
      viewer->glDrawElements(GL_TRIANGLES, GLsizei(m_ibo->size()/sizeof(GLuint)), GL_UNSIGNED_INT, 0);
    else
      viewer->glDrawElements(GL_LINES_ADJACENCY, GLsizei(m_ibo->size()/sizeof(GLuint)), GL_UNSIGNED_INT, 0);
    vao[0].release();
  }
  rendering_program.release();
  QOpenGLShaderProgram* line_program;
  vao[1].bind();
  if(!viewer->isOpenGL_4_3())
  {
    item->attribBuffers(viewer,Scene_item::PROGRAM_NO_SELECTION);
    line_program = item->getShaderProgram(Scene_item::PROGRAM_NO_SELECTION);
    line_program->bind();
  }
  else
  {
    item->attribBuffers(viewer,Scene_item::PROGRAM_SOLID_WIREFRAME);
    line_program = item->getShaderProgram(Scene_item::PROGRAM_SOLID_WIREFRAME);
    line_program->bind();
    QVector2D vp(viewer->width(), viewer->height());
    line_program->setUniformValue("viewport", vp);
    line_program->setUniformValue("near",(GLfloat)viewer->camera()->zNear());
    line_program->setUniformValue("far",(GLfloat)viewer->camera()->zFar());
    line_program->setUniformValue("width", 3.0f);
  }
  line_program->setAttributeValue("colors", QColor(Qt::black));
  viewer->glDepthRangef(0.00001f, 0.99999f);
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(v_box->size()/3));
  viewer->glDepthRangef(0.0f, 1.0f);
  vao[1].release();
  line_program->release();
}

GLint
Scene_image_item_priv::ibo_size() const
{
      m_ibo->bind();
    GLint nb_elts = m_ibo->size();
    m_ibo->release();

    return GLint(nb_elts/sizeof(GLuint));

  return 0;
}

void Scene_image_item::changed()
{
    d->initializeBuffers();
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

void Scene_image_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{ d->draw_gl(viewer); }

bool Scene_image_item::isGray() { return d->is_hidden;}

void Scene_image_item::invalidateOpenGLBuffers()
{
  d->v_box->clear();
  changed();
}
