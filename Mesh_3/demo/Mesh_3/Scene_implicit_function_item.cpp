#include "config.h"

#include "Scene_implicit_function_item.h"
#include <QColor>
#include <map>
#include <CGAL/gl.h>
#include <CGAL/Simple_cartesian.h>

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

#include "Color_ramp.h"


Scene_implicit_function_item::
Scene_implicit_function_item(Implicit_function_interface* f)
  : function_(f)
  , frame_(new ManipulatedFrame())
  , grid_size_(SCENE_IMPLICIT_GRID_SIZE)
  , max_value_(0.)
  , min_value_(0.)
  , blue_color_ramp_()
  , red_color_ramp_()
{
  compile_shaders();
  texture = new Texture(grid_size_,grid_size_);
  blue_color_ramp_.build_blue();
  red_color_ramp_.build_red();
  compute_min_max();
  compute_function_grid();
  connect(frame_, SIGNAL(modified()), this, SLOT(compute_function_grid()));
  are_buffers_initialized = false;
  texture_initialized = false;
}


Scene_implicit_function_item::~Scene_implicit_function_item()
{
    for(int i=0; i<vboSize; i++)
        buffers[i].destroy();
    for(int i=0; i<vaoSize; i++)
        vao[i].destroy();
    v_cube.clear();
    v_plan.clear();
    texture_map.clear();
  delete frame_;
}


Scene_implicit_function_item::Bbox
Scene_implicit_function_item::bbox() const
{
  return function_->bbox();
}


/**************************************************
****************SHADER FUNCTIONS******************/

void Scene_implicit_function_item::compile_shaders()
{

    for(int i=0; i< vboSize; i++)
        buffers[i].create();
    for(int i=0; i< vaoSize; i++)
        vao[i].create();

    //Vertex source code
    const char vertex_source[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"

        "uniform highp mat4 mvp_matrix;\n"
        "void main(void)\n"
        "{\n"
        "   gl_Position =  mvp_matrix * vertex; \n"
        "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
        "#version 120 \n"
        "uniform vec4 color; \n"
        "void main(void) { \n"
         "gl_FragColor = color; \n"
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

    //Vertex source code
    const char tex_vertex_source[] =
    {
         "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec2 tex_coord; \n"
        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 f_matrix;\n"
        "varying highp vec2 texc;\n"
        "void main(void)\n"
        "{\n"
        "   gl_Position = mvp_matrix * f_matrix * vertex;\n"
        "    texc = tex_coord;\n"
        "}"
    };
    //Fragment source code
    const char tex_fragment_source[] =
    {
        "#version 120 \n"
        "uniform sampler2D texture;\n"
        "varying highp vec2 texc;\n"
        "void main(void) { \n"
        "gl_FragColor = texture2D(texture, texc.st);\n"
        "} \n"
        "\n"
    };
    QOpenGLShader *tex_vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!tex_vertex_shader->compileSourceCode(tex_vertex_source))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader *tex_fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!tex_fragment_shader->compileSourceCode(tex_fragment_source))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!tex_rendering_program.addShader(tex_vertex_shader))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!tex_rendering_program.addShader(tex_fragment_shader))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!tex_rendering_program.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
    tex_rendering_program.bind();
}

void Scene_implicit_function_item::compute_elements()
{
    v_cube.resize(0);
    v_plan.resize(0);
    texture_map.resize(0);
    // The Quad
       {


           //A
         v_plan.push_back(bbox().xmin);
         v_plan.push_back(bbox().ymin);
         v_plan.push_back(0);


           //B
         v_plan.push_back(bbox().xmin);
         v_plan.push_back(bbox().ymax);
         v_plan.push_back(0);


           //C
         v_plan.push_back(bbox().xmax);
         v_plan.push_back(bbox().ymin);
         v_plan.push_back(0);



           //C
         v_plan.push_back(bbox().xmax);
         v_plan.push_back(bbox().ymin);
         v_plan.push_back(0);


           //B
         v_plan.push_back(bbox().xmin);
         v_plan.push_back(bbox().ymax);
         v_plan.push_back(0);


           //D
         v_plan.push_back(bbox().xmax);
         v_plan.push_back(bbox().ymax);
         v_plan.push_back(0);


           //UV Mapping
           texture_map.push_back(0.0);
           texture_map.push_back(0.0);

           texture_map.push_back(0.0);
           texture_map.push_back(1.0);

           texture_map.push_back(1.0);
           texture_map.push_back(0.0);

           texture_map.push_back(1.0);
           texture_map.push_back(0.0);

           texture_map.push_back(0.0);
           texture_map.push_back(1.0);

           texture_map.push_back(1.0);
           texture_map.push_back(1.0);
       }
    //The box
    {

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmin);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmin);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymax);
        v_cube.push_back(bbox().zmax);

        v_cube.push_back(bbox().xmax);
        v_cube.push_back(bbox().ymin);
        v_cube.push_back(bbox().zmax);

    }
    //The texture
   for(int i=0; i < texture->getWidth(); i++)
       for( int j=0 ; j < texture->getHeight() ; j++)
       {
          compute_texture(i,j);
       }
}

void Scene_implicit_function_item::initialize_buffers(Viewer* viewer) const
{

    rendering_program.bind();

    vao[0].bind();
    buffers[0].bind();
    buffers[0].allocate(v_cube.data(), static_cast<int>(v_cube.size()*sizeof(float)));
    vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(vertexLocation[0]);
    rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
    buffers[0].release();
    vao[0].release();
    rendering_program.release();
    tex_rendering_program.bind();
    //cutting plane
    vao[1].bind();
    buffers[1].bind();
    buffers[1].allocate(v_plan.data(), static_cast<int>(v_plan.size()*sizeof(float)));
    vertexLocation[1] = tex_rendering_program.attributeLocation("vertex");
    tex_rendering_program.bind();
    tex_rendering_program.setAttributeBuffer(vertexLocation[1],GL_FLOAT,0,3);
    tex_rendering_program.enableAttributeArray(vertexLocation[1]);
    buffers[1].release();
    tex_rendering_program.release();

    buffers[2].bind();
    buffers[2].allocate(texture_map.data(), static_cast<int>(texture_map.size()*sizeof(float)));
    tex_Location = tex_rendering_program.attributeLocation("tex_coord");
    tex_rendering_program.bind();
    tex_rendering_program.setAttributeBuffer(tex_Location,GL_FLOAT,0,2);
    tex_rendering_program.enableAttributeArray(tex_Location);
    buffers[2].release();
    tex_rendering_program.release();

    viewer->glBindTexture(GL_TEXTURE_2D, textureId);
    viewer->glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGB,
                 texture->getWidth(),
                 texture->getHeight(),
                 0,
                 GL_RGB,
                 GL_UNSIGNED_BYTE,
                 texture->getData());
    viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE );
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE );
    vao[1].release();
    tex_rendering_program.release();
    are_buffers_initialized = true;
}

void Scene_implicit_function_item::attrib_buffers(Viewer* viewer) const
{
    QMatrix4x4 mvpMatrix;
    QMatrix4x4 fMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvpMatrix.data()[i] = (float)mat[i];
    }
    frame_->getMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        fMatrix.data()[i] = (float)mat[i];
    }

    rendering_program.bind();
    colorLocation[0] = rendering_program.uniformLocation("color");
    mvpLocation[0] = rendering_program.uniformLocation("mvp_matrix");
    rendering_program.setUniformValue(mvpLocation[0], mvpMatrix);;
    rendering_program.release();

    tex_rendering_program.bind();
    colorLocation[1] = tex_rendering_program.uniformLocation("color");
    f_Location = tex_rendering_program.uniformLocation("f_matrix");
    mvpLocation[1] = tex_rendering_program.uniformLocation("mvp_matrix");
    tex_rendering_program.setUniformValue(mvpLocation[1], mvpMatrix);;
    tex_rendering_program.setUniformValue(f_Location, fMatrix);;
    tex_rendering_program.release();

}

void
Scene_implicit_function_item::draw(Viewer* viewer) const
{
  if(!texture_initialized)
  {
    viewer->glGenTextures(1, &textureId);
    texture_initialized = true;
  }
  if(!are_buffers_initialized)
    initialize_buffers(viewer);
  QColor color;
  vao[0].bind();
  attrib_buffers(viewer);
  rendering_program.bind();
  color.setRgbF(0.0,0.0,0.0);
  rendering_program.setUniformValue(colorLocation[0], color);
  glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(v_cube.size()/3));
  rendering_program.release();
  vao[0].release();

  viewer->glActiveTexture(GL_TEXTURE0);
  viewer->glBindTexture(GL_TEXTURE_2D, textureId);

  vao[1].bind();
  tex_rendering_program.bind();
  glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(v_plan.size()/3));
  tex_rendering_program.release();
  vao[1].release();

}

QString
Scene_implicit_function_item::toolTip() const
{
  return tr("<p>Function <b>%1</b>")
    .arg(this->name());
}

bool
Scene_implicit_function_item::supportsRenderingMode(RenderingMode m) const
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
Scene_implicit_function_item::
compute_function_grid()
{
  typedef CGAL::Simple_cartesian<double>  K;
  typedef K::Aff_transformation_3         Aff_transformation;
  typedef K::Point_3                      Point_3;
  
  // Get transformation
  const ::GLdouble* m = frame_->matrix();
  
  // OpenGL matrices are row-major matrices
  Aff_transformation t (m[0], m[4], m[8], m[12],
                        m[1], m[5], m[9], m[13],
                        m[2], m[6], m[10], m[14]);
  
  double diag = bbox().diagonal_length() * .6;
  
  const double dx = diag;
  const double dy = diag;
  const double z (0);
  
  for(int i=0 ; i<grid_size_ ; ++i)
  {
    double x = -diag/2. + double(i)/double(grid_size_) * dx;
    
    for(int j=0 ; j<grid_size_ ; ++j)
    {
      double y = -diag/2. + double(j)/double(grid_size_) * dy;
      
      Point_3 query = t( Point_3(x,y,z) );
      double v = function_->operator()(query.x(), query.y(), query.z());
      
      implicit_grid_[i][j] = Point_value(Point(query.x(),query.y(),query.z()),v);
    }
  }
  
  // Update display list
  this->changed();
}

void
Scene_implicit_function_item::
compute_min_max()
{
  max_value_ = 0;
  min_value_ = 0;
  
  double probes_nb = double(grid_size_) / 2;
  
  // Probe bounding box
  const Bbox& b = bbox();
  
  for ( int i = 0 ; i <= probes_nb ; ++i )
  {
    double x = b.xmin + double(i) * (b.xmax - b.xmin) / probes_nb;
    
    for ( int j = 0 ; j <= probes_nb ; ++j )
    {
      double y = b.ymin + double(j) * (b.ymax - b.ymin) / probes_nb;
      
      for ( int k = 0 ; k <= probes_nb ; ++k )
      {
        double z = b.zmin + double(k) * (b.zmax - b.zmin) / probes_nb;
        
        double v = (*function_)(x,y,z);
        max_value_ = (std::max)(v,max_value_);
        min_value_ = (std::min)(v,min_value_);
      }
    }
  }
}

void Scene_implicit_function_item::changed()
{
    compute_elements();
    are_buffers_initialized = false;
}

void Scene_implicit_function_item::compute_texture(int i, int j)
{

  double v =(implicit_grid_[i][j]).second;
       // determines grey level
       if ( v > 0 )
       {
           v = v/max_value_;
           GLdouble r = red_color_ramp_.r(v), g = red_color_ramp_.g(v), b = red_color_ramp_.b(v);
           texture->setData(i,j,255*r,255*g,255*b);
       }
       else
       {
           v = v/min_value_;
           GLdouble r = blue_color_ramp_.r(v), g = blue_color_ramp_.g(v), b = blue_color_ramp_.b(v);
           texture->setData(i,j,255*r,255*g,255*b);
       }
}


