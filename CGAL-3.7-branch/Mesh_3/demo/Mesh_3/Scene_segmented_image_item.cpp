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

namespace {
  
  unsigned char image_data(const Image& im,
                           unsigned int i, unsigned int j, unsigned int k)
  {
    if ( i<im.xdim() && j<im.ydim() && k<im.zdim() )
      return CGAL::IMAGEIO::static_evaluate<unsigned char>(im.image(),i,j,k);
    else
      return 0;
  }
}


Scene_segmented_image_item::Scene_segmented_image_item(Image* im)
  : m_image(im)
  , m_initialized(false)
  , m_draw_edges(true)
{
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  if(gl_vbo_available()) {
    ::glGenBuffers(3,m_vbo);
    ::glGenBuffers(1,&m_ibo);
  }
#endif // SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE

  initialize_buffers();
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

  const unsigned int& xdim = m_image->xdim();
  const unsigned int& ydim = m_image->ydim();
  const unsigned int& zdim = m_image->zdim();
  
  // -----------------------------------
  // Scan colors
  // -----------------------------------
  typedef std::map<unsigned char, QColor> ValueMap;
  ValueMap values;
  for(unsigned int i=0;i<xdim;i+=5)
  { 
    for(unsigned int j=0;j<ydim;j+=5)
    { 
      for(unsigned int k=0;k<zdim;k+=5)
      {
        unsigned char value = image_data(*m_image,i,j,k);
        
        if ( 0 < value )
        {
          values.insert(std::make_pair(value,QColor()));
        }
      }
    }
  }
  
  int i=0;
  const double starting_hue = 300./360.; // magenta
  for ( ValueMap::iterator it = values.begin(), end = values.end() ;
       it != end ; ++it, ++i )
  {
    double hue =  starting_hue + 1./values.size() * i;
    if ( hue > 1. ) { hue -= 1.; }
    it->second = QColor::fromHsvF(hue, 1., 0.8);
  }
  
  
  // -----------------------------------
  // Get data
  // -----------------------------------
  
  // Get step size
  int steps_max = 75;
  int size_max = (std::max)((std::max)(xdim, ydim), zdim);
  int d_max = (std::max)(1, size_max/steps_max);
  int dx=d_max, dy=d_max, dz=d_max;
  
  // Normals
  float a = std::sqrt(1.f/3.f);
  GLfloat normal_array[24] = { -a,-a,-a, -a,-a,a, -a,a,a, -a,a,-a,
                                a,-a,-a, a,-a,a, a,a,a, a,a,-a };
  
  // Cube faces
  GLuint indice_array[24] = { 0,1,2,3, 0,1,5,4, 0,4,7,3,
                              1,5,6,2, 2,3,7,6, 4,5,6,7 };
  
  // Cube vertices
  int cube_array[24] = { 0,0,0, 0,0,1, 0,1,1, 0,1,0,
                         1,0,0, 1,0,1, 1,1,1, 1,1,0 };
  
  // Tester (avoids drawing of interior cubes)
  int tester[18] = { -1,0,0, 0,-1,0, 0,0,-1, 1,0,0, 0,1,0, 0,0,1 };
                     // 0,1,1, 1,0,1, 1,1,0, 2,1,1, 1,2,1, 1,1,2 };
  
  // Stores Gl elements
  std::vector<float> vertices, normals, colors;
  std::vector<unsigned int> indices;
  
  unsigned int delta = 0;
  for(unsigned int i=0;i<xdim;i+=dx)
  { 
    for(unsigned int j=0;j<ydim;j+=dy)
    { 
      for(unsigned int k=0;k<zdim;k+=dz)
      {
        unsigned char value = image_data(*m_image,i,j,k);
        
        if(value > 0)
        {
          // Don't draw interior cubes
          unsigned int t=0;
          while ( t < (sizeof(tester)/sizeof(int)) 
                  && image_data(*m_image,
                                i + dx*tester[t],
                                j + dy*tester[t+1],
                                k + dz*tester[t+2] ) > 0 )
          {
            t+=3;
          }
          
          if ( sizeof(tester)/sizeof(int) == t )
            continue;
          
          // Vertices
          float x = m_image->vx() * i;
          float y = m_image->vy() * j;
          float z = m_image->vz() * k;
          
          float xdx = dx * m_image->vx();
          float ydy = dy * m_image->vy();
          float zdz = dz * m_image->vz();
          
          for ( int c=0; c<24; c+=3 )
          { 
            vertices.push_back(x + cube_array[c] * xdx );
            vertices.push_back(y + cube_array[c+1] * ydy );
            vertices.push_back(z + cube_array[c+2] * zdz );
          }
          
          // Colors
          for ( int c=0; c<24; c+=3)
          {
            unsigned char data = image_data(*m_image,
                                            i+dx*cube_array[c],
                                            j+dy*cube_array[c+1],
                                            k+dz*cube_array[c+2]);
            
            QColor color = (data > 0) ? values[data] : values[value];
            //color = color.darker(150);
            
            colors.push_back(color.red()/255.f);
            colors.push_back(color.green()/255.f);
            colors.push_back(color.blue()/255.f);
          }
          
          // Normals
          for ( int c=0; c<24; ++c )
            normals.push_back(normal_array[c]);
          
          // Indices
          for ( int c=0; c<24; ++c )
            indices.push_back(indice_array[c]+delta);
          
          delta += 8;
        }
      }
    }
  }
  
  // -----------------------------------
  // Create GL buffers
  // -----------------------------------
  CGAL_assertion( vertices.size() == normals.size() 
                 && vertices.size() == colors.size() );
  
  std::size_t vertices_nb = vertices.size();
  
  GLfloat* vertices_array = new GLfloat[vertices_nb];
  GLfloat* normals_array = new GLfloat[vertices_nb];
  GLfloat* colors_array = new GLfloat[vertices_nb];
  
  i = 0;
  for ( std::vector<float>::iterator 
       vit = vertices.begin(), vend = vertices.end(),
       nit = normals.begin(),  nend = normals.end(),
       cit = colors.begin(),   cend = colors.end() ;
       vit != vend && nit != nend && cit != cend ;
       ++vit, ++nit, ++cit, ++i )
  {
    vertices_array[i] = *vit;
    normals_array[i] = *nit;
    colors_array[i] = *cit;
  }
  
  std::size_t indices_size = indices.size();
  GLuint* indices_array = new GLuint[indices_size];
  
  i = 0;
  for ( std::vector<unsigned int>::iterator it = indices.begin(),
       end = indices.end() ; it != end ; ++it, ++i )
  {
    indices_array[i] = *it;
  }
  
  // Vertex buffer
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[0]);
  ::glBufferData(GL_ARRAY_BUFFER, vertices_nb*sizeof(GLfloat), vertices_array, GL_STATIC_DRAW);
  
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[1]);
  ::glBufferData(GL_ARRAY_BUFFER, vertices_nb*sizeof(GLfloat), normals_array, GL_STATIC_DRAW);
  
  ::glBindBuffer(GL_ARRAY_BUFFER, m_vbo[2]);
  ::glBufferData(GL_ARRAY_BUFFER, vertices_nb*sizeof(GLfloat), colors_array, GL_STATIC_DRAW);
  
  // Indices buffer
  ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
  ::glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_size*sizeof(GLuint), indices_array, GL_STATIC_DRAW);
  
  // Close buffers
  ::glBindBuffer(GL_ARRAY_BUFFER, 0);
  ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  
  // Cleanup
  delete vertices_array;
  delete normals_array;
  delete colors_array;
  delete indices_array;

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

