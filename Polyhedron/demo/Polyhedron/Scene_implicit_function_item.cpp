#include "Scene_implicit_function_item.h"
#include <QColor>
#include <QApplication>
#include <map>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Qt/manipulatedFrame.h>

#include "Color_ramp.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/double.h>
using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Triangle_container Tc;
typedef Edge_container Ec;

inline
bool is_nan(double d)
{
    return !CGAL::Is_valid<double>()( d );
}

struct Scene_implicit_function_item_priv
{
  Scene_implicit_function_item_priv(Implicit_function_interface* f, Scene_implicit_function_item* parent)
    : function_(f)
    , frame_(new ManipulatedFrame())
    , need_update_(true)
    , grid_size_(SCENE_IMPLICIT_GRID_SIZE)
    , max_value_(0.)
    , min_value_(0.)
    , blue_color_ramp_()
    , red_color_ramp_()
  {
    init=false;
    item = parent;
    blue_color_ramp_.build_blue();
    red_color_ramp_.build_red();
  }
  ~Scene_implicit_function_item_priv()
  {
    delete frame_;
  }
  typedef CGAL::qglviewer::Vec                  Point;
  typedef std::pair <Point,double>        Point_value;
  typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;
  void compute_min_max();
  void initialize_buffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_vertices_and_texmap(void);
  void compute_texture(int, int);
  Implicit_function_interface* function_;
  ManipulatedFrame* frame_;
  mutable bool need_update_;
  int grid_size_;
  double max_value_;
  double min_value_;
  bool init;
  mutable Point_value implicit_grid_[SCENE_IMPLICIT_GRID_SIZE][SCENE_IMPLICIT_GRID_SIZE];
  Color_ramp blue_color_ramp_;
  Color_ramp red_color_ramp_;
  mutable std::vector<float> positions_cube;
  mutable std::vector<float> positions_grid;
  mutable std::vector<float> positions_tex_quad;
  mutable std::vector<float> texture_map;
  std::size_t nb_quad;
  std::size_t nb_cube;
  std::size_t nb_grid;
  GLuint vao;
  GLuint buffer[4];
  mutable QOpenGLShaderProgram *program;
  Scene_implicit_function_item* item;
};

void Scene_implicit_function_item_priv::initialize_buffers(CGAL::Three::Viewer_interface *viewer = 0) const
{
    item->getTriangleContainer(0)->initializeBuffers(viewer);
    item->getTriangleContainer(0)->setFlatDataSize(nb_quad);

    item->getEdgeContainer(0)->initializeBuffers(viewer);
    item->getEdgeContainer(0)->setFlatDataSize(nb_cube);

    item->getEdgeContainer(1)->initializeBuffers(viewer);
    item->getEdgeContainer(1)->setFlatDataSize(nb_grid);

    positions_tex_quad.clear();
    positions_tex_quad.shrink_to_fit();
    texture_map.clear();
    texture_map.shrink_to_fit();
    positions_cube.clear();
    positions_cube.shrink_to_fit();
    positions_grid.clear();
    positions_grid.shrink_to_fit();
}

void Scene_implicit_function_item_priv::compute_vertices_and_texmap(void)
{
    positions_tex_quad.resize(0);
    positions_cube.resize(0);
    positions_grid.resize(0);
    texture_map.resize(0);

    const CGAL::Three::Scene_item::Bbox& b = item->bbox();
    float x,y,z;
    z = static_cast<float>(b.zmax()+b.zmin())/2.0f;
    x = static_cast<float>(b.xmax()+b.xmin())/2.0f;
    y = static_cast<float>(b.ymax()+b.ymin())/2.0f;
    // The Quad
    {


        //A
        positions_tex_quad.push_back(b.xmin()-x);
        positions_tex_quad.push_back(b.ymin()-z);
        positions_tex_quad.push_back(0);


        //B
        positions_tex_quad.push_back(b.xmin()-x);
        positions_tex_quad.push_back(b.ymax()-y);
        positions_tex_quad.push_back(0);


        //C
        positions_tex_quad.push_back(b.xmax()-x);
        positions_tex_quad.push_back(b.ymax()-y);
        positions_tex_quad.push_back(0);



        //A
        positions_tex_quad.push_back(b.xmin()-x);
        positions_tex_quad.push_back(b.ymin()-y);
        positions_tex_quad.push_back(0);


        //C
        positions_tex_quad.push_back(b.xmax()-x);
        positions_tex_quad.push_back(b.ymax()-y);
        positions_tex_quad.push_back(0);


        //D
        positions_tex_quad.push_back(b.xmax()-x);
        positions_tex_quad.push_back(b.ymin()-y);
        positions_tex_quad.push_back(0);



        texture_map.push_back(0.0);
        texture_map.push_back(0.0);

        texture_map.push_back(0.0);
        texture_map.push_back(1.0);

        texture_map.push_back(1.0);
        texture_map.push_back(1.0);

        texture_map.push_back(0.0);
        texture_map.push_back(0.0);

        texture_map.push_back(1.0);
        texture_map.push_back(1.0);

        texture_map.push_back(1.0);
        texture_map.push_back(0.0);



    }
    //The grid
    {
      double dx((b.xmax()-b.xmin())/10.0), dy((b.ymax()-b.ymin())/10.0);
        for(int u = 0; u < 11; u++)
        {

            positions_grid.push_back(b.xmin()-x + dx* u);
            positions_grid.push_back(b.ymin()-y);
            positions_grid.push_back(0);

            positions_grid.push_back(b.xmin()-x + dx* u);
            positions_grid.push_back(b.ymax()-y);
            positions_grid.push_back(0);
        }
        for(int v=0; v<11; v++)
        {

            positions_grid.push_back(b.xmin()-x);
            positions_grid.push_back(b.ymin()-y + v * dy);
            positions_grid.push_back(0);

            positions_grid.push_back(b.xmax()-x);
            positions_grid.push_back(b.ymin()-y + v * dy);
            positions_grid.push_back(0);
        }

    }
    //the Box
    {
      const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();

        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);

        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmin()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmin()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymax()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);


        positions_cube.push_back(b.xmax()+offset.x);
        positions_cube.push_back(b.ymin()+offset.y);
        positions_cube.push_back(b.zmax()+offset.z);

    }

    //The texture
    for( int i=0 ; i < item->getTriangleContainer(0)->getTextureSize().width() ; i++ )
    {
        for( int j=0 ; j < item->getTriangleContainer(0)->getTextureSize().height() ; j++)
        {
            compute_texture(i,j);
        }
    }
}

Scene_implicit_function_item::
Scene_implicit_function_item(Implicit_function_interface* f)
{
  d = new Scene_implicit_function_item_priv(f, this);
  //Generates an integer which will be used as ID for each buffer
  d->compute_min_max();
  setTriangleContainer(0, new Tc(Vi::PROGRAM_WITH_TEXTURE, false));
  setEdgeContainer(1, new Ec(Vi::PROGRAM_NO_SELECTION, false));
  setEdgeContainer(0, new Ec(Vi::PROGRAM_NO_SELECTION, false));
  getTriangleContainer(0)->setTextureSize(QSize(d->grid_size_-1,d->grid_size_-1));
  connect(d->frame_, SIGNAL(modified()), this, SLOT(plane_was_moved()));
  plane_was_moved();
  invalidateOpenGLBuffers();
}


Scene_implicit_function_item::~Scene_implicit_function_item()
{
  delete d;
}


void
Scene_implicit_function_item::compute_bbox() const
{
    setBbox(d->function_->bbox());
}

void
Scene_implicit_function_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  if(!d->init)
  {
    compute_function_grid();
    double offset_x = (bbox().xmin() + bbox().xmax()) / 2;
    double offset_y = (bbox().ymin() + bbox().ymax()) / 2;
    double offset_z = (bbox().zmin() + bbox().zmax()) / 2;
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    d->frame_->setPosition(offset_x+offset.x, offset_y+offset.y, offset_z+offset.z);
    d->frame_->setOrientation(1., 0, 0, 0);
    d->init = true;
  }

  if(d->frame_->isManipulated()) {
      if(d->need_update_) {
          compute_function_grid();
          d->need_update_ = false;
      }
  }

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

    QMatrix4x4 f_mat;
    GLdouble d_mat[16];
    d->frame_->getMatrix(d_mat);
    //Convert the GLdoubles matrices in GLfloats
    for (int i=0; i<16; ++i){
        f_mat.data()[i] = GLfloat(d_mat[i]);
    }
    getTriangleContainer(0)->setFrameMatrix(f_mat);
    getTriangleContainer(0)->draw(viewer, true);
}

void
Scene_implicit_function_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
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

    getEdgeContainer(0)->setColor(QColor(0,0,0));
    getEdgeContainer(0)->draw(viewer, true);
    QMatrix4x4 f_mat;
    GLdouble d_mat[16];
    d->frame_->getMatrix(d_mat);
    //Convert the GLdoubles matrices in GLfloats
    for (int i=0; i<16; ++i){
        f_mat.data()[i] = double(d_mat[i]);
    }
    getEdgeContainer(1)->setFrameMatrix(f_mat);
    getEdgeContainer(1)->setColor(QColor(155, 155, 155));
    getEdgeContainer(1)->draw(viewer, true);
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

void Scene_implicit_function_item_priv::compute_texture(int i, int j)
{
    double v = (implicit_grid_[i][j]).second;

    if(is_nan(v)) {
        item->getTriangleContainer(0)->setTextureData(i,j,51,51,51);
    } else
        // determines grey level
        if ( v > 0 )
        {
            v = v/max_value_;
            GLdouble r = red_color_ramp_.r(v), g = red_color_ramp_.g(v), b = red_color_ramp_.b(v);
            item->getTriangleContainer(0)->setTextureData(i,j,255*r,255*g,255*b);
        }
        else
        {
            v = v/min_value_;
            GLdouble r = blue_color_ramp_.r(v), g = blue_color_ramp_.g(v), b = blue_color_ramp_.b(v);
            item->getTriangleContainer(0)->setTextureData(i,j,255*r,255*g,255*b);
        }
}


void
Scene_implicit_function_item::
compute_function_grid() const
{
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    typedef CGAL::Simple_cartesian<double>  K;
    typedef K::Aff_transformation_3         Aff_transformation;
    typedef K::Point_3                      Point_3;

    // Get transformation
    const GLdouble* m = d->frame_->matrix();

    // OpenGL matrices are row-major matrices
    Aff_transformation t (m[0], m[4], m[8], m[12]-offset.x,
            m[1], m[5], m[9], m[13]-offset.y,
            m[2], m[6], m[10], m[14]-offset.z);

    double diag = (CGAL::sqrt((bbox().xmax()-bbox().xmin())*(bbox().xmax()-bbox().xmin()) + (bbox().ymax()-bbox().ymin())*(bbox().ymax()-bbox().ymin())  + (bbox().zmax()-bbox().zmin())*(bbox().zmax()-bbox().zmin()) )) * .6;

    const double dx = diag;
    const double dy = diag;
    const double z (0);

    int nb_quad = d->grid_size_ - 1;

    for(int i=0 ; i<d->grid_size_ ; ++i)
    {
        double x = -diag/2. + double(i)/double(nb_quad) * dx;

        for(int j=0 ; j<d->grid_size_ ; ++j)
        {
            double y = -diag/2. + double(j)/double(nb_quad) * dy;

            Point_3 query = t( Point_3(x, y, z) );
            double v = d->function_->operator()(query.x(), query.y(), query.z());

            d->implicit_grid_[i][j] = Scene_implicit_function_item_priv::Point_value(Scene_implicit_function_item_priv::Point(query.x(),query.y(),query.z()),v);
        }
    }

    // Update
    const_cast<Scene_implicit_function_item*>(this)->invalidateOpenGLBuffers();

}

void
Scene_implicit_function_item_priv::
compute_min_max()
{
    if(function_->get_min_max(min_value_, max_value_))
        return;

    double probes_nb = double(grid_size_) / 2;

    // Probe bounding box
    const CGAL::Three::Scene_item::Bbox& b = item->bbox();

    for ( int i = 0 ; i <= probes_nb ; ++i )
    {
        double x = b.xmin() + double(i) * (b.xmax() - b.xmin()) / probes_nb;

        for ( int j = 0 ; j <= probes_nb ; ++j )
        {
            double y = b.ymin() + double(j) * (b.ymax() - b.ymin()) / probes_nb;

            for ( int k = 0 ; k <= probes_nb ; ++k )
            {
                double z = b.zmin() + double(k) * (b.zmax() - b.zmin()) / probes_nb;

                double v = (*function_)(x,y,z);
                if(is_nan(v)) continue;
                max_value_ = (std::max)(v,max_value_);
                min_value_ = (std::min)(v,min_value_);
            }
        }
    }
}

void
Scene_implicit_function_item::invalidateOpenGLBuffers()
{
  compute_bbox();
  setBuffersFilled(false);
  getTriangleContainer(0)->reset_vbos(ALL);
  getEdgeContainer(0)->reset_vbos(ALL);
  getEdgeContainer(1)->reset_vbos(ALL);
}


void Scene_implicit_function_item::updateCutPlane()
{ // just handle deformation - paint like selection is handled in eventFilter()
  if(d->need_update_) {
    compute_function_grid();
    d->need_update_= false;
  }
}

Implicit_function_interface* Scene_implicit_function_item::function() const { return d->function_; }
Scene_implicit_function_item::ManipulatedFrame* Scene_implicit_function_item::manipulatedFrame() { return d->frame_; }
void Scene_implicit_function_item::plane_was_moved() { d->need_update_ = true; QTimer::singleShot(0, this, SLOT(updateCutPlane()));}

void Scene_implicit_function_item::initializeBuffers(Viewer_interface *v) const
{
  d->initialize_buffers(v);
}

void Scene_implicit_function_item::computeElements() const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  d->compute_vertices_and_texmap();
  getTriangleContainer(0)->allocate(
        Tc::Flat_vertices,
        d->positions_tex_quad.data(),
        static_cast<int>(d->positions_tex_quad.size()*sizeof(float)));
  d->nb_quad = d->positions_tex_quad.size();

  getTriangleContainer(0)->allocate(
        Tc::Texture_map,
        d->texture_map.data(),
        static_cast<int>(d->texture_map.size()*sizeof(float)));

  getEdgeContainer(0)->allocate(
        Ec::Vertices,
        d->positions_cube.data(),
        static_cast<int>(d->positions_cube.size()*sizeof(float)));
  d->nb_cube= d->positions_cube.size();

  getEdgeContainer(1)->allocate(
        Ec::Vertices,
        d->positions_grid.data(),
        static_cast<int>(d->positions_grid.size()*sizeof(float)));
  d->nb_grid= d->positions_grid.size();
  setBuffersFilled(true);
  QApplication::restoreOverrideCursor();
}
