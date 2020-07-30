#include "config.h"

#include <boost/config.hpp>
#if defined(BOOST_MSVC)
#  pragma warning( disable : 4503)
#endif

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include "Scene_c3t3_item.h"

#include <fstream>
#include <map>
#include <set>
#include <cmath>

#include <QFileInfo>
#include <QFileDialog>
#include <QAction>
#include <QMainWindow>
#include <QColor>
#include <QPixmap>
#include <QBitmap>
#include "ui_Rib_dialog.h"


// Constants
#define CGAL_RIB_NON_TRANSPARENT_MATERIAL_ALPHA 1
using namespace CGAL::Three;

class C3t3_rib_exporter_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "c3t3_rib_exporter_plugin.json")

public:
  C3t3_rib_exporter_plugin();
  virtual ~C3t3_rib_exporter_plugin() {}

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionCreateRib;
  }
  bool applicable(QAction*)const{
    Scene_c3t3_item* item = qobject_cast<Scene_c3t3_item*>(scene->item(scene->mainSelectionIndex()));
    return item && item->is_valid();}

public Q_SLOTS:
  void create_rib();
  void height_changed(int i);
  void width_changed(int i);

private:
  typedef Tr::Bare_point                    Bare_point;
  typedef Tr::Weighted_point                Weighted_point;

  typedef Geom_traits::Vector_3             Vector_3;
  typedef Geom_traits::Plane_3              Plane;
  typedef Geom_traits::FT                   FT;
  typedef Geom_traits::Aff_transformation_3 Aff_transformation_3;

  typedef CGAL::qglviewer::Vec qglVec;

  enum Rib_exporter_mode { CUT=0, MESH, TRIANGULATION };

  struct Rib_exporter_parameters
  {
    // Materials
    double sphere_radius;
    double cylinder_radius;

    // Lights
    bool ambientOn;
    double ambientIntensity;
    bool shadowOn;
    double shadowIntensity;

    // Picture
    int width;
    int height;
    Rib_exporter_mode mode;
    bool is_preview;
  };

private:
  void update_mask();

  bool get_parameters_from_dialog();

  QStringList nameFilters() const;
  bool save(const Scene_c3t3_item&, const QFileInfo& fileinfo);
  void init_maps(const C3t3& c3t3, const QColor& color);
  void init_point_radius(const C3t3& c3t3);
  void init_parameters();

  Bare_point camera_coordinates(const Bare_point& p);
  void fill_points_and_edges_map(const C3t3& c3t3);

  void add_edge(const Bare_point& p, const Bare_point& q, const QColor& color);
  void add_vertex(const Bare_point& p, const QColor& color);

  void write_header(const std::string& filename, std::ofstream& out);

  void write_lights(std::ofstream& out);
  void write_turn_background_light(bool turn_on, std::ofstream& out);

  void write_facets(const C3t3& c3t3, std::ofstream& out);
  void write_facets(const C3t3& c3t3, const Plane& plane, std::ofstream& out);
  void write_surface_cells(const C3t3& c3t3, const Plane& plane, std::ofstream& out);
  void write_cells_intersecting_a_plane(const C3t3& c3t3, const Plane& plane, std::ofstream& out);
  void write_cells_on_the_positive_side_of_a_plane(const C3t3& c3t3, const Plane& plane, std::ofstream& out);

  void write_triangle(const Bare_point& p, const Bare_point& q, const Bare_point& r,
                      const QColor& color, const QColor& edge_color, std::ofstream& out);
  void write_tetrahedron (const Bare_point& p, const Bare_point& q, const Bare_point& r, const Bare_point& s,
                          const QColor& color, const QColor& edge_color, std::ofstream& out);

  void write_point(const Bare_point& p, std::ofstream& out);
  void write_point_sphere(const Bare_point& p, std::ofstream& out);

  void write_edge_cylinder(const Bare_point& p, const Bare_point& q, std::ofstream& out);

  // Writes data which has been stored during triangle drawing
  void write_edges_flat(std::ofstream& out);
  void write_edges_volumic(std::ofstream& out);
  void write_vertices_volumic(std::ofstream& out);

  void write_color(const QColor& color, bool use_transparency, std::ofstream& out);
  void write_opacity(const double alpha, std::ofstream& out);

  // Background
  void write_background(const QColor& color, std::ofstream& out);

private:
  QAction* actionCreateRib;

  // Viewer

  typedef std::map<C3t3::Surface_patch_index, QColor> Surface_map;
  typedef std::map<C3t3::Subdomain_index, QColor> Subdomain_map;

  Surface_map surface_map_;
  Subdomain_map subdomain_map_;

  typedef std::map<std::pair<Bare_point, Bare_point>,QColor> Edge_map;
  typedef std::map<Bare_point, QColor> Vertex_map;

  Edge_map edges_;
  Vertex_map vertices_;

  double zmax_;
  double diag_;

  // Cache data to avoid writing too much lines in rib file
  QColor prev_color_;
  double prev_alpha_;
  const Scene_c3t3_item* prev_c3t3_;

  Rib_exporter_parameters parameters_;
};


C3t3_rib_exporter_plugin::
C3t3_rib_exporter_plugin()
  : actionCreateRib(NULL)
  , zmax_(0)
  , diag_(0)
  , prev_color_(0,0,0)
  , prev_alpha_(1)
  , prev_c3t3_(NULL)
{

}


void
C3t3_rib_exporter_plugin::
init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
{
  this->scene = scene_interface;
  this->mw = mainWindow;

  actionCreateRib = new QAction("Export C3t3 to RIB", mw);
  if( NULL != actionCreateRib )
  {
    actionCreateRib->setProperty("subMenuName", "Tetrahedral Mesh Generation");
    connect(actionCreateRib, SIGNAL(triggered()), this, SLOT(create_rib()));
  }
  init_parameters();
}


void
C3t3_rib_exporter_plugin::create_rib()
{
  if ( NULL == Three::activeViewer() )
  {
    std::cerr << "Can't find viewer" << std::endl;
    return;
  }

  // Get Scene_c3t3_item
  Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_c3t3_item* c3t3_item =
    qobject_cast<Scene_c3t3_item*>(scene->item(index));

  if ( NULL == c3t3_item )
  {
    return;
  }

  // Init data
  //if ( c3t3_item != prev_c3t3_ ) // Commented because it was causing problems
                                   // when changing the color of the c3t3
  {
    init_maps(c3t3_item->c3t3(), c3t3_item->color());
    init_point_radius(c3t3_item->c3t3());
    init_parameters();

    prev_c3t3_ = c3t3_item;
  }

  // Get parameters from user dialog
  if ( !get_parameters_from_dialog() )
  {

    QBitmap bitmap;
    bitmap.clear();
    Three::activeViewer()->setMask(bitmap);
    return;
  }

  // Disable Mask
  QBitmap bitmap;
  bitmap.clear();
  Three::activeViewer()->setMask(bitmap);

  // Save dialog
  QStringList filters;
  filters << nameFilters();
  filters << tr("All files (*)");

  QString filename = QFileDialog::getSaveFileName(mw,
                                                  tr("Save to File..."),
                                                  QString(),
                                                  filters.join(";;"));

  QFileInfo fileinfo(filename);

  // Save rib file
  save(*c3t3_item,fileinfo);

  std::cout << "Rib file created successfully" << std::endl;
}


void
C3t3_rib_exporter_plugin::
height_changed(int i)
{
  parameters_.height = i;
  update_mask();
}

void
C3t3_rib_exporter_plugin::
width_changed(int i)
{
  parameters_.width = i;
  update_mask();
}


void
C3t3_rib_exporter_plugin::
update_mask()
{
  double ratio = double(parameters_.width) / double(parameters_.height);

  if ( NULL == Three::activeViewer() )
  {
    std::cerr << "Can't find viewer..." << std::endl;
    return;
  }
  QBitmap bitmap;
  bitmap.setDevicePixelRatio(ratio);
  bitmap.fill();
  Three::activeViewer()->setMask(bitmap);
}



bool
C3t3_rib_exporter_plugin::
get_parameters_from_dialog()
{
  QDialog dialog(mw);
  Ui::Rib_dialog ui;
  ui.setupUi(&dialog);
  ui.ambientIntensity->setMaximum(1.0);
  ui.shadowIntensity->setMaximum(1.0);

  connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));
  connect(ui.resWidth, SIGNAL(valueChanged(int)), this, SLOT(width_changed(int)));
  connect(ui.resHeight, SIGNAL(valueChanged(int)), this, SLOT(height_changed(int)));

  // -----------------------------------
  // Set data
  // -----------------------------------

  // Materials
  ui.sphereRadius->setValue(parameters_.sphere_radius);
  ui.cylinderRadius->setValue(parameters_.cylinder_radius);

  // Lights
  ui.isAmbientOn->setChecked(parameters_.ambientOn);
  ui.ambientIntensity->setValue(parameters_.ambientIntensity);
  ui.isShadowOn->setChecked(parameters_.shadowOn);
  ui.shadowIntensity->setValue(parameters_.shadowIntensity);

  // Picture
  QStringList mode_list;
  mode_list << "Export Cut (draws current cut view)"
            << "Export Mesh (draws all surface facets)"
            << "Export Triangulation (draws all points and edges)";

  ui.exportMode->insertItems(0,mode_list);

  ui.resWidth->setValue(parameters_.width);
  ui.resHeight->setValue(parameters_.height);
  ui.exportMode->setCurrentIndex(static_cast<int>(parameters_.mode));
  ui.isPreview->setChecked(parameters_.is_preview);

  // Update mask

  update_mask();

  // -----------------------------------
  // Get data
  // -----------------------------------
  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return false;

  // Materials
  parameters_.sphere_radius = ui.sphereRadius->value();
  parameters_.cylinder_radius = ui.cylinderRadius->value();

  // Lights
  parameters_.ambientOn = ui.isAmbientOn->isChecked();
  parameters_.ambientIntensity = ui.ambientIntensity->value();
  parameters_.shadowOn = ui.isShadowOn->isChecked();
  parameters_.shadowIntensity = ui.shadowIntensity->value();

  // Picture
  parameters_.width = ui.resWidth->value();
  parameters_.height = ui.resHeight->value();
  parameters_.mode = static_cast<Rib_exporter_mode>(ui.exportMode->currentIndex());
  parameters_.is_preview = ui.isPreview->isChecked();

  return true;
}


QStringList
C3t3_rib_exporter_plugin::nameFilters() const
{
  return QStringList() << "RenderMan file (*.rib)";
}


bool
C3t3_rib_exporter_plugin::
save(const Scene_c3t3_item& c3t3_item, const QFileInfo& fileInfo)
{
  QString path = fileInfo.absoluteFilePath();
  std::ofstream rib_file (qPrintable(path));
  rib_file.precision(8);

  // Header
  QString basename = fileInfo.baseName();
  write_header(qPrintable(basename), rib_file);

  // Lights
  write_lights(rib_file);

  // Triangles
  switch ( parameters_.mode )
  {
    case CUT:
      rib_file << "Surface \"plastic\" \"Ka\" 0.65 \"Kd\" 0.85 \"Ks\" 0.25 \"roughness\" 0.1" << std::endl;
      write_facets(c3t3_item.c3t3(), c3t3_item.plane(), rib_file);

      rib_file << "Surface \"plastic\" \"Ka\" 0.65 \"Kd\" 0.65 \"Ks\" 0.35 \"roughness\" 0.2" << std::endl;
      //write_cells_intersecting_a_plane(c3t3_item.c3t3(), c3t3_item.plane(), rib_file);
      write_cells_on_the_positive_side_of_a_plane(c3t3_item.c3t3(), c3t3_item.plane(), rib_file);
      break;

    case MESH:
      rib_file << "Surface \"plastic\" \"Ka\" 0.65 \"Kd\" 0.85 \"Ks\" 0.25 \"roughness\" 0.1" << std::endl;
      //write_facets(c3t3_item.c3t3(), rib_file);
      write_surface_cells(c3t3_item.c3t3(), c3t3_item.plane(), rib_file);
      break;

    case TRIANGULATION:
      fill_points_and_edges_map(c3t3_item.c3t3());
      break;

    default:
      std::cerr << "Unexpected mode found" << std::endl;
      return false;
      break;
  }

  // Edges and vertices
  rib_file << "Surface \"plastic\" \"Ka\" 0.65 \"Kd\" 0.85 \"Ks\" 0.25 \"roughness\" 0.1" << std::endl;
  write_edges_volumic(rib_file);
  write_vertices_volumic(rib_file);

  // Background
  write_background(QColor(255,255,255), rib_file);

  rib_file << "WorldEnd" << std::endl;

  return true;
}

void
C3t3_rib_exporter_plugin::init_maps(const C3t3& c3t3, const QColor& color)
{
  surface_map_.clear();
  subdomain_map_.clear();
  edges_.clear();
  vertices_.clear();

  // Fill maps with 0 as value
  for ( C3t3::Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin(),
       fend = c3t3.facets_in_complex_end() ; fit != fend ; ++fit )
  {
    surface_map_.insert(std::make_pair(c3t3.surface_patch_index(*fit),QColor(0,0,0)));
  }

  for ( C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin(),
       cend = c3t3.cells_in_complex_end() ; cit != cend ; ++cit )
  {
    subdomain_map_.insert(std::make_pair(c3t3.subdomain_index(cit),QColor(0,0,0)));
  }

  // Fill value of maps
  size_t nb_colors = subdomain_map_.size(); // + surface_map_.size();

  // Starting hue
  double c = color.hueF();
  int i = 0;
  for ( Subdomain_map::iterator it = subdomain_map_.begin(), end = subdomain_map_.end();
       it != end ; ++it, ++i )
  {
    double hue = c + 1./nb_colors * i;
    if ( hue > 1 ) { hue -= 1.; }
    it->second = QColor::fromHsvF(hue, color.saturationF(), color.valueF());
  }
}


void
C3t3_rib_exporter_plugin::
init_point_radius(const C3t3& c3t3)
{
  const CGAL::Bbox_3 bbox = c3t3.bbox();

  const double xdelta = bbox.xmax() - bbox.xmin();
  const double ydelta = bbox.ymax() - bbox.ymin();
  const double zdelta = bbox.zmax() - bbox.zmin();

  diag_ = std::sqrt(xdelta*xdelta + ydelta*ydelta + zdelta*zdelta);

  parameters_.sphere_radius =  diag_ * 0.0015;
  parameters_.cylinder_radius =  diag_ * 0.00065;
}


void
C3t3_rib_exporter_plugin::
init_parameters()
{
  // Lights
  parameters_.ambientOn = true;
  parameters_.ambientIntensity = 0.20;
  parameters_.shadowOn = true;
  parameters_.shadowIntensity = 0.85;

  // Picture
  parameters_.width = 800;
  parameters_.height = 800;
  parameters_.mode = CUT;
  parameters_.is_preview = false;
}


C3t3_rib_exporter_plugin::Bare_point
C3t3_rib_exporter_plugin::
camera_coordinates(const Bare_point& p)
{
  qglVec p_vec ( p.x(), p.y(), p.z() );
  qglVec p_cam = Three::activeViewer()->camera()->cameraCoordinatesOf(p_vec);

  // Store maximal depth
  zmax_ = (std::max)(zmax_, double(-p_cam[2]));

  return Bare_point(p_cam[0],p_cam[1],p_cam[2]);
}


void
C3t3_rib_exporter_plugin::
fill_points_and_edges_map(const C3t3& c3t3)
{
  Geom_traits::Construct_point_3 wp2p
    = c3t3.triangulation().geom_traits().construct_point_3_object();

  for ( C3t3::Cells_in_complex_iterator it = c3t3.cells_in_complex_begin(),
       end = c3t3.cells_in_complex_end() ; it != end ; ++it )
  {
    const Bare_point& p1 = wp2p(it->vertex(0)->point());
    const Bare_point& p2 = wp2p(it->vertex(1)->point());
    const Bare_point& p3 = wp2p(it->vertex(2)->point());
    const Bare_point& p4 = wp2p(it->vertex(3)->point());

    const QColor& edge_color = subdomain_map_[c3t3.subdomain_index(it)];

    add_edge(p1,p2,edge_color);
    add_edge(p1,p3,edge_color);
    add_edge(p1,p4,edge_color);
    add_edge(p2,p3,edge_color);
    add_edge(p2,p4,edge_color);
    add_edge(p3,p4,edge_color);

    add_vertex(p1,edge_color);
    add_vertex(p2,edge_color);
    add_vertex(p3,edge_color);
    add_vertex(p4,edge_color);
  }
}


void
C3t3_rib_exporter_plugin::
add_edge(const Bare_point& p, const Bare_point& q, const QColor& color)
{
  if ( p < q )
  {
    edges_.insert(std::make_pair(std::make_pair(p,q),color));
  }
  else
  {
    edges_.insert(std::make_pair(std::make_pair(q,p),color));
  }
}


void
C3t3_rib_exporter_plugin::
add_vertex(const Bare_point& p, const QColor& color)
{
  vertices_.insert(std::make_pair(p,color));
}


void
C3t3_rib_exporter_plugin::
write_header(const std::string& filename, std::ofstream& out)
{
  out << "Option \"limits\" \"numthreads\" [16]" << std::endl
      << "Option \"searchpath\" \"shader\" \".:./shaders:%PIXIE_SHADERS%:%PIXIEHOME%/shaders\"" << std::endl;

  if ( ! parameters_.is_preview )
  {
    out << "Attribute \"visibility\" \"specular\" 1" << std::endl
        << "Attribute \"visibility\" \"transmission\" 1" << std::endl << std::endl;
  }

  out << "Display \""<< filename << ".tif\" \"file\" \"rgb\"" << std::endl;

  if ( ! parameters_.is_preview )
  {
    out << "Format " << parameters_.width << " " << parameters_.height << " 1" << std::endl;
  }
  else
  {
    double ratio = double(parameters_.height) / double(parameters_.width);

    int width = (ratio < 1.)  ? 300               : int(300. / ratio);
    int height = (ratio < 1.) ? int(ratio * 300.) : 300;

    out << "Format " << width << " " << height << " 1" << std::endl;
  }


  if ( parameters_.width > parameters_.height )
  {
    double ratio = double(parameters_.height) / double(parameters_.width);
    out << "ScreenWindow -1 1 " << -ratio << " " << ratio << std::endl;
  }
  else if ( parameters_.height > parameters_.width )
  {
    double ratio = double(parameters_.width) / double(parameters_.height);
    out << "ScreenWindow " << -ratio << " " << ratio << " -1 1" << std::endl;
  }

  out << "Projection \"perspective\" \"fov\" 45" << std::endl
      << "PixelSamples 4 4" << std::endl
      << "PixelFilter \"catmull-rom\" 3 3" << std::endl
      << "ShadingInterpolation \"smooth\"" << std::endl
      << "Rotate 180 0 0 1" << std::endl
      << "WorldBegin" << std::endl;
}


void
C3t3_rib_exporter_plugin::
write_lights(std::ofstream& out)
{
  if ( ! parameters_.is_preview )
  {
    // ShadowLight
    out << "LightSource \"shadowdistant\" 1 \"from\" [0 0 0] \"to\" [0 0 1]"
        << " \"shadowname\" \"raytrace\" \"intensity\" " << parameters_.shadowIntensity << std::endl;

    // Ambient light
    out << "LightSource \"ambientlight\" 2 \"intensity\" " << parameters_.ambientIntensity << std::endl;
  }
  else
  {
    out << "LightSource \"distantLight\" 1 \"from\" [0 0 0] \"to\" [0 0 1]"
    << " \"intensity\" 0.85" << std::endl;
  }

  // Background light
  out << "LightSource \"ambientlight\" 99 \"intensity\" 1" << std::endl;
  write_turn_background_light(false,out);
}


void
C3t3_rib_exporter_plugin::
write_turn_background_light(bool turn_on, std::ofstream& out)
{
  if (!turn_on)
  {
    out << "Illuminate 1 1" << std::endl;
    if ( ! parameters_.is_preview ) { out << "Illuminate 2 1" << std::endl; }
    out << "Illuminate 99 0" << std::endl;
  }
  else
  {
    out << "Illuminate 1 0" << std::endl;
    if ( ! parameters_.is_preview ) { out << "Illuminate 2 0" << std::endl; }
    out << "Illuminate 99 1" << std::endl;
  }
}


void
C3t3_rib_exporter_plugin::
write_facets(const C3t3& c3t3, std::ofstream& out)
{
  Geom_traits::Construct_point_3 wp2p
    = c3t3.triangulation().geom_traits().construct_point_3_object();

  for ( C3t3::Facets_in_complex_iterator it = c3t3.facets_in_complex_begin(),
       end = c3t3.facets_in_complex_end() ; it != end ; ++it )
  {
    const C3t3::Cell_handle& c = it->first;
    const int& k = it->second;

    const Bare_point& p1 = wp2p(c->vertex((k+1)&3)->point());
    const Bare_point& p2 = wp2p(c->vertex((k+2)&3)->point());
    const Bare_point& p3 = wp2p(c->vertex((k+3)&3)->point());

    QColor color = c3t3.is_in_complex(c) ? subdomain_map_[c3t3.subdomain_index(c)]
    : subdomain_map_[c3t3.subdomain_index(c->neighbor(k))];

    write_triangle(p1, p2, p3, color, color.darker(125), out );
  }
}


void
C3t3_rib_exporter_plugin::
write_facets(const C3t3& c3t3, const Plane& plane, std::ofstream& out)
{
  typedef EPICK::Oriented_side Side;

  Geom_traits::Construct_point_3 wp2p
    = c3t3.triangulation().geom_traits().construct_point_3_object();

  for ( C3t3::Facets_in_complex_iterator it = c3t3.facets_in_complex_begin(),
       end = c3t3.facets_in_complex_end() ; it != end ; ++it )
  {
    const C3t3::Cell_handle& c = it->first;
    const int& k = it->second;

    const Bare_point& p1 = wp2p(c->vertex((k+1)&3)->point());
    const Bare_point& p2 = wp2p(c->vertex((k+2)&3)->point());
    const Bare_point& p3 = wp2p(c->vertex((k+3)&3)->point());

    const Side s1 = plane.oriented_side(p1);
    const Side s2 = plane.oriented_side(p2);
    const Side s3 = plane.oriented_side(p3);

    if(   s1 == CGAL::ON_NEGATIVE_SIDE && s2 == CGAL::ON_NEGATIVE_SIDE
       && s3 == CGAL::ON_NEGATIVE_SIDE )
    {
      QColor color = c3t3.is_in_complex(c) ? subdomain_map_[c3t3.subdomain_index(c)]
                                           : subdomain_map_[c3t3.subdomain_index(c->neighbor(k))];

      write_triangle(p1, p2, p3, color, color.darker(125), out );
    }
  }
}

void
C3t3_rib_exporter_plugin::
write_surface_cells(const C3t3& c3t3, const Plane& /* plane */, std::ofstream& out)
{
  const Geom_traits& gt = c3t3.triangulation().geom_traits();
  Geom_traits::Construct_point_3 wp2p = gt.construct_point_3_object();

  for ( C3t3::Cells_in_complex_iterator it_cell = c3t3.cells_in_complex_begin(),
       end = c3t3.cells_in_complex_end() ; it_cell != end ; ++it_cell )
  {
    C3t3::Cell_handle c = it_cell;

    // Compute the number of surface vertices in c
    // Keep one of them in memory
    int num_2D_vertices = 0;
    int last_2D_vertex_index = -1;
    for (int i = 0 ; i < 4 ; ++i)
    {
      if (c3t3.in_dimension(c->vertex(i)) == 2)
      {
        ++num_2D_vertices;
        last_2D_vertex_index = i;
      }
    }

    //const int TRANSPARENCY_ALPHA_VALUE = 100;
    CGAL::Bbox_3 bbox = c3t3.bbox();
    float relPos = static_cast<float>((c->weighted_circumcenter(gt).x() - bbox.xmin())
                                      / (bbox.xmax() - bbox.xmin()));
    float TRANSPARENCY_ALPHA_VALUE =
      1.f -
      (relPos < 0.25f ?
        0.0f :
        (relPos > 0.75f ?
          1.0f :
          (relPos - 0.25f)*1.0f/0.5f
        )
      );


    /*
    // ONLY SURFACE FACETS ARE TRANSPARENT
    if (num_2D_vertices >= 3)
    {
      QColor basecolor = subdomain_map_[c3t3.subdomain_index(c)];
      QColor facecolor = basecolor.darker(150);
      QColor edgecolor = facecolor.darker(150);

      for (int i = 0 ; i < 4 ; ++i)
      {
        if (c3t3.in_dimension(c->vertex((i+1)%4)) == 2
            && c3t3.in_dimension(c->vertex((i+2)%4)) == 2
            && c3t3.in_dimension(c->vertex((i+3)%4)) == 2)
        {
          edgecolor.setAlpha(TRANSPARENCY_ALPHA_VALUE);
        }
        else
        {
          edgecolor.setAlpha(255);
        }

        write_triangle(wp2p(c->vertex((i+1)%4)->point()),
                       wp2p(c->vertex((i+2)%4)->point()),
                       wp2p(c->vertex((i+3)%4)->point()),
                       facecolor, edgecolor, out );
      }
    }*/


    // SURFACE CELLS ARE TRANSPARENT
    if (num_2D_vertices >= 2)
    {
      QColor basecolor = subdomain_map_[c3t3.subdomain_index(c)];
      QColor facecolor = basecolor.darker(150);
      QColor edgecolor = facecolor.darker(150);

      /*
      typedef Kernel::Oriented_side Side;

      // Transparency on the negative side of the plane
      const Side s0 = plane.oriented_side(wp2p(c->vertex(0)->point()));
      const Side s1 = plane.oriented_side(wp2p(c->vertex(1)->point()));
      const Side s2 = plane.oriented_side(wp2p(c->vertex(2)->point()));
      const Side s3 = plane.oriented_side(wp2p(c->vertex(3)->point()));
      if(   s0 == CGAL::ON_NEGATIVE_SIDE && s1 == CGAL::ON_NEGATIVE_SIDE
         && s2 == CGAL::ON_NEGATIVE_SIDE && s3 == CGAL::ON_NEGATIVE_SIDE )
      {
        edgecolor.setAlpha(TRANSPARENCY_ALPHA_VALUE);
      }
      else
      {
        edgecolor.setAlpha(255);
      }*/

      edgecolor.setAlphaF(TRANSPARENCY_ALPHA_VALUE);

      for (int i = 0 ; i < 4 ; ++i)
      {
        write_triangle(wp2p(c->vertex((i+1)%4)->point()),
                       wp2p(c->vertex((i+2)%4)->point()),
                       wp2p(c->vertex((i+3)%4)->point()),
                       facecolor, edgecolor, out );
      }
    }
    else if (num_2D_vertices == 1)
    {
      QColor basecolor = subdomain_map_[c3t3.subdomain_index(c)];
      QColor facecolor = basecolor.darker(150);
      QColor edgecolor = facecolor.darker(150);

      for (int i = 0 ; i < 4 ; ++i)
      {
        if (i == last_2D_vertex_index)
        {
          edgecolor.setAlphaF(TRANSPARENCY_ALPHA_VALUE);
        }
        else
        {
          edgecolor.setAlphaF(TRANSPARENCY_ALPHA_VALUE);
        }

        write_triangle(wp2p(c->vertex((i+1)%4)->point()),
                       wp2p(c->vertex((i+2)%4)->point()),
                       wp2p(c->vertex((i+3)%4)->point()),
                       facecolor, edgecolor, out );
      }
    }
  }
}


void
C3t3_rib_exporter_plugin::
write_cells_intersecting_a_plane(const C3t3& c3t3, const Plane& plane, std::ofstream& out)
{
  typedef EPICK::Oriented_side Side;

  Geom_traits::Construct_point_3 wp2p
    = c3t3.triangulation().geom_traits().construct_point_3_object();

  for ( C3t3::Cells_in_complex_iterator it = c3t3.cells_in_complex_begin(),
       end = c3t3.cells_in_complex_end() ; it != end ; ++it )
  {
    const Bare_point& p1 = wp2p(it->vertex(0)->point());
    const Bare_point& p2 = wp2p(it->vertex(1)->point());
    const Bare_point& p3 = wp2p(it->vertex(2)->point());
    const Bare_point& p4 = wp2p(it->vertex(3)->point());

    const Side s1 = plane.oriented_side(p1);
    const Side s2 = plane.oriented_side(p2);
    const Side s3 = plane.oriented_side(p3);
    const Side s4 = plane.oriented_side(p4);

    if(   s1 == CGAL::ON_ORIENTED_BOUNDARY || s2 == CGAL::ON_ORIENTED_BOUNDARY
       || s3 == CGAL::ON_ORIENTED_BOUNDARY || s4 == CGAL::ON_ORIENTED_BOUNDARY
       || s2 != s1 || s3 != s1 || s4 != s1 )
    {
      QColor basecolor = subdomain_map_[c3t3.subdomain_index(it)];
      QColor facecolor = basecolor.darker(150);
      QColor edgecolor = facecolor.darker(150);

      edgecolor.setAlpha(20);

      // Don't write facet twice
      if ( s1 != CGAL::ON_NEGATIVE_SIDE || s2 != CGAL::ON_NEGATIVE_SIDE || s3 != CGAL::ON_NEGATIVE_SIDE )
        write_triangle(p1, p2, p3, facecolor, edgecolor, out );

      if ( s1 != CGAL::ON_NEGATIVE_SIDE || s2 != CGAL::ON_NEGATIVE_SIDE || s4 != CGAL::ON_NEGATIVE_SIDE )
        write_triangle(p1, p2, p4, facecolor, edgecolor, out );

      if ( s1 != CGAL::ON_NEGATIVE_SIDE || s3 != CGAL::ON_NEGATIVE_SIDE || s4 != CGAL::ON_NEGATIVE_SIDE )
        write_triangle(p1, p3, p4, facecolor, edgecolor, out );

      if ( s2 != CGAL::ON_NEGATIVE_SIDE || s3 != CGAL::ON_NEGATIVE_SIDE || s4 != CGAL::ON_NEGATIVE_SIDE )
        write_triangle(p2, p3, p4, facecolor, edgecolor, out );
    }
  }
}

void
C3t3_rib_exporter_plugin::
write_cells_on_the_positive_side_of_a_plane(const C3t3& c3t3, const Plane& plane, std::ofstream& out)
{
  typedef EPICK::Oriented_side Side;

  Geom_traits::Construct_point_3 wp2p
    = c3t3.triangulation().geom_traits().construct_point_3_object();

  for ( C3t3::Cells_in_complex_iterator it = c3t3.cells_in_complex_begin(),
       end = c3t3.cells_in_complex_end() ; it != end ; ++it )
  {
    const Bare_point& p1 = wp2p(it->vertex(0)->point());
    const Bare_point& p2 = wp2p(it->vertex(1)->point());
    const Bare_point& p3 = wp2p(it->vertex(2)->point());
    const Bare_point& p4 = wp2p(it->vertex(3)->point());

    const Side s1 = plane.oriented_side(p1);
    const Side s2 = plane.oriented_side(p2);
    const Side s3 = plane.oriented_side(p3);
    const Side s4 = plane.oriented_side(p4);

    if( (   s1 == CGAL::ON_POSITIVE_SIDE || s2 == CGAL::ON_POSITIVE_SIDE
         || s3 == CGAL::ON_POSITIVE_SIDE || s4 == CGAL::ON_POSITIVE_SIDE )
        /*&&
        (   c3t3.surface_patch_index(it, 0) != C3t3::Surface_patch_index()
         || c3t3.surface_patch_index(it, 1) != C3t3::Surface_patch_index()
         || c3t3.surface_patch_index(it, 2) != C3t3::Surface_patch_index()
         || c3t3.surface_patch_index(it, 3) != C3t3::Surface_patch_index() )*/
    )
    {
      QColor basecolor = subdomain_map_[c3t3.subdomain_index(it)];
      QColor facecolor = basecolor.darker(150);
      QColor edgecolor = facecolor.darker(150);

      edgecolor.setAlpha(10);

      // Don't write facet twice
      if ( s1 != CGAL::ON_NEGATIVE_SIDE || s2 != CGAL::ON_NEGATIVE_SIDE || s3 != CGAL::ON_NEGATIVE_SIDE )
        write_triangle(p1, p2, p3, facecolor, edgecolor, out );

      if ( s1 != CGAL::ON_NEGATIVE_SIDE || s2 != CGAL::ON_NEGATIVE_SIDE || s4 != CGAL::ON_NEGATIVE_SIDE )
        write_triangle(p1, p2, p4, facecolor, edgecolor, out );

      if ( s1 != CGAL::ON_NEGATIVE_SIDE || s3 != CGAL::ON_NEGATIVE_SIDE || s4 != CGAL::ON_NEGATIVE_SIDE )
        write_triangle(p1, p3, p4, facecolor, edgecolor, out );

      if ( s2 != CGAL::ON_NEGATIVE_SIDE || s3 != CGAL::ON_NEGATIVE_SIDE || s4 != CGAL::ON_NEGATIVE_SIDE )
        write_triangle(p2, p3, p4, facecolor, edgecolor, out );
    }
  }
}


void
C3t3_rib_exporter_plugin::
write_triangle (const Bare_point& p, const Bare_point& q, const Bare_point& r,
                const QColor& color, const QColor& edge_color, std::ofstream& out)
{
  // Color
  write_color(color, true, out);

  // Triangle
  out << "Polygon \"P\" [";
  write_point(p,out);
  write_point(q,out);
  write_point(r,out);
  out << "]" << std::endl;

  // Edges (will be drawn later on)
  add_edge(p,q,edge_color);
  add_edge(p,r,edge_color);
  add_edge(q,r,edge_color);

  // Vertices (will be drawn later on)
  add_vertex(p,edge_color);
  add_vertex(q,edge_color);
  add_vertex(r,edge_color);
}

void
C3t3_rib_exporter_plugin::
write_tetrahedron (const Bare_point& p, const Bare_point& q, const Bare_point& r, const Bare_point& s,
                   const QColor& color, const QColor& edge_color, std::ofstream& out)
{
  // Color
  write_color(color, true, out);

  // Triangle
  out << "Polygon \"P\" [";
  write_point(p,out);
  write_point(q,out);
  write_point(r,out);
  write_point(s,out);
  out << "]" << std::endl;

  // Edges (will be drawn later on)
  add_edge(p,q,edge_color);
  add_edge(p,r,edge_color);
  add_edge(q,r,edge_color);
  add_edge(s,p,edge_color);
  add_edge(s,q,edge_color);
  add_edge(s,r,edge_color);

  // Vertices (will be drawn later on)
  add_vertex(p,edge_color);
  add_vertex(q,edge_color);
  add_vertex(r,edge_color);
  add_vertex(s,edge_color);
}

void
C3t3_rib_exporter_plugin::
write_point (const Bare_point& p, std::ofstream& out)
{
  // Transform point in camera coordinates
  const Bare_point& p_cam = camera_coordinates(p);

  // Write it
  out << " " << -p_cam.x() << " " << -p_cam.y() << " " << -p_cam.z() << " ";
}


void
C3t3_rib_exporter_plugin::
write_point_sphere(const Bare_point& p, std::ofstream& out)
{
  // Transform point in camera coordinates
  const Bare_point& p_cam = camera_coordinates(p);

  // radius
  const double& r = parameters_.sphere_radius;

  out << "Translate " << -p_cam.x() << " " << -p_cam.y() << " " << -p_cam.z() << std::endl;

  // Sphere radius zmin zmax thetamax
  out << "Sphere " << r << " " << -r << " " << r << " 360" << std::endl;
  out << "Identity" << std::endl;
}


void
C3t3_rib_exporter_plugin::
write_edge_cylinder(const Bare_point& p, const Bare_point& q, std::ofstream& out)
{
  // Transform point in camera coordinates
  const Bare_point& p_cam = camera_coordinates(p);
  const Bare_point& q_cam = camera_coordinates(q);

  double pq = CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(p_cam,q_cam)));

  Aff_transformation_3 t (CGAL::Translation(), Vector_3(p_cam,CGAL::ORIGIN));
  const Bare_point& q_cam_t = q_cam.transform(t);

  Vector_3 Oq (CGAL::ORIGIN,q_cam_t);
  Vector_3 Oz (FT(0),FT(0),FT(1));

  Vector_3 r_axis = CGAL::cross_product(Oq,Oz);
  double cos_angle = CGAL::to_double((Oq*Oz)/CGAL::sqrt(Oq.squared_length()));
  double angle = std::acos(cos_angle) * 180. / CGAL_PI;

  // radius
  const double& r = parameters_.cylinder_radius;

  out << "Translate " << -p_cam.x() << " " << -p_cam.y() << " " << -p_cam.z() << std::endl;
  out << "Rotate " << (angle+180.) << " " << -r_axis.x() << " " << -r_axis.y() << " " << -r_axis.z() << std::endl;

  // Cylinder radius zmin zmax thetamax
  out << "Cylinder " << r << " 0 " << pq << " 360" << std::endl;
  out << "Identity" << std::endl;
}


void
C3t3_rib_exporter_plugin::
write_edges_flat(std::ofstream& out)
{
  // Lights
  write_turn_background_light(true,out);

  out << "Surface \"constant\"" << std::endl;
  write_opacity(CGAL_RIB_NON_TRANSPARENT_MATERIAL_ALPHA, out);

  // Translation
  out << "Translate 0 0 -0.1" << std::endl;

  for ( Edge_map::iterator it = edges_.begin(), end = edges_.end() ;
       it != end ; ++it )
  {
    // Color
    write_color(it->second, true, out);

    // Edge
    out << "Curves \"linear\" [2] \"nonperiodic\" \"P\" [";
    write_point(it->first.first,out);
    write_point(it->first.second,out);
    out << "] \"constantwidth\" [0.15]" << std::endl;
  }
}


void
C3t3_rib_exporter_plugin::
write_edges_volumic(std::ofstream& out)
{
  // Material
  write_opacity(CGAL_RIB_NON_TRANSPARENT_MATERIAL_ALPHA, out);

  for ( Edge_map::iterator it = edges_.begin(), end = edges_.end() ;
       it != end ; ++it )
  {
    // Color
    write_color(it->second, true, out);
    // Edge
    write_edge_cylinder(it->first.first, it->first.second, out);
  }
}

void
C3t3_rib_exporter_plugin::
write_vertices_volumic(std::ofstream& out)
{
  // Material
  write_opacity(CGAL_RIB_NON_TRANSPARENT_MATERIAL_ALPHA, out);

  for ( Vertex_map::iterator it = vertices_.begin(), end = vertices_.end() ;
       it != end ; ++it )
  {
    // Color
    write_color(it->second, false, out);
    // Vertex
    write_point_sphere(it->first, out);
  }
}


void
C3t3_rib_exporter_plugin::
write_color(const QColor& color, bool use_transparency, std::ofstream& out)
{
  if ( prev_color_ == color )
  {
    return;
  }

  // Cache data
  prev_color_ = color;

  // Write opacity data
  if (use_transparency)
  {
    write_opacity(color.alphaF(),out);
  }

  // Write color data
  out << "Color [ " << color.redF() << " " << color.greenF() << " "
      << color.blueF() <<  " ]" << std::endl;
}


void
C3t3_rib_exporter_plugin::
write_opacity(const double alpha, std::ofstream& out)
{
  if ( alpha == prev_alpha_ )
  {
    return;
  }

  // Cache data
  prev_alpha_ = alpha;

  // Write opacity data
  out << "Opacity " << alpha << " " << alpha << " " << alpha << std::endl;
}


void
C3t3_rib_exporter_plugin::
write_background(const QColor& color, std::ofstream& out)
{
  write_turn_background_light(false,out);

  out << "Surface \"constant\"" << std::endl;
  write_color(color,false,out);

  double corner = zmax_ * 2.;
  double depth_pos = zmax_ * 1.6;

  out << "Polygon \"P\" [";
  out << " " << -corner << " " << -corner << " " << depth_pos << " ";
  out << " " <<  corner << " " << -corner << " " << depth_pos << " ";
  out << " " <<  corner << " " <<  corner << " " << depth_pos << " ";
  out << " " << -corner << " " <<  corner << " " << depth_pos << " ";
  out << "]" << std::endl;
}

#include "C3t3_rib_exporter_plugin.moc"
