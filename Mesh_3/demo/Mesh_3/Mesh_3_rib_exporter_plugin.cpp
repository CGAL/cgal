#include "Viewer.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_plugin_helper.h"

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


class Mesh_3_rib_exporter_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);

public:
  Mesh_3_rib_exporter_plugin();
  virtual ~Mesh_3_rib_exporter_plugin() {}
  
  virtual void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionCreateRib;
  }
  
public slots:
  void create_rib();

private:
  typedef Kernel::Point_3   Point_3;
  typedef Kernel::Vector_3  Vector_3;
  typedef Kernel::Plane_3   Plane;
  typedef Kernel::FT        FT;
  typedef Kernel::Aff_transformation_3 Aff_transformation_3;
  
  typedef qglviewer::Vec qglVec;
  
private:
  QStringList nameFilters() const;
  bool save(const Scene_item*, QFileInfo fileinfo);  
  void init_maps(const C3t3& c3t3, const QColor& color);
  
  Point_3 camera_coordinates(const Point_3& p) const;
  
  void write_facets(const C3t3& c3t3, const Plane& plane, std::ofstream& out);
  void write_cells(const C3t3& c3t3, const Plane& plane, std::ofstream& out);
  
  void write_triangle(const Point_3& p, const Point_3& q, const Point_3& r, 
                      const QColor& color, const QColor& edge_color, std::ofstream& out);
  
  void write_point (const Point_3& p, std::ofstream& out);
  void write_point_sphere(const QColor& color, const Point_3& p, std::ofstream& out);
  
  void write_edge_cylinder(const QColor& color, const Point_3& p, const Point_3& q, 
                           std::ofstream& out);
  
  // Writes data which has been stored during triangle drawing
  void write_edges_flat(std::ofstream& out);
  void write_edges_volumic(std::ofstream& out);
  void write_vertices_volumic(std::ofstream& out);
  
  // Utilities
  void write_color(const QColor& color, double use_transparency, std::ofstream& out);
  void write_background(const QColor& color, std::ofstream& out);
  
private:
  QAction* actionCreateRib;
  // Viewer should be in a separate lib to be used here
  QGLViewer* viewer_;
  
  typedef std::map<C3t3::Surface_index, QColor> Surface_map;
  typedef std::map<C3t3::Subdomain_index, QColor> Subdomain_map;
  
  Surface_map surface_map_;
  Subdomain_map subdomain_map_;
  
  typedef std::map<std::pair<Point_3,Point_3>,QColor> Edge_map;
  typedef std::map<Point_3,QColor> Vertex_map;
  
  Edge_map edges_;
  Vertex_map vertices_;
  
  double zmax_;
};


Mesh_3_rib_exporter_plugin::
Mesh_3_rib_exporter_plugin()
  : actionCreateRib(NULL)
  , viewer_(NULL)
  , zmax_(0)
{
  
}


void
Mesh_3_rib_exporter_plugin::
init(QMainWindow* mainWindow, Scene_interface* scene_interface)
{
  this->scene = scene_interface;
  this->mw = mainWindow;
  
  actionCreateRib = new QAction("Export C3t3 to RIB", mw);
  if( NULL != actionCreateRib )
  {
    connect(actionCreateRib, SIGNAL(triggered()), this, SLOT(create_rib()));
  }
  
  viewer_ = mw->findChild<QGLViewer*>("viewer");
  if ( NULL == viewer_ )
  {
    std::cerr << "Can't get QGLViewer" << std::endl;
  }
}


void
Mesh_3_rib_exporter_plugin::create_rib()
{
  if ( NULL == viewer_ )
  {
    std::cerr << "Can't find viewer" << std::endl;
    return;
  }
  
  Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  QStringList filters;
  filters << nameFilters();
  filters << tr("All files (*)");
  
  QString filename = QFileDialog::getSaveFileName(mw,
                                                  tr("Save to File..."),
                                                  QString(),
                                                  filters.join(";;"));
  
  QFileInfo fileinfo(filename);
  
  save(scene->item(index),fileinfo);
}
  

QStringList
Mesh_3_rib_exporter_plugin::nameFilters() const
{ 
  return QStringList() << "RenderMan file (*.rib)";
}


bool
Mesh_3_rib_exporter_plugin::save(const Scene_item* item, QFileInfo fileInfo)
{
  const Scene_c3t3_item* c3t3_item = qobject_cast<const Scene_c3t3_item*>(item);
  if ( NULL == c3t3_item )
  {
    return false;
  }
  
  init_maps(c3t3_item->c3t3(), c3t3_item->color());
  
  QString path = fileInfo.absoluteFilePath();
  std::ofstream obj_file (qPrintable(path));
  
  QString basename = fileInfo.baseName();
  
  obj_file << "Option \"limits\" \"numthreads\" [16]" << std::endl
           << "Option \"searchpath\" \"shader\" \".:./shaders:%PIXIE_SHADERS%:%PIXIE_HOME%/shaders\"" << std::endl
           << "Attribute \"visibility\" \"specular\" 1" << std::endl
           << "Attribute \"visibility\" \"transmission\" 1" << std::endl << std::endl;
  
  obj_file << "Display \""<< qPrintable(basename) << ".tif\" \"file\" \"rgb\"" << std::endl
           << "Format 600 600 1" << std::endl
           << "Projection \"perspective\" \"fov\" 48" << std::endl
           << "PixelSamples 4 4" << std::endl
           << "PixelFilter \"catmull-rom\" 3 3" << std::endl
           << "Rotate 180 0 0 1" << std::endl
           << "WorldBegin" << std::endl
           << "LightSource \"shadowdistant\" 1 \"from\" [0 0 0] \"to\" [0 0 1] \"shadowname\" \"raytrace\" \"intensity\" 0.7" << std::endl
           << "LightSource \"shadowdistant\" 2 \"from\" [1 -1 0] \"to\" [0 1 0] \"shadowname\" \"raytrace\" \"intensity\" 0.6" << std::endl
           << "LightSource \"ambientlight\" 3 \"intensity\" 0.1" << std::endl
           << "LightSource \"ambientlight\" 4 \"intensity\" 1" << std::endl
           << "Illuminate 4 0" << std::endl
           << "ShadingInterpolation \"smooth\"" << std::endl
           << "Surface \"plastic\" \"Ka\" 0.8 \"Kd\" 0.65 \"Ks\" 0.35 \"roughness\" 0.1" << std::endl;
           //<< "Surface \"matte\" \"Ka\" 0.8 \"Kd\" 0.65" << std::endl;
  
  write_facets(c3t3_item->c3t3(), c3t3_item->plane(), obj_file);
  write_cells(c3t3_item->c3t3(), c3t3_item->plane(), obj_file);
  write_edges_volumic(obj_file);
  write_vertices_volumic(obj_file);
  
  write_background(QColor(255,255,255), obj_file);
  
  obj_file << "WorldEnd" << std::endl;
  
  return true;
}

void
Mesh_3_rib_exporter_plugin::init_maps(const C3t3& c3t3, const QColor& color)
{
  surface_map_.clear();
  subdomain_map_.clear();
  edges_.clear();
  vertices_.clear();
  
  // Fill maps with 0 as value
  for ( C3t3::Facet_iterator fit = c3t3.facets_begin(), fend = c3t3.facets_end();
       fit != fend ; ++fit )
  {
    surface_map_.insert(std::make_pair(c3t3.surface_index(*fit),QColor(0,0,0)));
  }
  
  for ( C3t3::Cell_iterator cit = c3t3.cells_begin(), cend = c3t3.cells_end();
       cit != cend ; ++cit )
  {
    subdomain_map_.insert(std::make_pair(c3t3.subdomain_index(cit),QColor(0,0,0)));
  }
  
  // Fill value of maps
  int nb_colors = subdomain_map_.size(); // + surface_map_.size();
  
  // Starting hue
  double c = color.hueF();
  int i = 0;
//  for ( Surface_map::iterator sit = surface_map_.begin(), send = surface_map_.end();
//       sit != send ; ++sit, ++i )
//  {
//    double hue = c + 1./nb_colors * i;
//    if ( hue > 1 ) { hue -= 1.; }
//    sit->second = QColor::fromHsvF(hue, 1., 0.8);
//  }

  for ( Subdomain_map::iterator it = subdomain_map_.begin(), end = subdomain_map_.end();
       it != end ; ++it, ++i )
  {
    double hue = c + 1./nb_colors * i;
    if ( hue > 1 ) { hue -= 1.; }
    it->second = QColor::fromHsvF(hue, color.saturationF(), color.valueF());
  }
}


Mesh_3_rib_exporter_plugin::Point_3 
Mesh_3_rib_exporter_plugin::
camera_coordinates(const Point_3& p) const
{
  qglVec p_vec ( p.x(), p.y(), p.z() );
  qglVec p_cam = viewer_->camera()->cameraCoordinatesOf(p_vec);
  
  return Point_3(p_cam[0],p_cam[1],p_cam[2]);
}


void
Mesh_3_rib_exporter_plugin::
write_facets(const C3t3& c3t3, const Plane& plane, std::ofstream& out)
{
  typedef Kernel::Oriented_side Side;
  
  for ( C3t3::Facet_iterator it = c3t3.facets_begin(), end = c3t3.facets_end();
       it != end ; ++it )
  {
    const C3t3::Cell_handle& c = it->first;
    const int& k = it->second;
   
    const Point_3& p1 = c->vertex((k+1)&3)->point();
    const Point_3& p2 = c->vertex((k+2)&3)->point();
    const Point_3& p3 = c->vertex((k+3)&3)->point();

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
Mesh_3_rib_exporter_plugin::
write_cells(const C3t3& c3t3, const Plane& plane, std::ofstream& out)
{
  typedef Kernel::Oriented_side Side;
  
  for ( C3t3::Cell_iterator it = c3t3.cells_begin(), end = c3t3.cells_end();
       it != end ; ++it )
  {
    const Point_3& p1 = it->vertex(0)->point();
    const Point_3& p2 = it->vertex(1)->point();
    const Point_3& p3 = it->vertex(2)->point();
    const Point_3& p4 = it->vertex(3)->point();
    
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
Mesh_3_rib_exporter_plugin::
write_triangle (const Point_3& p, const Point_3& q, const Point_3& r,
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
  edges_.insert(std::make_pair(std::make_pair(p,q),edge_color));
  edges_.insert(std::make_pair(std::make_pair(p,r),edge_color));
  edges_.insert(std::make_pair(std::make_pair(q,r),edge_color));
  
  vertices_.insert(std::make_pair(p,edge_color));
  vertices_.insert(std::make_pair(q,edge_color));
}


void
Mesh_3_rib_exporter_plugin::
write_point (const Point_3& p, std::ofstream& out)
{
  // Transform point in camera coordinates
  Point_3 p_cam = camera_coordinates(p);
  
  // Write it
  out << " " << -p_cam.x() << " " << -p_cam.y() << " " << -p_cam.z() << " ";
  
  // Store maximal depth
  zmax_ = (std::max)(zmax_, double(-p_cam.z()));
}


void
Mesh_3_rib_exporter_plugin::
write_point_sphere(const QColor& color, const Point_3& p, std::ofstream& out)
{
  // Color
  write_color(color, false, out);
  
  // Transform point in camera coordinates
  Point_3 p_cam = camera_coordinates(p);
  
  // radius
  const double r = 0.6;
  
  out << "Translate " << -p_cam.x() << " " << -p_cam.y() << " " << -p_cam.z() << std::endl;
  
  // Sphere radius zmin zmax thetamax
  out << "Sphere " << r << " " << -r << " " << r << " 360" << std::endl;
  out << "Identity" << std::endl;
}


void
Mesh_3_rib_exporter_plugin::
write_edge_cylinder(const QColor& color, const Point_3& p, const Point_3& q, 
                    std::ofstream& out)
{
  // Color
  write_color(color, false, out);
  
  // Transform point in camera coordinates
  Point_3 p_cam = camera_coordinates(p);
  Point_3 q_cam = camera_coordinates(q);
  
  double pq = CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(p_cam,q_cam)));
  
  Aff_transformation_3 t (CGAL::Translation(), Vector_3(p_cam,CGAL::ORIGIN));
  Point_3 q_cam_t = q_cam.transform(t);
  
  Vector_3 Oq (CGAL::ORIGIN,q_cam_t);
  Vector_3 Oz (FT(0),FT(0),FT(1));
  
  Vector_3 r_axis = CGAL::cross_product(Oq,Oz);
  double cos_angle = CGAL::to_double((Oq*Oz)/CGAL::sqrt(Oq.squared_length()));
  double angle = std::acos(cos_angle) * 180. / CGAL_PI;
  
  // radius
  const double r = 0.3;
  out << "Translate " << -p_cam.x() << " " << -p_cam.y() << " " << -p_cam.z() << std::endl;
  out << "Rotate " << (angle+180.) << " " << -r_axis.x() << " " << -r_axis.y() << " " << -r_axis.z() << std::endl; 
  
  // Cylinder radius zmin zmax thetamax
  out << "Cylinder " << r << " 0 " << pq << " 360" << std::endl;
  out << "Identity" << std::endl;
}


void 
Mesh_3_rib_exporter_plugin::
write_edges_flat(std::ofstream& out)
{
  // Lights
  out << "Illuminate 1 0" << std::endl;
  out << "Illuminate 2 0" << std::endl;
  out << "Illuminate 3 0" << std::endl;
  out << "Illuminate 4 1" << std::endl;
  out << "Surface \"constant\"" << std::endl;
  out << "Opacity 1 1 1" << std::endl;
  
  // Translation
  out << "Translate 0 0 -0.1" << std::endl;
  
  for ( Edge_map::iterator it = edges_.begin(), end = edges_.end() ;
       it != end ; ++it )
  {
    // Color
    write_color(it->second, false, out);
    
    // Edge
    out << "Curves \"linear\" [2] \"nonperiodic\" \"P\" [";
    write_point(it->first.first,out);
    write_point(it->first.second,out);
    out << "] \"constantwidth\" [0.15]" << std::endl;
  }
}


void 
Mesh_3_rib_exporter_plugin::
write_edges_volumic(std::ofstream& out)
{
  // Material
  out << "Opacity 1 1 1" << std::endl;
  
  for ( Edge_map::iterator it = edges_.begin(), end = edges_.end() ;
       it != end ; ++it )
  {
    write_edge_cylinder(it->second, it->first.first, it->first.second, out);
    
    write_point_sphere(it->second, it->first.first, out);
    write_point_sphere(it->second, it->first.second, out);
  }
}

void 
Mesh_3_rib_exporter_plugin::
write_vertices_volumic(std::ofstream& out)
{
  // Material
  out << "Opacity 1 1 1" << std::endl;
  
  for ( Vertex_map::iterator it = vertices_.begin(), end = vertices_.end() ;
       it != end ; ++it )
  {
    write_point_sphere(it->second, it->first, out);
  }
}


void 
Mesh_3_rib_exporter_plugin::
write_color(const QColor& color, double use_transparency, std::ofstream& out)
{
  if (use_transparency)
  {
    double alpha = color.alphaF();
    out << "Opacity " << alpha << " " << alpha << " " << alpha << std::endl;
  }
  
  out << "Color [ " << color.redF() << " " << color.greenF() << " " 
      << color.blueF() <<  " ]" << std::endl;
}


void 
Mesh_3_rib_exporter_plugin::
write_background(const QColor& color, std::ofstream& out)
{
  out << "Illuminate 1 0" << std::endl;
  out << "Illuminate 2 0" << std::endl;
  out << "Illuminate 3 0" << std::endl;
  out << "Illuminate 4 1" << std::endl;
  
  out << "Surface \"constant\"" << std::endl;

//  out << "Surface \"matte\" \"Ka\" 0.9 \"Kd\" 0.65" << std::endl;
  write_color(color,false, out);
  
  double corner = zmax_ * 2.;
  double depth_pos = zmax_ * 1.6;
  
  out << "Polygon \"P\" [";
  out << " " << -corner << " " << -corner << " " << depth_pos << " ";
  out << " " <<  corner << " " << -corner << " " << depth_pos << " ";
  out << " " <<  corner << " " <<  corner << " " << depth_pos << " ";
  out << " " << -corner << " " <<  corner << " " << depth_pos << " ";
  out << "]" << std::endl;
}


#include <QtPlugin>
Q_EXPORT_PLUGIN2(Mesh_3_rib_exporter_plugin, Mesh_3_rib_exporter_plugin);
#include "Mesh_3_rib_exporter_plugin.moc"
