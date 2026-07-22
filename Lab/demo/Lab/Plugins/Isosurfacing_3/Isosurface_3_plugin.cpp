#include <QElapsedTimer>
#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include "Scene_polygon_soup_item.h"
#include "Scene_image_item.h"
#include "SMesh_type.h"
#include "Image_type.h"

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

using namespace CGAL::Three;
class CGAL_Lab_isosurface_3_plugin :
  public QObject,
  public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")
public:
  void init(QMainWindow*mw,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
    scene = scene_interface;
    this->mw = mw;
    QAction *actionDualLabel = new QAction("Dual contour label image", mw);
    actionDualLabel->setProperty("subMenuName","Isosurface");
    connect(actionDualLabel, SIGNAL(triggered()), this, SLOT(on_actionDualLabel_triggered()));
    _actions << actionDualLabel;
  }

  QList<QAction*> actions()const {return _actions;}

  bool applicable(QAction*) const {
    return
      qobject_cast<Scene_image_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionDualLabel_triggered();
private:
  QList<QAction*> _actions;
  Scene_interface* scene;
  QMainWindow* mw;
}; // end CGAL_Lab_convex_hull_plugin

class CGAL_Lab_isosurface_3_plugin_helper {
  typedef float Word_type;
  typedef std::size_t vertex_descriptor;
  typedef std::pair<Word_type, Word_type> surface_descriptor;

  struct Mesh_quad {
    typedef vertex_descriptor* iterator;
    typedef const vertex_descriptor* const_iterator;

    vertex_descriptor t[4];
    surface_descriptor surface;

    Mesh_quad() {}

    Mesh_quad(const Mesh_quad & other) {
      t[0] = other.t[0];
      t[1] = other.t[1];
      t[2] = other.t[2];
      t[3] = other.t[3];
      surface = other.surface;
    }

    Mesh_quad(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, vertex_descriptor v3, Word_type val1, Word_type val2) {
      if (val1 < val2) {
        t[0] = v3;
        t[1] = v2;
        t[2] = v1;
        t[3] = v0;
        surface = std::make_pair(val1, val2);
      }
      else {
        t[0] = v0;
        t[1] = v1;
        t[2] = v2;
        t[3] = v3;
        surface = std::make_pair(val2, val1);
      }
    }

    std::size_t size() const {
      return 4;
    }

    vertex_descriptor operator[](int i) const { return t[i]; }

    iterator begin() { return &t[0]; }
    const_iterator begin() const { return &t[0]; }
    iterator end() { return &t[4]; }
    const_iterator end() const { return &t[4]; }
  };

  struct Mesh
  {
    using Polygon = Mesh_quad;

    std::vector<Point_3> vertices;
    std::vector<Polygon> polygons;

    vertex_descriptor add_vertex(EPICK::Point_3 p) {
      vertices.push_back(p);
      return vertices.size()-1;
    }

    void add_face(const Polygon & poly) {
      polygons.push_back(Polygon(poly));
    }

    void add_face(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, vertex_descriptor v3, Word_type val1, Word_type val2)
    {
      polygons.push_back(Polygon(v0, v1, v2, v3, val1, val2));
    }
  };

public:
  CGAL_Lab_isosurface_3_plugin_helper(const Image* image)
  {
    image_ = image;
  }

  void dual_contouring_label_image()
  {
    int i,j,k;

    const int dx_ = 1;
    const int dy_ = 1;
    const int dz_ = 1;
    double vx_ = image_->vx();
    double vy_ = image_->vy();
    double vz_ = image_->vz();
    double tx_ = image_->tx();
    double ty_ = image_->ty();
    double tz_ = image_->tz();

    int xdim_ = static_cast<int>(image_->xdim());
    int ydim_ = static_cast<int>(image_->ydim());
    int zdim_ = static_cast<int>(image_->zdim());
    int xdim_p = xdim_+dx_;
    int ydim_p = ydim_+dy_;
    int xydim_p = ydim_p*xdim_p;

    std::map<std::size_t, vertex_descriptor> vertices;

    for ( k = -dz_ ; k < zdim_ ; k+=dz_ )
    {
      for ( j = -dy_ ; j < ydim_ ; j+=dy_ )
      {
        for ( i = -dx_ ; i < xdim_ ; i+=dx_ )
        {
          const Word_type & v0 = image_data(i, j, k);

          const Word_type & v1 = image_data(i+dx_, j, k);
          const Word_type & v2 = image_data(i, j+dy_, k);
          const Word_type & v3 = image_data(i, j, k+dz_);

          const Word_type & v4 = image_data(i+dx_, j+dy_, k);
          const Word_type & v5 = image_data(i, j+dy_, k+dz_);
          const Word_type & v6 = image_data(i+dx_, j, k+dz_);

          const Word_type & v7 = image_data(i+dx_, j+dy_, k+dz_);

          // if one is different
          if (v0 != v1 || v0 != v2 || v0 != v3 || v0 != v4 || v0 != v5 || v0 != v6 || v0 != v7) {
            unsigned int ip = i+dx_;
            unsigned int jp = j+dy_;
            unsigned int kp = k+dz_;
            unsigned int index = kp*xydim_p + jp*xdim_p + ip;
            double di = i+0.5*dx_;
            double dj = j+0.5*dy_;
            double dk = k+0.5*dz_;
            if (di < 0)
              di = 0.0;
            if (dj < 0)
              dj = 0.0;
            if (dk < 0)
              dk = 0.0;
            if (di > xdim_-dx_)
              di = xdim_-1;
            if (dj > ydim_-dy_)
              dj = ydim_-1;
            if (dk > zdim_-dz_)
              dk = zdim_-1;
            di = di * vx_ + tx_;
            dj = dj * vy_ + ty_;
            dk = dk * vz_ + tz_;

            vertex_descriptor vertex_0 = mesh_.add_vertex(Point_3(di, dj, dk));
            vertices.insert(std::make_pair(index, vertex_0));

            // check x direction
            if (v0 != v1)
            {
              vertex_descriptor vertex_1 = vertices.find(kp*xydim_p + j *xdim_p + ip)->second;
              vertex_descriptor vertex_2 = vertices.find(k *xydim_p + j *xdim_p + ip)->second;
              vertex_descriptor vertex_3 = vertices.find(k *xydim_p + jp*xdim_p + ip)->second;
              mesh_.add_face(vertex_0, vertex_1, vertex_2, vertex_3, v0, v1);
            }
            // check y direction
            if (v0 != v2)
            {
              vertex_descriptor vertex_1 = vertices.find(k *xydim_p + jp*xdim_p + ip)->second;
              vertex_descriptor vertex_2 = vertices.find(k *xydim_p + jp*xdim_p + i )->second;
              vertex_descriptor vertex_3 = vertices.find(kp*xydim_p + jp*xdim_p + i )->second;
              mesh_.add_face(vertex_0, vertex_1, vertex_2, vertex_3, v0, v2);
            }
            // check z direction
            if (v0 != v3)
            {
              vertex_descriptor vertex_1 = vertices.find(kp*xydim_p + jp*xdim_p + i )->second;
              vertex_descriptor vertex_2 = vertices.find(kp*xydim_p + j *xdim_p + i )->second;
              vertex_descriptor vertex_3 = vertices.find(kp*xydim_p + j *xdim_p + ip)->second;
              mesh_.add_face(vertex_0, vertex_1, vertex_2, vertex_3, v0, v3);
            }
          }
        }
      }
    }
  }

  void generate_polygon_soup(std::vector<Point_3> & vertices,
                             std::vector<std::vector<std::size_t>> & polygons,
                             std::vector<CGAL::IO::Color> & fcolors,
                             std::vector<CGAL::IO::Color> & vcolors,
                             const QColor & c)
  {
    vertices.clear();
    polygons.clear();
    fcolors.clear();
    vcolors.clear();

    const std::vector<Point_3> & m_vertices = mesh_.vertices;
    const std::vector<Mesh_quad> & m_polygons = mesh_.polygons;
    std::size_t polygon_size = m_polygons.size();

    vertices.insert(vertices.begin(), m_vertices.begin(), m_vertices.end());

    std::map<surface_descriptor, QColor> surface_color_map;
    compute_color_map(surface_color_map, c);

    fcolors.resize(polygon_size);
    polygons.resize(polygon_size);
    for (std::size_t i = 0; i < polygon_size; i++)
    {
      const Mesh_quad & quad = m_polygons[i];

      std::vector<std::size_t> polygon;
      polygon.resize(4);
      polygon[0] = quad.t[0];
      polygon[1] = quad.t[1];
      polygon[2] = quad.t[2];
      polygon[3] = quad.t[3];
      polygons[i] = polygon;

      const QColor & color = surface_color_map[quad.surface];
      fcolors[i] = CGAL::IO::Color(color.red(), color.green(), color.blue());
    }
  }

private:
  Mesh mesh_;
  const Image* image_;

  Word_type image_data(int i, int j, int k)
  {
    if ( i>=0 && static_cast<std::size_t>(i)<image_->xdim() &&
         j>=0 && static_cast<std::size_t>(j)<image_->ydim() &&
         k>=0 && static_cast<std::size_t>(k)<image_->zdim() )
      return image_->value(i, j, k);
    else
      return 0;
  }

  void compute_color_map(std::map<surface_descriptor, QColor> & surface_color_map, const QColor & c)
  {
    // obtain and order all subdomains and surfaces
    std::map<Word_type, QColor> subdomain_color_map;
    for (const Mesh_quad& quad : mesh_.polygons) {
      const surface_descriptor & surface = quad.surface;
      subdomain_color_map[surface.first] = QColor();
      subdomain_color_map[surface.second] = QColor();
      surface_color_map[surface] = QColor();
    }

    // assign default color to each subdomains (same colors as images)
    double nb_domains = subdomain_color_map.size();
    double i = 0;
    for (std::map<Word_type, QColor>::iterator it = subdomain_color_map.begin(),
         end = subdomain_color_map.end(); it != end; ++it, i += 1.)
    {
      double hue = c.hueF() + 1. / nb_domains * i;
      if (hue > 1) { hue -= 1.; }
      it->second = QColor::fromHsvF(hue, c.saturationF(), c.valueF());
    }

    // assign color for each surfaces (use subdomain color if the surface is in contact with background)
    nb_domains = surface_color_map.size();
    i = 0;
    double patch_hsv_value = fmod(c.valueF() + .5, 1.);
    for (std::map<surface_descriptor, QColor>::iterator it = surface_color_map.begin(),
         end = surface_color_map.end(); it != end; ++it, i += 1.)
    {
      QColor surface_color;
      if (it->first.first == 0) {
        surface_color = subdomain_color_map[it->first.second];
      }
      else {
        double hue = c.hueF() + 1. / nb_domains * i;
        if (hue > 1) { hue -= 1.; }
        surface_color =  QColor::fromHsvF(hue, c.saturationF(), patch_hsv_value);
      }

      it->second = surface_color;
    }
  }
};

void CGAL_Lab_isosurface_3_plugin::on_actionDualLabel_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_image_item* img_item =
    qobject_cast<Scene_image_item*>(scene->item(index));

  if( img_item )
  {
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    QElapsedTimer time;
    time.start();
    QElapsedTimer time_per_op;
    time_per_op.start();
    std::cout << "Dual contour label image...";

    CGAL_Lab_isosurface_3_plugin_helper helper(img_item->image());
    helper.dual_contouring_label_image();

    std::vector<Point_3> vertices;
    std::vector<std::vector<std::size_t>> polygons;
    std::vector<CGAL::IO::Color> fcolors;
    std::vector<CGAL::IO::Color> vcolors;

    helper.generate_polygon_soup(vertices, polygons, fcolors, vcolors, img_item->color());

    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

    Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item();
    new_item->load(vertices, polygons, fcolors, vcolors);
    new_item->setName(tr("%1 (surface)").arg(scene->item(index)->name()));
    new_item->setColor(img_item->color());
    new_item->setRenderingMode(FlatPlusEdges);
    scene->addItem(new_item);


    img_item->setVisible(false);

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

#include "Isosurface_3_plugin.moc"
