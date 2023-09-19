#include <QElapsedTimer>
#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include "Scene_surface_mesh_item.h"
#include "Scene_image_item.h"
#include "SMesh_type.h"
#include "Image_type.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

using namespace CGAL::Three;
class Polyhedron_demo_isosurface_3_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
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
}; // end Polyhedron_demo_convex_hull_plugin

class Polyhedron_demo_isosurface_3_plugin_helper {
  typedef float Word_type;
  typedef std::size_t vertex_descriptor;
  typedef std::pair<Word_type, Word_type> surface_descriptor;

  struct Mesh_quad {
    typedef vertex_descriptor* iterator;
    typedef const vertex_descriptor* const_iterator;

    vertex_descriptor t[4];
    surface_descriptor surface;

    Mesh_quad(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, vertex_descriptor v3, Word_type val1, Word_type val2) {
        t[0] = v0;
        t[1] = v1;
        t[2] = v2;
        t[3] = v3;
        if (val1 < val2) {
            surface = std::make_pair(val1, val2);
        }
        else {
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

    std::vector<EPICK::Point_3> vertices;
    std::vector<Polygon> polygons;

    vertex_descriptor add_vertex(EPICK::Point_3 p) {
        vertices.push_back(p);
        return vertices.size()-1;
    }

    void add_face(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, vertex_descriptor v3, Word_type val1, Word_type val2)
    {
        polygons.push_back(Polygon(v0, v1, v2, v3, val1, val2));
    }
  };

public:
  Polyhedron_demo_isosurface_3_plugin_helper(const Image* image)
  {
    image_ = image;
  }

  void dual_contouring_label_image()
  {
    int i,j,k;

    float dx_ = 1.0;
    float dy_ = 1.0;
    float dz_ = 1.0;
    double vx_ = image_->vx();
    double vy_ = image_->vy();
    double vz_ = image_->vz();
    double tx_ = image_->tx();
    double ty_ = image_->ty();
    double tz_ = image_->tz();

    int xdim_ = image_->xdim();
    int ydim_ = image_->ydim();
    int zdim_ = image_->zdim();
    int xdim_p = xdim_+dx_;
    int ydim_p = ydim_+dy_;
    int zdim_p = zdim_+dz_;

    std::map<std::size_t, vertex_descriptor> vertices;

    for ( k = -dz_ ; k < zdim_ ; k+=dz_ )
    {
        for ( j = -dy_ ; j < ydim_ ; j+=dy_ )
        {
            for ( i = -dx_ ; i < xdim_ ; i+=dx_ )
            {
                //treat_vertex(i,j,k);
                Word_type v0 = image_data(i, j, k);

                Word_type v1 = image_data(i+dx_, j, k);
                Word_type v2 = image_data(i, j+dy_, k);
                Word_type v3 = image_data(i, j, k+dz_);

                Word_type v4 = image_data(i+dx_, j+dy_, k);
                Word_type v5 = image_data(i, j+dy_, k+dz_);
                Word_type v6 = image_data(i+dx_, j, k+dz_);

                Word_type v7 = image_data(i+dx_, j+dy_, k+dz_);

                // if one is different
                if (v0 != v1 || v0 != v2 || v0 != v3 || v0 != v4 || v0 != v5 || v0 != v6 || v0 != v7) {
                    int ip = i+dx_;
                    int jp = j+dy_;
                    int kp = k+dz_;
                    int index = kp*ydim_p*xdim_p + jp*xdim_p + ip;
                    double di = double(i+0.5*dx_),dj = double(j+0.5*dy_),dk = double(k+0.5*dz_);
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
                        vertex_descriptor vertex_1 = vertices.find(kp*ydim_p*xdim_p + j *xdim_p + ip)->second;
                        vertex_descriptor vertex_2 = vertices.find(k *ydim_p*xdim_p + j *xdim_p + ip)->second;
                        vertex_descriptor vertex_3 = vertices.find(k *ydim_p*xdim_p + jp*xdim_p + ip)->second;
                        mesh_.add_face(vertex_0, vertex_1, vertex_2, vertex_3, v0, v1);
                    }
                    // check y direction
                    if (v0 != v2)
                    {
                        vertex_descriptor vertex_1 = vertices.find(kp*ydim_p*xdim_p + jp*xdim_p + i )->second;
                        vertex_descriptor vertex_2 = vertices.find(k *ydim_p*xdim_p + jp*xdim_p + i )->second;
                        vertex_descriptor vertex_3 = vertices.find(k *ydim_p*xdim_p + jp*xdim_p + ip)->second;
                        mesh_.add_face(vertex_0, vertex_1, vertex_2, vertex_3, v0, v2);
                    }
                    // check z direction
                    if (v0 != v3)
                    {
                        vertex_descriptor vertex_1 = vertices.find(kp*ydim_p*xdim_p + jp*xdim_p + i )->second;
                        vertex_descriptor vertex_2 = vertices.find(kp*ydim_p*xdim_p + j *xdim_p + i )->second;
                        vertex_descriptor vertex_3 = vertices.find(kp*ydim_p*xdim_p + j *xdim_p + ip)->second;
                        mesh_.add_face(vertex_0, vertex_1, vertex_2, vertex_3, v0, v3);
                    }
                }
            }
        }
    }
  }

  void convert_to_smesh(SMesh& out, QColor c = QColor(255, 0, 0))
  {
    // compute colors
    typedef std::map<surface_descriptor, QColor> Surface_map_type;
    Surface_map_type surface_map;
    Surface_map_type::iterator surface_map_it;
    for (const Mesh_quad& quad : mesh_.polygons) {
        surface_map_it = surface_map.find(quad.surface);
        if (surface_map_it == surface_map.end()) {
            surface_map[quad.surface] = QColor();
        }
    }
    const double nb_domains = surface_map.size();
    double i = 0;
    for (Surface_map_type::iterator it = surface_map.begin(),
         end = surface_map.end(); it != end; ++it, i += 1.)
    {
        double hue = c.hueF() + 1. / nb_domains * i;
        if (hue > 1) { hue -= 1.; }
        surface_map[it->first] = QColor::fromHsvF(hue, c.saturationF(), c.valueF());
    }

    // orient polygons
    CGAL::Polygon_mesh_processing::orient_polygon_soup(mesh_.vertices, mesh_.polygons);

    // make polygon mesh
    typedef std::pair<int, SMesh::Face_index> Index_link;
    std::vector<Index_link> meshface_to_smeshface;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(mesh_.vertices, mesh_.polygons, out,
            CGAL::parameters::polygon_to_face_output_iterator(std::back_inserter(meshface_to_smeshface))
        );

    // set facet color property
    auto color_property = out.add_property_map<SMesh::Face_index>("f:color", CGAL::IO::white()).first;
    for (const Index_link & f : meshface_to_smeshface) {
        const int & mesh_index = f.first;
        const SMesh::Face_index & smesh_index = f.second;
        const Mesh_quad & quad = mesh_.polygons[mesh_index];

        const QColor & color = surface_map[quad.surface];

        put(color_property, smesh_index, CGAL::IO::Color(color.red(), color.green(), color.blue()));
    }
  }
private:
  Mesh mesh_;
  const Image* image_;
  int nb_surfaces;

  Word_type image_data(std::size_t i, std::size_t j, std::size_t k)
  {
    if ( i>=0 && i<image_->xdim() && j>=0 && j<image_->ydim() && k>=0 && k<image_->zdim() )
        return image_->value(i, j, k);
    else
        return 0;
  }
};

void Polyhedron_demo_isosurface_3_plugin::on_actionDualLabel_triggered()
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
    std::cout << "Dual contour label image...";

    SMesh mesh;
    Polyhedron_demo_isosurface_3_plugin_helper helper(img_item->image());
    helper.dual_contouring_label_image();
    helper.convert_to_smesh(mesh, img_item->color());

    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

    Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(mesh);
    new_item->setName(tr("%1 (surfaces)").arg(scene->item(index)->name()));
    new_item->setColor(img_item->color());
    new_item->setRenderingMode(FlatPlusEdges);
    scene->addItem(new_item);

    img_item->setVisible(false);

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

#include "Isosurface_3_plugin.moc"
