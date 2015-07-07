#include <QApplication>
#include "Scene_polyhedron_item.h"
#include "Scene_plane_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/centroid.h>
#include <CGAL/bounding_box.h>
#include <CGAL/linear_least_squares_fitting_3.h>


#include <CGAL/Make_quad_soup.h> // output for plane fitting
#include <CGAL/Make_bar.h> // output for line fitting

#include "Kernel_type.h"
typedef Kernel::Plane_3 Plane;
typedef Kernel::Iso_cuboid_3 Iso_cuboid;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Line_3 Line;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;


class Polyhedron_demo_pca_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionFitPlane"
                         << "actionFitLine";
  }

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }


public Q_SLOTS:
  void on_actionFitPlane_triggered();
  void on_actionFitLine_triggered();

}; // end Polyhedron_demo_pca_plugin

void Polyhedron_demo_pca_plugin::on_actionFitPlane_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    // get triangles from the mesh
    std::list<Triangle> triangles;
    Polyhedron::Facet_iterator f;
    for(f = pMesh->facets_begin();
      f != pMesh->facets_end();
      ++f)
    {
      const Point& a = f->halfedge()->vertex()->point();
      const Point& b = f->halfedge()->next()->vertex()->point();
      const Point& c = f->halfedge()->prev()->vertex()->point();
      triangles.push_back(Triangle(a,b,c));
    }

    // fit plane to triangles
    Plane plane;
    std::cout << "Fit plane...";
    CGAL::linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());
    std::cout << "ok" << std::endl;

    // compute centroid
    Point center_of_mass = CGAL::centroid(triangles.begin(),triangles.end());

//     // compute bounding box diagonal
//     Iso_cuboid bbox = CGAL::bounding_box(pMesh->points_begin(),pMesh->points_end());

//     // compute scale (for rendering) using diagonal of bbox
//     Point cmin = (bbox.min)();
//     Point cmax = (bbox.max)();
//     FT diag = std::sqrt(CGAL::squared_distance(cmin,cmax));
//     Vector u1 = plane.base1();
//     u1 = u1 / std::sqrt(u1*u1);
//     u1 = u1 * 0.7 * diag;
//     Vector u2 = plane.base2();
//     u2 = u2 / std::sqrt(u2*u2);
//     u2 = u2 * 0.7 * diag;
//     std::list<Point> points;
//     points.push_back(center_of_mass + u1);
//     points.push_back(center_of_mass + u2);
//     points.push_back(center_of_mass - u1);
//     points.push_back(center_of_mass - u2);

    // add best fit plane as new polyhedron
//     Polyhedron *pFit = new Polyhedron;
//     typedef std::list<Point>::iterator Iterator;
//     Make_quad_soup<Polyhedron,Kernel,Iterator> quad;
//     quad.run(points.begin(),points.end(),*pFit);

//     Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pFit);
    Scene_plane_item* new_item = new Scene_plane_item(this->scene);
    new_item->setPosition(center_of_mass.x(),
                          center_of_mass.y(),
                          center_of_mass.z());
    const Vector& normal = plane.orthogonal_vector();
    new_item->setNormal(normal.x(), normal.y(), normal.z());
    new_item->setName(tr("%1 (plane fit)").arg(item->name()));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode(item->renderingMode());
    scene->addItem(new_item);

    QApplication::restoreOverrideCursor();
  }
}

void Polyhedron_demo_pca_plugin::on_actionFitLine_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    // get triangles from the mesh
    std::list<Triangle> triangles;
    Polyhedron::Facet_iterator f;
    for(f = pMesh->facets_begin();
      f != pMesh->facets_end();
      ++f)
    {
      const Point& a = f->halfedge()->vertex()->point();
      const Point& b = f->halfedge()->next()->vertex()->point();
      const Point& c = f->halfedge()->prev()->vertex()->point();
      triangles.push_back(Triangle(a,b,c));
    }

    // fit line to triangles
    Line line;
    std::cout << "Fit line...";
    CGAL::linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,CGAL::Dimension_tag<2>());
    std::cout << "ok" << std::endl;

    // compute centroid
    Point center_of_mass = CGAL::centroid(triangles.begin(),triangles.end());

    // compute bounding box diagonal
    Iso_cuboid bbox = CGAL::bounding_box(pMesh->points_begin(),pMesh->points_end());

    // compute scale for rendering using diagonal of bbox
    Point cmin = (bbox.min)();
    Point cmax = (bbox.max)();
    FT diag = std::sqrt(CGAL::squared_distance(cmin,cmax));

    // construct a 3D bar
    Vector u = line.to_vector();
    u = u / std::sqrt(u*u);

    Point a = center_of_mass + u * diag;
    Point b = center_of_mass - u * diag;

    Plane plane_a = line.perpendicular_plane(a);

    Vector u1 = plane_a.base1();
    u1 = u1 / std::sqrt(u1*u1);
    u1 = u1 * 0.01 * diag;
    Vector u2 = plane_a.base2();
    u2 = u2 / std::sqrt(u2*u2);
    u2 = u2 * 0.01 * diag;

    Point points[8];

    points[0] = a + u1;
    points[1] = a + u2;
    points[2] = a - u1;
    points[3] = a - u2;

    points[4] = b + u1;
    points[5] = b + u2;
    points[6] = b - u1;
    points[7] = b - u2;

    // add best fit line as new polyhedron bar
    Polyhedron *pFit = new Polyhedron;
    Make_bar<Polyhedron,Kernel> bar;
    bar.run(points,*pFit);

    Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pFit);
    new_item->setName(tr("%1 (line fit)").arg(item->name()));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode(item->renderingMode());
    scene->addItem(new_item);

    QApplication::restoreOverrideCursor();
  }
}

#include "Polyhedron_demo_pca_plugin.moc"
