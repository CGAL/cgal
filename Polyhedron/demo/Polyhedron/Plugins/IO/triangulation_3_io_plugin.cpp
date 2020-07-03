#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include "T3_type.h"
#include <iostream>
#include <fstream>
#include "Scene_triangulation_3_item.h"

class Triangulation_3_io_plugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "triangulation_3_io_plugin.json")

public:
  bool isDefaultLoader(const CGAL::Three::Scene_item *item) const override
  {
    //TODO: set it for triangulation_3_item
    //if(qobject_cast<const Scene_lcc_item*>(item))
    //  return true;
    return false;
  }
  QString name() const override{ return "triangulation_3_io_plugin"; }
  //todo: probably the same than for c3t3_io_plugn
  QString nameFilters() const override{ return
        ""; }

  //todo: same
  QString saveNameFilters() const override{
    return ".tr_tst";
  }

  //todo
  bool canLoad(QFileInfo) const override{ return true; }
  QList<CGAL::Three::Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override{

    // Open file
    //std::ifstream ifs(fileinfo.filePath().toUtf8());
    //if(!ifs) {
    //  std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    //  ok = false;
    //  Scene_triangulation_3_item* new_item = new Scene_triangulation_3_item();
    //  return QList<CGAL::Three::Scene_item*>();
    //}

    //temp {
    T3* tr = new Tr();
    CGAL::Random rng;
    typedef T3::Point Point;
    typedef T3::Cell_handle Cell_handle;

    while (tr->number_of_vertices() < 20)
      tr->insert(Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));

    const typename T3::Geom_traits::Plane_3
        plane(Point_3(0, 0, 0), Point_3(0, 1, 0), Point_3(0, 0, 1));

    for (Cell_handle c : tr->finite_cell_handles())
    {
      if (plane.has_on_positive_side(
            CGAL::centroid(Point_3(c->vertex(0)->point()), Point_3(c->vertex(1)->point()),
                           Point_3(c->vertex(2)->point()),Point_3( c->vertex(3)->point()))))
        c->set_subdomain_index(1);
      else
        c->set_subdomain_index(2);
    }

    CGAL_assertion(tr->is_valid(true));
    //}


    QString ext = fileinfo.suffix();
    bool res = true;
    if(ext == "")
    {
    }
    else
    {
    }
    if(!res)
    {
      ok = false;
      return QList<CGAL::Three::Scene_item*>();
    }
    Scene_triangulation_3_item* new_item = new Scene_triangulation_3_item(tr);
    //new_item->setName(fileinfo.fileName());
    new_item->invalidateOpenGLBuffers();
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(new_item);
    ok = true;
    return QList<CGAL::Three::Scene_item*>()<<new_item;
  }


  bool canSave(const CGAL::Three::Scene_item*)override{return true;}
  bool save(QFileInfo, QList<CGAL::Three::Scene_item*>& )override{
    return false;
  }

};


#include "triangulation_3_io_plugin.moc"
