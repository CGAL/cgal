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
  QString nameFilters() const override{ return "T3 files(*.tr_test)"; }


  //todo
  bool canLoad(QFileInfo) const override{ return true; }
  QList<CGAL::Three::Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override{

    // Open file
    std::ifstream ifs(fileinfo.filePath().toUtf8());
    if(!ifs) {
      std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
      ok = false;
      return QList<CGAL::Three::Scene_item*>();
    }
    T3* tr = new Tr();
    ifs >> *tr;
    if(ifs.fail() || !tr->is_valid(false)) {
      std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
      ok = false;
      return QList<CGAL::Three::Scene_item*>();
    }

    Scene_triangulation_3_item* new_item = new Scene_triangulation_3_item(tr);
    new_item->setName(fileinfo.fileName());
    new_item->invalidateOpenGLBuffers();
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(new_item);
    ok = true;
    return QList<CGAL::Three::Scene_item*>()<<new_item;
  }


  bool canSave(const CGAL::Three::Scene_item* item)override
  {
    return qobject_cast<const Scene_triangulation_3_item*>(item);
  }

  bool save(QFileInfo fileinfo, QList<Scene_item *> &items)override{
    for(int id : CGAL::Three::Three::scene()->selectionIndices())
    {
      Scene_item* item = CGAL::Three::Three::scene()->item(id);
      const Scene_triangulation_3_item* t3_item = qobject_cast<const Scene_triangulation_3_item*>(item);
      if (!t3_item)
      {
        continue;
      }

      QString path = fileinfo.absoluteFilePath();

      if(path.endsWith(".tr_test"))
      {
        std::ofstream out(fileinfo.filePath().toUtf8(),
                          std::ios_base::out|std::ios_base::binary);

        out << t3_item->triangulation();
        if( out.fail())
          return false;
        else
        {
          items.pop_front();
          return true;
        }
      }
    }
    return false;
  }

};


#include "triangulation_3_io_plugin.moc"
