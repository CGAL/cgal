#include "Scene_polyhedron_selection_item.h"
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <fstream>
using namespace CGAL::Three;
class Polyhedron_demo_selection_io_plugin :
        public QObject,
        public Polyhedron_demo_io_plugin_interface
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
    QString name() const { return "selection_io_plugin"; }
    QString nameFilters() const { return "Selection files (*.selection.txt)"; }

    bool canLoad() const { return true; }
    CGAL::Three::Scene_item* load(QFileInfo fileinfo) {
        if(fileinfo.suffix().toLower() != "txt") return 0;
        // There will be no actual loading at this step.
        // Polyhedron_demo_selection_plugin will trigger load when item in new_item_created
        Scene_polyhedron_selection_item* item = new Scene_polyhedron_selection_item();
        if(!item->load(fileinfo.filePath().toStdString())) {
            delete item;
            return NULL;
        }
        return item;
    }

    bool canSave(const CGAL::Three::Scene_item* scene_item) {
        return qobject_cast<const Scene_polyhedron_selection_item*>(scene_item);
    }
    bool save(const CGAL::Three::Scene_item* scene_item, QFileInfo fileinfo) {
        const Scene_polyhedron_selection_item* item = qobject_cast<const Scene_polyhedron_selection_item*>(scene_item);
        if(item == NULL) { return false; }

        return item->save(fileinfo.filePath().toStdString());
    }
};

#include <QtPlugin>
//Q_EXPORT_PLUGIN2(Polyhedron_demo_selection_io_plugin, Polyhedron_demo_selection_io_plugin)
#include "Selection_io_plugin.moc"
