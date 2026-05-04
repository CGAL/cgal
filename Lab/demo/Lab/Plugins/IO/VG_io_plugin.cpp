#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/CGAL_Lab_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/IO/VG.h>
#include <fstream>
#include <CGAL/Random.h>
#include "Color_map.h"


using namespace CGAL::Three;
class CGAL_Lab_vg_plugin :
  public QObject,
  public CGAL::Three::CGAL_Lab_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.IOPluginInterface/1.90" FILE "vg_io_plugin.json")

public:
  QString name() const { return "vg_plugin"; }
  QString nameFilters() const { return "Vertex Group files (*.vg);;"; }
  bool canLoad(QFileInfo fileinfo) const;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& );
};

bool CGAL_Lab_vg_plugin::canLoad(QFileInfo ) const {
  return true;
}

QList<Scene_item*> CGAL_Lab_vg_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene) {
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios_base::binary);

  if(!in)
    std::cerr << "Error!\n";

  auto constructor = [](unsigned int& type, const std::string& params) {
    if (params.empty())
      return std::to_string(type);
    else
      return std::to_string(type) + "_" + params;
    };

  std::vector<std::pair<std::string, std::vector<std::size_t>>> regions;
  Point_set point_set;

  std::vector<CGAL::IO::Color> colors;
  std::vector<std::string> labels;

  ok = CGAL::IO::read_VG(in, point_set, std::back_inserter(regions), CGAL::parameters::constructor(constructor).colors(std::ref(colors)).labels(std::ref(labels)));
  if(ok)
  {
    QList<Scene_item*> list;
    CGAL::Random rand(static_cast<unsigned int>(time(nullptr)));
    Scene_group_item* group = new Scene_group_item(fileinfo.baseName());
    group->setRenderingMode(Points);
    CGAL::Three::Three::scene()->addItem(group);

    std::vector<QColor> distinct_colors;
    compute_deterministic_color_map(QColor(80, 250, 80), regions.size(), std::back_inserter(distinct_colors));

    std::size_t idx = 0;
    for (const auto &r : regions) {
      Scene_points_with_normal_item* point_item =
        new Scene_points_with_normal_item;

      for (auto& item : r.second)
        point_item->point_set()->insert(point_set, item);

      if (idx + 1 < colors.size())
        point_item->setRgbColor(colors[idx].red(), colors[idx].green(), colors[idx].blue());
      else
        point_item->setRgbColor(distinct_colors[idx].red(), distinct_colors[idx].green(), distinct_colors[idx].blue());

      if (idx + 1 < labels.size())
        point_item->setName(QString::fromStdString(std::to_string(idx) + "_" + labels[idx]));
      else
        point_item->setName(QString::fromStdString(std::to_string(idx) + "_" + r.first));

      if (point_set.has_normal_map())
        point_item->setRenderingMode(Points);
      else
        point_item->setRenderingMode(CGAL::Three::Three::defaultPointSetRenderingMode());

      point_item->invalidateOpenGLBuffers();
      CGAL::Three::Three::scene()->addItem(point_item);
      CGAL::Three::Three::scene()->changeGroup(point_item, group);

      list << point_item;
      idx++;
    }

    return list;
  }

  return QList<Scene_item*>();
}

bool CGAL_Lab_vg_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return false;//qobject_cast<const Scene_points_with_normal_item*>(item);
}

bool CGAL_Lab_vg_plugin::save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  return false;
/*
  Scene_item* item = items.front();
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "vg" && extension != "VG")
    return false;

  // This plugin supports point sets
  Scene_points_with_normal_item* point_set_item =
    qobject_cast<Scene_points_with_normal_item*>(item);
  if(!point_set_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8().data());
  items.pop_front();
  Q_ASSERT(point_set_item->point_set() != nullptr);

  point_set_item->point_set()->reset_indices();

  return CGAL::IO::write_VG(out, *(point_set_item->point_set()));*/
}


#include "VG_io_plugin.moc"
