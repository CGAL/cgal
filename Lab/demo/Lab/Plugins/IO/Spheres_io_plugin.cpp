#include "Scene_spheres_item.h"

#include <QMainWindow>
#include <CGAL/Three/CGAL_Lab_io_plugin_interface.h>
#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Three/CGAL_Lab_plugin_helper.h>

#include <CGAL/Three/Three.h>
#include <fstream>
#include <QVariant>
#include <QMessageBox>
#include <QInputDialog>

using namespace CGAL::Three;

class CGAL_Lab_spheres_io_plugin :
  public QObject,
  public CGAL_Lab_io_plugin_interface,
  public CGAL_Lab_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface CGAL::Three::CGAL_Lab_io_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0" FILE "spheres_io_plugin.json")
    Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.IOPluginInterface/1.90")

public:
    // To silent a warning -Woverloaded-virtual
    // See https://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings

    using CGAL_Lab_io_plugin_interface::init;
    //! Configures the widget
    void init(QMainWindow* mainWindow,
              CGAL::Three::Scene_interface* scene_interface,
              Messages_interface*) override{
      //get the references
      this->scene = scene_interface;
      this->mw = mainWindow;

    }
  QString name() const override{ return "spheres_io_plugin"; }
  QString nameFilters() const override{ return "Spheres files (*.spheres.txt)"; }
  bool canLoad(QFileInfo fileinfo) const override;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override;

  bool canSave(const CGAL::Three::Scene_item*) override;
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>&) override;
  bool applicable(QAction* ) const override {
      return false;
  }
  QList<QAction*> actions() const override{

    return QList<QAction*>();
  }

  bool isDefaultLoader(const Scene_item* item) const override{
    if(qobject_cast<const Scene_spheres_item*>(item))
      return true;
    return false;
  }
};

bool CGAL_Lab_spheres_io_plugin::canLoad(QFileInfo fileinfo) const{
  if(!fileinfo.suffix().contains("cgal"))
    return true;
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    return false;
  }
  int first;
  if(!(in >> first)
     || first <= 0)
    return false;
  return true;
}


QList<Scene_item*>
CGAL_Lab_spheres_io_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene){

  // Open file
  std::ifstream ifs(fileinfo.filePath().toUtf8());
  if(!ifs) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    ok = false;
    return QList<Scene_item*>();
  }

  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    Scene_spheres_item* item = new Scene_spheres_item(nullptr, 0, false, false);
    item->setName(fileinfo.completeBaseName());
    ok = true;
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(item);
    return QList<Scene_item*>()<<item;
  }

  double x, y, z, r;
  std::vector<std::array<double, 4> > spheres;
  while (ifs >> x && ifs >> y && ifs >> z && ifs >> r)
    spheres.push_back({x, y, z, r});

  std::vector<QColor> colors;
  if (true) {
    colors.resize(spheres.size());
    for (QColor &c : colors)
      c = generate_random_color();
  }
  else {
    colors.reserve(spheres.size());
    compute_deterministic_color_map(QColor(180, 120, 130, 255), spheres.size(), std::back_inserter(colors));
  }

  Scene_spheres_item* item = new Scene_spheres_item(nullptr, spheres.size(), false, true);
  item->setName(fileinfo.completeBaseName());

  for (std::size_t i = 0;i<spheres.size();i++) {
    const std::array<double, 4>& s = spheres[i];
    item->add_sphere(CGAL::Epick::Sphere_3(CGAL::Epick::Point_3(s[0], s[1], s[2]), s[3] * s[3]), i, CGAL::IO::Color(colors[i].red(), colors[i].green(), colors[i].blue()));
  }

  std::cerr << "Number of spheres loaded: " << spheres.size() << std::endl;

  item->invalidateOpenGLBuffers();

  item->setRenderingMode(Gouraud);
  CGAL::Three::Three::scene()->addItem(item);
  item->computeElements();

  return QList<Scene_item*>()<<item;
}

bool CGAL_Lab_spheres_io_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_spheres_item*>(item) != 0;
}

bool CGAL_Lab_spheres_io_plugin::
save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  const Scene_spheres_item* sphere_item =
    qobject_cast<const Scene_spheres_item*>(item);

  if(!sphere_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());

  out.precision (std::numeric_limits<double>::digits10 + 2);

  if(!out) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return false;
  }

  return false;
}

#include "Spheres_io_plugin.moc"
