#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QInputDialog>
#include "Messages_interface.h"
#include "Scene_plane_item.h"
#include <CGAL/Three/Scene_group_item.h>

//This plugin creates an action in Operations that creates a plane item and adds it to a group.
class GroupItemPlugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  //This plugin is always applicable .
  bool applicable(QAction*) const
  {
    return true;
  }
  //the list of the actions of the plugin.
  QList<QAction*> actions() const
  {
    return _actions;
  }
  //this acts like a constructor for the plugin. It gets the references to the mainwindow and the scene, and connects the action.
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* sc  )
  {
    //gets the references
    this->scene = sc;
    this->mw = mainWindow;

    //creates the action
    QAction *actionMakeGroup= new QAction(QString("Create Group"), mw);
    //specifies the subMenu
    actionMakeGroup->setProperty("submenuName", "Basic");
    //links the action
    if(actionMakeGroup) {
      connect(actionMakeGroup, SIGNAL(triggered()),
              this, SLOT(makeGroup()));
      _actions << actionMakeGroup;
    }
  }

  void init(QMainWindow* mw, CGAL::Three::Scene_interface* sc, Messages_interface* mi)
  {
    //gets the reference to the message interface, to display text in the console widget
    this->messageInterface = mi;
    init(mw, sc);
  }
private Q_SLOTS:


  void makeGroup()
  {
    //creates an item
    Scene_plane_item *child = new Scene_plane_item(scene);
    child->setName("Trivial Plane");
    child->setColor(Qt::blue);
    child->setNormal(0.0,0.0,1.0);
    scene->addItem(child);
    //clears the selection to avoid adding unwanted items to the group.
    scene->setSelectedItem(-1);
    //Creates a new group
    Scene_group_item *group = new Scene_group_item("New group");
    //Then gives it its children
    scene->changeGroup(child, group);
    //adds it to the scene
    scene->addGroup(group);

  }

private:
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  //The reference to the scene
  CGAL::Three::Scene_interface* scene;
  //The reference to the main window
  QMainWindow* mw;
};

#include "Group_item_plugin.moc"
