#include "ui_Basic_dock_widget.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include "Messages_interface.h"

//! [dock]
class DockWidget :
    public QDockWidget,
    public Ui::BasicDockWidget
{
public:
  DockWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
   setupUi(this);
  }
};
//! [dock]
//This plugin crates an action in Operations that displays "Hello World" in the 'console' dockwidet.
class BasicPlugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  //decides if the plugin's actions will be displayed or not.
  bool applicable(QAction*) const
  {
    return true;
  }
  //the list of the actions of the plugin.
  QList<QAction*> actions() const
  {
    return _actions;
  }
  //! [init]
  //this acts like a constructor for the plugin. It gets the references to the mainwindow and the scene, and connects the action.
  void init(QMainWindow* mw, CGAL::Three::Scene_interface* sc, Messages_interface* mi)
  {
    //gets the reference to the message interface, to display text in the console widget
    this->messageInterface = mi;
    //get the references
    this->scene = sc;
    this->mw = mainWindow;

    //creates the action
    QAction *actionHelloWorld= new QAction(QString("Hello World"), mw);
    //specifies the subMenu
    actionHelloWorld->setProperty("submenuName", "Basic");
    //links the action
    if(actionHelloWorld) {
      connect(actionHelloWorld, SIGNAL(triggered()),
              this, SLOT(helloWorld()));
      _actions << actionHelloWorld;
    }

    dock_widget = new DockWidget("Print a number", mw);
    dock_widget->setVisible(false); // do not show at the beginning

    addDockWidget(dock_widget);

    connect(dock_widget->pushButton, SIGNAL(clicked(bool)),
            this, SLOT(on_dock_button_clicked()));
  }
  //! [init]

private Q_SLOTS:
//! [action]
  void helloWorld()
  {
    // dock widget should be instancied in init()
    if(dock_widget->isVisible()) { dock_widget->hide(); }
    else                         { dock_widget->show(); }
  }

  void on_dock_button_clicked()
  {
      messageInterface->information(QString("Here is your number :%1").arg(dock_widget->spinBox->value()));

  }
  //! [action]
  //! [closure]
  void closure()
  {
    dock_widget->hide();
  }
  //! [closure]
private:
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  DockWidget* dock_widget;
};

#include "Dock_widget_plugin.moc"



