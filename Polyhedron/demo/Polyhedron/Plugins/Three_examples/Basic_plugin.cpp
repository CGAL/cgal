/*
Change the value of EXAMPLE_COMPLEXITY in the first line to change the behavior :
  - 0 : prints "Hello World!" in the Info and console widgets
  - 1 : pops-up a simple dialog asking to enter an integer , then prints it in the Info and console widgets
  - 2 : pops-up a little more elaborated dialog asking to enter an integer , then prints it in the Info and console widgets if it was indeed an integer, else pops-up an error message box.
  */
#define EXAMPLE_COMPLEXITY 0
#include "ui_Basic_dialog.h"
//! [headers_plugin]
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QInputDialog>
#include <QMessageBox>
#include "CGAL/Three/Three.h"
//! [headers_plugin]
//! [dialog_plugin]
class ComplexDialog :
    public QDialog,
    public Ui::BasicDialog
{
  Q_OBJECT
public:
  ComplexDialog(QWidget* =0)
  {
    setupUi(this);

  }

};
//! [dialog_plugin]
//! [opening_plugin]

//This plugin creates an action in Operations depending on EXAMPLE_COMPLEXITY.
class BasicPlugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  //decides if the plugin's actions will be displayed or not.
  bool applicable(QAction*) const Q_DECL_OVERRIDE
  {
    return true;
  }
  //the list of the actions of the plugin.
  QList<QAction*> actions() const Q_DECL_OVERRIDE
  {
    return _actions;
  }
  //this acts like a constructor for the plugin. It gets the references to the main window and the scene, and connects the action.
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* sc, Messages_interface*) Q_DECL_OVERRIDE
  {
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
  }
private Q_SLOTS:

  //! [opening_plugin]
#if EXAMPLE_COMPLEXITY == 0
  //! [basic_plugin]
  void helloWorld()
  {
    CGAL::Three::Three::information(QString("Hello World!"));
  }
  //! [basic_plugin]
#elif EXAMPLE_COMPLEXITY == 1
  //! [basic_dialog_plugin]
  void helloWorld()
  {
    bool ok = false;
    const unsigned int parameter =
        QInputDialog::getInt((QWidget*)mw,
                             tr("Hello World"), // dialog title
                             tr("Hello dear user! What integer would you want me to display for you ? "), // field label
                             10, // default value = fast
                             0, // min
                             100, // max
                             1, // step
                             &ok);
    if(!ok) return;
    messageInterface->information(QString("You asked me to display %1, so here it is : %1").arg(parameter));
  }
  //! [basic_dialog_plugin]
#elif EXAMPLE_COMPLEXITY == 2
  //! [complex_dialog_plugin]
  void helloWorld()
  {
    //creates a new dialog
    ComplexDialog *dialog = new ComplexDialog();
    //opens the dialog
    if(!dialog->exec())
      return;
    //! [warningbox]
    QString result = dialog->lineEdit->text();
    bool ok = false;
    int int_res = result.toInt(&ok);
    if(!ok)
    {
      QMessageBox::warning(mw,
                           "ERROR",
                           tr("This is not an integer !")
                           );
      return;
    }
    //! [warningbox]
    messageInterface->information(QString("You asked me to display %1, so here it is : %1").arg(int_res));
  }
  //! [complex_dialog_plugin]
#endif
  //! [ending_plugin]
private:
  QList<QAction*> _actions;
  //The reference to the scene
  CGAL::Three::Scene_interface* scene;
  //The reference to the main window
  QMainWindow* mw;
};

#include "Basic_plugin.moc"
//! [ending_plugin]

