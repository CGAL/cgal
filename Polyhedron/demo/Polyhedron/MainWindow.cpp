#include <cmath>

#include "config.h"
#include "MainWindow.h"
#include "Scene.h"
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/TextRenderer.h>
#include <CGAL/Three/exceptions.h>
#include <CGAL/Qt/debug.h>

#include <QtDebug>
#include <QFileDialog>
#include <QFileInfo>
#include <QSettings>
#include <QHeaderView>
#include <QMenu>
#include <QMenuBar>
#include <QChar>
#include <QAction>
#include <QShortcut>
#include <QKeySequence>
#include <QLibrary>
#include <QPluginLoader>
#include <QMessageBox>
#include <QScrollBar>
#include <QColor>
#include <QColorDialog>
#include <QClipboard>
#include <QCloseEvent>
#include <QInputDialog>
#include <QTreeView>
#include <QSortFilterProxyModel>
#include <QMap>
#include <QSet>
#include <QStandardItemModel>
#include <QStandardItem>
#include <stdexcept>
#ifdef QT_SCRIPT_LIB
#  include <QScriptValue>
#  ifdef QT_SCRIPTTOOLS_LIB
#    include <QScriptEngineDebugger>
#  endif
#endif

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Scene_item_with_properties.h>
#include "ui_MainWindow.h"
#include "ui_Preferences.h"
#include "ui_Statistics_on_item_dialog.h"
#include "Show_point_dialog.h"
#include "File_loader_dialog.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <QGLViewer/manipulatedFrame.h>

#ifdef QT_SCRIPT_LIB
#  include <QScriptEngine>
#  include <QScriptValue>
#include "Color_map.h"
using namespace CGAL::Three;
QScriptValue 
myScene_itemToScriptValue(QScriptEngine *engine, 
                          CGAL::Three::Scene_item* const &in)
{ 
  return engine->newQObject(in); 
}

void myScene_itemFromScriptValue(const QScriptValue &object, 
                                 CGAL::Three::Scene_item* &out)
{
  out = qobject_cast<CGAL::Three::Scene_item*>(object.toQObject());
}
#endif // QT_SCRIPT_LIB

#ifdef QT_SCRIPT_LIB
#  ifdef QT_SCRIPTTOOLS_LIB

const QScriptEngineDebugger::DebuggerWidget debug_widgets[9] = {
  QScriptEngineDebugger::ConsoleWidget,
  QScriptEngineDebugger::StackWidget,
  QScriptEngineDebugger::ScriptsWidget,
  QScriptEngineDebugger::LocalsWidget,
  QScriptEngineDebugger::CodeWidget,
  QScriptEngineDebugger::CodeFinderWidget,
  QScriptEngineDebugger::BreakpointsWidget,
  QScriptEngineDebugger::DebugOutputWidget,
  QScriptEngineDebugger::ErrorLogWidget
};
const QString debug_widgets_names[9] = {
  "Script console",
  "Stack",
  "Scripts",
  "Locals",
  "Code",
  "CodeFinder",
  "Breakpoints",
  "DebugOutput",
  "ErrorLog"
};

#  endif
#endif

QScriptValue myPrintFunction(QScriptContext *context, QScriptEngine *engine)
{
  MainWindow* mw = qobject_cast<MainWindow*>(engine->parent());
  QString result;
  for (int i = 0; i < context->argumentCount(); ++i) {
    if (i > 0)
      result.append(" ");
    result.append(context->argument(i).toString());
  }

  if(mw) mw->message(QString("QtScript: ") + result, "");
  QTextStream (stdout) << (QString("QtScript: ") + result) << "\n";

  return engine->undefinedValue();
}

MainWindow::~MainWindow()
{
  delete ui;
  delete statistics_ui;
}
MainWindow::MainWindow(QWidget* parent)
  : CGAL::Qt::DemosMainWindow(parent)
{
  ui = new Ui::MainWindow;
  ui->setupUi(this);
  menu_map[ui->menuOperations->title()] = ui->menuOperations;
  // remove the Load Script menu entry, when the demo has not been compiled with QT_SCRIPT_LIB
#if !defined(QT_SCRIPT_LIB)
  ui->menuBar->removeAction(ui->actionLoadScript);
#endif
  // Save some pointers from ui, for latter use.
  sceneView = ui->sceneView;
  viewer = ui->viewer;
  // do not save the state of the viewer (anoying)
  viewer->setStateFileName(QString::null);

  // setup scene
  scene = new Scene(this);
  viewer->textRenderer()->setScene(scene);
  viewer->setScene(scene);
  ui->actionMaxTextItemsDisplayed->setText(QString("Set Maximum Text Items Displayed : %1").arg(viewer->textRenderer()->getMax_textItems()));
  {
    QShortcut* shortcut = new QShortcut(QKeySequence(Qt::ALT+Qt::Key_Q), this);
    connect(shortcut, SIGNAL(activated()),
            this, SLOT(setFocusToQuickSearch()));
  }

  proxyModel = new QSortFilterProxyModel(this);
  proxyModel->setSourceModel(scene);
  SceneDelegate *delegate = new SceneDelegate(this);
  delegate->setProxy(proxyModel);
  delegate->setScene(scene);


  connect(ui->searchEdit, SIGNAL(textChanged(QString)),
          proxyModel, SLOT(setFilterFixedString(QString)));
  sceneView->setModel(proxyModel);

  // setup the sceneview: delegation and columns sizing...
  sceneView->setItemDelegate(delegate);
  resetHeader();

  // setup connections
  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          this, SLOT(updateInfo()));
  
  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          this, SLOT(updateDisplayInfo()));

  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          viewer, SLOT(update()));

  connect(scene, SIGNAL(updated()),
          viewer, SLOT(update()));

  connect(scene, SIGNAL(updated()),
          this, SLOT(selectionChanged()));

  connect(scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)),
          this, SLOT(removeManipulatedFrame(CGAL::Three::Scene_item*)));

  connect(scene, SIGNAL(updated_bbox(bool)),
          this, SLOT(updateViewerBBox(bool)));

  connect(scene, SIGNAL(selectionChanged(int)),
          this, SLOT(selectSceneItem(int)));

  connect(scene, SIGNAL(itemPicked(const QModelIndex &)),
          this, SLOT(recenterSceneView(const QModelIndex &)));

  connect(sceneView->selectionModel(), 
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(updateInfo()));

  connect(sceneView->selectionModel(), 
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(updateDisplayInfo()));

  connect(sceneView->selectionModel(), 
          SIGNAL(selectionChanged ( const QItemSelection & , const QItemSelection & ) ),
          this, SLOT(selectionChanged()));

  sceneView->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(sceneView, SIGNAL(customContextMenuRequested(const QPoint & )),
          this, SLOT(showSceneContextMenu(const QPoint &)));

  connect(sceneView, SIGNAL(expanded(QModelIndex)),
          this, SLOT(setExpanded(QModelIndex)));

  connect(sceneView, SIGNAL(collapsed(QModelIndex)),
          this, SLOT(setCollapsed(QModelIndex)));
  connect(this, SIGNAL(collapsed(QModelIndex)),
          scene, SLOT(setCollapsed(QModelIndex)));
  connect(this, SIGNAL(expanded(QModelIndex)),
          scene, SLOT(setExpanded(QModelIndex)));

  connect(scene, SIGNAL(restoreCollapsedState()),
          this, SLOT(restoreCollapseState()));

  connect(viewer, SIGNAL(selected(int)),
          this, SLOT(selectSceneItem(int)));
  connect(viewer, SIGNAL(selectedPoint(double, double, double)),
          this, SLOT(showSelectedPoint(double, double, double)));

  connect(viewer, SIGNAL(selectionRay(double, double, double, 
                                      double, double, double)),
          scene, SIGNAL(selectionRay(double, double, double,
                                     double, double, double)));

  connect(viewer, SIGNAL(requestContextMenu(QPoint)),
          this, SLOT(contextMenuRequested(QPoint)));
  connect(viewer, SIGNAL(sendMessage(QString)),
          this, SLOT(information(QString)));

  // The contextMenuPolicy of infoLabel is now the default one, so that one
  // can easily copy-paste its text.
  // connect(ui->infoLabel, SIGNAL(customContextMenuRequested(const QPoint & )),
  //         this, SLOT(showSceneContextMenu(const QPoint &)));
  connect(ui->actionRecenterScene, SIGNAL(triggered()),
          viewer, SLOT(update()));
  connect(ui->actionAntiAliasing, SIGNAL(toggled(bool)),
          viewer, SLOT(setAntiAliasing(bool)));

  connect(ui->actionDrawTwoSides, SIGNAL(toggled(bool)),
          viewer, SLOT(setTwoSides(bool)));
  connect(ui->actionQuickCameraMode, SIGNAL(toggled(bool)),
          viewer, SLOT(setFastDrawing(bool)));
  connect(ui->actionSwitchProjection, SIGNAL(toggled(bool)),
          viewer, SLOT(SetOrthoProjection(bool)));

  // add the "About CGAL..." and "About demo..." entries
  this->addAboutCGAL();
  this->addAboutDemo(":/cgal/Polyhedron_3/about.html");

  // Connect the button "addButton" with actionLoad
  ui->addButton->setDefaultAction(ui->actionLoad);
  // Same with "removeButton" and "duplicateButton"
  ui->removeButton->setDefaultAction(ui->actionErase);
  ui->duplicateButton->setDefaultAction(ui->actionDuplicate);

  // Connect actionQuit (Ctrl+Q) and qApp->quit()
  connect(ui->actionQuit, SIGNAL(triggered()),
          this, SLOT(quit()));
  // Connect "Select all items"
  connect(ui->actionSelectAllItems, SIGNAL(triggered()),
          this, SLOT(selectAll()));

  connect(ui->actionColorItems, SIGNAL(triggered()),
          this, SLOT(colorItems()));

  // Recent files menu
  this->addRecentFiles(ui->menuFile, ui->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));

  // Reset the "Operation menu"
  clearMenu(ui->menuOperations);

#ifdef QT_SCRIPT_LIB
  std::cerr << "Enable scripts.\n";
  script_engine = new QScriptEngine(this);
  qScriptRegisterMetaType<CGAL::Three::Scene_item*>(script_engine,
                                       myScene_itemToScriptValue,
                                       myScene_itemFromScriptValue);
#  ifdef QT_SCRIPTTOOLS_LIB
  QScriptEngineDebugger* debugger = new QScriptEngineDebugger(this);
  debugger->setObjectName("qt script debugger");
  QAction* debuggerMenuAction = 
    menuBar()->addMenu(debugger->createStandardMenu());
  debuggerMenuAction->setText(tr("Qt Script &Debug"));
  for(unsigned int i = 0; i < 9; ++i)
  {
    QDockWidget* dock = new QDockWidget(debug_widgets_names[i], this);
    dock->setObjectName(debug_widgets_names[i]);
    dock->setWidget(debugger->widget(debug_widgets[i]));
    this->addDockWidget(Qt::BottomDockWidgetArea, dock);
    dock->hide();
  }
  debugger->setAutoShowStandardWindow(false);
  debugger->attachTo(script_engine);
#  endif // QT_SCRIPTTOOLS_LIB
  QScriptValue fun = script_engine->newFunction(myPrintFunction);
  script_engine->globalObject().setProperty("print", fun);
  
  //  evaluate_script("print('hello', 'world', 'from QtScript!')");
  QScriptValue mainWindowObjectValue = script_engine->newQObject(this);
  script_engine->globalObject().setProperty("main_window", mainWindowObjectValue);

  QScriptValue sceneObjectValue = script_engine->newQObject(scene);
  mainWindowObjectValue.setProperty("scene", sceneObjectValue);
  script_engine->globalObject().setProperty("scene", sceneObjectValue);

  QScriptValue viewerObjectValue = script_engine->newQObject(viewer);
  mainWindowObjectValue.setProperty("viewer", viewerObjectValue);
  script_engine->globalObject().setProperty("viewer", viewerObjectValue);

  QScriptValue cameraObjectValue = script_engine->newQObject(viewer->camera());
  viewerObjectValue.setProperty("camera", cameraObjectValue);
  script_engine->globalObject().setProperty("camera", cameraObjectValue);

  evaluate_script("var plugins = new Array();");
#  ifdef QT_SCRIPTTOOLS_LIB
  QScriptValue debuggerObjectValue = script_engine->newQObject(debugger);
  script_engine->globalObject().setProperty("debugger", debuggerObjectValue);
#  endif
#endif

  readSettings(); // Among other things, the column widths are stored.

  // Load plugins, and re-enable actions that need it.
  loadPlugins();

  // Setup the submenu of the View menu that can toggle the dockwidgets
  Q_FOREACH(QDockWidget* widget, findChildren<QDockWidget*>()) {
    ui->menuDockWindows->addAction(widget->toggleViewAction());
  }
  ui->menuDockWindows->removeAction(ui->dummyAction);


  this->readState("MainWindow", Size|State);

  //Manages the group_item creation
  actionAddToGroup= new QAction("Add New Group", this);

  if(actionAddToGroup) {
    connect(actionAddToGroup, SIGNAL(triggered()),
            this, SLOT(makeNewGroup()));
  }

  QMenu* menuFile = findChild<QMenu*>("menuFile");
  insertActionBeforeLoadPlugin(menuFile, actionAddToGroup);
  statistics_dlg = NULL;
  statistics_ui = new Ui::Statistics_on_item_dialog();

  actionResetDefaultLoaders = new QAction("Reset Default Loaders",this);

#ifdef QT_SCRIPT_LIB
  // evaluate_script("print(plugins);");
  Q_FOREACH(QAction* action, findChildren<QAction*>()) {
    if(action->objectName() != "") {
      QScriptValue objectValue = script_engine->newQObject(action);
      script_engine->globalObject().setProperty(action->objectName(),
                                                objectValue);
    }
  }
  // debugger->action(QScriptEngineDebugger::InterruptAction)->trigger();
#endif

  // setup menu filtering
  connect(ui->menuOperations, SIGNAL(aboutToShow()), this, SLOT(filterOperations()));
}

//Recursive function that do a pass over a menu and its sub-menus(etc.) and hide them when they are empty
void filterMenuOperations(QMenu* menu)
{
    Q_FOREACH(QAction* action, menu->actions()) {
        if(QMenu* menu = action->menu())
        {
            filterMenuOperations(menu);
            action->setVisible(!(menu->isEmpty()));
        }
    }

}

void MainWindow::filterOperations()
{
  Q_FOREACH(const PluginNamePair& p, plugins) {
    Q_FOREACH(QAction* action, p.first->actions()) {
        action->setVisible( p.first->applicable(action) );
    }
  }
  // do a pass over all menus in Operations and their sub-menus(etc.) and hide them when they are empty
  filterMenuOperations(ui->menuOperations);
}

#include <CGAL/Three/exceptions.h>

void MainWindow::evaluate_script(QString script,
                                 const QString& filename,
                                 const bool quiet) {
  QScriptContext* context = script_engine->currentContext();
  QScriptValue object = context->activationObject();
  QScriptValue former_current_filename = object.property("current_filename");;
  object.setProperty("current_filename", filename);

  QScriptValue value = script_engine->evaluate(script, filename);
  if(script_engine->hasUncaughtException()) {
    QScriptValue js_exception = script_engine->uncaughtException();
    QScriptValue js_bt =js_exception.property("backtrace");
    QStringList bt = script_engine->uncaughtExceptionBacktrace();
    if(js_bt.isValid()) {
      QStringList other_bt;
      qScriptValueToSequence(js_bt, other_bt);
      if(!other_bt.isEmpty()) bt = other_bt;
    }
    if(!quiet) {
      QTextStream err(stderr);
      err << "Qt Script exception:\n"
          << js_exception.toString()
          << "\nBacktrace:\n";
      Q_FOREACH(QString line, bt) {
        err << "  " << line << "\n";
      }
    }
    throw CGAL::Three::Script_exception
       (script_engine->uncaughtException().toString(), bt);
  }
  else if(!quiet && !value.isNull() && !value.isUndefined()) {
    QTextStream(stderr) << "Qt Script evaluated to \""
                        << value.toString() << "\"\n";
  }

  object.setProperty("current_filename", former_current_filename);
}

void MainWindow::evaluate_script_quiet(QString script,
                                       const QString& filename)
{
  evaluate_script(script, filename, true);
}

void MainWindow::enableScriptDebugger(bool b /* = true */)
{
  Q_UNUSED(b);
#ifdef QT_SCRIPT_LIB
#  ifdef QT_SCRIPTTOOLS_LIB
  QScriptEngineDebugger* debugger =
    findChild<QScriptEngineDebugger*>("qt script debugger");
  if(debugger) {
    if(b) {
      debugger->action(QScriptEngineDebugger::InterruptAction)->trigger();
    }
    else {
      std::cerr << "Detach the script debugger\n";
      debugger->detach();
    }
  }
  return;
#  endif
#endif
  // If we are here, then the debugger is not available
  this->error(tr("Your version of Qt is too old, and for that reason "
                 "the Qt Script Debugger is not available."));
}

namespace {
bool actionsByName(QAction* x, QAction* y) {
  return x->text() < y->text();
}
}

//Recursively creates all subMenus containing an action.
// In the current implementation, there is a bug if a menu
// and a submenu have the same name (cf map menu_map).
void MainWindow::setMenus(QString name, QString parentName, QAction* a )
{
  QString menuName, subMenuName;
  if (name.isNull())
    return;
  int slash_index = name.indexOf('/');

  if(slash_index==-1)
    menuName= name; // no extra sub-menu
  else
  {
    int l = name.length();
    menuName=name.mid(0,slash_index);
    subMenuName=name.mid(slash_index+1,l-slash_index-1);
    // recursively create sub-menus
    setMenus(subMenuName, menuName, a);
  }

  //Create the menu if it does not already exist
  if(!menu_map.contains(menuName))
    menu_map[menuName] = new QMenu(menuName, this);

  //Create the parent menu if it does not already exist
  if(!menu_map.contains(parentName))
    menu_map[parentName] = new QMenu(parentName, this);
  // add the submenu in the menu
  menu_map[parentName]->addMenu(menu_map[menuName]);

  // only add the action in the last submenu
  if(slash_index==-1)
  {
    ui->menuOperations->removeAction(a);
    menu_map[menuName]->addAction(a);
  }
}

void MainWindow::load_plugin(QString fileName, bool blacklisted)
{
    if(fileName.contains("plugin") && QLibrary::isLibrary(fileName)) {
      //set plugin name
      QFileInfo fileinfo(fileName);
      //set plugin name
      QString name = fileinfo.fileName();
      name.remove(QRegExp("^lib"));
      name.remove(QRegExp("\\..*"));
      //do not load it if it is in the blacklist
      if(blacklisted)
      {
        if ( plugin_blacklist.contains(name) ){
          qDebug("### Ignoring plugin \"%s\".", qPrintable(fileName));
          return;
        }
      }
      QDebug qdebug = qDebug();
      qdebug << "### Loading \"" << fileName.toUtf8().data() << "\"... ";
      QPluginLoader loader;
      loader.setFileName(fileinfo.absoluteFilePath());
      QObject *obj = loader.instance();
      if(obj) {
        obj->setObjectName(name);
        bool init1 = initPlugin(obj);
        bool init2 = initIOPlugin(obj);
        if (!init1 && !init2)
          qdebug << "not for this program";
        else
          qdebug << "success";
      }
      else {
        qdebug << "error: " << qPrintable(loader.errorString());
      }
    }
}

void MainWindow::loadPlugins()
{
  Q_FOREACH(QObject *obj, QPluginLoader::staticInstances())
  {
    initPlugin(obj);
    initIOPlugin(obj);
  }
  QList<QDir> plugins_directories;
  QString dirPath = qApp->applicationDirPath();
  plugins_directories<<dirPath;

  QDir msvc_dir(dirPath);
  QString build_dir_name = msvc_dir.dirName();//Debug or Release for msvc
  msvc_dir.cdUp();

  QFileInfoList filist = QDir(dirPath).entryInfoList();
  filist << msvc_dir.entryInfoList();

  Q_FOREACH(QFileInfo fileinfo, filist)
  {
      //checks if the path leads to a directory
      if(fileinfo.baseName().contains("Plugins"))
      {
        QString plugins_dir = fileinfo.absolutePath();
        plugins_dir.append("/").append(fileinfo.baseName());

        Q_FOREACH(QString package_dir,
                  QDir(plugins_dir).entryList(QDir::Dirs))
        {
          QString package_dir_path(plugins_dir);
          package_dir_path.append("/").append(package_dir);

          QString libdir_path(package_dir_path);
          libdir_path.append("/").append(build_dir_name);

          if (QDir(libdir_path).exists())
            plugins_directories << QDir(libdir_path);
          else
            plugins_directories << QDir(package_dir_path);
        }
      }
  }
  QString env_path = qgetenv("POLYHEDRON_DEMO_PLUGINS_PATH");
  if(!env_path.isEmpty()) {
    Q_FOREACH (QString pluginsDir, 
               env_path.split(":", QString::SkipEmptyParts)) {
      QDir dir(pluginsDir);
      if(dir.isReadable())
        plugins_directories << dir;
    }
  }

  QSet<QString> loaded;
  Q_FOREACH (QDir pluginsDir, plugins_directories) {
    qDebug("# Looking for plugins in directory \"%s\"...",
           qPrintable(pluginsDir.absolutePath()));
    Q_FOREACH(QString fileName, pluginsDir.entryList(QDir::Files))
    {
      QString abs_name = pluginsDir.absoluteFilePath(fileName);
      if(loaded.find(abs_name) == loaded.end())
      {
        load_plugin(abs_name, true);
        loaded.insert(abs_name);
      }
    }
  }
  updateMenus();
}
  //Creates sub-Menus for operations.
void MainWindow::updateMenus()
{
  QList<QAction*> as = ui->menuOperations->actions();
  Q_FOREACH(QAction* a, as)
  {
    QString menuPath = a->property("subMenuName").toString();
    setMenus(menuPath, ui->menuOperations->title(), a);
  }
  // sort the operations menu by name
  as = ui->menuOperations->actions();
  qSort(as.begin(), as.end(), actionsByName);
  ui->menuOperations->clear();
  ui->menuOperations->addActions(as);
}

bool MainWindow::hasPlugin(const QString& pluginName) const
{
  Q_FOREACH(const PluginNamePair& p, plugins) {
    if(p.second == pluginName) return true;
  }
  return false;
}

bool MainWindow::initPlugin(QObject* obj)
{
  QObjectList childs = this->children();
  CGAL::Three::Polyhedron_demo_plugin_interface* plugin =
    qobject_cast<CGAL::Three::Polyhedron_demo_plugin_interface*>(obj);
  if(plugin) {
    // Call plugin's init() method
    obj->setParent(this);
    plugin->init(this, this->scene, this);
    plugins << qMakePair(plugin, obj->objectName());
#ifdef QT_SCRIPT_LIB
    QScriptValue objectValue =
      script_engine->newQObject(obj);
    script_engine->globalObject().setProperty(obj->objectName(), objectValue);
    evaluate_script_quiet(QString("plugins.push(%1);").arg(obj->objectName()));
#endif

    Q_FOREACH(QAction* action, plugin->actions()) {
      // If action does not belong to the menus, add it to "Operations" menu
      if(!childs.contains(action)) {
        ui->menuOperations->addAction(action);
      }
      // Show and enable menu item
      addAction(action);
    }
    return true;
  }
  else 
    return false;
}

bool MainWindow::initIOPlugin(QObject* obj)
{
  CGAL::Three::Polyhedron_demo_io_plugin_interface* plugin =
    qobject_cast<CGAL::Three::Polyhedron_demo_io_plugin_interface*>(obj);
  if(plugin) {
    io_plugins << plugin;
    return true;
  }
  else 
    return false;
}

void MainWindow::clearMenu(QMenu* menu)
{
  Q_FOREACH(QAction* action, menu->actions())
  {
    QMenu* menu = action->menu();
    if(menu) {
      clearMenu(menu);
    }
    action->setVisible(false);
  }
  menu->menuAction()->setEnabled(false);
}

void MainWindow::addAction(QAction* action)
{
  if(!action) return;

  action->setVisible(true);
  action->setEnabled(true);
  Q_FOREACH(QWidget* widget, action->associatedWidgets())
  {
//     qDebug() << QString("%1 (%2)\n")
//       .arg(widget->objectName())
//       .arg(widget->metaObject()->className());
    QMenu* menu = qobject_cast<QMenu*>(widget);
    if(menu)
    {
      addAction(menu->menuAction());
    }
  }
}

void MainWindow::addAction(QString actionName,
                           QString actionText,
                           QString menuName) {
  QMenu* menu = 0;
  Q_FOREACH(QAction* action, findChildren<QAction*>()) {
    if(!action->menu()) continue;
    QString menuText = action->menu()->title();
    if(menuText != menuName) continue;
    menu = action->menu();
  }
  if(menu == 0) {
    menu = new QMenu(menuName, this);
    menuBar()->insertMenu(ui->menuView->menuAction(), menu);
  }
  QAction* action = new QAction(actionText, this);
  action->setObjectName(actionName);
  menu->addAction(action);
#ifdef QT_SCRIPT_LIB
  QScriptValue objectValue = script_engine->newQObject(action);
  script_engine->globalObject().setProperty(action->objectName(),
                                            objectValue);
#endif
}

void MainWindow::viewerShow(float xmin,
                            float ymin,
                            float zmin,
                            float xmax,
                            float ymax,
                            float zmax)
{
  qglviewer::Vec
    min_(xmin, ymin, zmin),
    max_(xmax, ymax, zmax);

  if(min_ == max_) return viewerShow(xmin, ymin, zmin);

#if QGLVIEWER_VERSION >= 0x020502
  viewer->camera()->setPivotPoint((min_+max_)*0.5);
#else
  viewer->camera()->setRevolveAroundPoint((min_+max_)*0.5);
#endif

  qglviewer::ManipulatedCameraFrame backup_frame(*viewer->camera()->frame());
  viewer->camera()->fitBoundingBox(min_, max_);
  qglviewer::ManipulatedCameraFrame new_frame(*viewer->camera()->frame());
  *viewer->camera()->frame() = backup_frame;
  viewer->camera()->interpolateTo(new_frame, 1.f);
  viewer->setVisualHintsMask(1);
}

void MainWindow::viewerShow(float x, float y, float z) {
  // viewer->camera()->lookAt(qglviewer::Vec(x, y, z));

  qglviewer::ManipulatedCameraFrame backup_frame(*viewer->camera()->frame());
  viewer->camera()->fitSphere(qglviewer::Vec(x, y, z),
                              viewer->camera()->sceneRadius()/100);
  qglviewer::ManipulatedCameraFrame new_frame(*viewer->camera()->frame());
  *viewer->camera()->frame() = backup_frame;
  viewer->camera()->interpolateTo(new_frame, 1.f);
  viewer->setVisualHintsMask(1);

#if QGLVIEWER_VERSION >= 0x020502
  viewer->camera()->setPivotPoint(qglviewer::Vec(x, y, z));
#else
  viewer->camera()->setRevolveAroundPoint(qglviewer::Vec(x, y, z));
#endif
}

void MainWindow::message(QString message, QString colorName, QString font) {
  if (message.endsWith('\n')) {
    message.remove(message.length()-1, 1);
  }
  statusBar()->showMessage(message, 5000);
  message = "<font color=\"" + colorName + "\" style=\"font-style: " + font + ";\" >" +
    message + "</font><br>";
  message = "[" + QTime::currentTime().toString() + "] " + message;
  ui->consoleTextEdit->append(message);
  ui->consoleTextEdit->verticalScrollBar()->setValue(ui->consoleTextEdit->verticalScrollBar()->maximum());
}

void MainWindow::information(QString text) {
  this->message("INFO: " + text, "");
}

void MainWindow::warning(QString text) {
  this->message("WARNING: " + text, "blue");
}

void MainWindow::error(QString text) {
  this->message("ERROR: " + text, "red");
}

void MainWindow::updateViewerBBox(bool recenter = true)
{
  const Scene::Bbox bbox = scene->bbox();
#if QGLVIEWER_VERSION >= 0x020502
    qglviewer::Vec center = viewer->camera()->pivotPoint();
#else
    qglviewer::Vec center = viewer->camera()->revolveAroundPoint();
#endif
  const double xmin = bbox.xmin();
  const double ymin = bbox.ymin();
  const double zmin = bbox.zmin();
  const double xmax = bbox.xmax();
  const double ymax = bbox.ymax();
  const double zmax = bbox.zmax();


  qglviewer::Vec 
    vec_min(xmin, ymin, zmin),
    vec_max(xmax, ymax, zmax),
    bbox_center((xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2);
  qglviewer::Vec offset(0,0,0);
  double l_dist = (std::max)((std::abs)(bbox_center.x - viewer->offset().x),
                      (std::max)((std::abs)(bbox_center.y - viewer->offset().y),
                          (std::abs)(bbox_center.z - viewer->offset().z)));
  if((std::log2)(l_dist) > 13.0 )
    for(int i=0; i<3; ++i)
    {
      offset[i] = -bbox_center[i];

    }
  if(offset != viewer->offset())
  {
    viewer->setOffset(offset);
    for(int i=0; i<scene->numberOfEntries(); ++i)
    {
      scene->item(i)->invalidateOpenGLBuffers();
      scene->item(i)->itemChanged();
    }
  }


  viewer->setSceneBoundingBox(vec_min,
                              vec_max);
  if(recenter)
  {
    viewer->camera()->showEntireScene();
  }
  else
  {
#if QGLVIEWER_VERSION >= 0x020502
    viewer->camera()->setPivotPoint(center);
#else
    viewer->camera()->setRevolveAroundPoint(center);
#endif
  }
}

void MainWindow::reloadItem() {
  QAction* sender_action = qobject_cast<QAction*>(sender());
  if(!sender_action) return;
  
  Scene_item* item = (Scene_item*)sender_action->data().value<void*>();
  if(!item) {
    std::cerr << "Cannot reload item: "
              << "the reload action has not item attached\n";
    return;
  }
  if(!item) {
    std::cerr << "Cannot reload item: "
              << "the reload action has a QObject* pointer attached\n"
              << "that is not a Scene_item*\n";
    return;
  }
  QString filename = item->property("source filename").toString();
  QString loader_name = item->property("loader_name").toString();
  if(filename.isEmpty() || loader_name.isEmpty()) {
    std::cerr << "Cannot reload item: "
              << "the item has no \"source filename\" or no \"loader_name\" attached\n";
    return;
  }

  CGAL::Three::Polyhedron_demo_io_plugin_interface* fileloader = findLoader(loader_name);
  QFileInfo fileinfo(filename);

  CGAL::Three::Scene_item* new_item = loadItem(fileinfo, fileloader);

  new_item->setName(item->name());
  new_item->setColor(item->color());
  new_item->setRenderingMode(item->renderingMode());
  new_item->setVisible(item->visible());
  Scene_item_with_properties *property_item = dynamic_cast<Scene_item_with_properties*>(new_item);
  if(property_item)
    property_item->copyProperties(item);
  scene->replaceItem(scene->item_id(item), new_item, true);
  new_item->invalidateOpenGLBuffers();
  item->deleteLater();
}

CGAL::Three::Polyhedron_demo_io_plugin_interface* MainWindow::findLoader(const QString& loader_name) const {
  Q_FOREACH(CGAL::Three::Polyhedron_demo_io_plugin_interface* io_plugin,
            io_plugins) {
    if(io_plugin->name() == loader_name) {
      return io_plugin;
    }
  }
  throw std::invalid_argument(QString("No loader found with the name %1 available")
                              .arg(loader_name).toStdString()) ;
}

bool MainWindow::file_matches_filter(const QString& filters,
                                     const QString& filename )
{
  QFileInfo fileinfo(filename);
  QString filename_striped=fileinfo.fileName();

  //match all filters between ()
  QRegExp all_filters_rx("\\((.*)\\)");

  QStringList split_filters = filters.split(";;");
  Q_FOREACH(const QString& filter, split_filters) {
    //extract filters
    if ( all_filters_rx.indexIn(filter)!=-1 ){
      Q_FOREACH(const QString& pattern,all_filters_rx.cap(1).split(' ')){
        QRegExp rx(pattern);
        rx.setPatternSyntax(QRegExp::Wildcard);
        if ( rx.exactMatch(filename_striped) ){
          return true;
        }
      }
    }
  }
  return false;
}

void MainWindow::open(QString filename)
{
  QFileInfo fileinfo(filename);

#ifdef QT_SCRIPT_LIB
  // Handles the loading of script file from the command line arguments,
  // and the special command line arguments that start with "javascript:"
  // or "qtscript:"
  QString program;
  if(filename.startsWith("javascript:")) {
    program=filename.right(filename.size() - 11);
  }
  if(filename.startsWith("qtscript:")) {
    program=filename.right(filename.size() - 9);
  }
  if(filename.endsWith(".js")) {
    loadScript(fileinfo);
    return;
  }
  if(!program.isEmpty())
  {
    {
      QTextStream(stderr) << "Execution of script \"" 
                          << filename << "\"\n";
                          // << filename << "\", with following content:\n"
                          // << program;
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    evaluate_script(program, filename);
    QApplication::restoreOverrideCursor();
    return;
  }
#endif

  if ( !fileinfo.exists() ){
    QMessageBox::warning(this,
                         tr("Cannot open file"),
                         tr("File %1 does not exist.")
                         .arg(filename));
    return;
  }


  QStringList selected_items;
  QStringList all_items;

  QMap<QString,QString>::iterator dfs_it = 
    default_plugin_selection.find( fileinfo.completeSuffix() );
  
  if ( dfs_it==default_plugin_selection.end() )
  {
    // collect all io_plugins and offer them to load if the file extension match one name filter
    // also collect all available plugin in case of a no extension match
    Q_FOREACH(CGAL::Three::Polyhedron_demo_io_plugin_interface* io_plugin, io_plugins) {
      if ( !io_plugin->canLoad() ) continue;
      all_items << io_plugin->name();
      if ( file_matches_filter(io_plugin->loadNameFilters(), filename) )
        selected_items << io_plugin->name();
    }
  }
  else
    selected_items << *dfs_it;
  
  bool ok;
  std::pair<QString, bool> load_pair;
  
  switch( selected_items.size() )
  {
    case 1:
      load_pair = std::make_pair(selected_items.first(), false);
      ok=true;
      break;
    case 0:
      load_pair = File_loader_dialog::getItem(fileinfo.fileName(), all_items, &ok);
      break;
    default:
      load_pair = File_loader_dialog::getItem(fileinfo.fileName(), selected_items, &ok);
  }

  viewer->makeCurrent();
  if(!ok || load_pair.first.isEmpty()) { return; }
  
  if (load_pair.second)
  {
    connect(actionResetDefaultLoaders, SIGNAL(triggered()),
            this, SLOT(reset_default_loaders()));
    default_plugin_selection[fileinfo.completeSuffix()]=load_pair.first;
    insertActionBeforeLoadPlugin(ui->menuFile, actionResetDefaultLoaders);
  }
  
  
  QSettings settings;
  settings.setValue("OFF open directory",
                    fileinfo.absoluteDir().absolutePath());
  CGAL::Three::Scene_item* scene_item = loadItem(fileinfo, findLoader(load_pair.first));
  if(scene_item != 0) {
    this->addToRecentFiles(fileinfo.absoluteFilePath());
  }

  selectSceneItem(scene->addItem(scene_item));

  CGAL::Three::Scene_group_item* group =
          qobject_cast<CGAL::Three::Scene_group_item*>(scene_item);
  if(group)
    scene->redraw_model();
}

bool MainWindow::open(QString filename, QString loader_name) {
  QFileInfo fileinfo(filename); 
  boost::optional<CGAL::Three::Scene_item*> item_opt;
  CGAL::Three::Scene_item* item = 0;
  try {
    item_opt = wrap_a_call_to_cpp
      ([this, fileinfo, loader_name]()
       {
         return loadItem(fileinfo, findLoader(loader_name));
       },
       this, __FILE__, __LINE__
       );
    if(!item_opt) return false;
    else item = *item_opt;
  }
  catch(std::logic_error e) {
    std::cerr << e.what() << std::endl;
    return false;
  }
  selectSceneItem(scene->addItem(item));

  CGAL::Three::Scene_group_item* group =
          qobject_cast<CGAL::Three::Scene_group_item*>(item);
  if(group)
    scene->redraw_model();

  return true;
}


CGAL::Three::Scene_item* MainWindow::loadItem(QFileInfo fileinfo, CGAL::Three::Polyhedron_demo_io_plugin_interface* loader) {
  CGAL::Three::Scene_item* item = NULL;
  if(!fileinfo.isFile() || !fileinfo.isReadable()) {
    throw std::invalid_argument(QString("File %1 is not a readable file.")
                                .arg(fileinfo.absoluteFilePath()).toStdString());
  }

  QApplication::setOverrideCursor(Qt::WaitCursor);
  item = loader->load(fileinfo);
  QApplication::restoreOverrideCursor();
  if(!item) {
    throw std::logic_error(QString("Could not load item from file %1 using plugin %2")
                           .arg(fileinfo.absoluteFilePath()).arg(loader->name()).toStdString());
  }

  item->setProperty("source filename", fileinfo.absoluteFilePath());
  item->setProperty("loader_name", loader->name());
  return item;
}


void MainWindow::setFocusToQuickSearch()
{
  ui->searchEdit->setFocus(Qt::ShortcutFocusReason);
}

void MainWindow::selectSceneItem(int i)
{
  if(i < 0 || i >= scene->numberOfEntries()) {
    sceneView->selectionModel()->clearSelection();
    updateInfo();
    updateDisplayInfo();
  }
  else {
    QItemSelection s =
      proxyModel->mapSelectionFromSource(scene->createSelection(i));

    sceneView->selectionModel()->select(s,
                                        QItemSelectionModel::ClearAndSelect);
  }
}


void MainWindow::showSelectedPoint(double x, double y, double z)
{
  static double x_prev = 0;
  static double y_prev = 0;
  static double z_prev = 0;
  double dist = std::sqrt((x-x_prev)*(x-x_prev) + (y-y_prev)*(y-y_prev) + (z-z_prev)*(z-z_prev)); 
  information(QString("Selected point: (%1, %2, %3) distance to previous: %4").
              arg(x, 0, 'g', 10).
              arg(y, 0, 'g', 10).
              arg(z, 0, 'g', 10).
              arg(dist,0,'g',10));
  x_prev = x;
  y_prev = y;
  z_prev = z;
}

void MainWindow::unSelectSceneItem(int i)
{
  removeSceneItemFromSelection(i);
}

void MainWindow::addSceneItemInSelection(int i)
{
  QItemSelection s =
    proxyModel->mapSelectionFromSource(scene->createSelection(i));
  sceneView->selectionModel()->select(s, QItemSelectionModel::Select);
  scene->itemChanged(i);
}

void MainWindow::removeSceneItemFromSelection(int i)
{
  QItemSelection s =
    proxyModel->mapSelectionFromSource(scene->createSelection(i));
  sceneView->selectionModel()->select(s,
                                      QItemSelectionModel::Deselect);
  scene->itemChanged(i);
}

void MainWindow::selectAll()
{
  QItemSelection s =
    proxyModel->mapSelectionFromSource(scene->createSelectionAll());
  sceneView->selectionModel()->select(s, 
                                      QItemSelectionModel::ClearAndSelect);
}

int MainWindow::getSelectedSceneItemIndex() const
{
  QModelIndexList selectedRows = sceneView->selectionModel()->selectedIndexes();
  if(selectedRows.size() == 0)
    return -1;
  else {
    QModelIndex i = proxyModel->mapToSource(selectedRows.first());
    return scene->getIdFromModelIndex(i);
  }
}

QList<int> MainWindow::getSelectedSceneItemIndices() const
{
  QModelIndexList selectedIndices = sceneView->selectionModel()->selectedIndexes();
  QList<int> result;
  Q_FOREACH(QModelIndex index, selectedIndices) {
      int temp = scene->getIdFromModelIndex(proxyModel->mapToSource(index));
      if(!result.contains(temp))
          result<<temp;
  }
  return result;
}

void MainWindow::selectionChanged()
{
  scene->setSelectedItemIndex(getSelectedSceneItemIndex());
  scene->setSelectedItemsList(getSelectedSceneItemIndices());
  CGAL::Three::Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(item != NULL && item->manipulatable()) {
    viewer->setManipulatedFrame(item->manipulatedFrame());
  } else {
    viewer->setManipulatedFrame(0);
  }
  if(viewer->manipulatedFrame() == 0) {
    Q_FOREACH(CGAL::Three::Scene_item* item, scene->entries()) {
      if(item->manipulatable() && item->manipulatedFrame() != 0) {
        if(viewer->manipulatedFrame() != 0) {
          // there are at least two possible frames
          viewer->setManipulatedFrame(0);
          break;
        } else {
          viewer->setManipulatedFrame(item->manipulatedFrame());
        }
      }
    }
  }
  if(viewer->manipulatedFrame() != 0) {
    connect(viewer->manipulatedFrame(), SIGNAL(modified()),
            this, SLOT(updateInfo()));
  }
  viewer->update();
}

void MainWindow::contextMenuRequested(const QPoint& global_pos) {
  int index = scene->mainSelectionIndex();
  showSceneContextMenu(index, global_pos);
}

void MainWindow::showSceneContextMenu(int selectedItemIndex,
                                      const QPoint& global_pos)
{
  CGAL::Three::Scene_item* item = scene->item(selectedItemIndex);
  if(!item) return;

  const char* prop_name = "Menu modified by MainWindow.";

  QMenu* menu = item->contextMenu();
  if(menu) {
    bool menuChanged = menu->property(prop_name).toBool();
    if(!menuChanged) {
      if(item->has_stats())
      {
        QAction* actionStatistics =
            menu->addAction(tr("Statistics..."));
        actionStatistics->setObjectName("actionStatisticsOnPolyhedron");
        connect(actionStatistics, SIGNAL(triggered()),
                this, SLOT(statisticsOnItem()));
      }
      menu->addSeparator();
      if(!item->property("source filename").toString().isEmpty()) {
        QAction* reload = menu->addAction(tr("&Reload Item from File"));
        reload->setData(qVariantFromValue((void*)item));
        connect(reload, SIGNAL(triggered()),
                this, SLOT(reloadItem()));
      }
      QAction* saveas = menu->addAction(tr("&Save as..."));
      saveas->setData(qVariantFromValue((void*)item));
      connect(saveas,  SIGNAL(triggered()),
              this, SLOT(on_actionSaveAs_triggered()));
      QAction* showobject = menu->addAction(tr("&Zoom to this Object"));
      showobject->setData(qVariantFromValue((void*)item));
      connect(showobject, SIGNAL(triggered()),
              this, SLOT(viewerShowObject()));

      menu->setProperty(prop_name, true);
    }
  }
  menu->addMenu(ui->menuOperations);
  if(menu)
    menu->exec(global_pos);
}

void MainWindow::showSceneContextMenu(const QPoint& p) {
  QWidget* sender = qobject_cast<QWidget*>(this->sender());
  if(!sender) return;

  int index = -1;
  if(sender == sceneView) {
      QModelIndex modelIndex = sceneView->indexAt(p);
      if(!modelIndex.isValid())
      {
          const char* prop_name = "Menu modified by MainWindow.";

          QMenu* menu = ui->menuFile;
          if(menu) {
              bool menuChanged = menu->property(prop_name).toBool();
              if(!menuChanged) {
                  menu->setProperty(prop_name, true);
              }
          }
          if(menu)
              menu->exec(sender->mapToGlobal(p));
          return;
      }
      else
      {
          index = scene->getIdFromModelIndex(proxyModel->mapToSource(modelIndex));
          scene->setSelectedItemIndex(index);
      }
  }
  else {
    index = scene->mainSelectionIndex();
  }

  showSceneContextMenu(index, sender->mapToGlobal(p));
}

void MainWindow::removeManipulatedFrame(CGAL::Three::Scene_item* item)
{
  if(item->manipulatable() &&
     item->manipulatedFrame() == viewer->manipulatedFrame()) {
    viewer->setManipulatedFrame(0);
  }
}

void MainWindow::updateInfo() {
  CGAL::Three::Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(item) {
    QString item_text = item->toolTip();
    QString item_filename = item->property("source filename").toString();
    if(item->bbox()!=CGAL::Bbox_3())
      item_text += QString("<div>Bounding box: min (%1,%2,%3), max (%4,%5,%6)</div>")
          .arg(item->bbox().xmin())
          .arg(item->bbox().ymin())
          .arg(item->bbox().zmin())
          .arg(item->bbox().xmax())
          .arg(item->bbox().ymax())
          .arg(item->bbox().zmax());
    if(!item_filename.isEmpty()) {
      item_text += QString("<div>File:<i> %1</div>").arg(item_filename);
    }
    ui->infoLabel->setText(item_text);
  }
  else
    ui->infoLabel->clear();
}
void MainWindow::updateDisplayInfo() {
  CGAL::Three::Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(item)
    ui->displayLabel->setPixmap(item->graphicalToolTip());
  else 
    ui->displayLabel->clear();

}

void MainWindow::readSettings()
{
    QSettings settings;
    // enable anti-aliasing 
    ui->actionAntiAliasing->setChecked(settings.value("antialiasing", false).toBool());
    // read plugin blacklist
    QStringList blacklist=settings.value("plugin_blacklist",QStringList()).toStringList();
    Q_FOREACH(QString name,blacklist){ plugin_blacklist.insert(name); }
    set_facegraph_mode_adapter(settings.value("polyhedron_mode", true).toBool());
}

void MainWindow::writeSettings()
{
  this->writeState("MainWindow");
  {
    QSettings settings;
    settings.setValue("antialiasing", 
                      ui->actionAntiAliasing->isChecked());
    //setting plugin blacklist
    QStringList blacklist;
    Q_FOREACH(QString name,plugin_blacklist){ blacklist << name; }
    if ( !blacklist.isEmpty() ) settings.setValue("plugin_blacklist",blacklist);
    else settings.remove("plugin_blacklist");
    //setting polyhedron mode
    settings.setValue("polyhedron_mode", this->property("is_polyhedron_mode").toBool());
  }
  std::cerr << "Write setting... done.\n";
}

void MainWindow::quit()
{
  close();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    for(int i=0; i<plugins.size(); i++)
    {
      plugins[i].first->closure();
    }
  writeSettings();
  event->accept();
}

bool MainWindow::loadScript(QString filename)
{
  QFileInfo fileinfo(filename);
  boost::optional<bool> opt = wrap_a_call_to_cpp
    ([this, fileinfo] {
      return loadScript(fileinfo);
    }, this, __FILE__, __LINE__, CGAL::Three::PARENT_CONTEXT);
  if(!opt) return false;
  else return *opt;
}

bool MainWindow::loadScript(QFileInfo info)
{
#if defined(QT_SCRIPT_LIB)
  QString program;
  QString filename = info.absoluteFilePath();
  QFile script_file(filename);
  script_file.open(QIODevice::ReadOnly);
  if(!script_file.isReadable()) {
    throw std::ios_base::failure(script_file.errorString().toStdString());
  }
  program = script_file.readAll();
  if(!program.isEmpty())
  {
    QTextStream(stderr) 
      << "Execution of script \"" 
      << filename << "\"\n";
    evaluate_script(program, filename);
    return true;
  }
#endif
  return false;
}

void MainWindow::throw_exception() {
  wrap_a_call_to_cpp([]() {
      throw std::runtime_error("Exception thrown in "
                               "MainWindow::throw_exception()");
    }, this, __FILE__, __LINE__);
}

void MainWindow::on_actionLoadScript_triggered() 
{
#if defined(QT_SCRIPT_LIB)
  QString filename = QFileDialog::getOpenFileName(
    this,
    tr("Select a script to run..."),
    ".",
    "QTScripts (*.js);;All Files (*)");

  loadScript(QFileInfo(filename));
#endif
}

void MainWindow::on_actionLoad_triggered()
{
  QStringList filters;
  // we need to special case our way out of this
  filters << "All Files (*)";


  typedef QMap<QString, CGAL::Three::Polyhedron_demo_io_plugin_interface*> FilterPluginMap;
  FilterPluginMap filterPluginMap;
  
  Q_FOREACH(CGAL::Three::Polyhedron_demo_io_plugin_interface* plugin, io_plugins) {
    QStringList split_filters = plugin->loadNameFilters().split(";;");
    Q_FOREACH(const QString& filter, split_filters) {
      FilterPluginMap::iterator it = filterPluginMap.find(filter);
      if(it != filterPluginMap.end()) {
        qDebug() << "Duplicate Filter: " << it.value()->name();
        qDebug() << "This filter will not be available.";
      } else {
        filterPluginMap[filter] = plugin;
      }
      filters << filter;
    }
  }
  QSettings settings;
  QString directory = settings.value("OFF open directory",
                                     QDir::current().dirName()).toString();

  QFileDialog dialog(this);
  dialog.setDirectory(directory);
  dialog.setNameFilters(filters);
  dialog.setFileMode(QFileDialog::ExistingFiles);

  if(dialog.exec() != QDialog::Accepted) { return; }
  viewer->update();
  FilterPluginMap::iterator it = 
    filterPluginMap.find(dialog.selectedNameFilter());
  
  CGAL::Three::Polyhedron_demo_io_plugin_interface* selectedPlugin = NULL;

  if(it != filterPluginMap.end()) {
    selectedPlugin = it.value();
  }

  std::size_t nb_files = dialog.selectedFiles().size();
  std::vector<QColor> colors_;
  colors_.reserve(nb_files);
  compute_color_map(QColor(100, 100, 255),//Scene_item's default color
                    static_cast<unsigned>(nb_files),
                    std::back_inserter(colors_));
  std::size_t nb_item = -1;
  Q_FOREACH(const QString& filename, dialog.selectedFiles()) {
    CGAL::Three::Scene_item* item = NULL;
    if(selectedPlugin) {
      QFileInfo info(filename);
      item = loadItem(info, selectedPlugin);
      item->setColor(colors_[++nb_item]);
      Scene::Item_id index = scene->addItem(item);
      selectSceneItem(index);
      CGAL::Three::Scene_group_item* group =
              qobject_cast<CGAL::Three::Scene_group_item*>(item);
      if(group)
        scene->redraw_model();
      this->addToRecentFiles(filename);
    } else {
      open(filename);
      scene->item(scene->numberOfEntries()-1)->setColor(colors_[++nb_item]);
    }
  }
}

void MainWindow::on_actionSaveAs_triggered()
{
  Scene_item* item = NULL;
  QAction* sender_action = qobject_cast<QAction*>(sender());
  if(sender_action && !sender_action->data().isNull()) {
    item = (Scene_item*)sender_action->data().value<void*>();
  }

  if(!item)
  {
    item = scene->item(scene->mainSelectionIndex());
  }
  if(!item)
    return;

  QVector<CGAL::Three::Polyhedron_demo_io_plugin_interface*> canSavePlugins;
  QStringList filters;
  Q_FOREACH(CGAL::Three::Polyhedron_demo_io_plugin_interface* plugin, io_plugins) {
    if(plugin->canSave(item)) {
      canSavePlugins << plugin;
      filters += plugin->saveNameFilters();
    }
  }
  QString ext;
  if(!filters.isEmpty())
  {
    QRegExp extensions("\\(\\*\\..+\\)");
    extensions.indexIn(filters.first().split(";;").first());
    ext = extensions.cap();
    filters << tr("All files (*)");
  }
  if(canSavePlugins.isEmpty()) {
    QMessageBox::warning(this,
                         tr("Cannot save"),
                         tr("The selected object %1 cannot be saved.")
                         .arg(item->name()));
    return;
  }
  QString caption = tr("Save %1 to File...%2").arg(item->name()).arg(ext);
  //remove `)`
  ext.chop(1);
  //remove `(*.`
  ext = ext.right(ext.size()-3);
  QString filename = 
    QFileDialog::getSaveFileName(this,
                                 caption,
                                 QString("%1.%2").arg(item->name()).arg(ext),
                                 filters.join(";;"));
  if(filename.isEmpty())
    return;

  viewer->update();
  save(filename, item);
}

void MainWindow::save(QString filename, CGAL::Three::Scene_item* item) {
  QFileInfo fileinfo(filename);
  bool saved = false;
  Q_FOREACH(CGAL::Three::Polyhedron_demo_io_plugin_interface* plugin, io_plugins) {
    if(  plugin->canSave(item) &&
        file_matches_filter(plugin->saveNameFilters(),filename) )
    {
      if(plugin->save(item, fileinfo))
      {
        saved = true;
        break;
      }
    }
  }
  if(!saved)
    QMessageBox::warning(this,
                         tr("Cannot save"),
                         tr("The selected object %1 was not saved. (Maybe a wrong extension ?)")
                         .arg(item->name()));
}

void MainWindow::on_actionSaveSnapshot_triggered()
{
  viewer->saveSnapshot(false);
}

bool MainWindow::on_actionErase_triggered()
{
  int next_index = scene->erase(scene->selectionIndices());
  //Secure the case where erase triggers other items deletions
  if(scene->numberOfEntries()< next_index +1 )
    next_index = -1;
  selectSceneItem(next_index);
  return next_index >= 0;
}

void MainWindow::on_actionEraseAll_triggered()
{
  scene->setSelectedItem(0);
  while(on_actionErase_triggered()) {
  }
}

void MainWindow::on_actionDuplicate_triggered()
{
  int index = scene->duplicate(getSelectedSceneItemIndex());
  selectSceneItem(index);
}

void MainWindow::on_actionShowHide_triggered()
{
  Q_FOREACH(QModelIndex index, sceneView->selectionModel()->selectedRows())
  {
    int i = scene->getIdFromModelIndex(proxyModel->mapToSource(index));
    CGAL::Three::Scene_item* item = scene->item(i);
    item->setVisible(!item->visible());
    scene->itemChanged(i);
  }
}

void MainWindow::on_actionSetPolyhedronA_triggered()
{
  int i = getSelectedSceneItemIndex();
  scene->setItemA(i);
}

void MainWindow::on_actionSetPolyhedronB_triggered()
{
  int i = getSelectedSceneItemIndex();
  scene->setItemB(i);
}
void MainWindow::on_actionPreferences_triggered()
{
  QDialog dialog(this);
  Ui::PreferencesDialog prefdiag;
  prefdiag.setupUi(&dialog);
  if(this->property("is_polyhedron_mode").toBool())
    prefdiag.polyRadioButton->setChecked(true);
  else
    prefdiag.smRadioButton->setChecked(true);
  connect(prefdiag.polyRadioButton, &QRadioButton::toggled,
          this, &MainWindow::set_facegraph_mode_adapter);
  
  QStandardItemModel* iStandardModel = new QStandardItemModel(this);
  //add blacklisted plugins
  Q_FOREACH(QString name, plugin_blacklist)
  {
    QStandardItem* item =  new QStandardItem(name);
    item->setCheckable(true);
    item->setCheckState(Qt::Checked);
    iStandardModel->appendRow(item);
  }

  //add operations plugins
  Q_FOREACH(PluginNamePair pair,plugins){
    QStandardItem* item =  new QStandardItem(pair.second);
    item->setCheckable(true);
    iStandardModel->appendRow(item);
  }
  
  //add io-plugins
  Q_FOREACH(CGAL::Three::Polyhedron_demo_io_plugin_interface* plugin, io_plugins)
  {
    QStandardItem* item =  new QStandardItem(plugin->name());
    item->setCheckable(true);
    if ( plugin_blacklist.contains(plugin->name()) ) item->setCheckState(Qt::Checked);
    iStandardModel->appendRow(item);
  }

  //Setting the model
  prefdiag.listView->setModel(iStandardModel);
  
  dialog.exec();  
  
  if ( dialog.result() )
  {
    plugin_blacklist.clear();
    for (int k=0,k_end=iStandardModel->rowCount();k<k_end;++k)
    {
      QStandardItem* item=iStandardModel->item(k);
      if (item->checkState()==Qt::Checked)
        plugin_blacklist.insert(item->text());
    }
  }
  
  for (int k=0,k_end=iStandardModel->rowCount();k<k_end;++k) delete iStandardModel->item(k);
  delete iStandardModel;
}

void MainWindow::on_actionSetBackgroundColor_triggered()
{
  QColor c =  QColorDialog::getColor();
  if(c.isValid()) {
    viewer->setBackgroundColor(c);
  }
}

void MainWindow::on_actionLookAt_triggered()
{
  Show_point_dialog dialog(this);
  dialog.setWindowTitle(tr("Look at..."));
  int i = dialog.exec();
  if( i == QDialog::Accepted &&
      dialog.has_correct_coordinates() )
  {
    viewerShow((float)dialog.get_x(),
               (float)dialog.get_y(),
               (float)dialog.get_z());
  }
}

void MainWindow::viewerShowObject()
{
  Scene_item* item = NULL;
  QAction* sender_action = qobject_cast<QAction*>(sender());
  if(sender_action && !sender_action->data().isNull()) {
    item = (Scene_item*)sender_action->data().value<void*>();
  }
  if(item) {
    const Scene::Bbox bbox = item->bbox();
    viewerShow((float)bbox.xmin()+viewer->offset().x, (float)bbox.ymin()+viewer->offset().y, (float)bbox.zmin()+viewer->offset().z,
               (float)bbox.xmax()+viewer->offset().x, (float)bbox.ymax()+viewer->offset().y, (float)bbox.zmax()+viewer->offset().z);
  }
}

QString MainWindow::cameraString() const
{
  return viewer->dumpCameraCoordinates();
}

void MainWindow::on_actionDumpCamera_triggered()
{
  information(QString("Camera: %1")
              .arg(cameraString()));
}

void MainWindow::on_actionCopyCamera_triggered()
{
  qApp->clipboard()->setText(this->cameraString());
}

void MainWindow::on_actionPasteCamera_triggered()
{
  QString s = qApp->clipboard()->text();
  viewer->moveCameraToCoordinates(s, 0.5f);
}

void MainWindow::setAddKeyFrameKeyboardModifiers(::Qt::KeyboardModifiers m)
{
  viewer->setAddKeyFrameKeyboardModifiers(m);
}

void MainWindow::on_actionRecenterScene_triggered()
{
  updateViewerBBox();
  viewer->camera()->interpolateToFitScene();
}

void MainWindow::on_actionLoadPlugin_triggered()
{
    //pop a dialog of path selection, get the path and add it to plugins_directory

    QString filters("Library files (*.dll *.DLL *.so *.a *.sl *.dylib *.bundle);;"
                    "Any files (*)");

    QStringList paths = QFileDialog::getOpenFileNames(
                this,
                tr("Select the directory containing your plugins:"),
                ".",filters);
    Q_FOREACH(QString name, paths)
      load_plugin(name, false);

    updateMenus();
}

void MainWindow::recurseExpand(QModelIndex index)
{
    int row = index.row();
    if(index.child(0,0).isValid())
    {
        recurseExpand(index.child(0,0));
    }

    QString name = scene->item(scene->getIdFromModelIndex(index))->name();
        CGAL::Three::Scene_group_item* group =
                qobject_cast<CGAL::Three::Scene_group_item*>(scene->item(scene->getIdFromModelIndex(index)));
        if(group && group->isExpanded())
        {
            sceneView->setExpanded(proxyModel->mapFromSource(index), true);
        }
        else if (group && !group->isExpanded()){
            sceneView->setExpanded(proxyModel->mapFromSource(index), false);
        }

        if( index.sibling(row+1,0).isValid())
            recurseExpand(index.sibling(row+1,0));
}
void MainWindow::restoreCollapseState()
{
    QModelIndex modelIndex = scene->index(0,0,scene->invisibleRootItem()->index());
    if(modelIndex.isValid())
        recurseExpand(modelIndex);
    resetHeader();
}
void MainWindow::makeNewGroup()
{
    Scene_group_item * group = new Scene_group_item();
    scene->addItem(group);
}

void MainWindow::on_upButton_pressed()
{
    scene->moveRowUp();
}

void MainWindow::on_downButton_pressed()
{
    scene->moveRowDown();
}

void MainWindow::recenterSceneView(const QModelIndex &id)
{
    if(id.isValid())
    {
        // mapFromSource is necessary to convert the QModelIndex received
        // from the Scene into a valid QModelIndex in the view, beacuse of
        // the proxymodel
        sceneView->scrollTo(proxyModel->mapFromSource(id));
    }
}

void MainWindow::statisticsOnItem()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  if (statistics_dlg == NULL)
  {
    statistics_dlg = new QDialog(this);
    statistics_ui->setupUi(statistics_dlg);
    connect(statistics_ui->okButtonBox, SIGNAL(accepted()),
            statistics_dlg, SLOT(accept()));
    connect(statistics_ui->updateButton, SIGNAL(clicked()),
            this, SLOT(statisticsOnItem()));
  }
  statistics_ui->label_htmltab->setText(get_item_stats());

  statistics_dlg->show();
  statistics_dlg->raise();

  QApplication::restoreOverrideCursor();
}

/* Creates a string containing an html table. This string is constructed by appending each parts of each row, so that the data can
  depend on the number of selected items. This String is then returned.*/
QString MainWindow::get_item_stats()
{
  //1st step : get all classnames of the selected items
  QList<QString> classnames;
  Q_FOREACH(int id, getSelectedSceneItemIndices())
  {
    QString classname = scene->item(id)->metaObject()->className();
    if(!classnames.contains(classname))
      classnames << classname;
  }
  //2nd step : separate the selection in lists corresponding to their classname
  QVector< QList<Scene_item*> > items;
  items.resize(classnames.size());
  Q_FOREACH(int id, getSelectedSceneItemIndices())
  {
    Scene_item* s_item = scene->item(id);
    for(int i=0; i<items.size(); i++)
      if(classnames.at(i).contains(s_item->metaObject()->className()))
      {
        items[i] << s_item;
        break;
      }
  }
  //last step :: making tables for each type of item
  QString str;
  for(int i=0; i< classnames.size(); i++)
  {
    CGAL::Three::Scene_item::Header_data data = items[i].at(0)->header();
    int title = 0;
    int titles_limit =0;
    if(data.titles.size()>0)
    {
      //1st row : item names
      str.append("<html> <table border=1>""<tr><td colspan = 2></td>");
      Q_FOREACH(Scene_item* sit, items[i])
      {
        str.append(QString("<td>%1</td>").arg(sit->name()));
      }



      for(int j=0; j<data.categories.size(); j++)
      {
        str.append(QString("<tr><th rowspan=%1> %2 </th>")
                   .arg(QString::number(data.categories[j].second))
                   .arg(data.categories[j].first));
        titles_limit+=data.categories[j].second;
        str.append(QString("<td> %1 </td>").arg(data.titles.at(title)));
        Q_FOREACH(Scene_item* sit, items[i])
        {
          str.append(QString("<td>%1</td>").arg(sit->computeStats(title)));
        }
        title++;
        for(;title<titles_limit; title++)
        {
          str.append(QString("</tr><tr><td> %1 </td>").arg(data.titles.at(title)));
          Q_FOREACH(Scene_item* sit, items[i])
          {
            str.append(QString("<td>%1</td>").arg(sit->computeStats(title)));
          }
        }

        str.append("</tr>");
      }

      str.append(QString("</tr>""</table></html>"));
    }
  }
  return str;
}

void MainWindow::setCollapsed(QModelIndex index)
{
  Q_EMIT collapsed(proxyModel->mapToSource(index));
}

void MainWindow::setExpanded(QModelIndex index)
{
  Q_EMIT expanded(proxyModel->mapToSource(index));
}


void MainWindow::on_actionMaxTextItemsDisplayed_triggered()
{
  bool ok;
  bool valid;
  QString text = QInputDialog::getText(this, tr("Maximum Number of Text Items"),
                                       tr("Maximum Text Items Diplayed:"), QLineEdit::Normal,
                                       QString("%1").arg(viewer->textRenderer()->getMax_textItems()), &ok);
  text.toInt(&valid);
  if (ok && valid){
    viewer->textRenderer()->setMax(text.toInt());
    ui->actionMaxTextItemsDisplayed->setText(QString("Set Maximum Text Items Displayed : %1").arg(text.toInt()));
  }
}

void MainWindow::resetHeader()
{
  sceneView->header()->setStretchLastSection(false);
  scene->invisibleRootItem()->setColumnCount(5);
  sceneView->header()->setSectionResizeMode(Scene::NameColumn, QHeaderView::Stretch);
  sceneView->header()->setSectionResizeMode(Scene::ColorColumn, QHeaderView::Fixed);
  sceneView->header()->setSectionResizeMode(Scene::RenderingModeColumn, QHeaderView::ResizeToContents);
  sceneView->header()->setSectionResizeMode(Scene::ABColumn, QHeaderView::Fixed);
  sceneView->header()->setSectionResizeMode(Scene::VisibleColumn, QHeaderView::Fixed);
  sceneView->header()->resizeSection(Scene::ColorColumn, sceneView->header()->fontMetrics().width("_#_"));
  sceneView->resizeColumnToContents(Scene::RenderingModeColumn);
  sceneView->header()->resizeSection(Scene::ABColumn, sceneView->header()->fontMetrics().width(QString("_AB_")));
  sceneView->header()->resizeSection(Scene::VisibleColumn, sceneView->header()->fontMetrics().width(QString("_View_")));
}

void MainWindow::reset_default_loaders()
{
  default_plugin_selection.clear();

  const char* prop_name = "Menu modified by MainWindow.";
  QMenu* menu = ui->menuFile;
  if(!menu)
    return;
  bool menuChanged = menu->property(prop_name).toBool();
  if(!menuChanged) {
    menu->setProperty(prop_name, true);
  }
  QList<QAction*> menuActions = menu->actions();
  menu->removeAction(actionResetDefaultLoaders);
}

void MainWindow::insertActionBeforeLoadPlugin(QMenu* menu, QAction* actionToInsert)
{
  if(menu)
  {
    QList<QAction*> menuActions = menu->actions();
    if(!menuActions.contains(actionToInsert))
      menu->insertAction(ui->actionLoadPlugin, actionToInsert);
  }
}

void MainWindow::colorItems()
{
  std::size_t nb_files = scene->selectionIndices().size();
  std::vector<QColor> colors_;
  colors_.reserve(nb_files);
  compute_color_map(scene->item(scene->selectionIndices().last())->color(),
                    static_cast<unsigned>(nb_files),
                    std::back_inserter(colors_));
  std::size_t nb_item = -1;
  Q_FOREACH(int id, scene->selectionIndices())
  {
    scene->item(id)->setColor(colors_[++nb_item]);
  }
  viewer->update();
}
// Only used to make the doc clearer. Only the adapter is actueally used in the code,
// for signal/slots reasons.
void MainWindow::set_face_graph_default_type(Face_graph_mode m)
{
  this->setProperty("is_polyhedron_mode", m);
}

void MainWindow::set_facegraph_mode_adapter(bool is_polyhedron)
{
  if(is_polyhedron)
   set_face_graph_default_type(POLYHEDRON);
  else
    set_face_graph_default_type(SURFACE_MESH);
}
