#include "config.h"
#include "MainWindow.h"
#include "Scene.h"
#include "Scene_item.h"
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
#include <QMap>

#include <stdexcept>

#ifdef QT_SCRIPT_LIB
#  include <QScriptValue>
#  ifdef QT_SCRIPTTOOLS_LIB
#    include <QScriptEngineDebugger>
#  endif
#endif

#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_io_plugin_interface.h"

#include "ui_MainWindow.h"

#include "Show_point_dialog.h"

#ifdef QT_SCRIPT_LIB
#  include <QScriptEngine>
#  include <QScriptValue>

QScriptValue 
myScene_itemToScriptValue(QScriptEngine *engine, 
                          Scene_item* const &in)
{ 
  return engine->newQObject(in); 
}

void myScene_itemFromScriptValue(const QScriptValue &object, 
                                 Scene_item* &out)
{
  out = qobject_cast<Scene_item*>(object.toQObject()); 
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

  return engine->undefinedValue();
}

MainWindow::~MainWindow()
{
  delete ui;
}

MainWindow::MainWindow(QWidget* parent)
  : CGAL::Qt::DemosMainWindow(parent)
{
  ui = new Ui::MainWindow;
  ui->setupUi(this);

  // remove the Load Script menu entry, when the demo has not been compiled with QT_SCRIPT_LIB
#if !defined(QT_SCRIPT_LIB)
  ui->menuBar->removeAction(ui->actionLoad_Script);
#endif
  
  // Save some pointers from ui, for latter use.
  sceneView = ui->sceneView;
  viewer = ui->viewer;

  // do not save the state of the viewer (anoying)
  viewer->setStateFileName(QString::null);

  // setup scene
  scene = new Scene(this);
  viewer->setScene(scene);
  sceneView->setModel(scene);

  // setup the sceneview: delegation and columns sizing...
  sceneView->setItemDelegate(new SceneDelegate(this));

  sceneView->header()->setStretchLastSection(false);
  sceneView->header()->setResizeMode(Scene::NameColumn, QHeaderView::Stretch);
  sceneView->header()->setResizeMode(Scene::NameColumn, QHeaderView::Stretch);
  sceneView->header()->setResizeMode(Scene::ColorColumn, QHeaderView::ResizeToContents);
  sceneView->header()->setResizeMode(Scene::RenderingModeColumn, QHeaderView::Fixed);
  sceneView->header()->setResizeMode(Scene::ABColumn, QHeaderView::Fixed);
  sceneView->header()->setResizeMode(Scene::VisibleColumn, QHeaderView::Fixed);

  sceneView->resizeColumnToContents(Scene::ColorColumn);
  sceneView->resizeColumnToContents(Scene::RenderingModeColumn);
  sceneView->resizeColumnToContents(Scene::ABColumn);
  sceneView->resizeColumnToContents(Scene::VisibleColumn);

  // setup connections
  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          this, SLOT(updateInfo()));

  
  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          this, SLOT(updateDisplayInfo()));

  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          viewer, SLOT(updateGL()));

  connect(scene, SIGNAL(updated()),
          viewer, SLOT(update()));

  connect(scene, SIGNAL(updated()),
          this, SLOT(selectionChanged()));

  connect(scene, SIGNAL(itemAboutToBeDestroyed(Scene_item*)),
          this, SLOT(removeManipulatedFrame(Scene_item*)));

  connect(scene, SIGNAL(updated_bbox()),
          this, SLOT(updateViewerBBox()));

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
  connect(ui->infoLabel, SIGNAL(customContextMenuRequested(const QPoint & )),
          this, SLOT(showSceneContextMenu(const QPoint &)));

  connect(ui->actionRecenterScene, SIGNAL(triggered()),
          viewer->camera(), SLOT(interpolateToFitScene()));
  connect(ui->actionRecenterScene, SIGNAL(triggered()),
          viewer, SLOT(update()));

  connect(ui->actionAntiAliasing, SIGNAL(toggled(bool)),
          viewer, SLOT(setAntiAliasing(bool)));

  connect(ui->actionDraw_two_sides, SIGNAL(toggled(bool)),
          viewer, SLOT(setTwoSides(bool)));

  // enable anti-aliasing by default
  // ui->actionAntiAliasing->setChecked(true);

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
  connect(ui->actionSelect_all_items, SIGNAL(triggered()),
          this, SLOT(selectAll()));

  // Recent files menu
  this->addRecentFiles(ui->menuFile, ui->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(open(QString)));

  // Reset the "Operation menu"
  clearMenu(ui->menuOperations);

#ifdef QT_SCRIPT_LIB
  std::cerr << "Enable scripts.\n";
  script_engine = new QScriptEngine(this);
  qScriptRegisterMetaType<Scene_item*>(script_engine,
                                       myScene_itemToScriptValue,
                                       myScene_itemFromScriptValue);
#  ifdef QT_SCRIPTTOOLS_LIB
  QScriptEngineDebugger* debugger = new QScriptEngineDebugger(this);
  debugger->setObjectName("qt script debugger");
  QAction* debuggerMenuAction = 
    menuBar()->addMenu(debugger->createStandardMenu());
  debuggerMenuAction->setText(tr("Qt Script &debug"));
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

  // Load plugins, and re-enable actions that need it.
  loadPlugins();

  // Setup the submenu of the View menu that can toggle the dockwidgets
  Q_FOREACH(QDockWidget* widget, findChildren<QDockWidget*>()) {
    ui->menuDockWindows->addAction(widget->toggleViewAction());
  }
  ui->menuDockWindows->removeAction(ui->dummyAction);

  readSettings(); // Among other things, the column widths are stored.

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

void MainWindow::filterOperations()
{
  Q_FOREACH(const PluginNamePair& p, plugins) {
    if(p.first->applicable()) {
      Q_FOREACH(QAction* action, p.first->actions()) {
        action->setVisible(true);
      }
    } else {
      Q_FOREACH(QAction* action, p.first->actions()) {
        action->setVisible(false);
      }
    }
  }

  // do a pass over all menus in Operations and hide them when they are empty
  Q_FOREACH(QAction* action, ui->menuOperations->actions()) {
    if(QMenu* menu = action->menu()) {
      action->setVisible(!(menu->isEmpty()));
    }
  }
}


#ifdef QT_SCRIPT_LIB
void MainWindow::evaluate_script(QString script,
                                 const QString& /*filename*/,
                                 const bool quiet) {
  QScriptValue value = script_engine->evaluate(script);
  if(script_engine->hasUncaughtException()) {
    QTextStream err(stderr);
    err << "Qt Script exception:\n"
        << script_engine->uncaughtException().toString()
        << "\nBacktrace:\n";
    Q_FOREACH(QString line, script_engine->uncaughtExceptionBacktrace()) {
      err << "  " << line << "\n";
    }
  }
  else if(!quiet && !value.isNull() && !value.isUndefined()) {
    QTextStream(stderr) << "Qt Script evaluated to \""
                        << value.toString() << "\"\n";
  }
}

void MainWindow::evaluate_script_quiet(QString script,
                                       const QString& filename)
{
  evaluate_script(script, filename, true);
}
#endif

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
  }
  return;
#  endif
#endif
  // If we are here, then the debugger is not available
  this->error(tr("Your version of Qt is too old, and for that reason"
                 "the Qt Script Debugger is not available."));
}

namespace {
bool actionsByName(QAction* x, QAction* y) {
  return x->text() < y->text();
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
  plugins_directories << qApp->applicationDirPath();
  QString env_path = qgetenv("POLYHEDRON_DEMO_PLUGINS_PATH");
  if(!env_path.isEmpty()) {
    Q_FOREACH (QString pluginsDir, 
               env_path.split(":", QString::SkipEmptyParts)) {
      QDir dir(pluginsDir);
      if(dir.isReadable())
        plugins_directories << dir;
    }
  }
  Q_FOREACH (QDir pluginsDir, plugins_directories) {
    qDebug("# Looking for plugins in directory \"%s\"...",
           qPrintable(pluginsDir.absolutePath()));
    Q_FOREACH (QString fileName, pluginsDir.entryList(QDir::Files)) {
      if(fileName.contains("plugin") && QLibrary::isLibrary(fileName)) {
        qDebug("### Loading \"%s\"...", qPrintable(fileName));
        QPluginLoader loader;
        loader.setFileName(pluginsDir.absoluteFilePath(fileName));
        QObject *obj = loader.instance();
        if(obj) {
          QString name = fileName;
          name.remove(QRegExp("^lib"));
          name.remove(QRegExp("\\..*"));
          obj->setObjectName(name);
          initPlugin(obj);
          initIOPlugin(obj);
        }
        else {
          qDebug("Error loading \"%s\": %s",
                 qPrintable(fileName),
                 qPrintable(loader.errorString()));
        }
      }
    }
  }

  // sort the operations menu by name
  QList<QAction*> actions = ui->menuOperations->actions();
  qSort(actions.begin(), actions.end(), actionsByName);
  ui->menuOperations->clear();
  ui->menuOperations->addActions(actions);
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
  Polyhedron_demo_plugin_interface* plugin =
    qobject_cast<Polyhedron_demo_plugin_interface*>(obj);
  if(plugin) {
    // Call plugin's init() method
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
  Polyhedron_demo_io_plugin_interface* plugin =
    qobject_cast<Polyhedron_demo_io_plugin_interface*>(obj);
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

void MainWindow::viewerShow(float x, float y, float z) {
  viewer->camera()->setRevolveAroundPoint(qglviewer::Vec(x, y, z));
  // viewer->camera()->lookAt(qglviewer::Vec(x, y, z));

  qglviewer::ManipulatedCameraFrame backup_frame(*viewer->camera()->frame());
  viewer->camera()->fitSphere(qglviewer::Vec(x, y, z),
                              viewer->camera()->sceneRadius()/100);
  qglviewer::ManipulatedCameraFrame new_frame(*viewer->camera()->frame());
  *viewer->camera()->frame() = backup_frame;
  viewer->camera()->interpolateTo(new_frame, 1.f);
  viewer->setVisualHintsMask(1);
}

void MainWindow::message(QString message, QString colorName, QString font) {
  if (message.endsWith('\n')) {
    message.remove(message.length()-1, 1);
  }
  std::cerr << qPrintable(message) << std::endl;
  statusBar()->showMessage(message, 5000);
  message = "<font color=\"" + colorName + "\" style=\"font-style: " + font + ";\" >" +
    message + "</font><br>";
  message = "[" + QTime::currentTime().toString() + "] " + message;
  ui->consoleTextEdit->insertHtml(message);
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

void MainWindow::updateViewerBBox()
{
  const Scene::Bbox bbox = scene->bbox();
  const double xmin = bbox.xmin;
  const double ymin = bbox.ymin;
  const double zmin = bbox.zmin;
  const double xmax = bbox.xmax;
  const double ymax = bbox.ymax;
  const double zmax = bbox.zmax;
  // qDebug() << QString("Bounding box: (%1, %2, %3) - (%4, %5, %6)\n")
  // .arg(xmin).arg(ymin).arg(zmin).arg(xmax).arg(ymax).arg(zmax);
  qglviewer::Vec 
    vec_min(xmin, ymin, zmin),
    vec_max(xmax, ymax, zmax);
  viewer->setSceneBoundingBox(vec_min,
                              vec_max);
  viewer->camera()->showEntireScene();
}

void MainWindow::reload_item() {
  QAction* sender_action = qobject_cast<QAction*>(sender());
  if(!sender_action) return;
  
  bool ok;
  int item_index = sender_action->data().toInt(&ok);
  QObject* item_object = scene->item(item_index);
  if(!ok || !item_object || sender_action->data().type() != QVariant::Int) {
    std::cerr << "Cannot reload item: "
              << "the reload action has not item attached\n";
    return;
  }
  Scene_item* item = qobject_cast<Scene_item*>(item_object);
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
              << "the item has no \"source filename\" or no \"loader_name\" attached attached\n";
    return;
  }

  Polyhedron_demo_io_plugin_interface* fileloader = find_loader(loader_name);
  QFileInfo fileinfo(filename);

  Scene_item* new_item = load_item(fileinfo, fileloader);

  new_item->setName(item->name());
  new_item->setColor(item->color());
  new_item->setRenderingMode(item->renderingMode());
  new_item->setVisible(item->visible());
  new_item->changed();
  scene->replaceItem(item_index, new_item);
  delete item;
}

Polyhedron_demo_io_plugin_interface* MainWindow::find_loader(const QString& loader_name) const {
  Q_FOREACH(Polyhedron_demo_io_plugin_interface* io_plugin, 
            io_plugins) {
    if(io_plugin->name() == loader_name) {
      return io_plugin;
    }
  }
  throw std::invalid_argument(QString("No loader found with the name %1 available")
                              .arg(loader_name).toStdString()) ;
}

void MainWindow::open(QString filename)
{
  QFileInfo fileinfo(filename);
  QString filename_striped=fileinfo.fileName();

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
    load_script(fileinfo);
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
    evaluate_script(program, filename);
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

  //match all filters between ()
  QRegExp all_filters_rx("\\((.*)\\)");
  
  // collect all io_plugins and offer them to load if the file extension match one name filter
  // also collect all available plugin in case of a no extension match
  QStringList selected_items;
  QStringList all_items;
  Q_FOREACH(Polyhedron_demo_io_plugin_interface* io_plugin, io_plugins) {
    all_items << io_plugin->name();
    QStringList split_filters = io_plugin->nameFilters().split(";;");
    bool stop=false;
    Q_FOREACH(const QString& filter, split_filters) {
      //extract filters
      if ( all_filters_rx.indexIn(filter)!=-1 ){
        Q_FOREACH(const QString& pattern,all_filters_rx.cap(1).split(' ')){
          QRegExp rx(pattern);
          rx.setPatternSyntax(QRegExp::Wildcard);
          if ( rx.exactMatch(filename_striped) ){
            selected_items << io_plugin->name();
            stop=true;
            break;
          }
        }
        if (stop) break;
      }
    }
  }
  
  bool ok;
  QString loader_name;

  switch( selected_items.size() )
  {
    case 1:
      loader_name=selected_items.first();
      ok=true;
      break;
    case 0:
      loader_name=QInputDialog::getItem(this, tr("Select a loader"), tr("Available loaders for %1 :").arg(fileinfo.fileName()), all_items, 0, false, &ok);
      break;
    default:
      loader_name=QInputDialog::getItem(this, tr("Select a loader"), tr("Available loaders for %1 :").arg(fileinfo.fileName()), selected_items, 0, false, &ok);
  }
  
  if(!ok || loader_name.isEmpty()) { return; }
  
  Scene_item* scene_item = load_item(fileinfo, find_loader(loader_name));
  selectSceneItem(scene->addItem(scene_item));
}

bool MainWindow::open(QString filename, QString loader_name) {
  QFileInfo fileinfo(filename); 
  Scene_item* item;
  try {
    item = load_item(fileinfo, find_loader(loader_name));
  }
  catch(std::logic_error e) {
    std::cerr << e.what() << std::endl;
    return false;
  }
  selectSceneItem(scene->addItem(item));
  return true;
}


Scene_item* MainWindow::load_item(QFileInfo fileinfo, Polyhedron_demo_io_plugin_interface* loader) {
  Scene_item* item = NULL;
  if(!fileinfo.isFile() || !fileinfo.isReadable()) {
    throw std::invalid_argument(QString("File %1 is not a readable file.")
                                .arg(fileinfo.absoluteFilePath()).toStdString());
  }

  item = loader->load(fileinfo);
  if(!item) {
    throw std::logic_error(QString("Could not load item from file %1 using plugin %2")
                           .arg(fileinfo.absoluteFilePath()).arg(loader->name()).toStdString());
  }

  item->setProperty("source filename", fileinfo.absoluteFilePath());
  item->setProperty("loader_name", loader->name());
  return item;
}

void MainWindow::selectSceneItem(int i)
{
  if(i < 0 || i >= scene->numberOfEntries()) {
    sceneView->selectionModel()->clearSelection();
    updateInfo();
    updateDisplayInfo();
  }
  else {
    sceneView->selectionModel()->select(scene->createSelection(i),
                                       QItemSelectionModel::ClearAndSelect);
  }
}


void MainWindow::showSelectedPoint(double x, double y, double z)
{
  information(QString("Selected point: (%1, %2, %3)").
              arg(x, 0, 'g', 10).
              arg(y, 0, 'g', 10).
              arg(z, 0, 'g', 10));
}

void MainWindow::unSelectSceneItem(int i)
{
  removeSceneItemFromSelection(i);
}

void MainWindow::addSceneItemInSelection(int i)
{
  sceneView->selectionModel()->select(scene->createSelection(i),
                                     QItemSelectionModel::Select);
  scene->itemChanged(i);
}

void MainWindow::removeSceneItemFromSelection(int i)
{
  sceneView->selectionModel()->select(scene->createSelection(i),
                                     QItemSelectionModel::Deselect);
  scene->itemChanged(i);
}

void MainWindow::selectAll()
{
  sceneView->selectionModel()->select(scene->createSelectionAll(), 
                                     QItemSelectionModel::ClearAndSelect);
}

int MainWindow::getSelectedSceneItemIndex() const
{
  QModelIndexList selectedRows = sceneView->selectionModel()->selectedRows();
  if(selectedRows.size() != 1)
    return -1;
  else
    return selectedRows.first().row();
}

QList<int> MainWindow::getSelectedSceneItemIndices() const
{
  QModelIndexList selectedRows = sceneView->selectionModel()->selectedRows();
  QList<int> result;
  Q_FOREACH(QModelIndex index, selectedRows) {
    result << index.row();
  }
  return result;
}

void MainWindow::selectionChanged()
{
  scene->setSelectedItem(getSelectedSceneItemIndex());
  scene->setSelectedItemsList(getSelectedSceneItemIndices());
  Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(item != NULL && item->manipulatable()) {
    viewer->setManipulatedFrame(item->manipulatedFrame());
  } else {
    viewer->setManipulatedFrame(0);
  }
  if(viewer->manipulatedFrame() == 0) {
    Q_FOREACH(Scene_item* item, scene->entries()) {
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
  viewer->updateGL();
}

void MainWindow::contextMenuRequested(const QPoint& global_pos) {
  int index = scene->mainSelectionIndex();
  showSceneContextMenu(index, global_pos);
}

void MainWindow::showSceneContextMenu(int selectedItemIndex,
                                      const QPoint& global_pos)
{
  Scene_item* item = scene->item(selectedItemIndex);
  if(!item) return;

  const char* prop_name = "Menu modified by MainWindow.";

  QMenu* menu = item->contextMenu();
  if(menu && !item->property("source filename").toString().isEmpty()) {
    bool menuChanged = menu->property(prop_name).toBool();
    if(!menuChanged) {
      menu->addSeparator();
      QAction* reload = menu->addAction(tr("Reload item from file"));
      reload->setData(qVariantFromValue(selectedItemIndex));
      connect(reload, SIGNAL(triggered()),
              this, SLOT(reload_item()));
      menu->setProperty(prop_name, true);
    }
  }
  if(menu)
    menu->exec(global_pos);
}

void MainWindow::showSceneContextMenu(const QPoint& p) {
  QWidget* sender = qobject_cast<QWidget*>(this->sender());
  if(!sender) return;

  int index = -1;
  if(sender == sceneView) {
    QModelIndex modelIndex = sceneView->indexAt(p);
    if(!modelIndex.isValid()) return;

    index = modelIndex.row();
  }
  else {
    index = scene->mainSelectionIndex();
  }

  showSceneContextMenu(index, sender->mapToGlobal(p));
}

void MainWindow::removeManipulatedFrame(Scene_item* item)
{
  if(item->manipulatable() &&
     item->manipulatedFrame() == viewer->manipulatedFrame()) {
    viewer->setManipulatedFrame(0);
  }
}

void MainWindow::updateInfo() {
  Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(item) {
    QString item_text = item->toolTip();
    QString item_filename = item->property("source filename").toString();
    if(!item_filename.isEmpty()) {
      item_text += QString("<br /><i>File: %1").arg(item_filename);
    }
    ui->infoLabel->setText(item_text);
  }
  else 
    ui->infoLabel->clear();
}

void MainWindow::updateDisplayInfo() {
  Scene_item* item = scene->item(getSelectedSceneItemIndex());
  if(item)
    ui->displayLabel->setPixmap(item->graphicalToolTip());
  else 
    ui->displayLabel->clear();
}

void MainWindow::readSettings()
{
  this->readState("MainWindow", Size|State);
}

void MainWindow::writeSettings()
{
  this->writeState("MainWindow");
  std::cerr << "Write setting... done.\n";
}

void MainWindow::quit()
{
  close();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
  writeSettings();
  event->accept();
}

bool MainWindow::load_script(QString filename)
{
  QFileInfo fileinfo(filename);
  return load_script(fileinfo);
}

bool MainWindow::load_script(QFileInfo info)
{
#if defined(QT_SCRIPT_LIB)
  QString program;
  QString filename = info.absoluteFilePath();
  QFile script_file(filename);
  script_file.open(QIODevice::ReadOnly);
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

void MainWindow::on_actionLoad_Script_triggered() 
{
#if defined(QT_SCRIPT_LIB)
  QString filename = QFileDialog::getOpenFileName(
    this,
    tr("Select a script to run..."),
    ".",
    "QTScripts (*.js);;All Files (*)");

  load_script(QFileInfo(filename));
#endif
}

void MainWindow::on_actionLoad_triggered()
{
  QStringList filters;
  // we need to special case our way out of this
  filters << "All Files (*)";

  QStringList extensions;

  typedef QMap<QString, Polyhedron_demo_io_plugin_interface*> FilterPluginMap;
  FilterPluginMap filterPluginMap;
  
  Q_FOREACH(Polyhedron_demo_io_plugin_interface* plugin, io_plugins) {
    QStringList split_filters = plugin->nameFilters().split(";;");
    Q_FOREACH(const QString& filter, split_filters) {
      FilterPluginMap::iterator it = filterPluginMap.find(filter);
      if(it != filterPluginMap.end()) {
        qDebug() << "Duplicate Filter: " << it.value();
        qDebug() << "This filter will not be available.";
      } else {
        filterPluginMap[filter] = plugin;
      }
      filters << filter;
    }
  }


  QFileDialog dialog(this);
  dialog.setNameFilters(filters);
  if(dialog.exec() != QDialog::Accepted) { return; }
  
  FilterPluginMap::iterator it = 
    filterPluginMap.find(dialog.selectedNameFilter());
  
  Polyhedron_demo_io_plugin_interface* selectedPlugin = NULL;

  if(it != filterPluginMap.end()) {
    selectedPlugin = it.value();
  }

  Q_FOREACH(const QString& filename, dialog.selectedFiles()) {
    Scene_item* item = NULL;
    if(selectedPlugin) {
      QFileInfo info(filename);
      item = load_item(info, selectedPlugin);
      Scene::Item_id index = scene->addItem(item);
      selectSceneItem(index);
      this->addToRecentFiles(filename);
    } else {
      open(filename);
    }
  }
}

void MainWindow::on_actionSaveAs_triggered()
{
  QModelIndexList selectedRows = sceneView->selectionModel()->selectedRows();
  if(selectedRows.size() != 1)
    return;
  Scene_item* item = scene->item(getSelectedSceneItemIndex());

  if(!item)
    return;

  QVector<Polyhedron_demo_io_plugin_interface*> canSavePlugins;
  QStringList filters;
  Q_FOREACH(Polyhedron_demo_io_plugin_interface* plugin, io_plugins) {
    if(plugin->canSave(item)) {
      canSavePlugins << plugin;
      filters += plugin->nameFilters();
    }
  }
  filters << tr("All files (*)");

  if(canSavePlugins.isEmpty()) {
    QMessageBox::warning(this,
                         tr("Cannot save"),
                         tr("The selected object %1 cannot be saved.")
                         .arg(item->name()));
    return;
  }

  QString filename = 
    QFileDialog::getSaveFileName(this,
                                 tr("Save to File..."),
                                 QString(),
                                 filters.join(";;"));
  save(filename, item);
}

void MainWindow::save(QString filename, Scene_item* item) {
  QFileInfo fileinfo(filename);

  Q_FOREACH(Polyhedron_demo_io_plugin_interface* plugin, io_plugins) {
    if(plugin->canSave(item)) {
      if(plugin->save(item, fileinfo))
        break;
    }
  }
}

bool MainWindow::on_actionErase_triggered()
{
  int next_index = scene->erase(scene->selectionIndices());
  selectSceneItem(next_index);
  return next_index >= 0;
}

void MainWindow::on_actionEraseAll_triggered()
{
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
    int i = index.row();
    Scene_item* item = scene->item(i);
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

void MainWindow::on_actionSetBackgroundColor_triggered()
{
  QColor c =  QColorDialog::getColor();
  if(c.isValid()) {
    viewer->setBackgroundColor(c);
  }
}

void MainWindow::on_action_Look_at_triggered()
{
  Show_point_dialog dialog(this);
  int i = dialog.exec();
  if( i == QDialog::Accepted &&
      dialog.has_correct_coordinates() )
  {
    viewerShow((float)dialog.get_x(),
               (float)dialog.get_y(),
               (float)dialog.get_z());
  }
}

QString MainWindow::camera_string() const
{
  return viewer->dumpCameraCoordinates();
}

void MainWindow::on_actionDumpCamera_triggered()
{
  information(QString("Camera: %1")
              .arg(camera_string()));
}

void MainWindow::on_action_Copy_camera_triggered()
{
  qApp->clipboard()->setText(this->camera_string());
}

void MainWindow::on_action_Paste_camera_triggered()
{
  QString s = qApp->clipboard()->text();
  viewer->moveCameraToCoordinates(s, 0.5f);
}

void MainWindow::setAddKeyFrameKeyboardModifiers(::Qt::KeyboardModifiers m)
{
  viewer->setAddKeyFrameKeyboardModifiers(m);
}
