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
#include <QStandardItemModel>
#include <QStandardItem>

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
#include "ui_Preferences.h"

#include "Show_point_dialog.h"
#include "File_loader_dialog.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <QGLViewer/manipulatedFrame.h>

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

  {
    QShortcut* shortcut = new QShortcut(QKeySequence(Qt::ALT+Qt::Key_Q), this);
    connect(shortcut, SIGNAL(activated()),
            this, SLOT(setFocusToQuickSearch()));
  }

  proxyModel = new QSortFilterProxyModel(this);
  proxyModel->setSourceModel(scene);

  connect(ui->searchEdit, SIGNAL(textChanged(QString)),
          proxyModel, SLOT(setFilterFixedString(QString)));
  sceneView->setModel(proxyModel);

  // setup the sceneview: delegation and columns sizing...
  sceneView->setItemDelegate(new SceneDelegate(this));

  sceneView->header()->setStretchLastSection(false);
  sceneView->header()->setSectionResizeMode(Scene::NameColumn, QHeaderView::Stretch);
  sceneView->header()->setSectionResizeMode(Scene::NameColumn, QHeaderView::Stretch);
  sceneView->header()->setSectionResizeMode(Scene::ColorColumn, QHeaderView::ResizeToContents);
  sceneView->header()->setSectionResizeMode(Scene::RenderingModeColumn, QHeaderView::Fixed);
  sceneView->header()->setSectionResizeMode(Scene::ABColumn, QHeaderView::Fixed);
  sceneView->header()->setSectionResizeMode(Scene::VisibleColumn, QHeaderView::Fixed);
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
          viewer, SLOT(updateGL()));

  connect(scene, SIGNAL(updated()),
          this, SLOT(selectionChanged()));

  connect(scene, SIGNAL(itemAboutToBeDestroyed(Scene_item*)),
          this, SLOT(removeManipulatedFrame(Scene_item*)));

  connect(scene, SIGNAL(updated_bbox()),
          this, SLOT(updateViewerBBox()));

  connect(scene, SIGNAL(selectionChanged(int)),
          this, SLOT(selectSceneItem(int)));

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

  // The contextMenuPolicy of infoLabel is now the default one, so that one
  // can easily copy-paste its text.
  // connect(ui->infoLabel, SIGNAL(customContextMenuRequested(const QPoint & )),
  //         this, SLOT(showSceneContextMenu(const QPoint &)));

  connect(ui->actionRecenterScene, SIGNAL(triggered()),
          viewer, SLOT(update()));

  connect(ui->actionAntiAliasing, SIGNAL(toggled(bool)),
          viewer, SLOT(setAntiAliasing(bool)));

  connect(ui->actionDraw_two_sides, SIGNAL(toggled(bool)),
          viewer, SLOT(setTwoSides(bool)));

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

  readSettings(); // Among other things, the column widths are stored.

  // Load plugins, and re-enable actions that need it.
  loadPlugins();

  // Setup the submenu of the View menu that can toggle the dockwidgets
  Q_FOREACH(QDockWidget* widget, findChildren<QDockWidget*>()) {
    ui->menuDockWindows->addAction(widget->toggleViewAction());
  }
  ui->menuDockWindows->removeAction(ui->dummyAction);


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
    Q_FOREACH(QAction* action, p.first->actions()) {
        action->setVisible( p.first->applicable(action) );
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
                                 const QString& filename,
                                 const bool quiet) {
  QScriptValue value = script_engine->evaluate(script, filename);
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
  this->error(tr("Your version of Qt is too old, and for that reason "
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
        //set plugin name
        QString name = fileName;
        name.remove(QRegExp("^lib"));
        name.remove(QRegExp("\\..*"));
        //do not load it if it is in the blacklist
        if ( plugin_blacklist.contains(name) ){
          qDebug("### Ignoring plugin \"%s\".", qPrintable(fileName));
          continue;
        }
        QDebug qdebug = qDebug();
        qdebug << "### Loading \"" << fileName.toUtf8().data() << "\"... ";
        QPluginLoader loader;
        loader.setFileName(pluginsDir.absoluteFilePath(fileName));
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
#if QGLVIEWER_VERSION >= 0x020502
  viewer->camera()->setPivotPoint(qglviewer::Vec(x, y, z));
#else
  viewer->camera()->setRevolveAroundPoint(qglviewer::Vec(x, y, z));
#endif
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
  new_item->invalidate_buffers();
  scene->replaceItem(item_index, new_item, true);
  item->deleteLater();
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


  QStringList selected_items;
  QStringList all_items;

  QMap<QString,QString>::iterator dfs_it = 
    default_plugin_selection.find( fileinfo.completeSuffix() );
  
  if ( dfs_it==default_plugin_selection.end() )
  {
    // collect all io_plugins and offer them to load if the file extension match one name filter
    // also collect all available plugin in case of a no extension match
    Q_FOREACH(Polyhedron_demo_io_plugin_interface* io_plugin, io_plugins) {
      if ( !io_plugin->canLoad() ) continue;
      all_items << io_plugin->name();
      if ( file_matches_filter(io_plugin->nameFilters(), filename) )
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
  
  if(!ok || load_pair.first.isEmpty()) { return; }
  
  if (load_pair.second)
     default_plugin_selection[fileinfo.completeSuffix()]=load_pair.first;
  
  
  QSettings settings;
  settings.setValue("OFF open directory",
                    fileinfo.absoluteDir().absolutePath());

  Scene_item* scene_item = load_item(fileinfo, find_loader(load_pair.first));
  if(scene_item != 0) {
    this->addToRecentFiles(fileinfo.absoluteFilePath());
  }
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
  QModelIndexList selectedRows = sceneView->selectionModel()->selectedRows();
  if(selectedRows.size() != 1)
    return -1;
  else {
    QModelIndex i = proxyModel->mapToSource(selectedRows.first());
    return i.row();
  }
}

QList<int> MainWindow::getSelectedSceneItemIndices() const
{
  QModelIndexList selectedRows = sceneView->selectionModel()->selectedRows();
  QList<int> result;
  Q_FOREACH(QModelIndex index, selectedRows) {
    result << proxyModel->mapToSource(index).row();
  }
  return result;
}

void MainWindow::selectionChanged()
{
  scene->setSelectedItemIndex(getSelectedSceneItemIndex());
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
  if(menu) {
    bool menuChanged = menu->property(prop_name).toBool();
    if(!menuChanged) {
      menu->addSeparator();
      if(!item->property("source filename").toString().isEmpty()) {
        QAction* reload = menu->addAction(tr("&Reload item from file"));
        reload->setData(qVariantFromValue(selectedItemIndex));
        connect(reload, SIGNAL(triggered()),
                this, SLOT(reload_item()));
      }
      QAction* saveas = menu->addAction(tr("&Save as..."));
      saveas->setData(qVariantFromValue(selectedItemIndex));
      connect(saveas,  SIGNAL(triggered()),
              this, SLOT(on_actionSaveAs_triggered()));
      QAction* showobject = menu->addAction(tr("&Zoom to this object"));
      showobject->setData(qVariantFromValue(selectedItemIndex));
      connect(showobject, SIGNAL(triggered()),
              this, SLOT(viewerShowObject()));

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

    index = proxyModel->mapToSource(modelIndex).row();
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
  {
    QSettings settings;
    // enable anti-aliasing 
    ui->actionAntiAliasing->setChecked(settings.value("antialiasing", false).toBool());
    // read plugin blacklist
    QStringList blacklist=settings.value("plugin_blacklist",QStringList()).toStringList();
    Q_FOREACH(QString name,blacklist){ plugin_blacklist.insert(name); }
  }
  this->readState("MainWindow", Size|State);
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
  }
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

  QSettings settings;
  QString directory = settings.value("OFF open directory",
                                     QDir::current().dirName()).toString();

  QFileDialog dialog(this);
  dialog.setDirectory(directory);
  dialog.setNameFilters(filters);
  dialog.setFileMode(QFileDialog::ExistingFiles);

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
  int index = -1;
  QAction* sender_action = qobject_cast<QAction*>(sender());
  if(sender_action && !sender_action->data().isNull()) {
    index = sender_action->data().toInt();
  }

  if(index < 0) {
    QModelIndexList selectedRows = sceneView->selectionModel()->selectedRows();
    if(selectedRows.size() != 1)
      return;
    index = getSelectedSceneItemIndex();
  }
  Scene_item* item = scene->item(index);

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

  QString caption = tr("Save %1 to File...").arg(item->name());
  QString filename = 
    QFileDialog::getSaveFileName(this,
                                 caption,
                                 QString(),
                                 filters.join(";;"));
  save(filename, item);
}

void MainWindow::save(QString filename, Scene_item* item) {
  QFileInfo fileinfo(filename);

  Q_FOREACH(Polyhedron_demo_io_plugin_interface* plugin, io_plugins) {
    if(  plugin->canSave(item) &&
        file_matches_filter(plugin->nameFilters(),filename) )
    {
      if(plugin->save(item, fileinfo))
        break;
    }
  }
}

void MainWindow::on_actionSaveSnapshot_triggered()
{
  viewer->saveSnapshot(false);
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
    int i = proxyModel->mapToSource(index).row();
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
void MainWindow::on_actionPreferences_triggered()
{
  QDialog dialog(this);
  Ui::PreferencesDialog prefdiag;
  prefdiag.setupUi(&dialog);
  
  
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
  Q_FOREACH(Polyhedron_demo_io_plugin_interface* plugin, io_plugins)
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

void MainWindow::viewerShowObject()
{
  int index = -1;
  QAction* sender_action = qobject_cast<QAction*>(sender());
  if(sender_action && !sender_action->data().isNull()) {
    index = sender_action->data().toInt();
  }
  if(index >= 0) {
    const Scene::Bbox bbox = scene->item(index)->bbox();
    viewerShow((float)bbox.xmin, (float)bbox.ymin, (float)bbox.zmin,
               (float)bbox.xmax, (float)bbox.ymax, (float)bbox.zmax);
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

void MainWindow::on_actionRecenterScene_triggered()
{
  updateViewerBBox();
  viewer->camera()->interpolateToFitScene();
}
