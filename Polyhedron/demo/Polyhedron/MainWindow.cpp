#ifdef CGAL_USE_SSH
#  include "CGAL/Use_ssh.h"
#endif
#include <cmath>

#include "config.h"
#include "MainWindow.h"
#include "Scene.h"
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/TextRenderer.h>
#include <CGAL/Three/exceptions.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Qt/debug.h>
#include <CGAL/double.h>

#include <QJsonArray>
#include <QtDebug>
#include <QFileDialog>
#include <QFileInfo>
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
#include <QStandardItemModel>
#include <QStandardItem>
#include <QTreeWidgetItem>
#include <QTreeWidget>
#include <QDockWidget>
#include <QSpinBox>
#include <stdexcept>
#include <fstream>
#include <QElapsedTimer>
#include <QWidgetAction>
#include <QJsonArray>
#include <QSequentialIterable>
#include <QDir>
#ifdef QT_SCRIPT_LIB
#  include <QScriptValue>
#  ifdef QT_SCRIPTTOOLS_LIB
#    include <QScriptEngineDebugger>
#  endif
#endif

#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Scene_item_with_properties.h>
#include "ui_SubViewer.h"
#include "ui_MainWindow.h"
#include "ui_Preferences.h"
#include "ui_Details.h"
#include "ui_Statistics_on_item_dialog.h"
#include "ui_SSH_dialog.h"
#include "Show_point_dialog.h"
#include "File_loader_dialog.h"
#include "Viewer.h"

#include <CGAL/Qt/manipulatedCameraFrame.h>
#include <CGAL/Qt/manipulatedFrame.h>

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
  searchAction->deleteLater();
  delete ui;
  delete statistics_ui;
}
MainWindow::MainWindow(const QStringList &keywords, bool verbose, QWidget* parent)
  : CGAL::Qt::DemosMainWindow(parent),
    accepted_keywords(keywords)
{
  bbox_need_update = true;
  ui = new Ui::MainWindow;
  ui->setupUi(this);
  menuBar()->setNativeMenuBar(false);
  searchAction = new QWidgetAction(0);
  CGAL::Three::Three::s_mainwindow = this;
  menu_map[ui->menuOperations->title()] = ui->menuOperations;
  this->verbose = verbose;
  // remove the Load Script menu entry, when the demo has not been compiled with QT_SCRIPT_LIB
#if !defined(QT_SCRIPT_LIB)
  ui->menuBar->removeAction(ui->actionLoadScript);
  ui->menuBar->removeAction(ui->on_actionLoad_a_Scene_from_a_Script_File);
#endif
  // Save some pointers from ui, for latter use.
  sceneView = ui->sceneView;
  viewer_window = new SubViewer(ui->mdiArea, this, nullptr);
  viewer = viewer_window->viewer;
  CGAL::Three::Three::s_mainviewer = viewer;
  viewer->setObjectName("mainViewer");
  viewer_window->showMaximized();
  viewer_window->setWindowFlags(
        Qt::SubWindow
        | Qt::CustomizeWindowHint
        | Qt::WindowMaximizeButtonHint
        | Qt::WindowSystemMenuHint
        | Qt::WindowTitleHint
        );
  viewer_window->setWindowTitle("Main Viewer");
  // setup scene
  scene = new Scene(this);
  CGAL::Three::Three::s_scene = scene;
  CGAL::Three::Three::s_connectable_scene = scene;
  {
    QShortcut* shortcut = new QShortcut(QKeySequence(Qt::ALT+Qt::Key_Q), this);
    connect(shortcut, SIGNAL(activated()),
            this, SLOT(setFocusToQuickSearch()));
    shortcut = new QShortcut(QKeySequence(Qt::Key_F5), this);
    connect(shortcut, SIGNAL(activated()),
            this, SLOT(reloadItem()));
    shortcut = new QShortcut(QKeySequence(Qt::Key_F11), this);
    connect(shortcut, SIGNAL(activated()),
            this, SLOT(toggleFullScreen()));
    shortcut = new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_R), this);
    connect(shortcut, &QShortcut::activated,
            this, &MainWindow::recenterScene);
    shortcut = new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_T), this);
    connect(shortcut, &QShortcut::activated,
            this,
            [](){
      Viewer* viewer = qobject_cast<Viewer*>(CGAL::Three::Three::activeViewer());
      bool b = viewer->property("draw_two_sides").toBool();
      viewer->setTwoSides(!b);
    }
            );
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

  connect(viewer, &Viewer::needNewContext,
          [this](){create();});


  connect(scene, SIGNAL(updated()),
          this, SLOT(selectionChanged()));

  connect(scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)),
          this, SLOT(removeManipulatedFrame(CGAL::Three::Scene_item*)));

  connect(scene, SIGNAL(updated_bbox(bool)),
          this, SLOT(invalidate_bbox(bool)));

  connect(scene, SIGNAL(selectionChanged(int)),
          this, SLOT(selectSceneItem(int)));
  connect(scene, SIGNAL(selectionChanged(QList<int>)),
          this, SLOT(selectSceneItems(QList<int>)));

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
  // setup menu filtering

  connect(sceneView->selectionModel(),
      QOverload<const QItemSelection & , const QItemSelection &>::of(&QItemSelectionModel::selectionChanged),
      this, [this](){filterOperations(false);});

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

  setupViewer(viewer, viewer_window);

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
    this->QMainWindow::addDockWidget(Qt::BottomDockWidgetArea, dock);
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
  operationSearchBar.setPlaceholderText("Filter...");
  searchAction->setDefaultWidget(&operationSearchBar);

  connect(&operationSearchBar, &QLineEdit::textChanged,
          this, [this](){filterOperations(true);});

  loadPlugins();
  accepted_keywords.clear();

  // Setup the submenu of the View menu that can toggle the dockwidgets
  Q_FOREACH(QDockWidget* widget, findChildren<QDockWidget*>()) {
    widget->setFeatures(QDockWidget::DockWidgetClosable | QDockWidget::DockWidgetFloatable);
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
  filterOperations(true);
  // debugger->action(QScriptEngineDebugger::InterruptAction)->trigger();
#endif
}

void addActionToMenu(QAction* action, QMenu* menu)
{
  auto actions = menu->actions();
  auto it = std::lower_bound(actions.begin(), actions.end(),
                             action->text().remove("&"),
                             [](QAction* a, QString text) {
                               return a->text().remove("&").compare(text) < 0;
                             });
  if(it == actions.end()) {
    menu->addAction(action);
  }
  else {
    menu->insertAction(*it, action);
  }
}

//Recursive function that do a pass over a menu and its sub-menus(etc.) and hide them when they are empty
void filterMenuOperations(QMenu* menu, QString filter, bool keep_from_here)
{
  QList<QAction*> buffer;
  Q_FOREACH(QAction* action, menu->actions())
    buffer.append(action);

  while(!buffer.isEmpty()){
    Q_FOREACH(QAction* action, buffer) {
      if(QMenu* submenu = action->menu())
      {
        bool keep = true;
        if(!keep_from_here){
          keep = submenu->menuAction()->text().contains(filter, Qt::CaseInsensitive);
          if(!keep)
          {
            Q_FOREACH(QAction* subaction, submenu->actions())
            {
              submenu->removeAction(subaction);
              buffer.append(subaction);
            }
          }
          else
          {
            addActionToMenu(submenu->menuAction(), menu);
          }
        }
        filterMenuOperations(submenu, filter, keep);
        action->setVisible(!(submenu->isEmpty()));

      }
      else if(action->text().remove("&").contains(filter, Qt::CaseInsensitive)){
        //menu->addAction(action);
        addActionToMenu(action, menu);
      }
      buffer.removeAll(action);
    }
  }
}

#ifdef Q_OS_WIN
void MainWindow::filterOperations(bool hide)
#else
void MainWindow::filterOperations(bool)
#endif
{
  //on some platforms editing an open menu slows everything like hell,
  //so we hide it for the time of the process.
#ifdef Q_OS_WIN
  if(hide)
    ui->menuOperations->hide();
#endif
  //return actions to their true menu
  Q_FOREACH(QMenu* menu, action_menu_map.values())
  {
    Q_FOREACH(QAction* action, menu->actions())
    {
      if(action != searchAction)
        menu->removeAction(action);
    }
  }

  Q_FOREACH(QAction* action, action_menu_map.keys())
  {
    QMenu* menu = action_menu_map[action];
    addActionToMenu(action, menu);
  }
  QString filter=operationSearchBar.text();
  Q_FOREACH(const PluginNamePair& p, plugins) {
    Q_FOREACH(QAction* action, p.first->actions()) {
      action->setVisible( p.first->applicable(action)
                          && (action->text().remove("&").contains(filter, Qt::CaseInsensitive)
                              || action->property("subMenuName")
                              .toString().contains(filter, Qt::CaseInsensitive)));
    }
  }
  // do a pass over all menus in Operations and their sub-menus(etc.) and hide them when they are empty
  filterMenuOperations(ui->menuOperations, filter, false);
#ifdef Q_OS_WIN
  if(hide)
    ui->menuOperations->show();
#endif
  operationSearchBar.setFocus();
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
  action_menu_map[menu_map[menuName]->menuAction()] = menu_map[parentName];

  // only add the action in the last submenu
  if(slash_index==-1)
  {
    ui->menuOperations->removeAction(a);
    menu_map[menuName]->addAction(a);
    action_menu_map[a] = menu_map[menuName];
  }
}

bool MainWindow::load_plugin(QString fileName, bool blacklisted)
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
        pluginsStatus_map[name] = QString("Blacklisted.");
        ignored_map[name] = true;
        //qDebug("### Ignoring plugin \"%s\".", qPrintable(fileName));
        PathNames_map[name].push_back(fileinfo.absoluteDir().absolutePath());
        return true;
      }
    }
    QDebug qdebug = qDebug();
    if(verbose)
      qdebug << "### Loading \"" << fileName.toUtf8().data() << "\"... ";
    QPluginLoader loader;
    loader.setFileName(fileinfo.absoluteFilePath());
    QJsonArray keywords = loader.metaData().value("MetaData").toObject().value("Keywords").toArray();
    QString date = loader.metaData().value("MetaData").toObject().value("ConfigDate").toString();
    QStringList s_keywords;
    for(int i = 0; i < keywords.size(); ++i)
    {
      s_keywords.append(keywords[i].toString());
    }
    plugin_metadata_map[name] = qMakePair(s_keywords, date);
    QObject *obj = loader.instance();
    bool do_load = accepted_keywords.empty();
    if(!do_load)
    {
      Q_FOREACH(QString k, s_keywords)
      {
        if(accepted_keywords.contains(k))
        {
          do_load = true;
          break;
        }
      }
    }
    if(do_load && obj) {
      obj->setObjectName(name);
      bool init1 = initPlugin(obj);
      bool init2 = initIOPlugin(obj);
      if (!init1 && !init2)
      {
        //qdebug << "not for this program";
        pluginsStatus_map[name] = QString("Not for this program.");
      }
      else
        //qdebug << "success";
        pluginsStatus_map[name] = QString("success");
    }
    else if(!do_load)
    {
      pluginsStatus_map[name]="Wrong Keywords.";
      ignored_map[name] = true;
    }
    else{
      //qdebug << "error: " << qPrintable(loader.errorString());
      pluginsStatus_map[name] = loader.errorString();

    }
    PathNames_map[name].push_back(fileinfo.absoluteDir().absolutePath());
    return true;
  }
  return false;
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
#if (QT_VERSION >= QT_VERSION_CHECK(5, 6, 0))
  QChar separator = QDir::listSeparator();
#else
#if defined(_WIN32)
  QChar separator = ';';
#else
  QChar separator = ':';
#endif
#endif
  if(!env_path.isEmpty()) {
#if defined(_WIN32)
    QString path = qgetenv("PATH");
    QByteArray new_path = path.append(env_path.prepend(separator)).toUtf8();
    qputenv("PATH", new_path);
#endif
    Q_FOREACH (QString pluginsDir,
               env_path.split(separator, CGAL_QT_SKIP_EMPTY_PARTS)) {
      QDir dir(pluginsDir);
      if(dir.isReadable())
        plugins_directories << dir;
    }
  }

  QSet<QString> loaded;
  Q_FOREACH (QDir pluginsDir, plugins_directories) {
    if(verbose)
      qDebug("# Looking for plugins in directory \"%s\"...",
             qPrintable(pluginsDir.absolutePath()));
    Q_FOREACH(QString fileName, pluginsDir.entryList(QDir::Files))
    {
      QString abs_name = pluginsDir.absoluteFilePath(fileName);
      if(loaded.find(abs_name) == loaded.end())
      {
        if(load_plugin(abs_name, true))
        {
          loaded.insert(abs_name);
        }
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
  std::sort(as.begin(), as.end(), actionsByName);
  ui->menuOperations->clear();
  ui->menuOperations->addAction(searchAction);
  ui->menuOperations->addActions(as);
  operationSearchBar.setFocus();
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
        action_menu_map[action] = ui->menuOperations;
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
    plugin->init();
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
  CGAL::qglviewer::Vec
      min_(xmin, ymin, zmin),
      max_(xmax, ymax, zmax);

  if(min_ == max_) return viewerShow(viewer, xmin, ymin, zmin);

  viewer->camera()->setPivotPoint((min_+max_)*0.5);

  CGAL::qglviewer::ManipulatedCameraFrame backup_frame(*viewer->camera()->frame());
  viewer->camera()->fitBoundingBox(min_, max_);
  CGAL::qglviewer::ManipulatedCameraFrame new_frame(*viewer->camera()->frame());
  *viewer->camera()->frame() = backup_frame;
  viewer->camera()->interpolateTo(new_frame, 1.f);
  viewer->setVisualHintsMask(1);
}

void MainWindow::viewerShow(Viewer_interface* vi, float x, float y, float z) {

  CGAL::qglviewer::ManipulatedCameraFrame backup_frame(*vi->camera()->frame());
  vi->camera()->fitSphere(CGAL::qglviewer::Vec(x, y, z),
                          vi->camera()->sceneRadius()/100);
  CGAL::qglviewer::ManipulatedCameraFrame new_frame(*vi->camera()->frame());
  *vi->camera()->frame() = backup_frame;
  vi->camera()->interpolateTo(new_frame, 1.f);
  vi->setVisualHintsMask(1);

  vi->camera()->setPivotPoint(CGAL::qglviewer::Vec(x, y, z));
}

void MainWindow::message(QString message, QString colorName, QString font) {
  if (message.endsWith('\n')) {
    message.remove(message.length()-1, 1);
  }
  statusBar()->showMessage(message, 5000);
  QTimer::singleShot(5000, [this]{this->statusBar()->setStyleSheet("");});
  message = "<font color=\"" + colorName + "\" style=\"font-style: " + font + ";\" >" +
      message + "</font><br>";
  message = "[" + QTime::currentTime().toString() + "] " + message;
  ui->consoleTextEdit->append(message);
  ui->consoleTextEdit->verticalScrollBar()->setValue(ui->consoleTextEdit->verticalScrollBar()->maximum());
}

void MainWindow::message_information(QString text) {
  statusBar()->setStyleSheet("color: blue");
  this->message("INFO: " + text, "blue");
}

void MainWindow::message_warning(QString text) {
  statusBar()->setStyleSheet("color: orange");
  this->message("WARNING: " + text, "orange");
}

void MainWindow::message_error(QString text) {
  statusBar()->setStyleSheet("color: red");
  this->message("ERROR: " + text, "red");
}

void MainWindow::updateViewersBboxes(bool recenter)
{
  if(bbox_need_update)
  {
  CGAL::qglviewer::Vec min, max;
  computeViewerBBox(min, max);
  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
  {
    if(v == NULL)
      continue;
    Viewer* vi = static_cast<Viewer*>(v);
    updateViewerBbox(vi, recenter, min, max);
  }
  bbox_need_update = false;
}

}

void MainWindow::computeViewerBBox(CGAL::qglviewer::Vec& vmin, CGAL::qglviewer::Vec& vmax)
{
  const Scene::Bbox bbox = scene->bbox();
  const double xmin = bbox.xmin();
  const double ymin = bbox.ymin();
  const double zmin = bbox.zmin();
  const double xmax = bbox.xmax();
  const double ymax = bbox.ymax();
  const double zmax = bbox.zmax();



  vmin = CGAL::qglviewer::Vec(xmin, ymin, zmin);
  vmax= CGAL::qglviewer::Vec(xmax, ymax, zmax);

  CGAL::qglviewer::Vec bbox_center((xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2);

  double bbox_diag = CGAL::approximate_sqrt(
        CGAL::square(xmax - xmin)
        + CGAL::square(ymax - ymin)
        + CGAL::square(zmax - zmin));

  CGAL::qglviewer::Vec offset(0,0,0);

  double l_dist = (std::max)((std::abs)(bbox_center.x - viewer->offset().x),
                             (std::max)((std::abs)(bbox_center.y - viewer->offset().y),
                                        (std::abs)(bbox_center.z - viewer->offset().z)));
  if((std::log2)(l_dist/bbox_diag) > 11.0 )
    for(int i=0; i<3; ++i)
    {
      offset[i] = -bbox_center[i];
    }
  if(offset != viewer->offset())
  {
    Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
    {
      if(v == NULL)
        continue;
      Viewer* vi = qobject_cast<Viewer*>(v);
      vi->setOffset(offset);
    }
    for(int i=0; i<scene->numberOfEntries(); ++i)
    {
      //      scene->item(i)->invalidate(Scene_item::GEOMETRY);
      scene->item(i)->invalidateOpenGLBuffers();
      scene->item(i)->itemChanged();
    }
  }
}

void MainWindow::reloadItem() {

  Scene_item* item = NULL;

  Q_FOREACH(Scene::Item_id id, scene->selectionIndices())
  {
    item = scene->item(id);
    if(!item)//secure items like selection items that get deleted when their "parent" item is reloaded.
      continue;
    QString filename = item->property("source filename").toString();
    QString loader_name = item->property("loader_name").toString();
    if(filename.isEmpty() || loader_name.isEmpty()) {
      this->warning(QString("Cannot reload item %1: "
                            "the item has no \"source filename\" or no \"loader_name\" attached\n").arg(item->name()));
      continue;
    }

    CGAL::Three::Polyhedron_demo_io_plugin_interface* fileloader = findLoader(loader_name);
    QFileInfo fileinfo(filename);
    bool ok;
    QList<Scene_item*> new_items = loadItem(fileinfo, fileloader, ok, false);
    if(!ok)
      return;
    QVariant varian = item->property("load_mates");
    if(!varian.isValid()) //typically when a soup is oriented, the soup_item is deleted and thus the varain points to an unexisting item.
    {
      Scene_item* new_item = new_items.front();
      new_item->setName(item->name());
      new_item->setColor(item->color());
      new_item->setRenderingMode(item->renderingMode());
      new_item->setVisible(item->visible());
      Scene_item_with_properties *property_item = dynamic_cast<Scene_item_with_properties*>(new_item);
      scene->replaceItem(scene->item_id(item), new_item, true);
      if(property_item)
        property_item->copyProperties(item);
      new_item->invalidateOpenGLBuffers();
      item->deleteLater();
      return;
    }
    QSequentialIterable iterable = varian.value<QSequentialIterable>();

       // Can use foreach:
    int mate_id = 0;
    Q_FOREACH(const QVariant &v, iterable)
    {
      Scene_item* mate = v.value<Scene_item*>();
      Scene_item* new_item = new_items[mate_id];
      new_item->setName(mate->name());
      new_item->setColor(mate->color());
      new_item->setRenderingMode(mate->renderingMode());
      new_item->setVisible(mate->visible());
      Scene_item_with_properties *property_item = dynamic_cast<Scene_item_with_properties*>(new_item);
      scene->replaceItem(scene->item_id(mate), new_item, true);
      if(property_item)
        property_item->copyProperties(mate);
      new_item->invalidateOpenGLBuffers();
      mate->deleteLater();
    }
  }
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
    for(CGAL::Three::Polyhedron_demo_io_plugin_interface* io_plugin : io_plugins) {
      if ( file_matches_filter(io_plugin->loadNameFilters(), filename.toLower()) )
      {
        if ( !io_plugin->canLoad(fileinfo) ) continue;
        all_items << io_plugin->name();
        if(io_plugin->isDefaultLoader(fileinfo.completeSuffix()))
          selected_items.prepend(io_plugin->name());
        else
          selected_items << io_plugin->name();
      }
    }
    //if no plugin is correct, offer them all.
    for(CGAL::Three::Polyhedron_demo_io_plugin_interface* io_plugin : io_plugins) {
        all_items << io_plugin->name();
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
  //viewer->makeCurrent();
  if(!ok || load_pair.first.isEmpty()) { return; }

  if (load_pair.second)
  {
    connect(actionResetDefaultLoaders, SIGNAL(triggered()),
            this, SLOT(reset_default_loaders()));
    default_plugin_selection[fileinfo.completeSuffix()]=load_pair.first;
    insertActionBeforeLoadPlugin(ui->menuFile, actionResetDefaultLoaders);
  }


  settings.setValue("OFF open directory",
                    fileinfo.absoluteDir().absolutePath());
  loadItem(fileinfo, findLoader(load_pair.first), ok);

  if(!ok)
    return;
  this->addToRecentFiles(fileinfo.absoluteFilePath());
  updateViewersBboxes(true);
}

bool MainWindow::open(QString filename, QString loader_name) {
  QFileInfo fileinfo(filename);
  boost::optional<bool> item_opt;
  try {
    item_opt = wrap_a_call_to_cpp
        ([this, fileinfo, loader_name]()
    {
      bool ok;
      loadItem(fileinfo, findLoader(loader_name), ok);
      return ok;
    },
    this, __FILE__, __LINE__
    );
    if(!item_opt) return false;
  }
  catch(std::logic_error& e) {
    std::cerr << e.what() << std::endl;
    return false;
  }
  return true;
}


QList<Scene_item*> MainWindow::loadItem(QFileInfo fileinfo,
                                        CGAL::Three::Polyhedron_demo_io_plugin_interface* loader,
                                        bool &ok,
                                        bool add_to_scene) {
  if(!fileinfo.isFile() || !fileinfo.isReadable()) {
    QMessageBox::warning(this, tr("Error"),
                         QString("File %1 is not a readable file.")
                         .arg(fileinfo.absoluteFilePath()));
  }
  QCursor tmp_cursor(Qt::WaitCursor);
  CGAL::Three::Three::CursorScopeGuard guard(tmp_cursor);
  QList<Scene_item*> result = loader->load(fileinfo, ok, add_to_scene);
  if(!ok)
  {
    QApplication::restoreOverrideCursor();
      QMessageBox::warning(this, tr("Error"),
                           QString("Could not load item from file %1 using plugin %2")
                                                      .arg(fileinfo.absoluteFilePath()).arg(loader->name()));
      return QList<Scene_item*>();
  }
  selectSceneItem(scene->item_id(result.back()));
  for(Scene_item* item : result)
  {
    CGAL::Three::Scene_group_item* group =
        qobject_cast<CGAL::Three::Scene_group_item*>(item);
    if(group)
      scene->redraw_model();
    item->setProperty("source filename", fileinfo.absoluteFilePath());
    item->setProperty("loader_name", loader->name());
    item->setProperty("load_mates",QVariant::fromValue(result));
  }
  return result;
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
    if(s.empty())
      return;
    QModelIndex mi = proxyModel->mapFromSource(scene->getModelIndexFromId(i).first());
    sceneView->setCurrentIndex(mi);
    sceneView->selectionModel()->select(s,
                                        QItemSelectionModel::ClearAndSelect);
    sceneView->scrollTo(s.indexes().first());
    sceneView->setCurrentIndex(sceneView->selectionModel()->selectedIndexes().first());
  }
}

void MainWindow::selectSceneItems(QList<int> is)
{
  if(is.first() < 0 || is.last() >= scene->numberOfEntries()) {
    sceneView->selectionModel()->clearSelection();
    updateInfo();
    updateDisplayInfo();
  }
  else {
    QItemSelection s =
      proxyModel->mapSelectionFromSource(scene->createSelection(is));

    QModelIndex i = proxyModel->mapFromSource(scene->getModelIndexFromId(is.first()).first());
    sceneView->setCurrentIndex(i);
    sceneView->selectionModel()->select(s,
                                        QItemSelectionModel::ClearAndSelect);
    if(!s.empty())
      sceneView->scrollTo(s.indexes().first());
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
  sceneView->selectAll();
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
  Q_FOREACH(CGAL::QGLViewer* vi, CGAL::QGLViewer::QGLViewerPool())
  {
    if(vi == NULL)
      continue;

    if(item != NULL && item->manipulatable()) {
      vi->setManipulatedFrame(item->manipulatedFrame());
    } else {
      vi->setManipulatedFrame(0);
    }
    if(vi->manipulatedFrame() == 0) {
      Q_FOREACH(CGAL::Three::Scene_item* item, scene->entries()) {
        if(item->manipulatable() && item->manipulatedFrame() != 0) {
          if(vi->manipulatedFrame() != 0) {
            // there are at least two possible frames
            vi->setManipulatedFrame(0);
            break;
          } else {
            vi->setManipulatedFrame(item->manipulatedFrame());
          }
        }
      }
    }
    if(vi->manipulatedFrame() != 0) {
      connect(vi->manipulatedFrame(), SIGNAL(modified()),
              this, SLOT(updateInfo()));
    }
    vi->update();
  }
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
        reload->setProperty("is_groupable", true);
        connect(reload, SIGNAL(triggered()),
                this, SLOT(reloadItem()));
      }
      QAction* saveas = menu->addAction(tr("&Save as..."));
      saveas->setData(QVariant::fromValue((void*)item));
      connect(saveas,  SIGNAL(triggered()),
              this, SLOT(on_actionSaveAs_triggered()));
      QAction* showobject = menu->addAction(tr("&Zoom to this Object"));
      showobject->setData(QVariant::fromValue((void*)item));
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
  if(scene->selectionIndices().isEmpty())return;
  int main_index = scene->selectionIndices().first();

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
      else if(scene->selectionIndices().size() > 1 )
      {
        QMap<QString, QAction*> menu_actions;
        QVector<QMenu*> slider_menus;
        bool has_stats = false;
        bool has_reload = false;
        Q_FOREACH(Scene::Item_id id, scene->selectionIndices())
        {
          if(!scene->item(id)->property("source filename").toString().isEmpty())
          {
            has_reload = true;
            break;
          }
        }
        Q_FOREACH(QAction* action, scene->item(main_index)->contextMenu()->actions())
        {
          if(action->property("is_groupable").toBool())
          {
            menu_actions[action->text()] = action;
            if(action->text() == QString("Alpha value"))
            {
              menu_actions["alpha slider"] = action->menu()->actions().last();
            }
            else if(action->text() == QString("Points Size"))
            {
              menu_actions["points slider"] = action->menu()->actions().last();
            }
            else if(action->text() == QString("Normals Length"))
            {
              menu_actions["normals slider"] = action->menu()->actions().last();
            }
            else if(action->text() == QString("Line Width"))
            {
              menu_actions["Line width"] = action->menu()->actions().last();
            }
          }

        }
        Q_FOREACH(Scene::Item_id index, scene->selectionIndices())
        {
          if(index == main_index)
            continue;

          CGAL::Three::Scene_item* item = scene->item(index);
          if(!item)
            continue;
          if(item->has_stats())
            has_stats = true;
        }
        QMenu menu;
        menu.addAction(actionAddToGroup);
        menu.insertSeparator(0);
        Q_FOREACH(QString name, menu_actions.keys())
        {
          if(name == QString("alpha slider")
             || name == QString("points slider")
             || name == QString("normals slider"))
            continue;
          if(name == QString("Alpha value"))
          {
            QWidgetAction* sliderAction = new QWidgetAction(&menu);
            QSlider* slider = new QSlider(&menu);
            slider->setMinimum(0);
            slider->setMaximum(255);
            slider->setValue(
                  qobject_cast<QSlider*>(
                    qobject_cast<QWidgetAction*>
                    (menu_actions["alpha slider"])->defaultWidget()
                  )->value());
            slider->setOrientation(Qt::Horizontal);
            sliderAction->setDefaultWidget(slider);

            connect(slider, &QSlider::valueChanged, [this, slider]()
            {
              Q_FOREACH(Scene::Item_id id, scene->selectionIndices())
              {
                Scene_item* item = scene->item(id);
                Q_FOREACH(QAction* action, item->contextMenu()->actions())
                {
                  if(action->text() == "Alpha value")
                  {
                    QWidgetAction* sliderAction = qobject_cast<QWidgetAction*>(action->menu()->actions().last());
                    QSlider* ac_slider = qobject_cast<QSlider*>(sliderAction->defaultWidget());
                    ac_slider->setValue(slider->value());
                    break;
                  }
                }
              }
            });
            QMenu* new_menu = new QMenu("Alpha value", &menu);
              new_menu->addAction(sliderAction);
              slider_menus.push_back(new_menu);
          }
          else if(name == QString("Points Size"))
          {
            QWidgetAction* sliderAction = new QWidgetAction(&menu);
            QSlider* slider = new QSlider(&menu);
            slider->setMinimum(1);
            slider->setMaximum(25);
            slider->setValue(
                  qobject_cast<QSlider*>(
                    qobject_cast<QWidgetAction*>
                    (menu_actions["points slider"])->defaultWidget()
                  )->value());
            slider->setOrientation(Qt::Horizontal);
            sliderAction->setDefaultWidget(slider);

            connect(slider, &QSlider::valueChanged, [this, slider]()
            {
              Q_FOREACH(Scene::Item_id id, scene->selectionIndices())
              {
                Scene_item* item = scene->item(id);
                Q_FOREACH(QAction* action, item->contextMenu()->actions())
                {
                  if(action->text() == "Points Size")
                  {
                    QWidgetAction* sliderAction = qobject_cast<QWidgetAction*>(action->menu()->actions().last());
                    QSlider* ac_slider = qobject_cast<QSlider*>(sliderAction->defaultWidget());
                    ac_slider->setValue(slider->value());
                    break;
                  }
                }
              }
            });
            QMenu* new_menu = new QMenu("Points Size", &menu);
              new_menu->addAction(sliderAction);
              slider_menus.push_back(new_menu);
          }
          else if(name == QString("Normals Length"))
          {
            QWidgetAction* sliderAction = new QWidgetAction(&menu);
            QSlider* slider = new QSlider(&menu);
            slider->setMinimum(0);
            slider->setMaximum(100);
            slider->setValue(
                  qobject_cast<QSlider*>(
                    qobject_cast<QWidgetAction*>
                    (menu_actions["normals slider"])->defaultWidget()
                  )->value());
            slider->setOrientation(Qt::Horizontal);
            sliderAction->setDefaultWidget(slider);

            connect(slider, &QSlider::valueChanged, [this, slider]()
            {
              Q_FOREACH(Scene::Item_id id, scene->selectionIndices())
              {
                Scene_item* item = scene->item(id);
                Q_FOREACH(QAction* action, item->contextMenu()->actions())
                {
                  if(action->text() == "Normals Length")
                  {
                    QWidgetAction* sliderAction = qobject_cast<QWidgetAction*>(action->menu()->actions().last());
                    QSlider* ac_slider = qobject_cast<QSlider*>(sliderAction->defaultWidget());
                    ac_slider->setValue(slider->value());
                    break;
                  }
                }
              }
            });
            QMenu* new_menu = new QMenu("Normals Length", &menu);
              new_menu->addAction(sliderAction);
              slider_menus.push_back(new_menu);
          }
          else if(name == QString("Line Width"))
          {
            QWidgetAction* sliderAction = new QWidgetAction(&menu);
            QSlider* slider = new QSlider(&menu);
            slider->setMinimum(1);
            float lineWidth[2];
            if(!viewer->isOpenGL_4_3())
              viewer->glGetFloatv(GL_LINE_WIDTH_RANGE, lineWidth);
            else
            {
              lineWidth[0] = 0;
              lineWidth[1] = 10;
            }
            slider->setMaximum(lineWidth[1]);
            slider->setValue(
                  qobject_cast<QSlider*>(
                    qobject_cast<QWidgetAction*>
                    (menu_actions["Line width"])->defaultWidget()
                  )->value());
            slider->setOrientation(Qt::Horizontal);
            sliderAction->setDefaultWidget(slider);

            connect(slider, &QSlider::valueChanged, [this, slider]()
            {
              Q_FOREACH(Scene::Item_id id, scene->selectionIndices())
              {
                Scene_item* item = scene->item(id);
                Q_FOREACH(QAction* action, item->contextMenu()->actions())
                {
                  if(action->text() == "Line Width")
                  {
                    QWidgetAction* sliderAction = qobject_cast<QWidgetAction*>(action->menu()->actions().last());
                    QSlider* ac_slider = qobject_cast<QSlider*>(sliderAction->defaultWidget());
                    ac_slider->setValue(slider->value());
                    break;
                  }
                }
              }
            });
            QMenu* new_menu = new QMenu("Line Width", &menu);
              new_menu->addAction(sliderAction);
              slider_menus.push_back(new_menu);
          }
          else
          {
            QAction* action = menu.addAction(name);
            connect(action, &QAction::triggered, this, &MainWindow::propagate_action);
          }
        }
        if(!slider_menus.empty())
        {
          Q_FOREACH(QMenu* m, slider_menus){
            menu.addMenu(m);
          }
          menu.insertSeparator(0);
        }
        if(has_stats)
        {
          QAction* actionStatistics =
              menu.addAction(tr("Statistics..."));
          actionStatistics->setObjectName("actionStatisticsOnPolyhedron");
          connect(actionStatistics, SIGNAL(triggered()),
                  this, SLOT(statisticsOnItem()));
        }
        if(has_reload)
        {
          QAction* reload = menu.addAction(tr("&Reload Item from File"));
          reload->setProperty("is_groupable", true);
          connect(reload, SIGNAL(triggered()),
                  this, SLOT(reloadItem()));
        }
        QAction* saveas = menu.addAction(tr("&Save as..."));
        connect(saveas,  SIGNAL(triggered()),
                this, SLOT(on_actionSaveAs_triggered()));
        menu.addMenu(ui->menuOperations);
        menu.exec(sender->mapToGlobal(p));
        return;
      }
  }
  showSceneContextMenu(main_index, sender->mapToGlobal(p));
  return;
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
    CGAL::Bbox_3 bbox = item->bbox();
    if(bbox !=CGAL::Bbox_3())
      item_text += QString("<div>Bounding box: min (%1,%2,%3), max (%4,%5,%6)</div>")
          .arg(bbox.xmin(),0, 'g', 17)
          .arg(bbox.ymin(),0, 'g', 17)
          .arg(bbox.zmin(),0, 'g', 17)
          .arg(bbox.xmax(),0, 'g', 17)
          .arg(bbox.ymax(),0, 'g', 17)
          .arg(bbox.zmax(),0, 'g', 17);
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
    viewer->setAntiAliasing(settings.value("antialiasing", false).toBool());
    viewer->setFastDrawing(settings.value("quick_camera_mode", true).toBool());
    scene->enableVisibilityRecentering(settings.value("offset_update", false).toBool());
    viewer->textRenderer()->setMax(settings.value("max_text_items", 10000).toInt());
    int val  = settings.value("transparency_pass_number", 4).toInt();
    if (val < 4 ) {
      val = 4;
      settings.setValue("transparency_pass_number", 4);
    }
    viewer->setTotalPass(val);
    CGAL::Three::Three::s_defaultSMRM = CGAL::Three::Three::modeFromName(
          settings.value("default_sm_rm", "flat+edges").toString());
    CGAL::Three::Three::s_defaultPSRM = CGAL::Three::Three::modeFromName(
          settings.value("default_ps_rm", "points").toString());
    // read plugin blacklist
    QStringList blacklist=settings.value("plugin_blacklist",QStringList()).toStringList();
    Q_FOREACH(QString name,blacklist){ plugin_blacklist.insert(name); }
    def_save_dir = settings.value("default_saveas_dir", QDir::homePath()).toString();
    this->default_point_size = settings.value("points_size").toInt();
    this->default_normal_length = settings.value("normals_length").toInt();
    this->default_lines_width = settings.value("lines_width").toInt();
    setProperty("ws_url", settings.value("ws_server_url").toString());
}

void MainWindow::writeSettings()
{
  this->writeState("MainWindow");
  {
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

#endif
}

void MainWindow::on_actionLoad_triggered()
{
  QStringList filters;
  // we need to special case our way out of this
  filters << "All Files (*)";


  typedef QMap<QString, CGAL::Three::Polyhedron_demo_io_plugin_interface*> FilterPluginMap;
  FilterPluginMap filterPluginMap;

  for(CGAL::Three::Polyhedron_demo_io_plugin_interface* plugin : io_plugins) {
    QStringList split_filters = plugin->loadNameFilters().split(";;");
    Q_FOREACH(const QString& filter, split_filters) {
      FilterPluginMap::iterator it = filterPluginMap.find(filter);
      if(it != filterPluginMap.end()) {
        if(verbose)
        {
          qDebug() << "Duplicate Filter: " << it.value()->name();
          qDebug() << "This filter will not be available.";
        }
      } else {
        filterPluginMap[filter] = plugin;
      }
      filters << filter;
    }
  }

  QString directory = settings.value("OFF open directory",
                                     QDir::current().dirName()).toString();

  QFileDialog dialog(this);
  dialog.setDirectory(directory);
  dialog.setNameFilters(filters);
  dialog.setFileMode(QFileDialog::ExistingFiles);

  if(dialog.exec() != QDialog::Accepted) { return; }
  for(auto v : CGAL::QGLViewer::QGLViewerPool())
    v->update();
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
      bool ok;
      QList<Scene_item*> result = loadItem(info, selectedPlugin, ok);
      if(!ok)
        continue;
      for(Scene_item* item : result)
      {
        if(!item->property("already_colored").toBool())
        {
          ++nb_item;
          item->setColor(colors_[nb_item]);
        }
        selectSceneItem(scene->item_id(item));
        CGAL::Three::Scene_group_item* group =
            qobject_cast<CGAL::Three::Scene_group_item*>(item);
        if(group)
          scene->redraw_model();
      }
      this->addToRecentFiles(filename);
    } else {
      int scene_size = scene->numberOfEntries();
      open(filename);
      item = scene->item(scene->numberOfEntries()-1);
      if(scene->numberOfEntries() != scene_size
         && !item->property("already_colored").toBool())
        item->setColor(colors_[++nb_item]);
    }
  }
}

void MainWindow::on_actionSaveAs_triggered()
{
  QList<Scene_item*> to_save;
  for(Scene::Item_id id : scene->selectionIndices())
  {
    Scene_item* item = scene->item(id);
    to_save.append(item);
  }
  while(!to_save.empty())
  {
    Scene_item* item = to_save.front();
    QVector<CGAL::Three::Polyhedron_demo_io_plugin_interface*> canSavePlugins;
    QStringList filters;
    QString sf;
    for(CGAL::Three::Polyhedron_demo_io_plugin_interface* plugin : io_plugins) {
      if(plugin->canSave(item)) {
        canSavePlugins << plugin;
        filters += plugin->saveNameFilters();
        if(plugin->isDefaultLoader(item))
          sf = plugin->saveNameFilters().split(";;").first();
      }
    }
    QRegExp extensions("\\(\\*\\..+\\)");
    QStringList filter_exts;
    if(filters.empty())
    {
      QMessageBox::warning(this,
                           tr("Cannot save"),
                           tr("The selected object %1 cannot be saved.")
                           .arg(item->name()));
          return;
    }
    Q_FOREACH(QString string, filters)
    {
      QStringList sl = string.split(";;");
      Q_FOREACH(QString s, sl){
        int pos = extensions.indexIn(s);
        if( pos >-1)
          filter_exts.append(extensions.capturedTexts());
      }
    }
    filters << tr("All files (*)");
    if(canSavePlugins.isEmpty()) {
      QMessageBox::warning(this,
                           tr("Cannot save"),
                           tr("The selected object %1 cannot be saved.")
                           .arg(item->name()));
      continue;
    }
    QString caption = tr("Save %1 to File...").arg(item->name());
    QString dir = item->property("source filename").toString();
    if(dir.isEmpty() &&
       !item->property("defaultSaveDir").toString().isEmpty())
    {
      dir = item->property("defaultSaveDir").toString();
    }
    else if(!last_saved_dir.isEmpty() && dir.isEmpty())
      dir = QString("%1/%2").arg(last_saved_dir).arg(item->defaultSaveName());
    else if(dir.isEmpty())
      dir = QString("%1/%2").arg(def_save_dir).arg(item->name());
    QString filename =
        QFileDialog::getSaveFileName(this,
                                     caption,
                                     dir,
                                     filters.join(";;"),
                                     &sf);

    if(filename.isEmpty())
      return;
    last_saved_dir = QFileInfo(filename).absoluteDir().path();
    extensions.indexIn(sf.split(";;").first());
    QString filter_ext, filename_ext;
    filter_ext = extensions.cap().split(" ").first();// in case of syntax like (*.a *.b)

    filter_ext.remove(")");
    filter_ext.remove("(");
    //remove *
    filter_ext=filter_ext.right(filter_ext.size()-1);

    QStringList filename_split = filename.split(".");
    filename_split.removeFirst();
    filename_ext = filename_split.join(".");
    filename_ext.push_front(".");

    QStringList final_extensions;
    Q_FOREACH(QString string, filter_exts)
    {
      Q_FOREACH(QString s, string.split(" ")){// in case of syntax like (*.a *.b)
        s.remove(")");
        s.remove("(");
        //remove *
        s=s.right(s.size()-1);
        final_extensions.append(s);
      }
    }
    bool ok = false;
    while(!ok)
    {
      if(final_extensions.contains(filename_ext))
      {
        ok = true;
      }
      else{
        QStringList shatterd_filename_ext = filename_ext.split(".");
        if(!shatterd_filename_ext.last().isEmpty())
        {
          shatterd_filename_ext.removeFirst();//removes ""
          shatterd_filename_ext.removeFirst();
          filename_ext = shatterd_filename_ext.join(".");
          filename_ext.push_front(".");
        }
        else
          break;
      }
    }
    if(!ok)
    {
      filename = filename.append(filter_ext);
    }
    for(auto v : CGAL::QGLViewer::QGLViewerPool())
      v->update();
    save(filename, to_save);
  }
}

void MainWindow::save(QString filename, QList<CGAL::Three::Scene_item*>& to_save) {
  QFileInfo fileinfo(filename);
  bool saved = false;
  for(CGAL::Three::Polyhedron_demo_io_plugin_interface* plugin : io_plugins) {
    if(  plugin->canSave(to_save.front()) &&
         file_matches_filter(plugin->saveNameFilters(),filename.toLower()) )
    {
      if(plugin->save(fileinfo, to_save))
      {
        saved = true;
        break;
      }
    }
  }
  if(!saved)
  {
    QMessageBox::warning(this,
                         tr("Cannot save"),
                         tr("The selected object %1 was not saved. (Maybe a wrong extension ?)")
                         .arg(to_save.front()->name()));
    to_save.pop_front();
  }
}

void MainWindow::on_actionSaveSnapshot_triggered()
{
  viewer->saveSnapshot();
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
  QList<int> all_ids;
  for(int i = 0; i < scene->numberOfEntries(); ++i)
    all_ids.push_back(i);
  scene->setSelectedItemsList(all_ids);
  on_actionErase_triggered();
}

void MainWindow::on_actionDuplicate_triggered()
{
  int index = scene->duplicate(getSelectedSceneItemIndex());
  selectSceneItem(index);
}

void MainWindow::on_actionShowHide_triggered()
{
  scene->setUpdatesEnabled(false);
  Q_FOREACH(QModelIndex index, sceneView->selectionModel()->selectedRows())
  {
    int i = scene->getIdFromModelIndex(proxyModel->mapToSource(index));
    CGAL::Three::Scene_item* item = scene->item(i);
    item->setVisible(!item->visible());
    item->redraw();
  }
  scene->setUpdatesEnabled(true);
//  updateViewersBboxes(false); //Not usable :when the scene changes scale, smaller items disappear.
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
  float lineWidth[2];
  if(!viewer->isOpenGL_4_3())
    viewer->glGetFloatv(GL_LINE_WIDTH_RANGE, lineWidth);
  else
  {
    lineWidth[0] = 0;
    lineWidth[1] = 10;
  }
  prefdiag.linesHorizontalSlider->setMinimum(lineWidth[0]);
  prefdiag.linesHorizontalSlider->setMaximum(lineWidth[1]);

  prefdiag.offset_updateCheckBox->setChecked(
        settings.value("offset_update", false).toBool());
  connect(prefdiag.offset_updateCheckBox, SIGNAL(toggled(bool)),
          scene, SLOT(enableVisibilityRecentering(bool)));

  prefdiag.antialiasingCheckBox->setChecked(settings.value("antialiasing", false).toBool());
  connect(prefdiag.antialiasingCheckBox, SIGNAL(toggled(bool)),
          viewer, SLOT(setAntiAliasing(bool)));

  prefdiag.quick_cameraCheckBox->setChecked(
        settings.value("quick_camera_mode", true).toBool());
  connect(prefdiag.quick_cameraCheckBox, SIGNAL(toggled(bool)),
          viewer, SLOT(setFastDrawing(bool)));
  prefdiag.max_itemsSpinBox->setValue(viewer->textRenderer()->getMax_textItems());

  connect(prefdiag.max_itemsSpinBox,static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
          this, [this](int i){
    setMaxTextItemsDisplayed(i);
  });
  prefdiag.transpSpinBox->setValue(viewer->total_pass());
  connect(prefdiag.transpSpinBox, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
              this, [this](int i)
  {
    setTransparencyPasses(i);
  });
  prefdiag.pointsHorizontalSlider->setValue(this->default_point_size);
  connect(prefdiag.pointsHorizontalSlider, &QSlider::valueChanged,
              this, [this](int i)
  {
    this->default_point_size = i;
  });
  prefdiag.normalsHorizontalSlider->setValue(this->default_normal_length);
  connect(prefdiag.normalsHorizontalSlider, &QSlider::valueChanged,
              this, [this](int i)
  {
    this->default_normal_length = i;
  });
  prefdiag.linesHorizontalSlider->setValue(this->default_lines_width);
  connect(prefdiag.linesHorizontalSlider, &QSlider::valueChanged,
              this, [this](int i)
  {
    this->default_lines_width = i;
  });
  connect(prefdiag.background_colorPushButton, &QPushButton::clicked,
          this, &MainWindow::setBackgroundColor);

  connect(prefdiag.default_save_asPushButton, &QPushButton::clicked,
          this, &MainWindow::setDefaultSaveDir);

  connect(prefdiag.lightingPushButton, &QPushButton::clicked,
          this, &MainWindow::setLighting_triggered);

  prefdiag.surface_meshComboBox->setCurrentText(CGAL::Three::Three::modeName(
                                                  CGAL::Three::Three::s_defaultSMRM));
  connect(prefdiag.surface_meshComboBox, &QComboBox::currentTextChanged,
          this, [this](const QString& text){
    this->s_defaultSMRM = CGAL::Three::Three::modeFromName(text);
  });

  prefdiag.point_setComboBox->setCurrentText(CGAL::Three::Three::modeName(
                                               CGAL::Three::Three::s_defaultPSRM));
  connect(prefdiag.point_setComboBox, &QComboBox::currentTextChanged,
          this, [this](const QString& text){
    this->s_defaultPSRM = CGAL::Three::Three::modeFromName(text);
  });

  connect(prefdiag.backFrontColor_pushButton, &QPushButton::clicked,
          this, [](){
    qobject_cast<Viewer*>(CGAL::Three::Three::activeViewer())->setBackFrontColors();
  });

  std::vector<QTreeWidgetItem*> items;
  QBrush successBrush(Qt::green),
      errorBrush(Qt::red),
      ignoredBrush(Qt::lightGray);

  //add blacklisted plugins
  Q_FOREACH (QString name, PathNames_map.keys())
  {
    QTreeWidgetItem *item = new QTreeWidgetItem(prefdiag.treeWidget);
    item->setText(1, name);
    if(plugin_blacklist.contains(name)){
      item->setCheckState(0, Qt::Unchecked);
    }
    else{
      item->setCheckState(0, Qt::Checked);
    }
    if(pluginsStatus_map[name] == QString("success"))
      item->setBackground(1, successBrush);
    else if(ignored_map[name]){
      item->setBackground(1, ignoredBrush);
    }
    else{
      item->setBackground(1, errorBrush);
    }
    items.push_back(item);
  }
  connect(prefdiag.detailsPushButton, &QPushButton::clicked,
          this, [this, prefdiag](){
    QStringList titles;
    titles << "Name" << "Keywords" << "ConfigDate";
    QDialog dialog(this);
    Ui::DetailsDialog detdiag;
    detdiag.setupUi(&dialog);
    QTreeWidgetItem *header = new QTreeWidgetItem(titles);
    detdiag.treeWidget->setHeaderItem(header);
    Q_FOREACH(QTreeWidgetItem* plugin_item, prefdiag.treeWidget->selectedItems())
    {
      QString name = plugin_item->text(1);
      QString keywords = plugin_metadata_map[name].first.join(", ");
      QString date = plugin_metadata_map[name].second;
      QStringList values;
      values << name << keywords << date;
      new QTreeWidgetItem(detdiag.treeWidget, values);
    }
    for(int i=0; i<3; ++i)
    {
      detdiag.treeWidget->resizeColumnToContents(i);
    }
    connect(detdiag.treeWidget, &QTreeWidget::clicked,
            this, [this, detdiag](){
      if(detdiag.treeWidget->selectedItems().isEmpty())
        detdiag.textBrowser->setText("");
      else {
        QString name = detdiag.treeWidget->selectedItems().first()->text(0);
        QString status = pluginsStatus_map[name];
        QString path = PathNames_map[name];
        detdiag.textBrowser->setText(QString("Path: %1 \nStatus: %2").arg(path).arg(status));
      }
    });
    dialog.exec();
  });
  connect(prefdiag.sshButton, &QPushButton::clicked,
          this, [this](){
    QDialog dialog(this);
    Ui::SSHDialog sshdiag;
    sshdiag.setupUi(&dialog);

#ifdef CGAL_USE_SSH
    sshdiag.userBox->setEnabled(true);
    sshdiag.serverBox->setEnabled(true);
    sshdiag.pkBox->setEnabled(true);
    sshdiag.privkBox->setEnabled(true);
    sshdiag.userEdit->setText(settings.value("ssh_user", QString()).toString());
    sshdiag.serverEdit->setText(settings.value("ssh_server", QString()).toString());
    sshdiag.publicEdit->setText(settings.value("ssh_public_key", QString()).toString());
    sshdiag.privkEdit->setText(settings.value("ssh_priv_key", QString()).toString());
    connect(sshdiag.pubButton, &QPushButton::clicked,
            this, [this, sshdiag](){
      QFileDialog diag(this,
                       "Public Key",
                       "",
                       "All Files (*)");
      diag.setFilter(QDir::Hidden|QDir::Files|QDir::Dirs|QDir::NoDotAndDotDot);
      if(!diag.exec())
        return;
      sshdiag.publicEdit->setText(diag.selectedFiles().front());
    });
    connect(sshdiag.privButton, &QPushButton::clicked,
            this, [this, sshdiag](){
      QFileDialog diag(this,
                       "Private Key",
                       "",
                       "All Files (*)");
      diag.setFilter(QDir::Hidden|QDir::Files|QDir::Dirs|QDir::NoDotAndDotDot);
      if(!diag.exec())
        return;
      sshdiag.privkEdit->setText(diag.selectedFiles().front());
    });
#else
    sshdiag.userBox->setEnabled(false);
    sshdiag.serverBox->setEnabled(false);
    sshdiag.pkBox->setEnabled(false);
    sshdiag.privkBox->setEnabled(false);
#endif
    sshdiag.wsEdit->setText(settings.value("ws_server_url", QString()).toString());

    dialog.exec();
    if ( dialog.result() )
    {
#ifdef CGAL_USE_SSH
      settings.setValue("ssh_user",
                        sshdiag.userEdit->text());
      settings.setValue("ssh_server",
                        sshdiag.serverEdit->text());
      settings.setValue("ssh_public_key",
                        sshdiag.publicEdit->text());
      settings.setValue("ssh_priv_key",
                        sshdiag.privkEdit->text());
#endif
      settings.setValue("ws_server_url",
                        sshdiag.wsEdit->text());
      setProperty("ws_url", sshdiag.wsEdit->text());
    }
  });

  dialog.exec();

  if ( dialog.result() )
  {
    plugin_blacklist.clear();

    for (std::size_t k=0; k<items.size(); ++k)
    {
     QTreeWidgetItem* item=items[k];
      if (item->checkState(0)==Qt::Unchecked)
        plugin_blacklist.insert(item->text(1));
    }

    //write settings
    settings.setValue("antialiasing",
                      prefdiag.antialiasingCheckBox->isChecked());
    settings.setValue("offset_update",
                      prefdiag.offset_updateCheckBox->isChecked());
    settings.setValue("quick_camera_mode",
                      prefdiag.quick_cameraCheckBox->isChecked());
    settings.setValue("transparency_pass_number",
                      viewer->total_pass());
    settings.setValue("max_text_items",
                      viewer->textRenderer()->getMax_textItems());
    settings.setValue("background_color",viewer->backgroundColor().name());
    settings.setValue("default_sm_rm", CGAL::Three::Three::modeName(
                        CGAL::Three::Three::defaultSurfaceMeshRenderingMode()));
    settings.setValue("default_ps_rm", CGAL::Three::Three::modeName(
                        CGAL::Three::Three::defaultPointSetRenderingMode()));
    settings.setValue("points_size", this->default_point_size);
    settings.setValue("normals_length", this->default_normal_length);
    settings.setValue("lines_width", this->default_lines_width);

  }
  else
  {
    readSettings();
  }
}

void MainWindow::setBackgroundColor()
{
  QColor c =  QColorDialog::getColor();
  if(c.isValid()) {
    Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
    {
      if(v == NULL)
        continue;
      v->setBackgroundColor(c);
      v->update();
    }
  }

}

void MainWindow::setLighting_triggered()
{
  qobject_cast<Viewer*>(CGAL::Three::Three::activeViewer())->setLighting();
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
    CGAL::qglviewer::Vec min((float)bbox.xmin()+viewer->offset().x, (float)bbox.ymin()+viewer->offset().y, (float)bbox.zmin()+viewer->offset().z),
        max((float)bbox.xmax()+viewer->offset().x, (float)bbox.ymax()+viewer->offset().y, (float)bbox.zmax()+viewer->offset().z);
    viewer->setSceneBoundingBox(min, max);
    viewerShow((float)min.x, (float)min.y, (float)min.z,
               (float)max.x, (float)max.y, (float)max.z);
  }
}
/* to check
QString MainWindow::cameraString() const
{
  const CGAL::qglviewer::Vec pos = viewer->camera()->position() - viewer->offset();
  const CGAL::qglviewer::Quaternion q = viewer->camera()->orientation();

  return QString("%1 %2 %3 %4 %5 %6 %7")
    .arg(pos[0])
    .arg(pos[1])
    .arg(pos[2])
    .arg(q[0])
    .arg(q[1])
    .arg(q[2])
    .arg(q[3]);
}*/
QString MainWindow::cameraString(CGAL::Three::Viewer_interface* v) const
{
  return v->dumpCameraCoordinates();
}

void MainWindow::setAddKeyFrameKeyboardModifiers(::Qt::KeyboardModifiers m)
{
  viewer->setAddKeyFrameKeyboardModifiers(m);
}

void MainWindow::recenterScene()
{
  //force the recomputaion of the bbox
  bbox_need_update = true;
  CGAL::qglviewer::Vec min, max;
  computeViewerBBox(min, max);
  updateViewerBbox(static_cast<Viewer*>(activeViewer()), true, min, max);
  activeViewer()->showEntireScene();
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
  if(scene->index(0,0,index).isValid())
  {
    recurseExpand(scene->index(0,0,index));
  }
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
  for(Scene::Item_id id : scene->selectionIndices())
  {
    scene->changeGroup(scene->item(id), group);
  }
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
    connect(statistics_ui->exportButton, &QPushButton::clicked,
            this, &MainWindow::exportStatistics);
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
  Q_FOREACH(int id, scene->selectionIndices())
  {
    Scene_item* item = scene->item(id);
    QString classname = item->property("classname").toString();
    if(classname.isEmpty())
       classname = item->metaObject()->className();
    if(!classnames.contains(classname))
      classnames << classname;
  }
  //2nd step : separate the selection in lists corresponding to their classname
  QVector< QList<Scene_item*> > items;
  items.resize(classnames.size());
  Q_FOREACH(int id, scene->selectionIndices())
  {
    Scene_item* s_item = scene->item(id);
    for(int i=0; i<items.size(); i++)
    {
      Scene_item* item = scene->item(id);
      QString classname = item->property("classname").toString();
      if(classname.isEmpty())
         classname = item->metaObject()->className();
      if(classnames.at(i).contains(classname))
      {
        items[i] << s_item;
        break;
      }
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


void MainWindow::setMaxTextItemsDisplayed(int val)
{
  for(auto v : CGAL::QGLViewer::QGLViewerPool())
    qobject_cast<CGAL::Three::Viewer_interface*>(v)->textRenderer()->setMax(val);
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
#if QT_VERSION >= QT_VERSION_CHECK(5, 11, 0)
  sceneView->header()->resizeSection(Scene::ColorColumn, sceneView->header()->fontMetrics().horizontalAdvance("_#_"));
  sceneView->resizeColumnToContents(Scene::RenderingModeColumn);
  sceneView->header()->resizeSection(Scene::ABColumn, sceneView->header()->fontMetrics().horizontalAdvance(QString("_AB_")));
  sceneView->header()->resizeSection(Scene::VisibleColumn, sceneView->header()->fontMetrics().horizontalAdvance(QString("_View_")));
#else
  sceneView->header()->resizeSection(Scene::ColorColumn, sceneView->header()->fontMetrics().width("_#_"));
  sceneView->resizeColumnToContents(Scene::RenderingModeColumn);
  sceneView->header()->resizeSection(Scene::ABColumn, sceneView->header()->fontMetrics().width(QString("_AB_")));
  sceneView->header()->resizeSection(Scene::VisibleColumn, sceneView->header()->fontMetrics().width(QString("_View_")));
#endif
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
  if(nb_files<2)
    return;
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
  for(auto v : CGAL::QGLViewer::QGLViewerPool())
    v->update();
}


void MainWindow::exportStatistics()
{
  std::vector<Scene_item*> items;
  Q_FOREACH(int id, getSelectedSceneItemIndices())
  {
    Scene_item* s_item = scene->item(id);
    items.push_back(s_item);
  }

  QString str;
  Q_FOREACH(Scene_item* sit, items)
  {
    CGAL::Three::Scene_item::Header_data data = sit->header();
    if(data.titles.size()>0)
    {
      int titles_limit =0;
      int title = 0;
      str.append(QString("%1: \n").arg(sit->name()));
      for(int j=0; j<data.categories.size(); j++)
      {
        str.append(QString("  %1: \n")
                   .arg(data.categories[j].first));
        titles_limit+=data.categories[j].second;
        for(;title<titles_limit; ++title)
        {
          str.append(QString("    %1: ").arg(data.titles.at(title)));
          str.append(QString("%1\n").arg(sit->computeStats(title)));
        }
      }
    }
  }

  QString filename =
      QFileDialog::getSaveFileName((QWidget*)sender(),
                                   "",
                                   QString("Statistics.txt"),
                                   "Text Files (*.txt)");
  if(filename.isEmpty())
    return;
  QFile output(filename);
  output.open(QIODevice::WriteOnly | QIODevice::Text);

  if(!output.isOpen()){
    qDebug() << "- Error, unable to open" << "outputFilename" << "for output";
  }
  QTextStream outStream(&output);
  outStream << str;
  output.close();
}

void MainWindow::propagate_action()
{
  QAction* sender = qobject_cast<QAction*>(this->sender());
  if(!sender) return;
  QString name = sender->text();
  Q_FOREACH(Scene::Item_id id, scene->selectionIndices())
  {
    Scene_item* item = scene->item(id);
    Q_FOREACH(QAction* action, item->contextMenu()->actions())
    {
      if(action->text() == name)
      {
        action->trigger();
        break;
      }
    }
  }
}

QString make_fullpath(const QString& filename, bool duplicate = false)
{
  QString fullpath = QString("%1/%2").arg(QDir::tempPath()).arg(filename);
  QString tmp_fullpath = fullpath;
  if(duplicate)
  {
    int i=0;
    while(QFileInfo(tmp_fullpath).exists())
    {
      QString basename = QFileInfo(tmp_fullpath).baseName();
      QString dir = QFileInfo(tmp_fullpath).dir().path();
      QString suffix= QFileInfo(fullpath).completeSuffix();
      tmp_fullpath=QString("%1/%2%3.%4").arg(dir).arg(basename).arg(++i).arg(suffix);
    }
  }
  return tmp_fullpath;
}
/*
 The two following functions allow to create files from string and strings from files.
 This is used as a workaround of the absence of stream management in our IO system.
 The whole to/from Base64 is used to avoid problems with binary formats. Everything is written
 as a base64 binary string, and converted back to what it was.
*/
QByteArray file_to_string(const char* filename)
{
  std::ifstream f(filename, std::ifstream::binary);
  // get size of file
  f.seekg (0,f.end);
  long size = f.tellg();
  f.seekg (0);
  std::ostringstream ss;
  // allocate memory for file content
    char* buffer = new char[size];

    // read content of infile
    f.read(buffer,size);

    // write to outfile
    ss.write(buffer,size);
    // release dynamically-allocated memory
    delete[] buffer;
  //ss.write( << f.rdbuf(); // reading data
  f.close();
  std::string st = ss.str();
  QByteArray ba(st.c_str(), static_cast<int>(st.size()));
  return ba;
}

QString MainWindow::write_string_to_file(const QString& str, const QString &filename)
{
  QString fullpath = make_fullpath(filename);
  std::ofstream f(fullpath.toStdString().c_str(), std::ofstream::binary);
  QByteArray compressed_item(str.toStdString().c_str());
  QByteArray item = qUncompress(QByteArray::fromBase64(compressed_item));
  QByteArray bb = item;
  f.write(bb.constData(),bb.toStdString().size());
  f.close();
  return fullpath;
}

void MainWindow::on_actionSa_ve_Scene_as_Script_triggered()
{
  if(scene->numberOfEntries() == 0)
    return;
  bool do_upload = false;
#ifdef CGAL_USE_SSH
  QString user = settings.value("ssh_user", QString()).toString();
  if(!user.isEmpty())
  {
    QMessageBox::StandardButton doyou =
        QMessageBox::question(this, tr("Upload ?"), tr("Do you wish to upload the scene"
                                                       " using the SSH preferences ?"));
    do_upload = (doyou == QMessageBox::Yes);
  }
#endif

  QString filename;

  if(do_upload){
    filename = QString("%1/save_scene.js").arg(QDir::tempPath());
  }else{
    filename = QFileDialog::getSaveFileName(this,
                                            "Save the Scene as a Script File",
                                            last_saved_dir,
                                            "Qt Script files (*.js)");
  }
  std::ofstream os(filename.toUtf8(), std::ofstream::binary);
  if(!os)
    return;
  CGAL::Three::Three::CursorScopeGuard cs(Qt::WaitCursor);
  std::vector<std::pair<QString, QString> > names;
  std::vector<std::pair<QString, QString> > loaders;
  std::vector<QColor> colors;
  std::vector<int> rendering_modes;
  QStringList not_saved;
  Polyhedron_demo_io_plugin_interface* camera_plugin = nullptr;
  for(int i = 0; i < scene->numberOfEntries(); ++i)
  {
    Scene_item* item = scene->item(i);
    QString loader;// = item->property("loader_name").toString();
    QString ext;
    for(Polyhedron_demo_io_plugin_interface* iop : io_plugins)
    {
      if(iop->name() == "camera_positions_plugin")
        camera_plugin = iop;
      if(iop->isDefaultLoader(item))
      {
        QString sf = iop->saveNameFilters().split(";;").first();
        //OFF Files (*.off)
        QRegularExpression re("\\(\\*\\.(.*)\\)");
        QRegularExpressionMatch rem = re.match(sf);
        if(!rem.hasMatch())
          continue;
        ext = rem.captured(1);
        QList<Scene_item*>to_save;
        to_save.append(item);
        QString savename(tr("%1.%2").arg(item->name()).arg(ext));
        QString fullpath = make_fullpath(savename, true);
        savename = QFileInfo(fullpath).fileName();
        iop->save(QFileInfo(fullpath), to_save);
        names.push_back(std::make_pair(savename, item->name()));
        loader=iop->name();
        break;
      }
    }
    if(loader.isEmpty())
    {
      QMessageBox::warning(this, "", tr("No plugin found for %1. Not saved.").arg(item->name()));
      continue;
    }
    loaders.push_back(std::make_pair(loader, ext));
    colors.push_back(item->color());
    rendering_modes.push_back(item->renderingMode());
  }
  bool has_camera_positions = false;
  if(camera_plugin)
  {
    QString fullpath = make_fullpath("camera_tmp.camera.txt");
    QList<Scene_item*> dummy;
    if(camera_plugin->save(QFileInfo(fullpath), dummy))
    {
      QByteArray item = file_to_string(fullpath.toStdString().c_str());
      os << "var camera_positions= [\'";
      os<<qCompress(item, 9).toBase64().toStdString().c_str();
      os << "\']\n" ;
      //delete temp file
      QFile tmp_file(fullpath);
      tmp_file.remove();
      has_camera_positions =true;
    }
  }
  if(loaders.empty())
    return;
  //path
  os << "var camera = \""<<viewer->dumpCameraCoordinates().toStdString()<<"\";\n";
  os << "var items = [";
  for(std::size_t i = 0; i< names.size() -1; ++i)
  {
    QString fullpath = make_fullpath(names[i].first);

    QByteArray item = file_to_string(fullpath.toStdString().c_str());
    os<<"[\'";
    os<<qCompress(item, 9).toBase64().toStdString().c_str();
    os << "\', \'"<<names[i].second.toStdString().c_str()<<"\']," ;
    //delete temp file
    QFile tmp_file(fullpath);
    tmp_file.remove();
  }
  QString fullpath = make_fullpath(names.back().first);
  QByteArray item = file_to_string(fullpath.toStdString().c_str());
  os<<"[\'";
  os<<qCompress(item, 9).toBase64().toStdString().c_str();
  os << "\', \'"<<names.back().second.toStdString().c_str()<<"\']];\n";
  //delete temp file
  QFile tmp_file(fullpath);
  tmp_file.remove();
  //plugin
  os << "var loaders = [";
  for(std::size_t i = 0; i< names.size() -1; ++i)
  {
    os << "[\'" << loaders[i].first.toStdString() << "\', \'"<<loaders[i].second.toStdString()<< "\'],";
  }
  os << "[\'" << loaders.back().first.toStdString() << "\', \'"<<loaders.back().second.toStdString()<< "\']];\n";

  //color
  os << "var colors = [";
  for(std::size_t i = 0; i< names.size() -1; ++i)
  {
    os << "[" << colors[i].red() <<", "<< colors[i].green() <<", "<< colors[i].blue() <<"], ";
  }
  os<<"[" << colors.back().red() <<", "<< colors.back().green() <<", "<< colors.back().blue() <<"]];\n";

  //rendering mode
  os << "var rendering_modes = [";
  for(std::size_t i = 0; i< names.size() -1; ++i)
  {
    os << rendering_modes[i] << ", ";
  }
  os << rendering_modes.back()<<"];\n";
  os << "items.forEach(function(item, index, array){\n";
  os<<"          var path=items[index][1];\n";
  os<<"          path+='.';\n";
  os<<"          path+=loaders[index][1];\n";
  os<<"          var fullpath = main_window.write_string_to_file(item[0], path);\n";
  os<<"          main_window.open(fullpath,loaders[index][0]);\n";
  os << "        var it = scene.item(scene.numberOfEntries-1);\n";
  os << "        var r = colors[index][0];\n";
  os << "        var g = colors[index][1];\n";
  os << "        var b = colors[index][2];\n";
  os << "        it.setRgbColor(r,g,b);\n";
  os << "        it.setRenderingMode(rendering_modes[index]);\n";
  os << "});\n";
  os << "viewer.moveCameraToCoordinates(camera, 0.05);\n";
  if(has_camera_positions)
  {
    os<<"  var path=\"cams.camera.txt\";\n";
    os<<"  var fullpath = main_window.write_string_to_file(camera_positions, path);\n";
    os<<"  main_window.open(fullpath,\'camera_positions_plugin\');\n";
  }
  os.close();
  if(!not_saved.empty())
    QMessageBox::warning(this,
                         "Items Not  Saved",
                         QString("The following items could not be saved: %1").arg(
                           not_saved.join(", ")));
#ifdef CGAL_USE_SSH
  using namespace CGAL::ssh_internal;
  if(do_upload)
  {
    QString server = settings.value("ssh_server", QString()).toString();
    QString pk = settings.value("ssh_public_key", QString()).toString();
    QString privK = settings.value("ssh_priv_key", QString()).toString();
    user = user.trimmed();
    server = server.trimmed();
    pk = pk.trimmed();
    privK=privK.trimmed();
    if(user.isEmpty()){
      return;
    }
    QString path;
    path = QInputDialog::getText(this,
                                 "",
                                 tr("Enter the name of your scene file."));
    if(path.isEmpty())
      return;
    if(!path.contains("Polyhedron_demo_"))
      path.prepend("Polyhedron_demo_");
    try{
      ssh_session session = NULL;
      bool res = establish_ssh_session_from_agent(session,
                                                  user.toStdString().c_str(),
                                                  server.toStdString().c_str(),
                                                  pk.toStdString().c_str());

      if(!res)
      {
        bool ok;
        QString pass;
        pass = QInputDialog::getText(this, "SSH Password",
                                     "Enter ssh key password:",
                                     QLineEdit::Password,
                                     tr(""),
                                     &ok);
        if(!ok)
        {
          ssh_free(session);
          return;
        }
        pass = pass.trimmed();
        res = establish_ssh_session(session,
                                    user.toStdString().c_str(),
                                    server.toStdString().c_str(),
                                    pk.toStdString().c_str(),
                                    privK.toStdString().c_str(),
                                    pass.toStdString().c_str());
      }

      if(!res)
      {
        QMessageBox::warning(this,
                             "Error",
                             "The SSH session could not be started.");
        ssh_free(session);
        return;
      }
      res = push_file(session,path.toStdString().c_str(), filename.toStdString().c_str());
      if(!res)
      {
        QMessageBox::warning(this,
                             "Error",
                             "The file could not be uploaded. Check your console for more information.");
        close_connection(session);
        ssh_free(session);
        return;
      }
      close_connection(session);
      ssh_free(session);
      QFile tmp_file(filename);
      tmp_file.remove();
    } catch( ssh::SshException e )
    {
      std::cout << "Error during connection : ";
      std::cout << e.getError() << std::endl;
    }
  }
#endif
}
void MainWindow::setTransparencyPasses(int val)
{
  viewer->setTotalPass(val);
  viewer->update();
}

void MainWindow::toggleFullScreen()
{
  QList<QDockWidget *> dockWidgets = findChildren<QDockWidget *>();
  if(visibleDockWidgets.isEmpty())
  {
    Q_FOREACH(QDockWidget * dock, dockWidgets)
    {
      if(dock->isVisible())
      {
        visibleDockWidgets.append(dock);
        dock->hide();
      }
    }
  }
  else
  {
    Q_FOREACH(QDockWidget * dock, visibleDockWidgets){
      dock->show();
    }
    visibleDockWidgets.clear();

  }
}

void MainWindow::setDefaultSaveDir()
{
  QString dirpath = QFileDialog::getExistingDirectory(this, "Set Default Save as Directory", def_save_dir);
  if(!dirpath.isEmpty())
    def_save_dir = dirpath;
  settings.setValue("default_saveas_dir", def_save_dir);
}


void MainWindow::setupViewer(Viewer* viewer, SubViewer* subviewer)
{
  // do not save the state of the viewer (anoying)
  viewer->setStateFileName(QString());
  viewer->textRenderer()->setScene(scene);
  viewer->setScene(scene);
  connect(scene, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex & )),
          viewer, SLOT(update()));
  connect(scene, SIGNAL(updated()),
          viewer, SLOT(update()));

  QAction* action = subviewer->findChild<QAction*>("actionRecenter");
  connect(action, SIGNAL(triggered()),
          viewer, SLOT(update()));
  connect(action, &QAction::triggered,
          subviewer, &SubViewer::recenter);
  action= subviewer->findChild<QAction*>("actionLookat");
  connect(action, SIGNAL(triggered()),
          subviewer, SLOT(lookat()));
  action= subviewer->findChild<QAction*>("actionColor");
  connect(action, &QAction::triggered,
          subviewer, &SubViewer::color);
  action= subviewer->findChild<QAction*>("actionDumpCamera");
  connect(action, &QAction::triggered,
          [this, viewer](){
    information(QString("Camera: %1")
                .arg(cameraString(viewer)));
  });
  action= subviewer->findChild<QAction*>("actionCopyCamera");
  connect(action, &QAction::triggered,
          [this, viewer](){
    qApp->clipboard()->setText(cameraString(viewer));
  });
  action= subviewer->findChild<QAction*>("actionPasteCamera");
  connect(action, &QAction::triggered,
          this, [viewer](){
    QString s = qApp->clipboard()->text();
    viewer->moveCameraToCoordinates(s, 0.5f);
  });
  action= subviewer->findChild<QAction*>("actionAntiAliasing");
  connect(action, SIGNAL(toggled(bool)),
          viewer, SLOT(setAntiAliasing(bool)));
  action= subviewer->findChild<QAction*>("actionDrawTwoSide");
  connect(action, SIGNAL(toggled(bool)),
          viewer, SLOT(setTwoSides(bool)));
  action= subviewer->findChild<QAction*>("actionQuick");
  connect(action, SIGNAL(toggled(bool)),
          viewer, SLOT(setFastDrawing(bool)));
  action= subviewer->findChild<QAction*>("actionOrtho");
  connect(action, SIGNAL(toggled(bool)),
          viewer, SLOT(SetOrthoProjection(bool)));
  action= subviewer->findChild<QAction*>("actionTotalPass");
  connect(action, &QAction::triggered,
          this, [this, viewer]() {
    bool ok;
    int nb = QInputDialog::getInt(this, "Set Maximum Number of Passes",
                                  "Enter number of transparency passes:",
                                  4, 4, 99, 1, &ok);
    if (!ok){
      return;
    }
    viewer->setTotalPass(nb);
  });

  action = subviewer->findChild<QAction*>("actionScaleScene");
  action->setCheckable(true);
  action->setChecked(false);
  connect(action, &QAction::triggered,
          viewer, &Viewer::scaleScene);

  action= subviewer->findChild<QAction*>("actionBackFrontShading");
  connect(action, SIGNAL(toggled(bool)),
          viewer, SLOT(setBackFrontShading(bool)));
  connect(viewer, SIGNAL(requestContextMenu(QPoint)),
          this, SLOT(contextMenuRequested(QPoint)));
  connect(viewer, SIGNAL(selected(int)),
          this, SLOT(selectSceneItem(int)));
  connect(viewer, SIGNAL(selectedPoint(double, double, double)),
          this, SLOT(showSelectedPoint(double, double, double)));

  connect(viewer, SIGNAL(selectionRay(double, double, double,
                                      double, double, double)),
          scene, SIGNAL(selectionRay(double, double, double,
                                     double, double, double)));

  connect(viewer, &Viewer::sendMessage,
          this, [](QString s){
    information(s);
  });

#ifdef CGAL_USE_WEBSOCKETS
  action= subviewer->viewer->findChild<QAction*>("actionShareCamera");
  connect(action, &QAction::toggled,
          this, [this, viewer](bool b)
  {
    if(!viewer){
      return;
    }
    QString session;
    if(b){
      bool ok;
      session = QInputDialog::getText(
            this,"Session",
            "Please enter the session name.\n"
            "Only the machines that enter the same session name will be connected.\n"
            "Several sessions can run simultaneously on a same server. ",
            QLineEdit::Normal, QString(), &ok);
      if(session.isEmpty() || !ok)
      {
        viewer->setShareCam(false, session);
        return;
      }
    }
    viewer->setShareCam(b, session);
  });
#endif

}

void MainWindow::on_actionAdd_Viewer_triggered()
{
  SubViewer* subviewer = new SubViewer(ui->mdiArea, this, viewer);
  Viewer* viewer2 = subviewer->viewer;
  viewer2->setManipulatedFrame(viewer->manipulatedFrame());
  CGAL::qglviewer::Vec min, max;
  computeViewerBBox(min, max);
  updateViewerBbox(viewer2, true, min, max);
  viewer2->setObjectName("viewer2");
  connect(viewer2, SIGNAL(doneInitGL(CGAL::Three::Viewer_interface*)),
          scene, SLOT(newViewer(CGAL::Three::Viewer_interface*)));
  connect(viewer2, &Viewer::contextIsDestroyed,
          this, [this, viewer2](){
    scene->removeViewer(viewer2);
    viewerDestroyed(viewer2);
  });

  setupViewer(viewer2, subviewer);
  viewer2->camera()->interpolateToFitScene();
  subviewer->show();
  ui->mdiArea->tileSubWindows();
  QPoint pos = viewer_window->pos();
  QSize size = viewer_window->size();
  viewer_window->move(subviewer->pos());
  viewer_window->resize(subviewer->size());
  subviewer->move(pos);
  subviewer->resize(size);
  newViewerCreated(viewer2);
}

void MainWindow::recenterViewer()
{
  scene->computeBbox();
  CGAL::qglviewer::Vec min, max;
  computeViewerBBox(min, max);
  Viewer* target = qobject_cast<Viewer*>(childAt(cursor().pos()));
  if(target)
  {
    scene->computeBbox();
    updateViewerBbox(target, true, min, max);
    target->camera()->interpolateToFitScene();
  }
}

void MainWindow::updateViewerBbox(Viewer *vi, bool recenter,
                                  CGAL::qglviewer::Vec min,
                                  CGAL::qglviewer::Vec max){
  CGAL::qglviewer::Vec center = viewer->camera()->pivotPoint();
  vi->setSceneBoundingBox(min,
                          max);
  if(recenter)
  {
    vi->resetFov();
    vi->camera()->showEntireScene();
  }
  else
  {
    vi->camera()->setPivotPoint(center);
  }
}

QObject* MainWindow::getDirectChild(QObject* widget)
{

  if(!widget->property("helpText").toString().isEmpty())
    return widget;
  return getDirectChild(widget->parent());
}

void MainWindow::on_action_Rearrange_Viewers_triggered()
{
  if(ui->mdiArea->subWindowList().size() == 1)
    ui->mdiArea->subWindowList().first()->showMaximized();
  else
  {
    ui->mdiArea->tileSubWindows();
    QMdiSubWindow* subviewer = qobject_cast<QMdiSubWindow*>(
          ui->mdiArea->childAt(ui->mdiArea->pos()));
    if(!subviewer)//should not happen but better safe than sorry
    {
      return;
    }
    QPoint pos = viewer_window->pos();
    QSize size = viewer_window->size();
    viewer_window->move(subviewer->pos());
    viewer_window->resize(subviewer->size());
    subviewer->move(pos);
    subviewer->resize(size);
  }
}

SubViewer::SubViewer(QWidget *parent, MainWindow* mw, Viewer* mainviewer)
  :QMdiSubWindow (parent),
    mw(mw),
    viewMenu(new QMenu(this)),
    is_main(false)
{
  if(mainviewer)
    viewer = new Viewer(this, mainviewer);
  else
  {
    viewer = new Viewer(this);
    is_main = true;
  }
  setWidget(viewer);
  QAction* actionRecenter = new QAction("Re&center Scene",this);
  actionRecenter->setObjectName("actionRecenter");
  viewMenu->addAction(actionRecenter);
  QAction* actionLookat = new QAction("&Look at...",this);
  actionLookat->setObjectName("actionLookat");
  viewMenu->addAction(actionLookat);
  QAction* actionColor = new QAction("Change &Background Color...",this);
  actionColor->setObjectName("actionColor");
  viewMenu->addAction(actionColor);
  QAction* actionDumpCamera = new QAction("&Dump Camera Coordinates",this);
  actionDumpCamera->setObjectName("actionDumpCamera");
  QAction* actionCopyCamera = new QAction("&Copy Camera",this);
  actionCopyCamera->setObjectName("actionCopyCamera");
  QAction* actionPasteCamera = new QAction("&Paste Camera",this);
  actionPasteCamera->setObjectName("actionPasteCamera");
  QMenu* cameraMenu = new QMenu("Ca&mera", mw);
  cameraMenu->addAction(actionDumpCamera);
  cameraMenu->addAction(actionCopyCamera);
  cameraMenu->addAction(actionPasteCamera);
  viewMenu->addMenu(cameraMenu);

  QAction* actionAntiAliasing = new QAction("&Antialiasing",this);
  actionAntiAliasing->setObjectName("actionAntiAliasing");
  actionAntiAliasing->setCheckable(true);
  actionAntiAliasing->setChecked(false);
  viewMenu->addAction(actionAntiAliasing);
  QAction* actionDrawTwoSide = new QAction("Draw &Two Sides",this);
  actionDrawTwoSide->setObjectName("actionDrawTwoSide");
  actionDrawTwoSide->setCheckable(true);
  actionDrawTwoSide->setChecked(false);
  viewMenu->addAction(actionDrawTwoSide);
  QAction* actionQuick = new QAction("&Quick Camera Mode",this);
  actionQuick->setObjectName("actionQuick");
  actionQuick->setCheckable(true);
  actionQuick->setChecked(true);
  viewMenu->addAction(actionQuick);
  QAction* actionOrtho = new QAction("&Orthographic Projection",this);
  actionOrtho->setObjectName("actionOrtho");
  actionOrtho->setCheckable(true);
  actionOrtho->setChecked(false);
  viewMenu->addAction(actionOrtho);
  QAction* actionTotalPass = new QAction("Set Transparency Pass &Number...",this);
  actionTotalPass->setObjectName("actionTotalPass");
  viewMenu->addAction(actionTotalPass);
#ifdef CGAL_USE_WEBSOCKETS
  QAction* actionShareCamera= new QAction("Join &WS Server",viewer);
  actionShareCamera->setObjectName("actionShareCamera");
  actionShareCamera->setCheckable(true);
  actionShareCamera->setChecked(false);
  viewMenu->addAction(actionShareCamera);
#endif

  QAction* actionBackFrontShading = new QAction("Activate Back/Front shading.",this);
  actionBackFrontShading->setObjectName("actionBackFrontShading");
  actionBackFrontShading->setCheckable(true);
  actionBackFrontShading->setChecked(false);
  viewMenu->addAction(actionBackFrontShading);

  QAction* actionScaleScene = new QAction("&Scale the Scene...",this);
  actionScaleScene->setObjectName("actionScaleScene");
  viewMenu->addAction(actionScaleScene);

  if(mainviewer)
    setAttribute(Qt::WA_DeleteOnClose);
  setWindowIcon(QIcon(":/cgal/icons/resources/menu.png"));
  setSystemMenu(viewMenu);
}

SubViewer::~SubViewer()
{
  viewer->deleteLater();
}

void SubViewer::recenter()
{
  CGAL::qglviewer::Vec min, max;
  mw->computeViewerBBox(min, max);
  mw->updateViewerBbox(viewer, true, min, max);
  viewer->camera()->interpolateToFitScene();
}

void SubViewer::lookat()
{
  Show_point_dialog dialog(mw);
  dialog.setWindowTitle(tr("Look at..."));
  int i = dialog.exec();
  if( i == QDialog::Accepted &&
      dialog.has_correct_coordinates() )
  {
    if (viewer->camera()->frame()->isSpinning())
      viewer->camera()->frame()->stopSpinning();
    mw->viewerShow(viewer,
                   (float)dialog.get_x() + viewer->offset().x,
                   (float)dialog.get_y() + viewer->offset().y,
                   (float)dialog.get_z() + viewer->offset().z);
  }
}

void SubViewer::color()
{
  QColor c =  QColorDialog::getColor();
  if(c.isValid()) {
    viewer->setBackgroundColor(c);
    viewer->update();
  }
}

void SubViewer::closeEvent(QCloseEvent *closeEvent)
{

  if(is_main)
  {
    QMessageBox::information(mw, "", "This is the main viewer. It cannot be closed.");
    closeEvent->ignore();
  }
  else
    QWidget::closeEvent(closeEvent);
}
void SubViewer::changeEvent(QEvent *event)
{
  QMdiSubWindow::changeEvent(event);
  if(event->type() == QEvent::WindowStateChange)
  {
    if(isMaximized())
    {
      QMenu* menu = mw->findChild<QMenu*>("menuView");
      Q_FOREACH(QAction* action, viewMenu->actions())
      {
        menu->addAction(action);
      }
      setWindowFlags(
              Qt::SubWindow
              | Qt::CustomizeWindowHint
              | Qt::WindowMaximizeButtonHint
              //| Qt::WindowSystemMenuHint
              | Qt::WindowTitleHint
              );
      QAction* action = mw->findChild<QAction*>("action_Rearrange_Viewers");
      action->setVisible(false);
      viewer->update();
    }
    else
    {
      QMenu* menu = mw->findChild<QMenu*>("menuView");
      Q_FOREACH(QAction* action, viewMenu->actions())
      {
        menu->removeAction(action);
      }
      setWindowFlags(
              Qt::SubWindow
              | Qt::CustomizeWindowHint
              | Qt::WindowMaximizeButtonHint
              | Qt::WindowSystemMenuHint
              | Qt::WindowTitleHint
              );
      QAction* action = mw->findChild<QAction*>("action_Rearrange_Viewers");
      action->setVisible(true);
      for(auto v : CGAL::QGLViewer::QGLViewerPool())
        v->update();
    }
  }
}

void MainWindow::invalidate_bbox(bool do_recenter)
{
  bbox_need_update = true;
  if(do_recenter)
    updateViewersBboxes(true);
}

void MainWindow::on_action_Save_triggered()
{
  if(QMessageBox::question(this, "Save", "Are you sure you want to override these files ?")
     == QMessageBox::No)
    return;
  QList<Scene_item*> to_save;

  for(Scene::Item_id id : scene->selectionIndices())
  {
    Scene_item* item = scene->item(id);
    if(!item->property("source filename").toString().isEmpty())
    {
      QString filename = item->property("source filename").toString();
      to_save.append(item);
      save(filename, to_save);
    }
  }
}

void MainWindow::on_actionLoad_a_Scene_from_a_Script_File_triggered()
{
  bool do_download = false;
  QString filename;

#ifdef CGAL_USE_SSH
  QString user = settings.value("ssh_user", QString()).toString();

  if(!user.isEmpty())
  {
    QMessageBox::StandardButton doyou =
        QMessageBox::question(this, tr("Download ?"), tr("Do you wish to download the scene"
                                                         " using the SSH preferences ?"));
    do_download= (doyou == QMessageBox::Yes);
  }
#endif

  if(do_download)
  {
    #ifdef CGAL_USE_SSH
    using namespace CGAL::ssh_internal;
    QString server = settings.value("ssh_server", QString()).toString();
    QString pk = settings.value("ssh_public_key", QString()).toString();
    QString privK = settings.value("ssh_priv_key", QString()).toString();
    user = user.trimmed();
    server = server.trimmed();
    pk = pk.trimmed();
    privK=privK.trimmed();

    try{
      ssh_session session = nullptr;
      bool res = establish_ssh_session_from_agent(session,
                                                  user.toStdString().c_str(),
                                                  server.toStdString().c_str(),
                                                  pk.toStdString().c_str());
      if(!res){
        bool ok;
        QString pass= QInputDialog::getText(this, "SSH Password",
                                     "Enter ssh key password:",
                                     QLineEdit::Password,
                                     tr(""),
                                     &ok);
        if(!ok)
        {
          ssh_free(session);
          return;
        }
        pass = pass.trimmed();
        res = establish_ssh_session(session,
                                    user.toStdString().c_str(),
                                    server.toStdString().c_str(),
                                    pk.toStdString().c_str(),
                                    privK.toStdString().c_str(),
                                    pass.toStdString().c_str());
      }
      if(!res)
      {
        QMessageBox::warning(this,
                             "Error",
                             "The SSH session could not be started.");
        ssh_free(session);
        return;
      }
      QStringList names;
      if(!CGAL::ssh_internal::explore_the_galaxy(session, names))
      {
        QMessageBox::warning(this,
                             "Error",
                             "Could not find remote directory.");
      }
      QString path;
      path = QInputDialog::getItem(this,
                                   "Choose a file",
                                   tr("Choose the scene file."),
                                   names);
      filename = QString("%1/load_scene.js").arg(QDir::tempPath());
      if(path.isEmpty())
      {
        ssh_free(session);
        return;
      }
      path.prepend("Polyhedron_demo_");
      path = tr("/tmp/%2").arg(path);
      res = pull_file(session,path.toStdString().c_str(), filename.toStdString().c_str());
      if(!res)
      {
        QMessageBox::warning(this,
                             "Error",
                             "The file could not be fetched. Check your console for more info.");
        close_connection(session);
        ssh_free(session);
        return;
      }
      close_connection(session);
      ssh_free(session);
    } catch( ssh::SshException e )
    {
      std::cout << "Error during connection : ";
      std::cout << e.getError() << std::endl;
    }
    #endif
  }
  else
  {
    filename =  QFileDialog::getOpenFileName(
          this,
          tr("Select a Whole Scene file..."),
          ".",
          "Whole Scene files (*.js)");
    if(filename.isEmpty())
      return;
  }
  loadScript(QFileInfo(filename));
  if(do_download){
    QFile tmp_file(filename);
    tmp_file.remove();
  }
}
