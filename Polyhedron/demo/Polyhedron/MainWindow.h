#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include "config.h"

#include <QtOpenGL/qgl.h>
#include <CGAL/Qt/DemosMainWindow.h>
#ifdef QT_SCRIPT_LIB
#  include  <QScriptEngine>
#endif

#include <QVector>
#include <QList>
#include <QFileInfo>
#include <QStringList>
#include <QSet>

class Scene;
class Viewer;
class QTreeView;
class QMenu;
namespace CGAL {
namespace Three{
class Polyhedron_demo_io_plugin_interface;
class Polyhedron_demo_plugin_interface;
}
}

class Scene_item;
class QSortFilterProxyModel;

namespace Ui {
  class MainWindow;
  class Add_polylines_dialog;
  class Add_point_set_dialog;
}

#include "Polyhedron_type_fwd.h"

#include "Messages_interface.h"

/*!
 * \brief The main window of the applicatipon.
 * It contains all the widgets, the menus and the elements of interface
 * of the application.*/

class MainWindow : 
  public CGAL::Qt::DemosMainWindow,
  public Messages_interface
{
  Q_OBJECT
  Q_INTERFACES(Messages_interface)
public:
  /*! \brief The constructor
   * It links the class with its UI file and sets it up.
   * It also saves pointers to viewer and the sceneView.
   * Then it creates and initializes the scene and do the
   * connexions with the UI. Finally it loads the plugins.*/

  MainWindow(QWidget* parent = 0);
  ~MainWindow();

  /*! Find an IO plugin.
   * throws std::invalid_argument if no loader with that argument can be found
   @returns the IO plugin associated with `loader_name`*/
  CGAL::Three::Polyhedron_demo_io_plugin_interface* find_loader(const QString& loader_name) const;
  
  /*! \brief Load an item with a given loader.
   * throws `std::logic_error` if loading does not succeed or
   * `std::invalid_argument` if `fileinfo` specifies an invalid file*/
  Scene_item* load_item(QFileInfo fileinfo, CGAL::Three::Polyhedron_demo_io_plugin_interface*);

public Q_SLOTS:
  void updateViewerBBox();
  void open(QString);

  //! given an extension file, returns true if `filename` matches the filter
  bool file_matches_filter(const QString& filters, const QString& filename);

  /*! Open a file with a given loader, and return true if it was successful.
   This slot is for use by scripts.*/
  bool open(QString filename, QString loader_name);

  /*! Reloads an item. Expects to be called by a QAction with the
   index of the item to be reloaded as data attached to the action.
   The index must identify a valid `Scene_item`.*/
  void reload_item();
  
  /*!
   * This is an overloaded function.
   * If QT_SCRIPT_LIB is defined, returns true if the script is valid.
   * If not, returns false.
   */
  bool load_script(QString filename);

  /*! If QT_SCRIPT_LIB is defined, returns true if the script is valid.
  * If not, returns false.
  */
  bool load_script(QFileInfo);

  /*!
   * Gives the keyboard input focus to the widget searchEdit.
   */
  void setFocusToQuickSearch();

  /*!
   * Clears the current selection and select the scene_item with
   * the index i in the sceneView. Calls itemChanged(i) from the scene.
   */
  void selectSceneItem(int i);
  /*!
   * Prints coordinates of a point and its distance to the last
   * position printed by this function.
   */
  void showSelectedPoint(double, double, double);
  /*!
   * Calls removeSceneItemFromSelection(i).
   */
  void unSelectSceneItem(int i);
  /*!
   * Clears the current selection and select all the scene_items
   * in the sceneView.
   */
  void selectAll();
  /*!
   * Adds the scene_item with the index i in the sceneView to the
   * current selection. Calls itemChanged(i) from the scene.
   */
  void addSceneItemInSelection(int i);

  /*!
   * Removes the scene_item with the index i in the sceneView to the
   * current selection. Calls itemChanged(i) from the scene.
   */
  void removeSceneItemFromSelection(int i);

  /*!
   * Calls setAddKeyFrameKeyboardModifiers(m) from the viewer.
   */
  void setAddKeyFrameKeyboardModifiers(Qt::KeyboardModifiers);
  /*!
   * Set the visibility of all actions from the targeted menu
   * and its submenus to false and disables them.
   */
  void clearMenu(QMenu*);
  /*!
   * Enables and sets the targeted action's visibility to true.
   * Does the same in all the menus of the widgets associated with
   * this action.
   */
  void addAction(QAction*);
  /*!
   * Creates an action with text actionText named actionName and
   * adds it to the menu menuName. If menuName does not exist,
   * it is created.
   */
  void addAction(QString actionName,
                 QString actionText,
                 QString menuName);
  /*!
   * Sets the scene center to the target position and makes the
   * scene slide to this new center. Also sets the pivotPoint of
   * the camera to this position.
   */
  void viewerShow(float, float, float);
  /*!
   * Sets the scene center to be the center of the target BBox.
   * Also sets the pivotPoint of the camera to this position.
   */
  void viewerShow(float, float, float, float, float, float);
  /*!
   * Centers the scene on the target object.
   */
  void viewerShowObject();
  /*!
   * Displays a text preceded by the mention "INFO :".
   */
  void information(QString);
  /*!
   * Displays a blue text preceded by the mention "WARNING :".
   */
  void warning(QString);
  /*!
   * Displays a red text preceded by the mention "ERROR :".
   */
  void error(QString);

    //!Displays a text in the chosen html color with the chosen html font.

  void message(QString, QString, QString = QString("normal"));

    //!Returns true if the target plugin is present. If not, returns false.
  bool hasPlugin(const QString&) const;
  /*!
   * If able, finds a script debugger and interrupts current action. Default
   * value for parameter is true.
   */
  void enableScriptDebugger(bool = true);

protected Q_SLOTS:

   //!Gets the new selected item(s) from the sceneView and updates the scene
   //!and viewer accordingly.
  /*!
   * Set the scene selected item or selected item list. Sets the manipulated
   * frame of the viewer to the new selected item's and calls updateGL().
   */
  void selectionChanged();
  /*!
   * Invoques a context menu for the currently selected item at the requested
   * position.
   */
  void contextMenuRequested(const QPoint& global_pos);
  /*!
   * Invoques a context menu for the requested item at the requested
   * position.
   */
  void showSceneContextMenu(int selectedItemIndex,
                            const QPoint& global_pos);

  //!This is an overloaded function. Invoques a context menu at the requested
  //!position.
  /*!
   * If the widget which received the request is not the sceneView, the index
   * chosen by default for the menu is the one of the currently selected item.
   * If it is the sceneView, then the index of the clicked item is collected.
   * If this index is valid, then it is used for the menu. If not, the function
   * returns.
   */
  void showSceneContextMenu(const QPoint& local_pos_of_treeview);

  //!Prints information about the currently selected item if able.
  void updateInfo();
  //!Prints graphical information about the currently selected item if able.
  void updateDisplayInfo();
  //!Sets the current manipulated frame to 0.
  void removeManipulatedFrame(Scene_item*);

  // settings
  //!Closes the main window.
  void quit();
  //!Reads the plugin_blacklist contents and apply them.
  void readSettings();
  //!Sets up the plugin_blacklist.
  void writeSettings();

  // load, erase, duplicate
  //!Loops on on_actionErase_triggered();
  void on_actionEraseAll_triggered();
  //!Opens a dialog to open one or several files.
  void on_actionLoad_triggered();
  //!Erases the selected items. Returns true if items remain in the sceneView.
  bool on_actionErase_triggered();
  //!Duplicates the selected item and selects the new item.
  void on_actionDuplicate_triggered();
  //!If QT_SCRIPT_LIB is defined, opens a dialog to choose a script.
  void on_actionLoad_Script_triggered();

  // Show/Hide
  //!Swap the visibility of the selected item(s).
  void on_actionShowHide_triggered();

  // Select A/B
  //!Sets the selected item as item_A.
  void on_actionSetPolyhedronA_triggered();
  //!Sets the selected item as Item_B.
  void on_actionSetPolyhedronB_triggered();

  //Preferences edition
  //!Opens the Preferences dialog.
  void on_actionPreferences_triggered();
  //!Opens a dialog to add polylines on the fly.
  void on_actionAdd_polylines_triggered();
  //!Opens a dialog to add a point set on the fly.
  void on_actionAdd_point_set_triggered();
  //!Adds a polyline
  void addPolylineButton_clicked();
  //!Adds a point set
  void addPointSetButton_clicked();
  //!Closes the dialog
  void closePolylinesButton_clicked();
  //!Closes the dialog
  void closePointSetButton_clicked();

  // save as...
  //!Opens a dialog to save selected item if able.
  void on_actionSaveAs_triggered(); 
  //!Calls the function save of the current plugin if able.
  void save(QString filename, Scene_item* item);
  //!Calls the function saveSnapShot of the viewer.
  void on_actionSaveSnapshot_triggered();
  //!Opens a Dialog to choose a color and make it the background color.
  void on_actionSetBackgroundColor_triggered();
  /*! Opens a Dialog to enter coordinates of the new center point and sets it
   * with viewerShow.
   *@see viewerShow(float, float, float, float, float, float)
   */
  void on_action_Look_at_triggered();
  //!Returns the position and orientation of the current camera frame.
  QString camera_string() const;
  /*! Prints the position and orientation of the current camera frame.
   * @see camera_string()
   */
  void on_actionDumpCamera_triggered();
  //!Sets the coordinates of the camera in the clipboard text.
  void on_action_Copy_camera_triggered();
  //!Gets coordinates from the clipboard and sets them as the current ones for
  //!the camera.
  void on_action_Paste_camera_triggered();
  //!Hides not available operations and show available operations in all the
  //!menus.
  void filterOperations();
  //!Updates the bounding box and moves the camera to fits the scene.
  void on_actionRecenterScene_triggered();
protected:
  QList<QAction*> createSubMenus(QList<QAction*>);
  /*! For each objects in the sceneView, loads the associated plugins.
   * Gets the property "submenuName" of all the actions and creates submenus.
   * Sorts the Operations menu by name.
   * @see initPlugin(QObject*);
   * @see initIOPlugin(QObject*);
   */
  void loadPlugins();
  /*!
   * \brief Initializes the plugins.
   * Makes pairs between plugins and object names and fills the Operations menu.
   * Called only once.
   */
  bool initPlugin(QObject*);
  //!Initializes the IO plugin for the target object.
  bool initIOPlugin(QObject*);
  /*!
   * Calls writeSettings() and set the flag accepted for the event.
   * @see writeSettings()
   */
  void closeEvent(QCloseEvent *event);
  /*! Returns the currently selected item in the sceneView. Returns -1
   * if none is selected.
   */
  int getSelectedSceneItemIndex() const;
  //! Returns a list of the selected items in the sceneView.
  QList<int> getSelectedSceneItemIndices() const;

private:
  QString strippedName(const QString &fullFileName);
  void setMenus(QString, QString, QAction *a);
  /// plugin black-list
  QSet<QString> plugin_blacklist;
  Scene* scene;
  Viewer* viewer;
  QSortFilterProxyModel* proxyModel;
  QTreeView* sceneView;
  Ui::MainWindow* ui;
  Ui::Add_polylines_dialog *add_polydiagui;
  Ui::Add_point_set_dialog *add_pointsetdiagui;
  QDialog *add_polydiag;
  QDialog *add_pointsetdiag;
  int nb_of_polylines;
  int nb_of_point_set;
  QVector<CGAL::Three::Polyhedron_demo_io_plugin_interface*> io_plugins;
  QMap<QString,QString> default_plugin_selection;
  // typedef to make Q_FOREACH work
  typedef QPair<CGAL::Three::Polyhedron_demo_plugin_interface*, QString> PluginNamePair;
  QVector<PluginNamePair > plugins;
#ifdef QT_SCRIPT_LIB
  QScriptEngine* script_engine;
public:
  /*! Evaluates a script and search for uncaught exceptions. If quiet is false, prints the
   *backtrace of the uncaught exceptions.
   */
  void evaluate_script(QString script, 
                       const QString & fileName = QString(),
                       const bool quiet = false);
  //! Calls evaluate_script(script, filename, true).
  void evaluate_script_quiet(QString script, 
                             const QString & fileName = QString());
#endif
};

#endif // ifndef MAINWINDOW_H
