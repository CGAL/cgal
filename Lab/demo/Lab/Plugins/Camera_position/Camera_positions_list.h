#ifndef CAMERA_POSITIONS_LIST_H
#define CAMERA_POSITIONS_LIST_H

#include <QDockWidget>
#include <QModelIndex>

namespace CGAL{namespace Three{class Viewer_interface;}}
class QListView;
class QStandardItemModel;

class Camera_positions_list : public QDockWidget {
  Q_OBJECT
public:
  Camera_positions_list(QWidget* parent);

public Q_SLOTS:
  void load(QString filename);
  bool save(QString filename);
protected Q_SLOTS:
  void on_plusButton_pressed();
  void on_minusButton_pressed();
  void on_upButton_pressed();
  void on_downButton_pressed();
  void on_openButton_pressed();
  void on_saveButton_pressed();
  void on_clearButton_pressed();
  void on_frontButton_pressed();
  void on_backButton_pressed();
  void on_topButton_pressed();
  void on_botButton_pressed();
  void on_leftButton_pressed();
  void on_rightButton_pressed();
  void activatedRow(QModelIndex index);

protected:
  void addItem(QString, QString);

private:
  int counter;
  QListView* m_listView;
  QStandardItemModel* m_model;
};

#endif

