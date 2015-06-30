#ifndef CAMERA_POSITIONS_LIST_H
#define CAMERA_POSITIONS_LIST_H

#include <QDockWidget>
#include <QModelIndex>

class Viewer_interface;
class QListView;
class QStandardItemModel;

class Camera_positions_list : public QDockWidget {
  Q_OBJECT
public:  
  Camera_positions_list(QWidget* parent);

  void setViewer(Viewer_interface*);

public Q_SLOTS:
  void load(QString filename);
protected Q_SLOTS:
  void on_plusButton_pressed();
  void on_minusButton_pressed();
  void on_upButton_pressed();
  void on_downButton_pressed();
  void on_openButton_pressed();
  void on_saveButton_pressed();
  void on_clearButton_pressed();

  void activatedRow(QModelIndex index);

protected:
  void addItem(QString, QString);

private:
  Viewer_interface* m_viewer;
  int counter;
  QListView* m_listView;
  QStandardItemModel* m_model;
};

#endif
