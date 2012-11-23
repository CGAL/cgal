#include "Camera_positions_list.h"

#include "ui_Camera_positions_list.h"
#include <QListView>
#include <QStandardItemModel>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>

#include "Viewer_interface.h"

#include <cassert>

Camera_positions_list::Camera_positions_list(QWidget* parent)
  : QDockWidget(parent), m_viewer(0), counter(0), m_model(new QStandardItemModel(this))
{
  Ui::Camera_positions_list ui;
  ui.setupUi(this);
  m_listView = ui.listView;
  m_listView->setModel(m_model);
  
  m_listView->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed);
  connect(m_listView, SIGNAL(activated(QModelIndex)),
          this, SLOT(activatedRow(QModelIndex)));
}

void Camera_positions_list::setViewer(Viewer_interface* viewer)
{
  m_viewer = viewer;
}

void Camera_positions_list::on_plusButton_pressed()
{
  if(!m_viewer) return;
  addItem(tr("Camera position #%1").arg(++counter),
          m_viewer->dumpCameraCoordinates());
}

void Camera_positions_list::addItem(QString text, QString data)
{
  QStandardItem* item = new QStandardItem(text);
  item->setData(data, Qt::UserRole);
  m_model->insertRow(m_model->rowCount(), item);
}

void Camera_positions_list::on_upButton_pressed()
{
  int row = m_listView->selectionModel()->currentIndex().row();
  m_model->insertRow(row-1, m_model->takeRow(row));
  m_listView->selectionModel()->setCurrentIndex(m_model->index(row-1, 0),
                                                QItemSelectionModel::Clear);
}

void Camera_positions_list::on_downButton_pressed()
{
  int row = m_listView->selectionModel()->currentIndex().row();
  m_model->insertRow(row+1, m_model->takeRow(row));
  m_listView->selectionModel()->setCurrentIndex(m_model->index(row+1, 0),
                                                QItemSelectionModel::Clear);
}

void Camera_positions_list::on_minusButton_pressed()
{
  Q_FOREACH(QModelIndex index, 
            m_listView->selectionModel()->selectedIndexes()) {
    m_model->removeRows(index.row(), 1);
  }
}

void Camera_positions_list::on_clearButton_pressed()
{
  m_model->clear();
}

// void Camera_positions_list::editItem(QListWidgetItem* item)
// {
//   std::cerr << "is_editable: " << m_listView->flags(item)QListWidget

//   m_listView->editItem(item);
// }

void Camera_positions_list::activatedRow(QModelIndex index)
{
  QString s = m_model->data(index, Qt::UserRole).toString();
  if(s.isNull()) return;
  m_viewer->moveCameraToCoordinates(s);
}

void Camera_positions_list::on_saveButton_pressed()
{
  QString filename =
    QFileDialog::getSaveFileName(this, 
                                 tr("Save camera coordinates to file"),
                                 QString(),
                                 tr("(*.camera.txt)"));
  QFile file(filename);
  file.open(QIODevice::WriteOnly);
  QTextStream out(&file);
  for(int i = 0; i < m_model->rowCount(); ++i)
  {
    QStandardItem* item = m_model->item(i);
    out << item->data(Qt::DisplayRole).toString()
        << "\n"
        << item->data(Qt::UserRole).toString()
        << "\n";
  }
  file.close();
}

void Camera_positions_list::on_openButton_pressed()
{
  QString filename =
    QFileDialog::getOpenFileName(this, 
                                 tr("Read camera coordinates from file"),
                                 QString(),
                                 tr("(*.camera.txt)"));
  load(filename);
}

void Camera_positions_list::load(QString filename) {
  QFile file(filename);
  std::clog << "Loading camera positions " << qPrintable(filename) << std::endl;
  file.open(QIODevice::ReadOnly);
  QTextStream input(&file);
  while(!input.atEnd()) {
    QString text = input.readLine(1000);
    QString coord = input.readLine(1000);
    if(text.isNull() || coord.isNull()) return;
    qglviewer::Frame frame;
    if(m_viewer->readFrame(coord, frame))
    {
      addItem(text,
              m_viewer->dumpFrame(frame));
    }
  }
}

#include "Camera_positions_list.moc"
