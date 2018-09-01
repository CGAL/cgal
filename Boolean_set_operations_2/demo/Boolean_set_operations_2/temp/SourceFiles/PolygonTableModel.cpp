#define CGAL_PolygonTableOptionsNumber 1

#include "PolygonTableModel.h"

PolygonTableModel::PolygonTableModel(QObject *parent)
  : QAbstractTableModel(parent)
{
}

int PolygonTableModel::rowCount(const QModelIndex &parent) const
{
  return m_polygonList.count();
}
int PolygonTableModel::columnCount(const QModelIndex &parent) const
{
  return CGAL_PolygonTableOptionsNumber;
}

QVariant PolygonTableModel::data(const QModelIndex &index, int role) const 
{
  if (!index.isValid()) {
    return QVariant();
  }

  if (index.row() >= m_polygonList.size()) {
    return QVariant();
  }

  if (role == Qt::DisplayRole) {
    PolygonListItem* polygonItem = m_polygonList.at(index.row());
    switch (index.column())
    {
    case 0:
      return polygonItem->getName();
    case 1:
      return polygonItem->getColor();
    case 2:
      return polygonItem->isVisible();
    default:
      return QVariant();
    }
  }
  else {
    return QVariant();
  }
}
QVariant PolygonTableModel::headerData(int section, Qt::Orientation orientation,
  int role) const 
{
  if (role != Qt::DisplayRole) {
    return QVariant();
  }

  if (orientation == Qt::Horizontal) {
    return QString("Column %1").arg(section);
  }
  else {
    return QString("Row %1").arg(section);
  }
}

Qt::ItemFlags PolygonTableModel::flags(const QModelIndex &index) const
{
  if (!index.isValid())
    return Qt::ItemIsEnabled;

  return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
}
bool PolygonTableModel::setData(const QModelIndex &index, const QVariant &value,
  int role)
{
  if (index.isValid() && role == Qt::EditRole) {
    emit dataChanged(index, index);
    return true;
  }

  return false;
}

bool PolygonTableModel::insertRows(int position, int rows, const QModelIndex &index)
{
  beginInsertRows(QModelIndex(), position, position + rows - 1);

  for (int row = 0; row < rows; ++row) {
    m_polygonList.insert(position, new PolygonListItem("Empty Polygon"));
  }

  endInsertRows();
  return true;
}
bool PolygonTableModel::removeRows(int position, int rows, const QModelIndex &index)
{
  beginRemoveRows(QModelIndex(), position, position + rows - 1);

  for (int row = 0; row < rows; ++row) {
    m_polygonList.removeAt(position);
  }

  endRemoveRows();
  return true;
}
bool PolygonTableModel::insertItem(int position, PolygonListItem* item) {
  bool result = insertRows(position, 1);
  m_polygonList[position] = item;
  return result;
}
bool PolygonTableModel::appendItem(PolygonListItem* item) {
  int position = m_polygonList.size();
  bool result = insertRows(position, 1);
  m_polygonList[position] = item;
  return result;
}

PolygonListItem* PolygonTableModel::getPolygonItem(int index) {
  return m_polygonList[index];
}




