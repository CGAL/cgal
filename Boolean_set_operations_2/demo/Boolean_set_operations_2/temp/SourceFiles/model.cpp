

/*
model.cpp

A simple model that uses a QStringList as its data source.
*/

#include "model.h"

/*!
Returns the number of items in the string list as the number of rows
in the model.
*/


int StringListModel::rowCount(const QModelIndex &parent) const
{
	return stringList.count();
}


/*!
Returns an appropriate value for the requested data.
If the view requests an invalid index, an invalid variant is returned.
Any valid index that corresponds to a string in the list causes that
string to be returned.
*/


QVariant StringListModel::data(const QModelIndex &index, int role) const
{
	if (!index.isValid())
		return QVariant();

	if (index.row() >= stringList.size())
		return QVariant();

	if (role == Qt::DisplayRole)
		return stringList.at(index.row());
	else
		return QVariant();
}


/*!
Returns the appropriate header string depending on the orientation of
the header and the section. If anything other than the display role is
requested, we return an invalid variant.
*/


QVariant StringListModel::headerData(int section, Qt::Orientation orientation,
	int role) const
{
	if (role != Qt::DisplayRole)
		return QVariant();

	if (orientation == Qt::Horizontal)
		return QString("Column %1").arg(section);
	else
		return QString("Row %1").arg(section);
}


/*!
Returns an appropriate value for the item's flags. Valid items are
enabled, selectable, and editable.
*/


Qt::ItemFlags StringListModel::flags(const QModelIndex &index) const
{
	if (!index.isValid())
		return Qt::ItemIsEnabled;

	return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
}


/*!
Changes an item in the string list, but only if the following conditions
are met:

* The index supplied is valid.
* The index corresponds to an item to be shown in a view.
* The role associated with editing text is specified.

The dataChanged() signal is emitted if the item is changed.
*/


bool StringListModel::setData(const QModelIndex &index,
	const QVariant &value, int role)
{
	if (index.isValid() && role == Qt::EditRole) {

		stringList.replace(index.row(), value.toString());
		emit dataChanged(index, index);
		return true;
	}

	return false;
}


/*!
Inserts a number of rows into the model at the specified position.
*/


bool StringListModel::insertRows(int position, int rows, const QModelIndex &parent)
{
	beginInsertRows(QModelIndex(), position, position + rows - 1);

	for (int row = 0; row < rows; ++row) {
		stringList.insert(position, "");
	}

	endInsertRows();
	return true;

}


/*!
Removes a number of rows from the model at the specified position.
*/


bool StringListModel::removeRows(int position, int rows, const QModelIndex &parent)
{
	beginRemoveRows(QModelIndex(), position, position + rows - 1);

	for (int row = 0; row < rows; ++row) {
		stringList.removeAt(position);
	}

	endRemoveRows();
	return true;

}






