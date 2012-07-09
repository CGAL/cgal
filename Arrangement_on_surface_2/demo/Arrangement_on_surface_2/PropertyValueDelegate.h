#ifndef PROPERTY_VALUE_DELEGATE_H
#define PROPERTY_VALUE_DELEGATE_H
#include <QtGui>

class PropertyValueDelegate : public QItemDelegate
{
Q_OBJECT

public:
    PropertyValueDelegate( QObject* parent = 0 );

public:
    QWidget* createEditor( QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index ) const;
    void setModelData( QWidget* editor, QAbstractItemModel* model, const QModelIndex& index ) const;
    bool eventFilter( QObject *object, QEvent *event );

public slots:
    void commit( );

};
#endif // PROPERTY_VALUE_DELEGATE_H
