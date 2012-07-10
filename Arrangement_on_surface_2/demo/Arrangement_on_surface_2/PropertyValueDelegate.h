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

class PositiveSpinBox : public QSpinBox
{
Q_OBJECT
Q_PROPERTY( unsigned int value READ value WRITE setValue USER true )

public:
    PositiveSpinBox( QWidget* parent );
    void setValue( unsigned int );
    unsigned int value( ) const;
};

#endif // PROPERTY_VALUE_DELEGATE_H
