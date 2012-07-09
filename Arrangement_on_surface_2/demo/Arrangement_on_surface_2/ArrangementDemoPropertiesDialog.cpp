#include "ArrangementDemoPropertiesDialog.h"
#include "ui_ArrangementDemoPropertiesDialog.h"
#include "ArrangementDemoWindow.h"
#include "PropertyValueDelegate.h"
#include "DeleteCurveMode.h"

ArrangementDemoPropertiesDialog::
ArrangementDemoPropertiesDialog( ArrangementDemoWindow* parent_, Qt::WindowFlags f ):
    QDialog( parent_, f ),
    parent( parent_ ),
    ui( new Ui::ArrangementDemoPropertiesDialog )
{
    this->setupUi( );
}

void
ArrangementDemoPropertiesDialog::
setupUi( )
{
    this->ui->setupUi( this );
    PropertyValueDelegate* myDelegate = new PropertyValueDelegate;
    this->ui->tableWidget->setItemDelegate( myDelegate );

    // populate the table widget with items
    QTableWidgetItem* edgeColorItem = new QTableWidgetItem;
    QTableWidgetItem* vertexColorItem = new QTableWidgetItem;
    QTableWidgetItem* edgeWidthItem = new QTableWidgetItem;
    QTableWidgetItem* vertexRadiusItem = new QTableWidgetItem;
    QTableWidgetItem* deleteCurveModeItem = new QTableWidgetItem;

    this->ui->tableWidget->setItem( 0, 0, edgeColorItem );
    this->ui->tableWidget->setItem( 1, 0, vertexColorItem );
    this->ui->tableWidget->setItem( 2, 0, edgeWidthItem );
    this->ui->tableWidget->setItem( 3, 0, vertexRadiusItem );
    this->ui->tableWidget->setItem( 4, 0, deleteCurveModeItem );

    // fill in the items with data
    this->updateUi( );
}

void
ArrangementDemoPropertiesDialog::
updateUi( )
{
    if ( this->parent == NULL )
    {
        return;
    }
    ArrangementDemoTabBase* currentTab = this->parent->getCurrentTab( );
    if ( currentTab == NULL )
    {
        return;
    }
    CGAL::Qt::ArrangementGraphicsItemBase* agi = currentTab->getArrangementGraphicsItem( );
    if ( agi == NULL )
    {
        return;
    }

    QPen vertexPen = agi->getVerticesPen( );
    QPen edgePen = agi->getEdgesPen( );
    QBrush vertexPenBrush = vertexPen.brush( );
    QBrush edgePenBrush = edgePen.brush( );
    QColor vertexColor = vertexPenBrush.color( );
    QColor edgeColor = edgePenBrush.color( );
    int edgeWidth = edgePen.width( );
    int vertexRadius = vertexPen.width( );

    QTableWidgetItem* edgeColorItem = this->ui->tableWidget->item( 0, 0 );
    QTableWidgetItem* edgeWidthItem = this->ui->tableWidget->item( 1, 0 );
    QTableWidgetItem* vertexColorItem = this->ui->tableWidget->item( 2, 0 );
    QTableWidgetItem* vertexRadiusItem = this->ui->tableWidget->item( 3, 0 );
    QTableWidgetItem* deleteCurveModeItem = this->ui->tableWidget->item( 4, 0 );

    edgeColorItem->setData( Qt::DisplayRole, edgeColor );
    edgeColorItem->setData( Qt::DecorationRole, edgeColor );
    edgeColorItem->setData( Qt::UserRole, QVariant::fromValue( edgeColor ) );

    edgeWidthItem->setData( Qt::DisplayRole, edgeWidth );

    vertexColorItem->setData( Qt::DisplayRole, vertexColor );
    vertexColorItem->setData( Qt::DecorationRole, vertexColor );
    vertexColorItem->setData( Qt::UserRole, QVariant::fromValue( vertexColor ) );

    vertexRadiusItem->setData( Qt::DisplayRole, vertexRadius );

    DeleteCurveMode deleteCurveMode;
    deleteCurveModeItem->setData( Qt::DisplayRole, DeleteCurveMode::ToString( deleteCurveMode ) );
    deleteCurveModeItem->setData( Qt::UserRole, QVariant::fromValue( deleteCurveMode ) );

#if 0
        ArrangementDemoPropertiesDialog->setWindowTitle(QApplication::translate("ArrangementDemoPropertiesDialog", "Dialog", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem = tableWidget->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QApplication::translate("ArrangementDemoPropertiesDialog", "Value", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem1 = tableWidget->verticalHeaderItem(0);
        ___qtablewidgetitem1->setText(QApplication::translate("ArrangementDemoPropertiesDialog", "Edge Color", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem2 = tableWidget->verticalHeaderItem(1);
        ___qtablewidgetitem2->setText(QApplication::translate("ArrangementDemoPropertiesDialog", "Edge Width", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem3 = tableWidget->verticalHeaderItem(2);
        ___qtablewidgetitem3->setText(QApplication::translate("ArrangementDemoPropertiesDialog", "Vertex Color", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem4 = tableWidget->verticalHeaderItem(3);
        ___qtablewidgetitem4->setText(QApplication::translate("ArrangementDemoPropertiesDialog", "Vertex Radius", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem5 = tableWidget->verticalHeaderItem(4);
        ___qtablewidgetitem5->setText(QApplication::translate("ArrangementDemoPropertiesDialog", "Delete Curve Mode", 0, QApplication::UnicodeUTF8));
#endif
}
