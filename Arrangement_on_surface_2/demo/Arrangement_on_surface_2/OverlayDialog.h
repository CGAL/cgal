#ifndef OVERLAY_DIALOG_H
#define OVERLAY_DIALOG_H
#include <QDialog>
#include "ui_OverlayDialog.h"
#include <vector>
#include <CGAL/Object.h>

class ArrangementDemoWindow;

class OverlayDialog : public QDialog
{
Q_OBJECT

public:
    typedef enum OverlayDialogRole {
        ARRANGEMENT = 32
    } OverlayDialogRole;

    OverlayDialog( ArrangementDemoWindow* parent, Qt::WindowFlags f = 0 );

    std::vector< CGAL::Object > selectedArrangements( ) const;

public slots:
    void on_pickPushButton_pressed( );
    void on_unpickPushButton_pressed( );

protected:
    void restrictSelection( QListWidgetItem* item );
    void unrestrictSelection( );

    Ui::OverlayDialog* ui;
};
#endif // OVERLAY_DIALOG_H
