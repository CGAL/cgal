#ifndef ARRANGEMENT_DEMO_PROPERTIES_DIALOG_H
#define ARRANGEMENT_DEMO_PROPERTIES_DIALOG_H
#include <QDialog>

class ArrangementDemoWindow;

namespace Ui
{
    class ArrangementDemoPropertiesDialog;
}

class ArrangementDemoPropertiesDialog : public QDialog
{
Q_OBJECT
public:
    enum PropertyKey {
        EDGE_COLOR_KEY,
        EDGE_WIDTH_KEY,
        VERTEX_COLOR_KEY,
        VERTEX_RADIUS_KEY,
        DELETE_CURVE_MODE_KEY
    };

    ArrangementDemoPropertiesDialog( ArrangementDemoWindow* parent_ = 0, Qt::WindowFlags f = 0 );
    QVariant property( int index );

protected:
    void setupUi( );
    void updateUi( );
    
    ArrangementDemoWindow* parent;
    Ui::ArrangementDemoPropertiesDialog* ui;
}; // class ArrangementDemoPropertiesDialog
#endif // ARRANGEMENT_DEMO_PROPERTIES_DIALOG_H
