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
    ArrangementDemoPropertiesDialog( ArrangementDemoWindow* parent_ = 0, Qt::WindowFlags f = 0 );

protected:
    void setupUi( );
    void updateUi( );
    
    ArrangementDemoWindow* parent;
    Ui::ArrangementDemoPropertiesDialog* ui;
}; // class ArrangementDemoPropertiesDialog
#endif // ARRANGEMENT_DEMO_PROPERTIES_DIALOG_H
