#ifndef NEW_TAB_DIALOG_H
#define NEW_TAB_DIALOG_H
#include <QDialog>
#include "ui_NewTabDialog.h"

class NewTabDialog : public QDialog
{
public:
    NewTabDialog( QWidget* parent = 0, Qt::WindowFlags f = 0 );
    int checkedId( ) const;

protected:
    Ui::NewTabDialog* ui;
    QButtonGroup* buttonGroup;
};
#endif // NEW_TAB_DIALOG_H
