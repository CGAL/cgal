// Note: this structure is inspired from the QInputDialog that allows
// the user to easily create a dialog to get one value (integer,
// double or string). This structure allows the user to easily create
// a form inside a QDialog to get as manu values as needed.

#ifndef CGAL_QMULTIPLEINPUTDIALOG_H
#define CGAL_QMULTIPLEINPUTDIALOG_H

#include <QDialog>
#include <QFormLayout>
#include <QDialogButtonBox>

class QMultipleInputDialog
{
  QDialog* dialog;
  QFormLayout* form;
public:
  QMultipleInputDialog (const char* name, QWidget* parent)
  {
    dialog = new QDialog (parent);
    dialog->setWindowTitle (name);
    form = new QFormLayout(dialog);
  }
  ~QMultipleInputDialog ()
  {
    delete dialog;
  }

  template <typename QObjectType>
  QObjectType* add (const char* name)
  {
    QObjectType* out = new QObjectType (dialog);
    form->addRow (QString(name), out);
    return out;
  }

  int exec()
  {
    QDialogButtonBox* oknotok = new QDialogButtonBox
      (QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
       Qt::Horizontal, dialog);
      
    form->addRow (oknotok);
    QObject::connect (oknotok, SIGNAL(accepted()), dialog, SLOT(accept()));
    QObject::connect (oknotok, SIGNAL(rejected()), dialog, SLOT(reject()));

    return dialog->exec();
  }
};



#endif // CGAL_QMULTIPLEINPUTDIALOG_H
