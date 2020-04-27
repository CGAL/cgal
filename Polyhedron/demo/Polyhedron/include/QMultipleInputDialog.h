// Note: this structure is inspired from the QInputDialog that allows
// the user to easily create a dialog to get one value (integer,
// double or string). This structure allows the user to easily create
// a form inside a QDialog to get as manu values as needed.

#ifndef CGAL_QMULTIPLEINPUTDIALOG_H
#define CGAL_QMULTIPLEINPUTDIALOG_H

#include <QDialog>
#include <QFormLayout>
#include <QDialogButtonBox>
#include <QRadioButton>

class QMultipleInputDialog
{
  QDialog* dialog;
  QFormLayout* form;
  std::map<std::string, QWidget*> map_widgets;


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
  QObjectType* add (const char* name, const char* key = NULL)
  {
    QObjectType* out = NULL;

    if (boost::is_same<QObjectType, QRadioButton>::value)
    {
      out = dynamic_cast<QObjectType*>(new QRadioButton (QString(name), dialog));
      form->addRow (out);
    }
    else
    {
      out = new QObjectType (dialog);
      form->addRow (QString(name), out);
    }

    if (key != NULL)
      map_widgets.insert (std::make_pair (key, out));

    return out;
  }

  template <typename QObjectType>
  const QObjectType* get (const char* key) const
  {
    typename std::map<std::string, QWidget*>::const_iterator
      found = map_widgets.find (key);
    if (found == map_widgets.end())
      return NULL;

    QWidget* widget_out = found->second;
    return qobject_cast<QObjectType*>(widget_out);
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

  void exec_no_cancel()
  {
    QDialogButtonBox* ok = new QDialogButtonBox
      (QDialogButtonBox::Ok, Qt::Horizontal, dialog);

    form->addRow (ok);
    QObject::connect (ok, SIGNAL(accepted()), dialog, SLOT(accept()));
    dialog->exec();
  }
};



#endif // CGAL_QMULTIPLEINPUTDIALOG_H
