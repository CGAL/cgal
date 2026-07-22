#include <QFileInfo>
#include <QLayout>
#include <QDialogButtonBox>
#include <QPushButton>

#include "ui_FileLoaderDialog.h"

class File_loader_dialog : public QDialog, private Ui::FileLoaderDialog
{
  Q_OBJECT
  public:
    File_loader_dialog(QWidget* parent = nullptr)
      : QDialog(parent)
    {
      setupUi(this);
    }
  static
    std::pair<QString,bool>
    getItem(QWidget* parent, QString filename, const QStringList& item_list, bool* ok)
    {
      File_loader_dialog dialog(parent);
      dialog.buttonBox->button(QDialogButtonBox::Ok)->setDefault(true);
      dialog.buttonBox->button(QDialogButtonBox::Ok)->setFocus();
      dialog.pluginBox->addItems(item_list);
      dialog.label->setText(tr("Available loaders for %1 :").arg(filename));
      QFileInfo fileinfo(filename);
      dialog.alwaysUse->setText(tr("Use for &all *.%1 files in this session").arg(fileinfo.completeSuffix()) );
      *ok = dialog.exec();
      return std::make_pair(
        dialog.pluginBox->currentText(),
        dialog.alwaysUse->isChecked()
      );
    }

};
