#include "Show_point_dialog.h"
#include "ui_Show_point_dialog.h"

#include <QClipboard>

Show_point_dialog::Show_point_dialog(QWidget* parent)
  : QDialog(parent)
  , ui(new Ui::Show_point_dialog)
  , m_has_correct_coordinates(false)
{
  ui->setupUi(this);

  QClipboard::Mode mode = QClipboard::Selection; 
  while(true) {
    QString clipboard_text = QApplication::clipboard()->text(mode);
    
    interprete_string(clipboard_text);
    if(m_has_correct_coordinates) {
      ui->lineEdit->setText(clipboard_text);
      ui->lineEdit->selectAll();
      break;
    } else {
      interprete_string(ui->lineEdit->text());
    }
    if(mode == QClipboard::Selection) mode = QClipboard::Clipboard;
    else break;
  }

  connect(ui->lineEdit, SIGNAL(textChanged(const QString&)),
          this, SLOT(interprete_string(const QString&)));
}

Show_point_dialog::~Show_point_dialog()
{
  delete ui;
}

void Show_point_dialog::interprete_string(const QString& string)
{
  QString double_re("([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)");
  QString not_double_char_re("[^0-9-+.eE]");
  QString full_re = QString("^") 
    + not_double_char_re + "*"
    + double_re
    + not_double_char_re + "+"
    + double_re
    + not_double_char_re + "+"
    + double_re 
    + "(" + not_double_char_re + "*"
    + "|" + not_double_char_re + "+" + double_re
    + ")"
    + "$";
  QRegExp re(full_re);
  if(re.exactMatch(string)) {
    // const double x = re.cap(1).toDouble();
    // const double y = re.cap(2).toDouble();
    // const double z = re.cap(3).toDouble();
    ui->coord_x->setText(QString(re.cap(1)));
    ui->coord_y->setText(QString(re.cap(2)));
    ui->coord_z->setText(QString(re.cap(3)));
    m_has_correct_coordinates = true;
  }
  else {
    ui->coord_x->setText(QString());
    ui->coord_y->setText(QString());
    ui->coord_z->setText(QString());
    m_has_correct_coordinates = false;
  }
}

double Show_point_dialog::get_x() const
{
  return ui->coord_x->text().toDouble();
}

double Show_point_dialog::get_y() const
{
  return ui->coord_y->text().toDouble();
}

double Show_point_dialog::get_z() const
{
  return ui->coord_z->text().toDouble();
}

bool Show_point_dialog::has_correct_coordinates() const
{
  return m_has_correct_coordinates;
}

#include "Show_point_dialog.moc"
