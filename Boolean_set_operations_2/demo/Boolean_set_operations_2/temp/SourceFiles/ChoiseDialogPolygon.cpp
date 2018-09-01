#include "ChoiseDialogPolygon.h"
#include "ui_ChoiseDialogPolygon.h"

ChoiseDialogPolygon::ChoiseDialogPolygon(PolygonTableModel* options,
                                         QWidget* parent) :
  QDialog(parent),
  ui(new Ui::ChoiseDialogPolygon)
{
  ui->setupUi(this);

  m_options = options;
  ui->listView->setModel(m_options);
}

ChoiseDialogPolygon::~ChoiseDialogPolygon()
{
    delete ui;
}

PolygonWithHoles* ChoiseDialogPolygon::getChoise()
{
    return m_choise;
}

void ChoiseDialogPolygon::on_buttonBox_accepted()
{
  m_choise = m_options->getPolygonItem(ui->listView->currentIndex().row())->
                        getPolygonItem()->polygon();
}
