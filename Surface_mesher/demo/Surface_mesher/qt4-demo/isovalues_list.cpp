#include "isovalues_list.h"
#include "ui_isovalues_list.h"
#include "colorlisteditor.h"

#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QTreeWidgetItemIterator>
#include <QHeaderView>
#include <QMetaProperty>
#include <QItemDelegate>
#include <QItemEditorFactory>
#include <QItemEditorCreatorBase>
#include <QPainter>
#include <QStringList>
#include <QString>
#include <QList>
#include <QVariant>
#include <QSettings>
#include <QUrl>
#include <QLineEdit>
#include <QDoubleValidator>

Isovalues_delegate::Isovalues_delegate(QWidget* parent) : QItemDelegate(parent) {}
void Isovalues_delegate::paint(QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index) const
{
  switch(index.column())
  {
  case Isovalues_list::Color: {
    painter->fillRect(option.rect, index.data().value<QColor>());
    drawFocus(painter, option, option.rect);
    break;
  }
  default:
    QItemDelegate::paint(painter, option, index);
  }
}

QWidget *Isovalues_delegate::createEditor(QWidget *parent,
                                          const QStyleOptionViewItem & option,
                                          const QModelIndex & index) const
{
  if(index.column() == Isovalues_list::Color)
  {
    return new ColorListEditor(parent);
  }
  else if(index.column() == Isovalues_list::Isovalue)
  {
    QLineEdit* lineedit = new QLineEdit(parent);
    lineedit->setAutoFillBackground(true);
    lineedit->setValidator(new QDoubleValidator(lineedit));
    return lineedit;
  }
  else return QItemDelegate::createEditor(parent, option, index);
}
  
void Isovalues_delegate::setEditorData(QWidget *editor,
                                       const QModelIndex &index) const
{
  if(index.column() == Isovalues_list::Color)
  {
    ColorListEditor* coloreditor = qobject_cast<ColorListEditor*>(editor);
    if(coloreditor)
      coloreditor->setColor(index.data().value<QColor>());
  }
  else if(index.column() == Isovalues_list::Isovalue)
  {
    QLineEdit* lineedit = qobject_cast<QLineEdit*>(editor);
    if(lineedit)
      lineedit->setText(index.data().toString());
  }
  else QItemDelegate::setEditorData(editor, index);
}
void Isovalues_delegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                      const QModelIndex &index) const
{
  if(index.column() == Isovalues_list::Color)
  {
    ColorListEditor* coloreditor = qobject_cast<ColorListEditor*>(editor);
    if(coloreditor)
    {
      model->setData(index, coloreditor->color());
      emit new_color(index);
    }
  }
  else if(index.column() == Isovalues_list::Isovalue)
  {
    QLineEdit* lineedit = qobject_cast<QLineEdit*>(editor);
    if(lineedit)
    {
      model->setData(index, lineedit->text().toDouble());
      emit new_isovalue(index);
    }
  }
  else QItemDelegate::setModelData(editor, model, index);
}

const double Isovalues_list::default_isovalue = 0.0;

Isovalues_list::Isovalues_list(QWidget* parent):
  QWidget(parent)
{
  Ui::Isovalues_list().setupUi(this);

  treeWidget = parent->findChild<QTreeWidget*>("treeWidget");
  Q_ASSERT_X(treeWidget, "Isovalues_list constructor", "cannot find widget \"treeWidget\"");

  treeWidget->sortByColumn(Isovalue, Qt::AscendingOrder);
  treeWidget->header()->setClickable(false);

  Isovalues_delegate* isovalues_delegate = new Isovalues_delegate(parent);

  treeWidget->setItemDelegate(isovalues_delegate);
  connect(isovalues_delegate, SIGNAL(new_isovalue(const QModelIndex&)),
          this, SIGNAL(isovalues_changed()));
  connect(isovalues_delegate, SIGNAL(new_color(const QModelIndex&)),
          this, SIGNAL(colors_changed()));
}

QColor Isovalues_list::color(const int i) const
{
  if(i < 0 || i > treeWidget->topLevelItemCount())
    return QColor();
  else
    return treeWidget->topLevelItem(i)->data(Color, Qt::DisplayRole).value<QColor>();
}

QColor Isovalues_list::color(const QTreeWidgetItem* item) const
{
    return item->data(Color, Qt::DisplayRole).value<QColor>();
}

int Isovalues_list::numberOfIsoValues() const
{
  return treeWidget->topLevelItemCount();
}

double Isovalues_list::isovalue(const int i) const
{
  if(i < 0 || i > numberOfIsoValues())
    return 0.;
  else
    return treeWidget->topLevelItem(i)->data(Isovalue, Qt::DisplayRole).toDouble();
}

QString Isovalues_list::name(const int i) const
{
  if(i < 0 || i > treeWidget->topLevelItemCount())
    return QString();
  else
    return treeWidget->topLevelItem(i)->data(Name, Qt::DisplayRole).toString();
}

bool Isovalues_list::enabled(const int i) const
{
  if(i < 0 || i > treeWidget->topLevelItemCount())
    return 0.;
  else
    return treeWidget->topLevelItem(i)->data(Isovalue, Qt::CheckStateRole).toDouble();
}

bool Isovalues_list::enabled(const QTreeWidgetItem* item) const
{
    return item->data(Isovalue, Qt::CheckStateRole).toDouble();
}

const QTreeWidgetItem* Isovalues_list::item(const int i) const
{
  if(i < 0 || i > treeWidget->topLevelItemCount())
    return 0;
  else
    return treeWidget->topLevelItem(i);
}

void Isovalues_list::save_values(QString filename) const
{
  QSettings settings;

  settings.beginGroup(QUrl::toPercentEncoding(filename));
  settings.beginWriteArray("isovalues");
  for (int i = 0; i < numberOfIsoValues(); ++i) {
    settings.setArrayIndex(i);
    settings.setValue("isovalue", isovalue(i));
    settings.setValue("color", color(i));
    settings.setValue("name", name(i));
    settings.setValue("enabled", enabled(i));
  }
 settings.endArray();
 settings.endGroup();
}

void Isovalues_list::load_values(QString filename) 
{
  QSettings settings;

  settings.beginGroup(QUrl::toPercentEncoding(filename));
  int nb = settings.beginReadArray("isovalues");
  for (int i = 0; i < nb; ++i) {
    settings.setArrayIndex(i);
    QTreeWidgetItem *newItem = new QTreeWidgetItem(treeWidget);
    newItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
    newItem->setData(Isovalue, Qt::CheckStateRole, settings.value("enabled").toBool() ? Qt::Checked : Qt::Unchecked);
    newItem->setData(Isovalue, Qt::DisplayRole, settings.value("isovalue").toDouble());
    newItem->setData(Color, Qt::DisplayRole, settings.value("color").value<QColor>());
    newItem->setData(Name, Qt::DisplayRole, settings.value("name").toString());
  }
  settings.endArray();
  settings.endGroup();
}

void Isovalues_list::on_minusButton_clicked()
{
  Q_FOREACH(QTreeWidgetItem* item, treeWidget->selectedItems())
  {
    treeWidget->invisibleRootItem()->removeChild(item);
    delete item;
  }
}

void Isovalues_list::on_plusButton_clicked()
{
  QTreeWidgetItem *newItem = new QTreeWidgetItem(treeWidget);
  newItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
  newItem->setData(Isovalue, Qt::CheckStateRole, Qt::Checked);
  newItem->setData(Isovalue, Qt::DisplayRole, default_isovalue);
  QStringList colors = QColor::colorNames();
  const int color_index = qrand() % colors.size();
  QColor color = QColor(colors[color_index]);
  newItem->setData(Color, Qt::DisplayRole, color);
  newItem->setData(Name, Qt::DisplayRole, "");
}

#include "isovalues_list.moc"
