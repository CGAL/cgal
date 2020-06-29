#include "values_list.h"
#include "ui_values_list.h"
#include "colorlisteditor.h"
#include <iostream>

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
#if QT_VERSION >= QT_VERSION_CHECK(5, 10, 0)
#include <QRandomGenerator>
#endif

Values_delegate::Values_delegate(QWidget* parent) : QItemDelegate(parent) {}
void Values_delegate::paint(QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index) const
{
  switch(index.column())
  {
  case Values_list::Color: {
    painter->fillRect(option.rect, index.data().value<QColor>());
    drawFocus(painter, option, option.rect);
    break;
  }
  default:
    QItemDelegate::paint(painter, option, index);
  }
}

QWidget *Values_delegate::createEditor(QWidget *parent,
                                          const QStyleOptionViewItem & option,
                                          const QModelIndex & index) const
{
  if(index.column() == Values_list::Color)
  {
    return new ColorListEditor(parent);
  }
  else if(index.column() == Values_list::Value)
  {
    QLineEdit* lineedit = new QLineEdit(parent);
    lineedit->setAutoFillBackground(true);
    lineedit->setValidator(new QDoubleValidator(lineedit));
    return lineedit;
  }
  else return QItemDelegate::createEditor(parent, option, index);
}

void Values_delegate::setEditorData(QWidget *editor,
                                       const QModelIndex &index) const
{
  if(index.column() == Values_list::Color)
  {
    ColorListEditor* coloreditor = qobject_cast<ColorListEditor*>(editor);
    if(coloreditor)
      coloreditor->setColor(index.data().value<QColor>());
  }
  else if(index.column() == Values_list::Value)
  {
    QLineEdit* lineedit = qobject_cast<QLineEdit*>(editor);
    if(lineedit)
      lineedit->setText(index.data().toString());
  }
  else QItemDelegate::setEditorData(editor, index);
}
void Values_delegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                      const QModelIndex &index) const
{
  if(index.column() == Values_list::Color)
  {
    ColorListEditor* coloreditor = qobject_cast<ColorListEditor*>(editor);
    if(coloreditor)
    {
      model->setData(index, coloreditor->color());
      Q_EMIT new_color(index);
    }
  }
  else if(index.column() == Values_list::Value)
  {
    QLineEdit* lineedit = qobject_cast<QLineEdit*>(editor);
    if(lineedit)
    {
      model->setData(index, lineedit->text().toDouble());
      Q_EMIT new_value(index);
    }
  }
  else QItemDelegate::setModelData(editor, model, index);
}

const double Values_list::default_value = 0.0;

Values_list::Values_list(QWidget* parent):
  QWidget(parent)
{
  Ui::Values_list().setupUi(this);

  treeWidget = parent->findChild<QTreeWidget*>("treeWidget");
  Q_ASSERT_X(treeWidget, "Values_list constructor", "cannot find widget \"treeWidget\"");

  treeWidget->sortByColumn(Value, Qt::AscendingOrder);

  treeWidget->header()->setSectionsClickable(false);


  Values_delegate* values_delegate = new Values_delegate(parent);

  treeWidget->setItemDelegate(values_delegate);
  connect(values_delegate, SIGNAL(new_value(const QModelIndex&)),
          this, SIGNAL(values_changed()));
  connect(values_delegate, SIGNAL(new_color(const QModelIndex&)),
          this, SIGNAL(colors_changed()));
  connect(this->treeWidget->model(),
          SIGNAL(dataChanged (const QModelIndex &, const QModelIndex &)),
          this, SIGNAL(changed()));

  connect(this, SIGNAL(changed()),
          this, SLOT(update_items_cache()));
}

QColor Values_list::color(const int i) const
{
  if(i < 0 || i > treeWidget->topLevelItemCount())
    return QColor();
  else
    return treeWidget->topLevelItem(i)->data(Color, Qt::DisplayRole).value<QColor>();
}

QColor Values_list::color(const QTreeWidgetItem* item) const
{
    return item->data(Color, Qt::DisplayRole).value<QColor>();
}

int Values_list::numberOfValues() const
{
  return treeWidget->topLevelItemCount();
}

double Values_list::value(const int i) const
{
  if(i < 0 || i > numberOfValues())
    return 0.;
  else
    return treeWidget->topLevelItem(i)->data(Value, Qt::DisplayRole).toDouble();
}

QString Values_list::name(const int i) const
{
  if(i < 0 || i > treeWidget->topLevelItemCount())
    return QString();
  else
    return treeWidget->topLevelItem(i)->data(Name, Qt::DisplayRole).toString();
}

bool Values_list::enabled(const int i) const
{
  if(i < 0 || i > treeWidget->topLevelItemCount())
    return false;
  else
    return treeWidget->topLevelItem(i)->data(Value, Qt::CheckStateRole).toInt() > 0;
}

bool Values_list::enabled(const QTreeWidgetItem* item) const
{
    return item->data(Value, Qt::CheckStateRole).toInt() > 0;
}

const QTreeWidgetItem* Values_list::item(const int i) const
{
  if(i < 0 || i > treeWidget->topLevelItemCount())
    return 0;
  else
    return treeWidget->topLevelItem(i);
}

void Values_list::save_values(QString filename) const
{
  QSettings settings;

  settings.beginGroup(QUrl::toPercentEncoding(filename));
  settings.beginWriteArray("values");
  for (int i = 0; i < numberOfValues(); ++i) {
    settings.setArrayIndex(i);
    settings.setValue("value", value(i));
    settings.setValue("color", color(i));
    settings.setValue("name", name(i));
    settings.setValue("enabled", enabled(i));
  }
 settings.endArray();
 settings.endGroup();
}

void Values_list::load_values(QString filename)
{
  QSettings settings;

  treeWidget->clear();
  settings.beginGroup(QUrl::toPercentEncoding(filename));
  int nb = settings.beginReadArray("values");
  for (int i = 0; i < nb; ++i) {
    settings.setArrayIndex(i);
    QTreeWidgetItem *newItem = new QTreeWidgetItem(treeWidget);
    newItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
    newItem->setData(Value, Qt::CheckStateRole, settings.value("enabled").toBool() ? Qt::Checked : Qt::Unchecked);
    newItem->setData(Value, Qt::DisplayRole, settings.value("value").toDouble());
    newItem->setData(Color, Qt::DisplayRole, settings.value("color").value<QColor>());
    newItem->setData(Name, Qt::DisplayRole, settings.value("name").toString());
  }
  settings.endArray();
  settings.endGroup();
  update_items_cache();
}

void Values_list::on_minusButton_clicked()
{
  Q_FOREACH(QTreeWidgetItem* item, treeWidget->selectedItems())
  {
    //   treeWidget->invisibleRootItem()->removeChild(item);
    delete item;
  }
  Q_EMIT values_changed();
}

void Values_list::on_plusButton_clicked()
{
  addValue();
}

void Values_list::addValue(const double i)
{
  QTreeWidgetItem *newItem = new QTreeWidgetItem(treeWidget);
  newItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
  newItem->setData(Value, Qt::CheckStateRole, Qt::Checked);
  newItem->setData(Value, Qt::DisplayRole, i);
  QStringList colors = QColor::colorNames();
#if QT_VERSION < QT_VERSION_CHECK(5, 10, 0)
  const int color_index = qrand() % colors.size();
#else
const int color_index = QRandomGenerator::global()->generate() % colors.size();
#endif

  QColor color = QColor(colors[color_index]);
  newItem->setData(Color, Qt::DisplayRole, color);
  newItem->setData(Name, Qt::DisplayRole, "");
  Q_EMIT values_changed();
}

void Values_list::update_items_cache() {
  items_cache.clear();
  for(int i = 0, nb = numberOfValues(); i < nb; ++i) {
    items_cache.insert(std::make_pair(value(i), item(i)));
  }
}

const QTreeWidgetItem* Values_list::search(const double value) const
{
  Items_cache::const_iterator it = items_cache.find(value);
  if(it != items_cache.end()) {
    return it->second;
  }
  else {
    return 0;
  }
}

void Values_list::setHeaderTitle(QString title)
{
  treeWidget->headerItem()->setText(0, title);
}

