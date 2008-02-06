#ifndef _ISOVALUES_LIST_H
#define _ISOVALUES_LIST_H

#include <QWidget>
#include <QString>
#include <QList>
#include <QModelIndex>
#include <QItemDelegate>

class QTreeWidget;
class QTreeWidgetItem;

class Isovalues_delegate : public QItemDelegate
{
  Q_OBJECT
public:
  Isovalues_delegate(QWidget* parent);

signals:
  void new_color(const QModelIndex&) const;
  void new_isovalue(const QModelIndex&) const;

protected:
  void paint(QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index) const;
  QWidget *createEditor(QWidget *parent,
                        const QStyleOptionViewItem & option,
                        const QModelIndex & index) const;
  void setEditorData(QWidget *editor,
                     const QModelIndex &index) const;
  void setModelData(QWidget *editor, QAbstractItemModel *model,
                    const QModelIndex &index) const;
};

class Isovalues_list : public QWidget
{
  Q_OBJECT 
public:
  enum Field { Isovalue = 0 , Color = 1, Name = 2};
  Isovalues_list(QWidget* parent);

  int numberOfIsoValues() const;
  QColor color(const int i) const;
  QColor color(const QTreeWidgetItem* i) const;
  double isovalue(const int i) const;
  QString name(const int i) const;
  bool enabled(const int i) const;
  bool enabled(const QTreeWidgetItem* i) const;
  const QTreeWidgetItem* item(const int i) const;

public slots:
  void save_values(QString) const;
  void load_values(QString);
  void on_plusButton_clicked();
  void on_minusButton_clicked();

signals:
  void colors_changed();
  void isovalues_changed();
private:
  QTreeWidget* treeWidget;

  static const double default_isovalue;
};

#endif // _ISOVALUES_LIST_H

