#ifndef _VALUES_LIST_H
#define _VALUES_LIST_H

#include <QWidget>
#include <QString>
#include <QList>
#include <QModelIndex>
#include <QItemDelegate>

class QTreeWidget;
class QTreeWidgetItem;

class Values_delegate : public QItemDelegate
{
  Q_OBJECT
public:
  Values_delegate(QWidget* parent);

signals:
  void new_color(const QModelIndex&) const;
  void new_value(const QModelIndex&) const;

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

class Values_list : public QWidget
{
  Q_OBJECT 
public:
  enum Field { Value = 0 , Color = 1, Name = 2};
  Values_list(QWidget* parent);
  
  // const accessors
  int numberOfValues() const;
  QColor color(const int i) const;
  QColor color(const QTreeWidgetItem* i) const;
  double value(const int i) const;
  QString name(const int i) const;
  bool enabled(const int i) const;
  bool enabled(const QTreeWidgetItem* i) const;
  const QTreeWidgetItem* item(const int i) const;
  const QTreeWidgetItem* search(const double value) const;

public slots:
  void save_values(QString) const;
  void load_values(QString);
  void on_plusButton_clicked();
  void on_minusButton_clicked();

  // setters
  void addValue(const double v = Values_list::default_value);

  void setHeaderTitle(QString);

private slots:
  void update_items_cache();

signals:
  void changed();
  void colors_changed();
  void values_changed();
private:
  QTreeWidget* treeWidget;
  typedef std::map<double, const QTreeWidgetItem*> Items_cache;
  Items_cache items_cache;

  static const double default_value;
};

#endif // _VALUES_LIST_H

