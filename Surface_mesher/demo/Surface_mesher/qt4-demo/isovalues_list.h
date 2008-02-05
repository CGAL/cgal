#ifndef _ISOVALUES_LIST_H
#define _ISOVALUES_LIST_H

#include <QWidget>
#include <QString>
#include <QList>

class QTreeWidget;
class QTreeWidgetItem;

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
  void color_changed(int);
  void isovalues_changed();
private:
  QTreeWidget* treeWidget;

  static const double default_isovalue;
};

#endif // _ISOVALUES_LIST_H

