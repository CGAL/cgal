#include "EventFilterManagerGroup.h"

EventFilterManagerGroup::EventFilterManagerGroup() :
  m_objectsMap(),
  m_filtersMap()
{
  QString m_currentFilter("");
  QString m_currentWatchedObject("");
}

EventFilterManagerGroup::~EventFilterManagerGroup()
{
}

void EventFilterManagerGroup::addObjectToWatch(QString name, QObject* object)
{
  EventFilterManager* eventFilterManager = new EventFilterManager();
  object->installEventFilter(eventFilterManager);
  m_objectsMap[name] = eventFilterManager;
}

void EventFilterManagerGroup::removeObjectFromWatch(QString name)
{
  EventFilterManager* eventFilterManager = m_objectsMap[name];
  eventFilterManager->clear();
  eventFilterManager->~EventFilterManager();
  m_objectsMap.remove(name);
  if (m_currentWatchedObject == name) {
    m_currentWatchedObject = "";
    m_currentFilter = "";
  }
  QList<QString>* filtersToRemove = new QList<QString>();
  QMap<QString, QString>::iterator imap;
  for (imap = m_filtersMap.begin(); imap != m_filtersMap.end(); ++imap) {
    if (imap.value() == name)
      filtersToRemove->append(imap.key());
  }
  QList<QString>::const_iterator ilist;
  for (ilist = filtersToRemove->cbegin(); ilist != filtersToRemove->cend(); ++ilist)
  {
    m_filtersMap.remove(ilist.i->t());
  }
}

void EventFilterManagerGroup::addFilterWidget(const QString targetName, const QString name, QObject* filter)
{
  m_filtersMap[name] = targetName;
  m_objectsMap[targetName]->addFilterWidget(name, filter);
}

QObject* EventFilterManagerGroup::removeFilterWidget(const QString name)
{
  QObject* filter = m_objectsMap[m_filtersMap[name]]->removeFilterWidget(name);
  m_filtersMap.remove(name);
  return filter;
}

QObject* EventFilterManagerGroup::getFilterWidget(const QString name)
{
  return m_objectsMap[m_filtersMap[name]]->getFilterWidget(name);
}

QString const EventFilterManagerGroup::getCurrentFilterName() const
{
  return m_currentFilter;
}

void EventFilterManagerGroup::setCurrentFilterName(const QString name)
{
  if (m_filtersMap.contains(name)) {
    QString newObject = m_filtersMap[name];
    if (m_currentWatchedObject != QString::null && newObject != m_currentWatchedObject)
      m_objectsMap[m_currentWatchedObject]->setCurrentFilterName("");
    m_objectsMap[newObject]->setCurrentFilterName(name);
    m_currentWatchedObject = newObject;
    m_currentFilter = name;
  }
  else if (name == "") {
    if (m_currentWatchedObject != "")
      m_objectsMap[m_currentWatchedObject]->setCurrentFilterName("");
    m_currentWatchedObject = "";
    m_currentFilter = "";
  }
  else {
    qDebug() << Q_FUNC_INFO << "The name \"" + name + "\" is not a name of any filter and is not an empty name";
  }
}

QList<QString> EventFilterManagerGroup::getFilterNames()
{
  auto x = m_filtersMap.keys();
  return std::move(x);
}

void EventFilterManagerGroup::clear()
{
  QList<QString>::const_iterator i;
  QList<QString> keys = m_objectsMap.keys();
  for (i = keys.cbegin(); i != keys.cend(); ++i) {
    m_objectsMap[*i]->clear();
    m_objectsMap[*i]->~EventFilterManager();
    m_objectsMap.remove(*i);
  }
  m_objectsMap.clear();
  m_filtersMap.clear();
}
