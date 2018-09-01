#include "EventFilterManager.h"

EventFilterManager::EventFilterManager(QObject* parent) :
  QObject(parent),
	m_filterMap()
{
	QString m_currentFilter("");
}

EventFilterManager::~EventFilterManager()
{
}

void EventFilterManager::addFilterWidget(QString name, QObject* filter)
{
	m_filterMap[name] = filter;
}

QObject* EventFilterManager::getFilterWidget(QString name)
{
	return m_filterMap[name];
}

QObject* EventFilterManager::removeFilterWidget(QString name)
{
	QObject* filter = m_filterMap[name];
	m_filterMap.remove(name);
	if (m_currentFilter == name) {
		m_currentFilter = "";
	}
	return filter;
}

QString const EventFilterManager::getCurrentFilterName() const
{
	return m_currentFilter;
}

void EventFilterManager::setCurrentFilterName(QString name)
{
	if (m_filterMap.contains(name) || name == "") {
		m_currentFilter = name;
	}
	else {
		qDebug() << Q_FUNC_INFO << "The name \"" + name + "\" is not a name of " +
                               "any filter and is not an empty name";
	}
}

QList<QString> EventFilterManager::getFilterNames()
{
	return m_filterMap.uniqueKeys();
}

void EventFilterManager::clear()
{
	m_currentFilter = "";
	m_filterMap.clear();
}

bool EventFilterManager::eventFilter(QObject* target, QEvent* event)
{
	if (m_filterMap.contains(m_currentFilter)) {
		return m_filterMap[m_currentFilter]->eventFilter(target, event);
	}
	return false;
}
