#include "Scene_group_item.h"
#include <QDebug>


Scene_group_item::Scene_group_item(QString name)
    :  Scene_item(0,0)
{
    this->name_ = name;
}

bool Scene_group_item::isFinite() const
{
    Q_FOREACH(Scene_item *item, children)
        if(!item->isFinite()){
            return false;
        }
    return true;
}

bool Scene_group_item::isEmpty() const {
    Q_FOREACH(Scene_item *item, children)
        if(!item->isEmpty()){
            return false;
        }
    return true;
}

Scene_group_item::Bbox Scene_group_item::bbox() const
{
    double xmax=0, ymax=0, zmax=0;
    double xmin=0, ymin=0, zmin=0;
    if(!children.isEmpty())
    {
        xmax = children.first()->bbox().xmax; ymax = children.first()->bbox().ymax; zmax = children.first()->bbox().zmax;
        xmin = children.first()->bbox().xmin; ymin = children.first()->bbox().ymin; zmin = children.first()->bbox().zmin;
    }

    Q_FOREACH(Scene_item* item, children)
    {
        xmax = std::max(xmax,item->bbox().xmax); ymax = std::max(ymax, item->bbox().ymax); zmax = std::max(zmax, item->bbox().zmax);
        xmin = std::min(xmin, item->bbox().xmin);  ymin = std::min(ymin,item->bbox().ymin); zmin = std::min(zmin,item->bbox().zmin);
    }
    return Bbox(xmin, ymin, zmin, xmax, ymax, zmax);
}


bool Scene_group_item::supportsRenderingMode(RenderingMode m) const {
    Q_FOREACH(Scene_item* item, children)
        if(!item->supportsRenderingMode(m))
            return false;
    return !children.isEmpty();

}

QString Scene_group_item::toolTip() const {
    QString str =
            QObject::tr( "<p>Number of children: %1<br />").arg(children.size());
    str+="</p>";
    str += QString("Bounding box: min (%1,%2,%3), max(%4,%5,%6)")
            .arg(bbox().xmin)
            .arg(bbox().ymin)
            .arg(bbox().zmin)
            .arg(bbox().xmax)
            .arg(bbox().ymax)
            .arg(bbox().zmax);
    return str;
}

void Scene_group_item::addChild(Scene_item* new_item)
{  
        children.append(new_item);
}
