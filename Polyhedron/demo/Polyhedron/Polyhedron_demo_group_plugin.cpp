#include "Polyhedron_demo_group_plugin.h"
#include "Polyhedron_demo_group_plugin.moc"
#include "Scene.h"

/**********
 * Models *
 **********/

QVariant GroupModel::data(const QModelIndex &index, int role) const
{
    if(index.row()<0 )
        return QVariant();

    if(role != Qt::DisplayRole )
        return QVariant();

    return groups()[index.row()]->name();
}
Qt::ItemFlags
GroupModel::flags ( const QModelIndex & index ) const
{
    if (index.isValid() && index.row() == 0) {
        return QAbstractListModel::flags(index) | ::Qt::ItemIsEditable;
    }
    else {
        return QAbstractListModel::flags(index);
    }
    return 0;
}

bool
GroupModel::setData(const QModelIndex &index,
               const QVariant &value,
               int role)
{

    if( role != ::Qt::EditRole || !index.isValid() )
      return false;

    if(index.row() != 0)
      return false;

    Scene_group_item* item = groups()[0];
    if(!item) return false;
    data_pool->setGroupName(value.toString());
    Q_EMIT dataChanged(index, index);
    return true;
}

QVariant ItemsModel::data(const QModelIndex &index, int role) const
{

    if(index.row()<0 || index.row() > items().size())
        return QVariant();

    if(role != Qt::DisplayRole)
        return QVariant();

    return items()[index.row()]->name();
}

/****************
 * Group Plugin *
 ****************/

void Polyhedron_demo_group_plugin::init(QMainWindow* mainWindow,
          Scene_interface* scene_interface,Messages_interface* m ) {
    //get the references
    trueScene = dynamic_cast<Scene*>(scene_interface);
    this->scene = scene_interface;
    this->mw = mainWindow;
    messages = m;
    group_model = new GroupModel(trueScene);
    items_model = new ItemsModel(trueScene);
    //creates and link the actions
    actionAddToGroup= new QAction("Add to group", mw);
    if(actionAddToGroup) {
        connect(actionAddToGroup, SIGNAL(triggered()),
                this, SLOT(GroupChoice()));
    }
}

void Polyhedron_demo_group_plugin::GroupChoice()
{
    Polyhedron_demo_group_plugin_dialog *dialog = new Polyhedron_demo_group_plugin_dialog(this);
    dialog->show();

}
/****************
 * Dialog Class *
 ****************/

Polyhedron_demo_group_plugin_dialog::Polyhedron_demo_group_plugin_dialog( Polyhedron_demo_group_plugin  *p_plugin)
{
    m_plugin= p_plugin;
    setupUi(this);
    itemsView->setModel(m_plugin->items_model);
    itemsView->setSelectionMode(QAbstractItemView::ExtendedSelection);
    groupsView->setModel(m_plugin->group_model);
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(selected_scene_items()));
    connect(buttonBox, SIGNAL(accepted()), m_plugin->trueScene, SLOT(group_added()));

}

void Polyhedron_demo_group_plugin_dialog::selected_scene_items()
{
    QModelIndexList selectedRows = itemsView->selectionModel()->selectedIndexes();
    QModelIndex selectedGroup= groupsView->selectionModel()->selectedIndexes().first();
    Q_FOREACH(QModelIndex index, selectedRows)
    {
        m_plugin->selected_items << m_plugin->items_model->items()[index.row()];
    }
    m_plugin->selected_group = m_plugin->group_model->groups()[selectedGroup.row()];
    m_plugin->add_to_group();

}
