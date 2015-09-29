#include "Polyhedron_demo_group_plugin.h"
#include "Polyhedron_demo_group_plugin.moc"
#include "Scene.h"

/**********
 * Models *
 **********/

QVariant GroupModel::data(const QModelIndex &index, int role) const
{


    if(index.row()<0 || index.row()>groups().size())
        return QVariant();

    if(role != Qt::DisplayRole)
        return QVariant();

    return groups()[index.row()]->name();
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
          Scene_interface* scene_interface) {
    //get the references
    Scene* trueScene = dynamic_cast<Scene*>(scene_interface);
    this->scene = scene_interface;
    this->mw = mainWindow;
    group_model = new GroupModel(trueScene);
    items_model = new ItemsModel(trueScene);

    //creates and link the actions
    actionAddToGroup= new QAction("Add to group", mw);
    if(actionAddToGroup) {
        connect(actionAddToGroup, SIGNAL(triggered()),
                this, SLOT(draw_triangle()));
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
    groupsView->setModel(m_plugin->group_model);
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(selected_scene_items()));
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
