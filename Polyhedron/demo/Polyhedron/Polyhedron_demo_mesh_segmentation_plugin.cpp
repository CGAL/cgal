#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include "ui_Mesh_segmentation_dialog.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_with_color_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Surface_mesh_segmentation.h>
#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QAction>
#include <QDebug>
#include <QObject>
//#include <QtConcurrentRun>
#include <boost/property_map/property_map.hpp>

class Polyhedron_demo_mesh_segmentation_plugin : 
    public QObject,
    public Polyhedron_demo_plugin_helper
{
    Q_OBJECT
        Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:

    QList<QAction*> actions() const {
        return QList<QAction*>() << actionSegmentation;
    }

    bool applicable() const {
      return 
        qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
    }    
    
    void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
        this->scene = scene_interface;
        this->mw = mainWindow;
        actionSegmentation = new QAction("Mesh Segmentation", mw);
        if(actionSegmentation) {
            connect(actionSegmentation, SIGNAL(triggered()),this, SLOT(on_actionSegmentation_triggered()));
        }
    }
    
    template <class FacetValuePropertyMap>
    void colorize( Polyhedron* mesh, FacetValuePropertyMap facet_value_pmap, std::vector<QColor>& color_vector, bool sdf);
    public slots:
        void on_actionSegmentation_triggered();
        void on_Partition_button_clicked();
        void on_SDF_button_clicked();
private:
    QAction* actionSegmentation;

    QDialog* dialog;
    Ui::Mesh_segmentation_dialog* ui_dialog;
};


void Polyhedron_demo_mesh_segmentation_plugin::on_actionSegmentation_triggered()
{   
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    Scene_polyhedron_with_color_item* item_colored = qobject_cast<Scene_polyhedron_with_color_item*>(scene->item(index));

    if(!item && !item_colored) { return; }
    
    dialog = new QDialog(mw);
    dialog->setAttribute(Qt::WA_DeleteOnClose, true);
    //dialog->setModal(true);
    
    ui_dialog = new Ui::Mesh_segmentation_dialog();
    ui_dialog->setupUi(dialog);
    connect(ui_dialog->Partition_button,  SIGNAL(clicked()), this, SLOT(on_Partition_button_clicked()));   
    connect(ui_dialog->SDF_button,  SIGNAL(clicked()), this, SLOT(on_SDF_button_clicked()));  
    dialog->show();    
}
//void calculata_sdf_values(CGAL::Surface_mesh_segmentation<Polyhedron>* segmentation)
//{
//    segmentation->calculate_sdf_values();
//}
void Polyhedron_demo_mesh_segmentation_plugin::on_SDF_button_clicked()
{
    Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    Scene_polyhedron_with_color_item* item_colored = qobject_cast<Scene_polyhedron_with_color_item*>(scene->item(index));

    if(!item && !item_colored) { return; }
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    int number_of_rays = ui_dialog->Number_of_rays_spin_box->value();
    double cone_angle = (ui_dialog->Cone_angle_spin_box->value()  / 180.0) * CGAL_PI;    
    
    CGAL::Surface_mesh_segmentation<Polyhedron>* segmentation = NULL;
    if(!item_colored)
    {    
        Polyhedron* pMesh = item->polyhedron();
        item->setVisible(false);

        Scene_polyhedron_with_color_item* new_item = new Scene_polyhedron_with_color_item(*pMesh);
       
        new_item->setName(tr("%1_segmented").arg(item->name()));
        std::vector<QColor> color_vector;

        CGAL::Surface_mesh_segmentation<Polyhedron>* segmentation
            = new CGAL::Surface_mesh_segmentation<Polyhedron>(*new_item->polyhedron());	
        segmentation->calculate_sdf_values(new_item->sdf_pmap, cone_angle, number_of_rays);
      
        colorize(new_item->polyhedron(), new_item->sdf_pmap, color_vector, true);
        
        new_item->set_color_vector(color_vector);
        new_item->segmentation = segmentation;
        
        Scene_interface::Item_id new_item_index = scene->addItem(new_item);
        scene->setSelectedItem(new_item_index);
        scene->itemChanged(new_item_index);  
    }
    else
    {     
        std::vector<QColor> color_vector;
        CGAL::Surface_mesh_segmentation<Polyhedron>* segmentation = item_colored->segmentation;

        segmentation->calculate_sdf_values(item_colored->sdf_pmap, cone_angle, number_of_rays);

        colorize(item_colored->polyhedron(), item_colored->sdf_pmap, color_vector, true);
        item_colored->set_color_vector(color_vector);
        scene->itemChanged(index);
        scene->setSelectedItem(index);                 
    }
    QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mesh_segmentation_plugin::on_Partition_button_clicked()
{    
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_polyhedron_with_color_item* item_colored = qobject_cast<Scene_polyhedron_with_color_item*>(scene->item(index));

    if(!item_colored) { return; }
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    int number_of_clusters = ui_dialog->Number_of_clusters_spin_box->value();   
    double smoothness = ui_dialog->Smoothness_spin_box->value();
    
    std::vector<QColor> color_vector;
    CGAL::Surface_mesh_segmentation<Polyhedron>* segmentation = item_colored->segmentation;
    
    typedef std::map<CGAL::Surface_mesh_segmentation<Polyhedron>::Facet_const_iterator, int> internal_map;
    internal_map facet_value_map_internal;
    boost::associative_property_map<internal_map> partition_pmap(facet_value_map_internal);
    
    segmentation->partition(item_colored->sdf_pmap, partition_pmap, number_of_clusters, smoothness);
    
    colorize(item_colored->polyhedron(), partition_pmap, color_vector, false);
    item_colored->set_color_vector(color_vector);
    scene->itemChanged(index);
    scene->setSelectedItem(index);
    
    QApplication::restoreOverrideCursor();
    //qApp->processEvents();
}

template <class FacetValuePropertyMap>
void Polyhedron_demo_mesh_segmentation_plugin::colorize(
     Polyhedron* mesh,
     FacetValuePropertyMap facet_value_pmap,
     std::vector<QColor>& color_vector, 
     bool sdf)
{
    QColor predefined_list[] = { Qt::red, Qt::darkRed, Qt::green, Qt::darkGreen, Qt::blue, Qt::darkBlue
        , Qt::cyan, Qt::darkCyan, Qt::magenta, Qt::darkMagenta, Qt::yellow, Qt::darkYellow, Qt::gray, Qt::darkGray, Qt::lightGray, Qt::white, Qt::black};
    
    for(Polyhedron::Facet_const_iterator facet_it = mesh->facets_begin(); 
        facet_it != mesh->facets_end(); ++facet_it)   
    {
        if(sdf)
        {
            int sdf_color = static_cast<int>(255.0 * boost::get(facet_value_pmap, facet_it));
            color_vector.push_back(QColor(sdf_color,sdf_color,sdf_color));            
        }
        else
        {
            color_vector.push_back(predefined_list[ static_cast<int>(boost::get(facet_value_pmap, facet_it)) % 17]);
        }        
    }
}


Q_EXPORT_PLUGIN2(Polyhedron_demo_mesh_segmentation_plugin, Polyhedron_demo_mesh_segmentation_plugin)

#include "Polyhedron_demo_mesh_segmentation_plugin.moc"
