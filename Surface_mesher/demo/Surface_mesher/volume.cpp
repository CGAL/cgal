#include <CGAL/basic.h>

#include "volume.h"

#include  <algorithm> // std::sort
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <boost/foreach.hpp>

#include <CGAL/Bbox_3.h>

#include "viewer.h"
#include "mainwindow.h"
#include "values_list.h"

#include "File_XT.h" // format XT from Total/ELF

#include <QApplication>
#include <QFileDialog>
#include <QAction>
#include <QStatusBar>
#include <QDoubleSpinBox>
#include <QMessageBox>
#include <QTreeWidgetItem>
#include <QTime>
#include <QColor>
#include <QColorDialog>
#include <QSettings>
#include <QUrl>
#include "Raw_image_dialog.h"

#include <CGAL/glu.h>

#include <CGAL/Surface_mesher/Standard_criteria.h>
// #include <CGAL/Surface_mesher/Image_surface_oracle_3.h>
#include <CGAL/Surface_mesher/Implicit_surface_oracle_3.h>
#include <CGAL/Surface_mesher/Vertices_on_the_same_psc_element_criterion.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <CGAL/make_surface_mesh.h>

struct Threshold : public std::unary_function<FT, unsigned char> {
  double isovalue;
  bool is_identity;

  Threshold(double isovalue) : isovalue(isovalue), is_identity(false) {}

  result_type operator()(FT value)
  {
    if(is_identity)
      return static_cast<unsigned char>(value);
    else if(value >=  isovalue)
      return 1;
    else
      return 0;
  }
};

class Classify_from_isovalue_list :
  public std::unary_function<FT, unsigned char> 
{
  typedef std::pair<FT, result_type> Isovalue;
  typedef std::vector<Isovalue> Isovalues;
  boost::shared_ptr<Isovalues> isovalues;
  bool is_identity;

  struct Sort_isovalues : std::binary_function<Isovalue, Isovalue, bool> 
  {
    bool operator()(const Isovalue& isoval1, const Isovalue& isoval2)
    {
      return isoval1.first < isoval2.first;
    }
  };
public:
  Classify_from_isovalue_list(Values_list * list)
    : is_identity(false)
  {
    isovalues = boost::shared_ptr<Isovalues>(new Isovalues(list->numberOfValues()));
    for(int i = 0, nbs = list->numberOfValues(); i < nbs; ++i )
      (*isovalues)[i] = std::make_pair(list->value(i), i);
    std::sort(isovalues->begin(), isovalues->end(), Sort_isovalues());
  }

  void set_identity(bool b) {
    is_identity = b;
  }

  result_type operator()(FT value)
  {
    if(is_identity) {
      return static_cast<unsigned char>(value);
    }
    result_type result = 0;
//     std::cerr << "isovalues: ";
    for(int i = 1, end = isovalues->size(); i <= end; ++i)
    {
//       std::cerr << (*isovalues)[i-1] << ", ";
      if(value >= (*isovalues)[i-1].first &&
         i >= result)
      {
        result = i;
      }
    }
//     if(result>1)
//       std::cerr << "result = "  << (int)result << "/" << list->numberOfValues() << std::endl;
//     else
//       std::cerr << std::endl;
    if(result>0)
      return (*isovalues)[result-1].second + 1;
    else
      return 0;
  }
};

class Generate_surface_identifiers :
  public std::binary_function<Classify_from_isovalue_list::result_type,
                              Classify_from_isovalue_list::result_type,
                              const QTreeWidgetItem*>
{
  Values_list* list;
  bool labellized;

public:
  Generate_surface_identifiers(Values_list* list)
    : list(list), labellized(false) {};

  void set_labellized_image(bool b)
  {
    labellized = b;
  }

  result_type operator()(const Classify_from_isovalue_list::result_type& a,
                         const Classify_from_isovalue_list::result_type& b)
  {
    if(labellized)
      return list->search((std::max)(a, b));
    else
      return list->item((std::min)(a, b));
  }
};

// class Classify_from_isovalue_list :
//   public std::unary_function<FT, const QTreeWidgetItem*> 
// {
//   typedef std::pair<FT, result_type> Isovalue;
//   typedef std::vector<Isovalue> Isovalues;
//   boost::shared_ptr<Isovalues> isovalues;

//   struct Sort_isovalues : std::binary_function<Isovalue, Isovalue, bool> 
//   {
//     bool operator()(const Isovalue& isoval1, const Isovalue& isoval2)
//     {
//       return isoval1.first < isoval2.first;
//     }
//   };
// public:
//   Classify_from_isovalue_list(Isovalues_list * list)
//   {
//     isovalues = boost::shared_ptr<Isovalues>(new Isovalues(list->numberOfIsoValues()));
//     for(int i = 0, nbs = list->numberOfIsoValues(); i < nbs; ++i )
//       (*isovalues)[i] = std::make_pair(list->isovalue(i), list->item(i));
//     std::sort(isovalues->begin(), isovalues->end(), Sort_isovalues());
//   }

//   result_type operator()(FT value)
//   {
//     int result = 0;
// //     std::cerr << "isovalues: ";
//     for(int i = 1, end = isovalues->size(); i <= end; ++i)
//     {
// //       std::cerr << (*isovalues)[i-1] << ", ";
//       if(value >= (*isovalues)[i-1].first &&
//          i >= result)
//       {
//         result = i;
//       }
//     }
//     if(result>1)
//       std::cerr << boost::format("result = %1%/%2%\n") % result % isovalues->size();
//     if(result>0)
//       return (*isovalues)[result-1].second;
//     else
//       return 0;
//   }
// };
Volume::Volume(MainWindow* mw) : 
  Surface(mw),
  m_sm_angle(30),
  m_sm_radius(0),
  m_sm_distance(0),
  m_relative_precision(0.000000001),
  m_view_surface(false),
  m_triangulation_color(QColor(Qt::green)),
  m_inverse_normals(false),
  two_sides(false),
  del(),
  c2t3(del),
  mw(mw),
  lists_draw_surface(),
  lists_draw_surface_is_valid(false),
#ifdef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
  list_draw_marching_cube(0),
  list_draw_marching_cube_is_valid(false),
  lists_draw_surface_mc(),
#endif // CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
  m_view_mc(false)
{
  spinBox_radius_bound = mw->findChild<QDoubleSpinBox*>("spinBox_radius_bound");
  spinBox_distance_bound = mw->findChild<QDoubleSpinBox*>("spinBox_distance_bound");
  Q_ASSERT_X(spinBox_radius_bound && spinBox_distance_bound,
             "Volume::Volume()", "Cannot find spinboxes!");

  values_list = mw->values;

  connect(spinBox_radius_bound, SIGNAL(valueChanged(double)),
          this, SLOT(set_radius_bound(double)));
  connect(spinBox_distance_bound, SIGNAL(valueChanged(double)),
          this, SLOT(set_distance_bound(double)));

  connect(mw->actionSurface_mesher, SIGNAL(triggered()),
          this, SLOT(display_surface_mesher_result()));

  connect(mw->actionInverse_normals, SIGNAL(toggled(bool)),
          this, SLOT(set_inverse_normals(bool)));
  m_inverse_normals = mw->actionInverse_normals->isChecked();

  connect(mw->actionDisplay_front_and_back, SIGNAL(toggled(bool)),
          this, SLOT(set_two_sides(bool)));
  two_sides = mw->actionDisplay_front_and_back->isChecked();

  connect(mw->actionDraw_triangles_edges, SIGNAL(toggled(bool)),
          this, SLOT(set_draw_triangles_edges(bool)));
  draw_triangles_edges = mw->actionDraw_triangles_edges->isChecked();

  connect(mw->actionUse_Gouraud_shading, SIGNAL(toggled(bool)),
          this, SLOT(set_use_gouraud(bool)));
  use_gouraud = mw->actionUse_Gouraud_shading->isChecked();

  connect(mw->actionShow_the_image_bounding_box, SIGNAL(toggled(bool)),
          this, SLOT(set_show_bbox(bool)));
  show_bbox = mw->actionShow_the_image_bounding_box->isChecked();

  connect(mw->actionShow_triangulation, SIGNAL(toggled(bool)),
          this, SLOT(set_draw_triangulation(bool)));
  m_draw_triangulation = mw->actionShow_triangulation->isChecked();

  connect(mw->actionTriangulation_edges_color, SIGNAL(triggered()),
          this, SLOT(set_triangulation_edges_color()));

  connect(this, SIGNAL(new_bounding_box(double, double, double, double, double, double)),
          mw->viewer, SLOT(interpolateToFitBoundingBox(double, double, double, double, double, double)));

  connect(values_list, SIGNAL(values_changed()),
          this, SLOT(changed_parameters()));
  connect(values_list, SIGNAL(changed()),
          mw->viewer, SLOT(updateGL()));
  connect(this, SIGNAL(changed()),
          this, SLOT(check_can_export_off()));

  connect(mw->labellizedRadioButton, SIGNAL(toggled(bool)),
	  this, SLOT(labellizedToogled(bool)));

  mw->actionExport_surface_mesh_to_OFF->setEnabled(false);
  connect(mw->actionExport_surface_mesh_to_OFF, SIGNAL(triggered()),
          this, SLOT(export_off()));

  connect(mw->actionSave, SIGNAL(triggered()),
          this, SLOT(save_image_to_inr()));

#ifdef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
  connect(mw->actionMarching_cubes, SIGNAL(triggered()),
          this, SLOT(display_marchin_cube()));
#endif
}

void Volume::set_inverse_normals(const bool b) {
  m_inverse_normals = b;

  list_draw_marching_cube = 0; // Invalidate the display list for the
                               // marching cube. See gl_draw_marchingcube()
                               // for an explanation.

  emit changed();
}

void Volume::set_two_sides(const bool b) {
  two_sides = b;
  emit changed();
}

void Volume::set_draw_triangles_edges(const bool b) {
  draw_triangles_edges = b;
  emit changed();
}

void Volume::set_draw_triangulation(const bool b) {
  m_draw_triangulation = b;
  emit changed();
}

void Volume::set_triangulation_edges_color() {
  const QColor color = QColorDialog::getColor(m_triangulation_color, mw);
  if (color.isValid()) {
    m_triangulation_color = color;
    emit changed();
  }
}

void Volume::set_use_gouraud(const bool b) {
  use_gouraud = b;
  emit changed();
}

void Volume::set_show_bbox(const bool b) {
  show_bbox = b;
  emit changed();
}

void Volume::only_in()
{
  mw->show_only("volume");
#ifndef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
  mw->actionMarching_cubes->setVisible(false);
#endif
}

#ifdef CGAL_USE_VTK
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageReader.h>
#include <vtkImageGaussianSmooth.h>

bool Volume::opendir(const QString& dirname) 
{
  bool result = true;
  if(!fileinfo.isReadable())
  {
    QMessageBox::warning(mw, mw->windowTitle(),
                         tr("Cannot read directory <tt>%1</tt>!").arg(dirname));
    status_message(tr("Opening of directory %1 failed!").arg(dirname));
    result = false;
  }
  else
  {
    vtkDICOMImageReader* dicom_reader = vtkDICOMImageReader::New();
    dicom_reader->SetDirectoryName(dirname.toUtf8());
    vtkImageGaussianSmooth* smoother = vtkImageGaussianSmooth::New();
    smoother->SetStandardDeviations(1., 1., 1.);
    smoother->SetInputConnection(dicom_reader->GetOutputPort());
    smoother->Update();
    vtkImageData* vtk_image = smoother->GetOutput();
    dicom_reader->SetReleaseDataFlag(false);
    vtk_image->SetReleaseDataFlag(false);
    vtk_image->Print(std::cerr);
    if(!m_image.read_vtk_image_data(vtk_image))
    {
      QMessageBox::warning(mw, mw->windowTitle(),
                           tr("Error with file <tt>%1/</tt>:\nunknown file format!").arg(dirname));
      status_message(tr("Opening of file %1/ failed!").arg(dirname));
      result = false;
    }
    else
    {
      status_message(tr("File %1/ successfully opened.").arg(dirname));
      finish_open();
      result = true;
    }
    dicom_reader->Delete();
    // smoother->Delete();
  }
  return result;
}

bool Volume::open_vtk(const QString& filename)
{
  only_in();

  fileinfo.setFile(filename);

  if(fileinfo.isDir())
  {
    return opendir(filename);
  }

  if(!fileinfo.isReadable())
  {
    QMessageBox::warning(mw, mw->windowTitle(),
                         tr("Cannot read file <tt>%1</tt>!").arg(filename));
    status_message(tr("Opening of file %1 failed!").arg(filename));
    return false;
  }
  else
  {
    vtkImageReader* vtk_reader = vtkImageReader::New();
    vtk_reader->SetFileName(filename.toUtf8());
    vtk_reader->SetDataScalarTypeToUnsignedChar();
    vtk_reader->SetDataExtent(0, 249, 0, 249, 0,  124);
    vtk_reader->SetDataSpacing(1., 1., 1.);
    vtk_reader->SetFileDimensionality(3);
    vtk_reader->Update();
    vtk_reader->Print(std::cerr);
    vtkImageData* vtk_image = vtk_reader->GetOutput();
    vtk_image->Print(std::cerr);
    if(!m_image.read_vtk_image_data(vtk_image))
    {
      QMessageBox::warning(mw, mw->windowTitle(),
                           tr("Error with file <tt>%1</tt>:\nunknown file format!").arg(filename));
      status_message(tr("Opening of file %1 failed!").arg(filename));
      return false;
    }
    else
    {
      status_message(tr("File %1 successfully opened.").arg(filename));
      finish_open();
      return true;
    }
  }
}

// Total 3D images (XT format, that is the old Inrimage format, 1994.
bool Volume::open_xt(const QString& filename)
{
  only_in();

  fileinfo.setFile(filename);

  if(fileinfo.isDir())
  {
    return false;
  }

  if(!fileinfo.isReadable())
  {
    QMessageBox::warning(mw, mw->windowTitle(),
                         tr("Cannot read file <tt>%1</tt>!").arg(filename));
    status_message(tr("Opening of file %1 failed!").arg(filename));
    return false;
  }
  else
  {
    long dimx, dimy, dimz;
    long word_dim;
    long header_size;
    const char* filename_stl = qPrintable(filename);
    CGAL::Total::lire_longueur_entete(filename_stl, &header_size);
    CGAL::Total::lire_nb_octet(filename_stl, &word_dim);
    CGAL::Total::lire_longueur_trace(filename_stl, &dimx);
    CGAL::Total::lire_nb_trace(filename_stl, &dimy);
    CGAL::Total::lire_nb_plan(filename_stl, &dimz);

    vtkImageReader* vtk_reader = vtkImageReader::New();
    vtk_reader->SetFileName(filename_stl);
    switch(word_dim) {
    case 8:
      vtk_reader->SetDataScalarTypeToUnsignedChar();
      break;
    case 16:
      vtk_reader->SetDataScalarTypeToUnsignedShort();
      break;
    default:
      return false;
    }
    vtk_reader->SetHeaderSize(header_size);
    vtk_reader->SetDataExtent(1, dimx, 1, dimy, 1,  dimz);
    vtk_reader->SetDataSpacing(1., 1., 1.);
    vtk_reader->SetFileDimensionality(3);
    vtk_reader->Update();
    vtk_reader->Print(std::cerr);
    vtkImageData* vtk_image = vtk_reader->GetOutput();
    vtk_image->Print(std::cerr);
    if(!m_image.read_vtk_image_data(vtk_image))
    {
      QMessageBox::warning(mw, mw->windowTitle(),
                           tr("Error with file <tt>%1</tt>:\nunknown file format!").arg(filename));
      status_message(tr("Opening of file %1 failed!").arg(filename));
      return false;
    }
    else
    {
      status_message(tr("File %1 successfully opened.").arg(filename));
      finish_open();
    }
    return true;
  }
}

#else // CGAL_USE_VTK
bool Volume::opendir(const QString&)
{
  return false;
}

bool Volume::open_xt(const QString&)
{
  return false;
}
#endif // CGAL_USE_VTK

bool Volume::open(const QString& filename)
{
  only_in();

  fileinfo.setFile(filename);

  if(fileinfo.isDir())
  {
    return opendir(filename);
  }

  if(!fileinfo.isReadable())
  {
    QMessageBox::warning(mw, mw->windowTitle(),
                         tr("Cannot read file <tt>%1</tt>!").arg(filename));
  }
  else
  {
    if(m_image.read(filename.toStdString().c_str()))
    {
      status_message(tr("File %1 successfully opened.").arg(filename));
      finish_open();
      return true;
    }
    else if(open_xt(filename)) {
      return true;
    }
    else 
    {
      QSettings settings;
      settings.beginGroup(QUrl::toPercentEncoding(fileinfo.absoluteFilePath()));
      if( settings.value("is_raw").toBool() &&
	  m_image.read_raw(filename.toStdString().c_str(),
			   settings.value("dim_x").toInt(),
			   settings.value("dim_y").toInt(),
			   settings.value("dim_z").toInt(),
			   settings.value("spacing_x").toDouble(),
			   settings.value("spacing_y").toDouble(),
			   settings.value("spacing_z").toDouble(),
			   settings.value("offset").toInt()) )
      {
	status_message(tr("File %1 successfully opened.").arg(filename));
	finish_open();
	return true;
      }
      else if(QMessageBox::warning(mw, mw->windowTitle(),
				   tr("Error with file <tt>%1</tt>:\n"
				      "unknown file format!\n"
				      "\n"
				      "Open it as a raw image?").arg(filename),
				   QMessageBox::Yes|QMessageBox::No) == QMessageBox::Yes) 
      {
	Raw_image_dialog raw_dialog;
	raw_dialog.label_file_size->setText(QString("%1 B").arg(fileinfo.size()));
	if( raw_dialog.exec() && 
	    m_image.read_raw(filename.toStdString().c_str(),
			     raw_dialog.dim_x->value(),
			     raw_dialog.dim_y->value(),
			     raw_dialog.dim_z->value(),
			     raw_dialog.spacing_x->value(),
			     raw_dialog.spacing_y->value(),
			     raw_dialog.spacing_z->value(),
			     raw_dialog.offset->value()) )
	{
	  status_message(tr("File %1 successfully opened.").arg(filename));
	  QSettings settings;
	  settings.beginGroup(QUrl::toPercentEncoding(fileinfo.absoluteFilePath()));
	  settings.setValue("is_raw", true);
	  settings.setValue("dim_x", raw_dialog.dim_x->value());
	  settings.setValue("dim_y", raw_dialog.dim_y->value());
	  settings.setValue("dim_z", raw_dialog.dim_z->value());
	  settings.setValue("spacing_x", raw_dialog.spacing_x->value());
	  settings.setValue("spacing_y", raw_dialog.spacing_y->value());
	  settings.setValue("spacing_z", raw_dialog.spacing_z->value());
	  settings.setValue("offset", raw_dialog.offset->value());
	  settings.endGroup();
	  finish_open();
	  return true;
	}
      }
    }
  }
  status_message(tr("Opening of file %1 failed!").arg(filename));
  return false;
}

void Volume::finish_open()
{
  m_image.finish_open();
  mw->viewer->camera()->setSceneBoundingBox(qglviewer::Vec(0, 0, 0),
                                            qglviewer::Vec(m_image.xmax(),
                                                           m_image.ymax(),
                                                           m_image.zmax()));

  mw->viewer->showEntireScene();
  values_list->load_values(fileinfo.absoluteFilePath());
  load_image_settings(fileinfo.absoluteFilePath());
  changed_parameters();
  emit changed();
}

void Volume::export_off()
{
  QFileDialog filedialog(mw, tr("Export surface to file"));
  filedialog.setFileMode(QFileDialog::AnyFile);
  filedialog.setFilter(tr("OFF files (*.off);;"
                          "All files (*)"));
  filedialog.setAcceptMode(QFileDialog::AcceptSave);
  filedialog.setDefaultSuffix("off");
  if(filedialog.exec())
  {
    const QString filename = filedialog.selectedFiles().front();
    std::cerr << "Saving to file \"" << filename.toLocal8Bit().data() << "\"...";
    std::ofstream out(filename.toUtf8());
    CGAL::output_surface_facets_to_off(out, c2t3);
    if(!out)
    {
      QMessageBox::warning(mw, mw->windowTitle(),
                           tr("Export to the OFF file <tt>%1</tt> failed!").arg(filename));
      status_message(tr("Export to the OFF file %1 failed!").arg(filename));
      std::cerr << " failed!\n";
    }
    else
    {
      std::cerr << " done.\n";
      status_message(tr("Successfull export to the OFF file %1.").arg(filename));
    }
  }
}

void Volume::save_image_to_inr()
{
  QFileDialog filedialog(mw, tr("Export image to Inrimage format"));
  filedialog.setFileMode(QFileDialog::AnyFile);
  filedialog.setFilter(tr("Inrimage files (*.inr);;"
                          "Compressed Inrimage files (*.inr.gz)"));
  filedialog.setAcceptMode(QFileDialog::AcceptSave);
  filedialog.setDefaultSuffix("inr.gz");
  if(filedialog.exec())
  {
    const QString filename = filedialog.selectedFiles().front();
    std::cerr << "Saving image to file \"" << filename.toLocal8Bit().data() << "\"...";
    const int result = ::_writeImage(m_image.image(), filename.toUtf8());
    if(result != ImageIO_NO_ERROR)
    {
      QMessageBox::warning(mw, mw->windowTitle(),
                           tr("Export to the Inrimage file <tt>%1</tt> failed!").arg(filename));
      status_message(tr("Export to the Inrimage file %1 failed!").arg(filename));
      std::cerr << " failed!\n";
    }
    else
    {
      std::cerr << " done.\n";
      status_message(tr("Successfull export to the Inrimage file %1.").arg(filename));
    }
  }
}

void Volume::check_can_export_off()
{
  mw->actionExport_surface_mesh_to_OFF->setEnabled(m_view_surface);// || m_view_mc);
}

void Volume::status_message(QString string)
{
  std::cerr << qPrintable(string) << std::endl;
  mw->statusBar()->showMessage(string);
}

void Volume::busy() const 
{
  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
}

void Volume::not_busy() const 
{
  QApplication::restoreOverrideCursor();
}

void Volume::display_marchin_cube()
{
#ifdef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
  if(m_surface_mc.empty())
  {
    QTime total_time;
    total_time.start();

    values_list->save_values(fileinfo.absoluteFilePath());

    unsigned int nx = m_image.xdim();
    unsigned int ny = m_image.ydim();
    unsigned int nz = m_image.zdim();
    if(nx * ny * nz == 0)
    {
      status_message("No volume loaded.");
      return;
    }

    mc_timer.reset();
    busy();
    status_message("Marching cubes...");

    mc_timer.start();
    m_surface_mc.clear();

    if(mc.ntrigs()!=0)
      mc.clean_all();
    mc.set_resolution(nx,ny,nz);
    mc.init_all();
    mc.set_ext_data(static_cast<unsigned char*>(m_image.image()->data));

    nbs_of_mc_triangles.resize(values_list->numberOfValues());

    for(int value_id = 0; 
        value_id < values_list->numberOfValues();
        ++value_id)
    {
      status_message(tr("Marching cubes, isovalue #%1...").arg(value_id));

      // set data
//       for(unsigned int i=0;i<nx;i++)
//         for(unsigned int j=0;j<ny;j++)
//           for(unsigned int k=0;k<nz;k++)
//           {
//             const float& value = m_image.value(i,j,k);
//             mc.set_data(value,i,j,k);
//           }
      // compute scaling ratio
      if(value_id > 0)
        mc.init_temps();
      mc.run(values_list->value(value_id),
             m_image.vx(),
             m_image.vy(),
             m_image.vz());
      mc.clean_temps();

      std::vector<double> facets;
      mc.get_facets(facets);

      mc_timer.stop();
      const unsigned int begin = value_id == 0 ? 0 : nbs_of_mc_triangles[value_id-1];
      const unsigned int nbt = facets.size() / 9;
      for(unsigned int i=begin;i<nbt;i++)
      {
        const Point a(facets[9*i],   facets[9*i+1], facets[9*i+2]);
        const Point b(facets[9*i+3], facets[9*i+4], facets[9*i+5]);
        const Point c(facets[9*i+6], facets[9*i+7], facets[9*i+8]);
        const Triangle_3 t(a,b,c);
        const Vector u = t[1] - t[0];
        const Vector v = t[2] - t[0];
        Vector n = CGAL::cross_product(u,v);
        n = n / std::sqrt(n*n);
        m_surface_mc.push_back(Facet(t,n,values_list->item(value_id)));
      }
      nbs_of_mc_triangles[value_id]=m_surface_mc.size();
      mc_timer.start();
    }
    mc_timer.stop();
    not_busy();
    mc_total_time = total_time.elapsed();

    // invalidate the display list
    lists_draw_surface_mc_is_valid = false;
    list_draw_marching_cube_is_valid = false;
  }
  CGAL::Bbox_3 bbox(0,0,0,0,0,0);
  for(std::vector<Facet>::const_iterator
        it = m_surface_mc.begin(), end = m_surface_mc.end();
      it != end; ++it)
  {
    bbox = bbox + it->get<0>().bbox();
  }

  m_view_mc = true;
  m_view_surface = false;
  emit changed();
  if(!m_surface_mc.empty())
  {
    emit new_bounding_box(bbox.xmin(),
                          bbox.ymin(),
                          bbox.zmin(),
                          bbox.xmax(),
                          bbox.ymax(),
                          bbox.zmax());
  }
  status_message(tr("Marching cubes done. %2 facets in %1s (CPU time), total time is %3s.")
                 .arg(mc_timer.time())
                 .arg(m_surface_mc.size())
                 .arg(mc_total_time/1000.));

  save_image_settings(fileinfo.absoluteFilePath());
#endif // CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
}

void Volume::display_surface_mesher_result()
{
  if(m_surface.empty() || // Either the surface is not computed.
     m_view_surface) // Or it is computed and displayed, and one want
                     // to recompute it.
  {
    QTime total_time;
    total_time.start();

    values_list->save_values(fileinfo.absoluteFilePath());

    unsigned int nx = m_image.xdim();
    unsigned int ny = m_image.ydim();
    unsigned int nz = m_image.zdim();
    if(nx * ny * nz == 0)
    {
      status_message("No volume loaded.");
      return;
    }

    m_surface.clear();
    sm_timer.reset();
    busy();

    status_message("Surface meshing...");

    sm_timer.start();

    c2t3.clear();
    del.clear();
    Sphere bounding_sphere(m_image.center(),m_image.radius()*m_image.radius());

    Classify_from_isovalue_list classify(values_list);
    Generate_surface_identifiers generate_ids(values_list);

    m_image.set_interpolation(mw->interpolationCheckBox->isChecked());
    if(mw->labellizedRadioButton->isChecked()) {
      std::cerr << "Labellized image\n";
    }
    m_image.set_labellized(mw->labellizedRadioButton->isChecked());
    classify.set_identity(mw->labellizedRadioButton->isChecked());
    generate_ids.set_labellized_image(mw->labellizedRadioButton->isChecked());

    // definition of the surface
    Surface_3 surface(m_image, bounding_sphere, m_relative_precision);
//     Threshold threshold(m_image.isovalue());

    // surface mesh traits class
    typedef CGAL::Surface_mesher::Implicit_surface_oracle_3<Kernel,
      //     typedef CGAL::Surface_mesher::Image_surface_oracle_3<Kernel,
      Surface_3, 
      Classify_from_isovalue_list,
      Generate_surface_identifiers> Oracle;
    Oracle oracle(classify, generate_ids);

    if(mw->searchSeedsCheckBox->isChecked())
    {
      typedef std::vector<std::pair<Point, double> > Seeds;
      Seeds seeds;
      {
	std::cerr << "Search seeds...\n";
	std::set<unsigned char> domains;
	search_for_connected_components(std::back_inserter(seeds),
					CGAL::inserter(domains),
					classify);
	std::cerr << "Found " << seeds.size() << " seed(s).\n";

	if(mw->labellizedRadioButton->isChecked() && 
	   values_list->numberOfValues() == 0) 
	{
	  Q_FOREACH(unsigned char label, domains) {
	    if(label != 0) {
	      values_list->addValue(label);
	    }
	  }
	}
      }
      std::ofstream seeds_out("seeds.off");
      std::ofstream segments_out("segments.txt");
      seeds_out.precision(18);
      seeds_out << "OFF\n" << seeds.size() << " 0 0\n";
      segments_out.precision(18);
      for(Seeds::const_iterator it = seeds.begin(), end = seeds.end();
	  it != end; ++it)
      {
        seeds_out << it->first << std::endl;
	CGAL::Random_points_on_sphere_3<Point> random_points_on_sphere_3(it->second);
	Oracle::Intersect_3 intersect = oracle.intersect_3_object();
	for(int i = 0; i < 20; ++i)
	{
	  const Point test = it->first + (*random_points_on_sphere_3++ - CGAL::ORIGIN);
	  CGAL::Object o = intersect(surface, Segment_3(it->first, test));
	  if (const Point* intersection = CGAL::object_cast<Point>(&o)) {
            segments_out << "2 " << it->first << " " << *intersection << std::endl;
	    del.insert(*intersection);
          }
	  else 
	  {
	    std::cerr << 
	      boost::format("Error. Segment (%1%, %2%) does not intersect the surface! values=(%3%, %4%)\n")
	      % it->first % test
	      % surface(it->first) % surface(test);
	  }
	}
      }
    }
    else {
      oracle.construct_initial_points_object()(surface, 
					       CGAL::inserter(c2t3.triangulation()),
					       20);
    }

    std::ofstream points_out("initial-points.off");
    points_out.precision(18);
    points_out << "OFF\n" << c2t3.triangulation().number_of_vertices() << " 0 0\n";
    BOOST_FOREACH(const Tr::Vertex& v,
                  std::make_pair(c2t3.triangulation().vertices_begin(),
                                 c2t3.triangulation().vertices_end()))
    {
      points_out << v.point() << std::endl;
    }

    std::cerr << boost::format("Number of initial points: %1%\n") % del.number_of_vertices();

    // defining meshing criteria
    typedef CGAL::Surface_mesher::Refine_criterion<Tr> Criterion;
    CGAL::Surface_mesher::Curvature_size_criterion<Tr>
      curvature_size_criterion (m_sm_distance);
    CGAL::Surface_mesher::Uniform_size_criterion<Tr>
      uniform_size_criterion (m_sm_radius);
    CGAL::Surface_mesher::Aspect_ratio_criterion<Tr>
      aspect_ratio_criterion (m_sm_angle);
    CGAL::Surface_mesher::Vertices_on_the_same_psc_element_criterion<Tr, Surface_3>
      vertices_on_the_same_psc_element_criterion(surface);
    
    std::vector<Criterion*> criterion_vector;
    criterion_vector.push_back(&aspect_ratio_criterion);
    criterion_vector.push_back(&uniform_size_criterion);
    criterion_vector.push_back(&curvature_size_criterion);
    if(mw->sameIndexCheckBox->isChecked()) {
      criterion_vector.push_back(&vertices_on_the_same_psc_element_criterion);
      std::cerr << "vertices_on_the_same_psc_element_criterion is activated.\n";
    }

    typedef CGAL::Surface_mesher::Standard_criteria<Criterion> Criteria;
    Criteria criteria(criterion_vector);
    std::cerr << "Surface_mesher... angle=" << m_sm_angle << ", radius= " << m_sm_radius
              << ", distance=" << m_sm_distance << "\n";

    typedef CGAL::Surface_mesher_generator<C2t3,
      Oracle,
      Criteria,
      CGAL::Manifold_tag,
      CGAL_SURFACE_MESHER_VERBOSITY
      >::type Surface_mesher_manifold;
      
    typedef CGAL::Surface_mesher_generator<C2t3,
      Oracle,
      Criteria,
      CGAL::Non_manifold_tag,
      CGAL_SURFACE_MESHER_VERBOSITY
      >::type Surface_mesher_non_manifold; 

    if(mw->manifoldCheckBox->isChecked()) {
      // meshing surface
      std::cerr << "manifold criteria is activated.\n";
//       make_surface_mesh(c2t3, surface, oracle, criteria,
// 			CGAL::Manifold_tag(), 0);
      Surface_mesher_manifold manifold_mesher(c2t3, surface, oracle, criteria);
      manifold_mesher.refine_mesh();
    }
    else {
//       m_view_surface = true;
      Surface_mesher_non_manifold non_manifold_mesher(c2t3, surface, oracle, criteria);
#if 0
      int nb_steps = 0;
//       direct_draw = true;
      non_manifold_mesher.init();
      while(!non_manifold_mesher.is_algorithm_done()) {
	CGAL::Null_mesh_visitor null_visitor;
	non_manifold_mesher.one_step(null_visitor);
	if(++nb_steps % 1000 == 0) {
	  CGAL::Timer timer;
	  std::cerr << "(process events...";
	  timer.start();
	  list_draw_marching_cube_is_valid = false;
	  lists_draw_surface_is_valid = false;
	  for(Tr::Finite_cells_iterator 
		cit = del.finite_cells_begin(),
		end = del.finite_cells_end();
	      cit != end; ++cit)
	  {
	    cit->info() = classify(surface(cit->circumcenter()));
	  }
// 	  emit changed();
	  qApp->processEvents();
	  timer.stop();
	  std::cerr << timer.time() << " secondes)\n";
	}
      }
#else
      non_manifold_mesher.refine_mesh();
#endif
    }
    sm_timer.stop();
    not_busy();
    direct_draw = false;

    for(Tr::Finite_cells_iterator 
	  cit = del.finite_cells_begin(),
	  end = del.finite_cells_end();
	cit != end; ++cit)
    {
      cit->info() = classify(surface(cit->circumcenter()));
    }
    // get output surface
    for(C2t3::Facet_iterator
          fit = c2t3.facets_begin(), end = c2t3.facets_end();
        fit != end; ++fit)
    {
      const Tr::Cell_handle& cell = fit->first;
      const int index = fit->second;

      // here "left" means nothing
      const Point left_circumcenter = cell->circumcenter();
      const Point right_circumcenter = cell->neighbor(index)->circumcenter();

      const Triangle_3 t = 
        Triangle_3(cell->vertex(del.vertex_triple_index(index, 0))->point(),
                   cell->vertex(del.vertex_triple_index(index, 1))->point(),
                   cell->vertex(del.vertex_triple_index(index, 2))->point());
      const Vector u = t[1] - t[0];
      const Vector v = t[2] - t[0];
      Vector n = CGAL::cross_product(u,v);
      n = n / std::sqrt(n*n);
      if(mw->labellizedRadioButton->isChecked()) 
      {
	m_surface.push_back(Facet(t,
				  n,
				  values_list->search((std::max)(surface(left_circumcenter), 
								 surface(right_circumcenter)))));
      }
      else {
	m_surface.push_back(Facet(t,n,cell->vertex(del.vertex_triple_index(index, 0))->point().element_index()));
      }
    }

    // invalidate the display list
    lists_draw_surface_is_valid = false;
    sm_total_time = total_time.elapsed();
  }

  CGAL::Bbox_3 bbox(0,0,0,0,0,0);
  for(std::vector<Facet>::const_iterator
        it = m_surface.begin(), end = m_surface.end();
      it != end; ++it)
  {
    bbox = bbox + it->get<0>().bbox();
  }

  // toggle visualization
  m_view_mc = false;
  m_view_surface = true;
  emit changed();
  if(!m_surface.empty())
  {
    emit new_bounding_box(bbox.xmin(),
                          bbox.ymin(),
                          bbox.zmin(),
                          bbox.xmax(),
                          bbox.ymax(),
                          bbox.zmax());
  }
  status_message(tr("Surface meshing done. %1 facets in %2s (CPU time), total time is %3s.")
                 .arg(m_surface.size())
                 .arg(sm_timer.time())
                 .arg(sm_total_time/1000.));
  save_image_settings(fileinfo.absoluteFilePath());
}

void Volume::draw()
{
  float	ambient[]  =   { 0.25f,
                         0.20725f,
                         0.20725f,
                         0.922f };
  float	diffuse[]  =   { 1.0f,
                         0.829f,
                         0.829f,
                         0.922f };

  float	specular[]  = {  0.296648f,
                         0.296648f,
                         0.296648f,
                         0.522f };

  float	emission[]  = {  0.3f,
                         0.3f,
                         0.3f,
                         1.0f };
  float shininess[] = {  11.264f };

  // apply
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission);

  ::glEnable(GL_LINE_SMOOTH);

  if(two_sides)
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  else
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE); // default

  // draw surface mesh
  if(m_view_surface)
  {
    ::glEnable(GL_LIGHTING);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    ::glColor3f(0.2f, 0.2f, 1.f);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(3.0f,-3.0f);
    gl_draw_surface();

    if(draw_triangles_edges)
    {
      ::glDisable(GL_LIGHTING);
      ::glLineWidth(1.);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      ::glColor3ub(0,0,0);
      ::glDisable(GL_POLYGON_OFFSET_FILL);
      gl_draw_surface();
    }
  }

#ifdef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
  // draw MC surface mesh
  if(m_view_mc)
  {
    ::glEnable(GL_LIGHTING);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    ::glColor3f(0.2f, 0.2f, 1.f);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(3.0f,-3.0f);
    gl_draw_surface_mc();

    if(draw_triangles_edges)
    {
      ::glDisable(GL_LIGHTING);
      ::glLineWidth(1.);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      ::glColor3ub(0,0,0);
      ::glDisable(GL_POLYGON_OFFSET_FILL);
      gl_draw_surface_mc();
    }
  }
#endif // CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE

  if(show_bbox) {
    ::glDisable(GL_LIGHTING);
    m_image.gl_draw_bbox(3.0f,0,0,0);
  }

  if(!m_view_mc && m_draw_triangulation)
  {
    // draw the triangualtion
    mw->viewer->qglColor(m_triangulation_color);
    ::glLineWidth(1.0);
    ::glBegin(GL_LINES);
    for(Tr::Finite_edges_iterator 
          eit = del.finite_edges_begin(), 
          end = del.finite_edges_end();
        eit != end; ++eit) 
    {
      const Point p1 = eit->first->vertex(eit->second)->point();
      const Point p2 = eit->first->vertex(eit->third)->point();
      ::glVertex3d(p1.x(),p1.y(),p1.z());
      ::glVertex3d(p2.x(),p2.y(),p2.z());
    }
    ::glEnd();
  }
}

void Volume::set_radius_bound(double d)
{ 
  m_sm_radius = FT(d);
  changed_parameters();
}

void Volume::set_distance_bound(double d)
{ 
  m_sm_distance = FT(d);
  changed_parameters();
}

#ifdef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
void Volume::gl_draw_surface_mc()
{
  if(use_gouraud)
  {
    gl_draw_marchingcube();
    return;
  }
  
  if(lists_draw_surface_mc_is_valid)
  {
    for(int i = 0, nbs = values_list->numberOfValues(); i < nbs; ++i )
    {
      if(values_list->enabled(i))
      {
        mw->viewer->qglColor(values_list->color(i));
        ::glCallList(lists_draw_surface_mc[i]);
      }
    }
  }
  else
  {
    lists_draw_surface_mc.resize(values_list->numberOfValues(), 0);
    for(int i = 0, nbs = values_list->numberOfValues(); i < nbs; ++i )
    {
      if(!lists_draw_surface_mc[i])
      {
        lists_draw_surface_mc[i] = ::glGenLists(1);
      }

      std::cerr << boost::format("(Re-)Generating list #%1% for marching cube surface #%2%"
                                 " in gl_draw_surface(), ()\n")
        % lists_draw_surface_mc[i]
        % i;

      mw->viewer->qglColor(values_list->color(i));

      if(lists_draw_surface_mc[i])             // If
        ::glNewList(lists_draw_surface_mc[i],  // lists_draw_surface[i]==0 then something
                    values_list->enabled(i) // got wrong in the list generation.
                    ? GL_COMPILE_AND_EXECUTE
                    : GL_COMPILE);


      gl_draw_surface(m_surface_mc.begin(),
                      m_surface_mc.end(),
                      values_list->item(i));
        
      if(lists_draw_surface_mc[i]) // If lists_draw_surface[i]==0 then
      {                            // something got wrong in the list
        ::glEndList();             // generation.
      }
    }
    lists_draw_surface_mc_is_valid = (::glGetError() == GL_NO_ERROR);
  }
}
#endif // CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE

void Volume::gl_draw_surface()
{
//   if(mw->labellizedRadioButton->isChecked()) {
//     mw->viewer->qglColor(::Qt::blue);
//     gl_draw_surface(m_surface.begin(),
// 		    m_surface.end(),
// 		    0);
//   }
//   else
  if(direct_draw) {
    ::glBegin(GL_TRIANGLES);
    unsigned int counter = 0;
    for(Tr::Finite_cells_iterator
	  cit = del.finite_cells_begin(), end = del.finite_cells_end();
	cit != end; ++cit)
    {
      for(int facet_index = 0; facet_index < 4; ++facet_index)
      {
	const Tr::Cell_handle& facet_cell = cit;
	if(c2t3.face_status(facet_cell, facet_index) == C2t3::NOT_IN_COMPLEX) {
	  continue;
	}
	const Point& a = facet_cell->vertex(del.vertex_triple_index(facet_index, 0))->point();
	const Point& b = facet_cell->vertex(del.vertex_triple_index(facet_index, 1))->point();
	const Point& c = facet_cell->vertex(del.vertex_triple_index(facet_index, 2))->point();
	Vector n = CGAL::cross_product(b-a,c-a);
	n = n / std::sqrt(n*n); // unit normal
	if(m_inverse_normals) {
	  ::glNormal3d(-n.x(),-n.y(),-n.z());
	} else {
	  ::glNormal3d(n.x(),n.y(),n.z());
	}
	mw->viewer->qglColor(values_list->color(values_list->search(facet_cell->info())));
	::glVertex3d(a.x(),a.y(),a.z());
	::glVertex3d(b.x(),b.y(),b.z());
	::glVertex3d(c.x(),c.y(),c.z());
	++counter;
      }
      ::glEnd();
    }
    return;
  }
  if(!direct_draw && lists_draw_surface_is_valid)
  {
    for(int i = 0, nbs = values_list->numberOfValues(); i < nbs; ++i )
    {
      if(values_list->enabled(i))
      {
        mw->viewer->qglColor(values_list->color(i));
        ::glCallList(lists_draw_surface[i]);
      }
    }
  }
  else
  {
    lists_draw_surface.resize(values_list->numberOfValues(), 0);
    for(int i = 0, nbs = values_list->numberOfValues(); i < nbs; ++i )
    {
      if(!lists_draw_surface[i])
      {
        lists_draw_surface[i] = ::glGenLists(1);
      }

      std::cerr << boost::format("(Re-)Generating list #%1% for surface #%2%"
                                 " in gl_draw_surface(), ()\n")
        % lists_draw_surface[i]
        % i;
        
      mw->viewer->qglColor(values_list->color(i));

      if(!direct_draw && lists_draw_surface[i]) // If
        ::glNewList(lists_draw_surface[i],      // lists_draw_surface[i]==0
                    values_list->enabled(i)     // then something got wrong
                    ? GL_COMPILE_AND_EXECUTE    // in the list generation.
                    : GL_COMPILE);

      if(!mw->labellizedRadioButton->isChecked()) 
      {
	gl_draw_surface(m_surface.begin(),
			m_surface.end(),
			values_list->item(i));
      }
      else 
      {
	const unsigned char volume_index = values_list->value(i);

	::glBegin(GL_TRIANGLES);
	unsigned int counter = 0;
	for(C2t3::Facet_iterator
	      fit = c2t3.facets_begin(), end = c2t3.facets_end();
	    fit != end; ++fit)
	{
	  Tr::Cell_handle facet_cell = fit->first;
	  int facet_index = fit->second;
	  Tr::Cell_handle opposite_cell = facet_cell->neighbor(facet_index);
	  int opposite_index = opposite_cell->index(facet_cell);

	  if( facet_cell->info() != volume_index ) {
	    if( opposite_cell->info() == volume_index ) {
	      std::swap(facet_cell, opposite_cell);
	      std::swap(facet_index, opposite_index);
	    }
	    else 
	      continue; // go to next facet
	  }
	  const Point& a = opposite_cell->vertex(del.vertex_triple_index(opposite_index, 0))->point();
	  const Point& b = opposite_cell->vertex(del.vertex_triple_index(opposite_index, 1))->point();
	  const Point& c = opposite_cell->vertex(del.vertex_triple_index(opposite_index, 2))->point();
	  Vector n = CGAL::cross_product(b-a,c-a);
	  n = n / std::sqrt(n*n); // unit normal
	  if(m_inverse_normals) {
	    ::glNormal3d(-n.x(),-n.y(),-n.z());
	  } else {
	    ::glNormal3d(n.x(),n.y(),n.z());
	  }
	  ::glVertex3d(a.x(),a.y(),a.z());
	  ::glVertex3d(b.x(),b.y(),b.z());
	  ::glVertex3d(c.x(),c.y(),c.z());
	  ++counter;
	}
	::glEnd();
	std::cerr << boost::format("(c2t3) number of facets: %1%\n")
	  % counter;
      }

      if(!direct_draw && lists_draw_surface[i])
	                        // If lists_draw_surface[i]==0 then
      {                         // something got wrong in the list
        ::glEndList();          // generation.
      }
    }
    lists_draw_surface_is_valid = (::glGetError() == GL_NO_ERROR);
  }
}

template <typename Iterator>
void Volume::gl_draw_surface(Iterator begin, Iterator end, const QTreeWidgetItem* i)
{
  ::glBegin(GL_TRIANGLES);
  unsigned int counter = 0;
  for(Iterator it = begin; it != end; ++it)
  {
    const Facet& f = *it;

    if(f.get<2>() != i) continue;

    const Vector& n = f.get<1>();

    if(m_inverse_normals)
      ::glNormal3d(-n.x(),-n.y(),-n.z());
    else
      ::glNormal3d(n.x(),n.y(),n.z());

    const Triangle_3& t = f.get<0>();
    const Point& a = t[0];
    const Point& b = t[1];
    const Point& c = t[2];

    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
    ::glVertex3d(c.x(),c.y(),c.z());
    ++counter;
  }
  ::glEnd();
  std::cerr << boost::format("number of facets: %1%\n")
    % counter;
}

void Volume::changed_parameters()
{
  m_surface.clear();
#ifdef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
  m_surface_mc.clear();
#endif
  list_draw_marching_cube_is_valid = false;
  lists_draw_surface_is_valid = false;
  c2t3.clear();
  del.clear();
  m_view_mc = m_view_surface = false;
  emit changed();
}

#ifdef CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE
void Volume::gl_draw_one_marching_cube_vertex(int i)
{
  if(!m_inverse_normals)
    glArrayElement(i);
  else
  {
    const Vertex* const vertex = mc.vert(i);
    ::glNormal3d(-vertex->nx, -vertex->ny, -vertex->nz);
    ::glVertex3d(vertex->x, vertex->y, vertex->z);
  }
}

void Volume::gl_draw_marchingcube()
{
  if(list_draw_marching_cube_is_valid)
    ::glCallList(list_draw_marching_cube);
  else
  {
    if(!list_draw_marching_cube)
      list_draw_marching_cube = ::glGenLists(1);
    std::cerr << boost::format("(Re-)Generating list #%1% for"
                               " gl_draw_marchingcube()\n")
      % list_draw_marching_cube;

    if(list_draw_marching_cube)          // If list_draw_marching_cube==0 then
    ::glNewList(list_draw_marching_cube, // something got wrong in the list
                GL_COMPILE_AND_EXECUTE); // generation.

    ::glVertexPointer(3, GL_DOUBLE, sizeof(Vertex), mc.vertices());
    ::glNormalPointer(GL_DOUBLE, sizeof(Vertex), &(mc.vertices()->nx));
    ::glEnableClientState(GL_VERTEX_ARRAY);

    // because of that conditionnal, the display list has to be
    // reconstructed each time m_inverse_normals is toggled.
    if(!m_inverse_normals)
      ::glEnableClientState(GL_NORMAL_ARRAY);

    for(int i = 0, nbs = values_list->numberOfValues(); i < nbs; ++i)
    {
      const int begin = i == 0 ? 0 : nbs_of_mc_triangles[i-1];
      const int end = nbs_of_mc_triangles[i];
      mw->viewer->qglColor(values_list->color(i));
      ::glBegin(GL_TRIANGLES);
      for(int i = begin; i < end; ++i)
      {
        const MC_Triangle* const trig = mc.trig(i);
        gl_draw_one_marching_cube_vertex(trig->v1);
        gl_draw_one_marching_cube_vertex(trig->v2);
        gl_draw_one_marching_cube_vertex(trig->v3);
      }
      ::glEnd();
    }
    if(list_draw_marching_cube > 0) // If list_draw_marching_cube==0 then
    {                               // something got wrong in the list
      ::glEndList();                // generation.
      list_draw_marching_cube_is_valid = (::glGetError() == GL_NO_ERROR);
    }
    if(!list_draw_marching_cube_is_valid)
      std::cerr << boost::format("OpenGL error: %1%\n") 
        % ::gluErrorString(::glGetError());
  }
}
#endif // CGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE

void Volume::save_image_settings(QString filename)
{
  QSettings settings;
  settings.beginGroup(QUrl::toPercentEncoding(filename));
  settings.setValue("labellized", mw->labellizedRadioButton->isChecked());
  settings.endGroup();
}

void Volume::load_image_settings(QString filename)
{
  QSettings settings;
  settings.beginGroup(QUrl::toPercentEncoding(filename));
  mw->labellizedRadioButton->setChecked(settings.value("labellized").toBool());
  settings.endGroup();
}

void Volume::labellizedToogled(bool toggled)
{
  if(toggled) {
    values_list->setHeaderTitle(tr("Label"));
  }
  else {
    values_list->setHeaderTitle(tr("Iso-Value"));
  }
}

#include "volume.moc"
