// PoissonDoc.cpp : implementation of the CPoissonDoc class
//

// This demo + Gyroviz
#include "stdafx.h"
#include "Poisson.h"
#include "DialogOptions.h"
#include "PoissonDoc.h"
#include "compute_normal.h"
#include "read_pwc_point_set.h"
#include "read_g23_point_set.h"
#include "remove_outliers_wrt_camera_cone_angle.h"
#include "normal_orientation_wrt_cameras.h"

// CGAL
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#define CGAL_C2T3_USE_FILE_WRITER_OFF
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// This package
#include <CGAL/IO/read_off_point_set.h>
#include <CGAL/IO/write_off_point_set.h>
#include <CGAL/IO/read_xyz_point_set.h>
#include <CGAL/IO/write_xyz_point_set.h>
#include <CGAL/IO/output_surface_facets_to_triangle_soup.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/merge_simplify_point_set.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/radial_orient_normals.h>
#include <CGAL/Peak_memory_sizer.h>
#include <CGAL/surface_reconstruction_points_assertions.h>

// STL
#include <deque>
#include <vector>
#include <iterator>
#include <fstream>
#include <cassert>
#include <math.h>
#ifndef M_PI
  #define M_PI       3.14159265358979323846
#endif

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CPoissonDoc

IMPLEMENT_DYNCREATE(CPoissonDoc, CDocument)

BEGIN_MESSAGE_MAP(CPoissonDoc, CDocument)
    ON_COMMAND(ID_FILE_SAVE_SURFACE, OnFileSaveSurface)
    ON_UPDATE_COMMAND_UI(ID_FILE_SAVE_SURFACE, OnUpdateFileSaveSurface)
    ON_COMMAND(ID_FILE_SAVE_AS, OnFileSaveAs)
    ON_UPDATE_COMMAND_UI(ID_FILE_SAVE_AS, OnUpdateFileSaveAs)
    ON_COMMAND(ID_EDIT_OPTIONS, OnEditOptions)
    ON_COMMAND(ID_EDIT_DELETE, OnEditDelete)
    ON_UPDATE_COMMAND_UI(ID_EDIT_DELETE, OnUpdateEditDelete)
    ON_COMMAND(ID_EDIT_RESET_SELECTION, OnEditResetSelection)
    ON_UPDATE_COMMAND_UI(ID_EDIT_RESET_SELECTION, OnUpdateEditResetSelection)
    ON_COMMAND(ID_ALGORITHMS_ESTIMATENORMALSBYPCA, OnAlgorithmsEstimateNormalsByPCA)
    ON_COMMAND(ID_ALGORITHMS_ESTIMATENORMALBYJETFITTING, OnAlgorithmsEstimateNormalsByJetFitting)
    ON_COMMAND(ID_ALGORITHMS_ORIENTNORMALSCAMERAS, OnAlgorithmsOrientNormalsWrtCameras)
    ON_COMMAND(ID_ALGORITHMS_ORIENTNORMALSMST, OnAlgorithmsOrientNormalsWithMST)
	ON_COMMAND(ID_ALGORITHMS_SMOOTHUSINGJETFITTING, OnAlgorithmsSmoothUsingJetFitting)
    ON_COMMAND(ID_MODE_POINT_SET, OnModePointSet)
    ON_UPDATE_COMMAND_UI(ID_MODE_POINT_SET, OnUpdateModePointSet)
    ON_COMMAND(ID_MODE_POISSON, OnModePoisson)
    ON_UPDATE_COMMAND_UI(ID_MODE_POISSON, OnUpdateModePoisson)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_SMOOTHUSINGJETFITTING, OnUpdateAlgorithmsSmoothUsingJetFitting)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_ESTIMATENORMALSBYPCA, OnUpdateAlgorithmsEstimateNormalsByPCA)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_ESTIMATENORMALBYJETFITTING, OnUpdateAlgorithmsEstimateNormalByJetFitting)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_ORIENTNORMALSMST, OnUpdateAlgorithmsOrientNormalsWithMST)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_ORIENTNORMALSCAMERAS, OnUpdateAlgorithmsOrientNormalsWrtCameras)
    ON_COMMAND(ID_ALGORITHMS_OUTLIER_REMOVAL_WRT_CAMERAS_CONE_ANGLE, OnAlgorithmsOutlierRemovalWrtCamerasConeAngle)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_OUTLIER_REMOVAL_WRT_CAMERAS_CONE_ANGLE, OnUpdateAlgorithmsOutlierRemovalWrtCamerasConeAngle)
    ON_COMMAND(ID_ALGORITHMS_OUTLIER_REMOVAL, OnOutlierRemoval)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_OUTLIER_REMOVAL, OnUpdateOutlierRemoval)
    ON_COMMAND(ID_ANALYSIS_AVERAGE_SPACING, OnAnalysisAverageSpacing)
    ON_UPDATE_COMMAND_UI(ID_ANALYSIS_AVERAGE_SPACING, OnUpdateAnalysisAverageSpacing)
	ON_COMMAND(ID_RECONSTRUCTION_ONE_STEP_NORMALIZED, OnOneStepPoissonReconstructionWithNormalizedDivergence)
    ON_COMMAND(ID_RECONSTRUCTION_APSS_RECONSTRUCTION, OnReconstructionApssReconstruction)
    ON_UPDATE_COMMAND_UI(ID_RECONSTRUCTION_APSS_RECONSTRUCTION, OnUpdateReconstructionApssReconstruction)
    ON_COMMAND(ID_MODE_APSS, OnModeAPSS)
    ON_UPDATE_COMMAND_UI(ID_MODE_APSS, OnUpdateModeAPSS)
    ON_COMMAND(ID_POINT_CLOUD_SIMPLIFICATION_BY_CLUSTERING, OnPointCloudSimplificationByClustering)
    ON_UPDATE_COMMAND_UI(ID_POINT_CLOUD_SIMPLIFICATION_BY_CLUSTERING, OnUpdatePointCloudSimplificationByClustering)
    ON_COMMAND(ID_POINT_CLOUD_SIMPLIFICATION_RANDOM, OnPointCloudSimplificationRandom)
    ON_UPDATE_COMMAND_UI(ID_POINT_CLOUD_SIMPLIFICATION_RANDOM, OnUpdatePointCloudSimplificationRandom)
    ON_COMMAND(ID_RADIAL_NORMAL_ORIENTATION, OnRadialNormalOrientation)
    ON_UPDATE_COMMAND_UI(ID_RADIAL_NORMAL_ORIENTATION, OnUpdateRadialNormalOrientation)
    ON_COMMAND(ID_FLIP_NORMALS, OnFlipNormals)
    ON_UPDATE_COMMAND_UI(ID_FLIP_NORMALS, OnUpdateFlipNormals)
    ON_UPDATE_COMMAND_UI(ID_RECONSTRUCTION_ONE_STEP_NORMALIZED, OnUpdateOneStepPoissonReconstructionWithNormalizedDivergence)
END_MESSAGE_MAP()


// CPoissonDoc construction/destruction

CPoissonDoc::CPoissonDoc()
: m_poisson_function(NULL),
  m_apss_function(NULL),
  m_surface_mesher_c2t3(m_surface_mesher_dt)
{
  m_edit_mode = NO_EDIT_MODE; // No points yet

  // Poisson options
  m_sm_angle_poisson = 20.0; // Min triangle angle (degrees). 20 = fast, 30 guaranties convergence.
  m_sm_radius_poisson = 0.1; // Max triangle radius w.r.t. point set radius. 0.1 is fine.
  m_sm_distance_poisson = 0.002; // Approximation error w.r.t. p.s.r. For Poisson: 0.01 = fast, 0.002 = smooth.

  // APSS options
  m_sm_angle_apss = 20.0; // Min triangle angle (degrees). 20 = fast, 30 guaranties convergence.
  m_sm_radius_apss = 0.1; // Max triangle radius w.r.t. point set radius. 0.1 is fine.
  m_sm_distance_apss = 0.003; // Approximation error w.r.t. p.s.r. (APSS). 0.015 = fast, 0.003 = smooth.
                              // Note: 1.5 * Poisson's distance gives roughly the same number of triangles.
  m_nb_neighbors_apss = 24; // #neighbors to compute APPS sphere fitting. 12 = fast, 24 = robust (GG).

  // Average Spacing options
  m_nb_neighbors_avg_spacing = 7; // K-nearest neighbors = 1 ring (average spacing)

  // Smoothing options
  m_nb_neighbors_smooth_jet_fitting = 0.1 /* % */; // K-nearest neighbors (smooth points by Jet Fitting)

  // Normals Computing options
  m_nb_neighbors_pca_normals = 0.15 /* % */; // K-nearest neighbors (estimate normals by PCA)
  m_nb_neighbors_jet_fitting_normals = 0.1 /* % */; // K-nearest neighbors (estimate normals by Jet Fitting)
  m_nb_neighbors_mst = 18; // K-nearest neighbors = 3 rings (orient normals by MST)

  // Outlier Removal options
  m_min_cameras_cone_angle = 0.5 /* degrees */; // min angle of camera's cone
  m_threshold_percent_avg_knn_sq_dst = 5.0 /* % */; // percentage of outliers to remove
  m_nb_neighbors_remove_outliers = 0.05 /* % */; // K-nearest neighbors (outlier removal)

  // Point Set Simplification options
  m_clustering_step = 0.004; // Grid's step for simplification by clustering
  m_random_simplification_percentage = 50.0 /* % */; // percentage of random points to remove
}

CPoissonDoc::~CPoissonDoc()
{
  // Clean up current mode
  CloseMode();
}

// CPoissonDoc diagnostics

#ifdef _DEBUG
void CPoissonDoc::AssertValid() const
{
  CDocument::AssertValid();
}

void CPoissonDoc::Dump(CDumpContext& dc) const
{
  CDocument::Dump(dc);
}
#endif //_DEBUG


// File >> Open implementation
BOOL CPoissonDoc::OnOpenDocument(LPCTSTR lpszPathName)
{
  CGAL::Timer task_timer; task_timer.start();

  if (!CDocument::OnOpenDocument(lpszPathName))
    return FALSE;

  status_message("Load point set %s...",lpszPathName);

  // get extension
  CString file = lpszPathName;
  CString extension = lpszPathName;
  extension = extension.Right(4);
  extension.MakeLower();

  // set current path
  CString path = lpszPathName;
  path = path.Left(path.ReverseFind('\\'));
  SetCurrentDirectory(path);

  // if .off extension
  if(extension.CompareNoCase(".off") == 0)
  {
    // Is this OFF file a mesh or a point cloud?
    std::ifstream header_stream(lpszPathName);
    CGAL::File_scanner_OFF header(header_stream, true /* verbose */);
    if(!header_stream || header.size_of_vertices() == 0)
    {
      prompt_message("Unable to read file");
      return FALSE;
    }
    bool is_mesh = (header.size_of_facets() > 0);
    header_stream.close();

    // Read OFF file as a mesh and compute normals from connectivity
    if (is_mesh)
    {
      // Read the mesh file in a polyhedron
      std::ifstream stream(lpszPathName);
      typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
      Polyhedron input_mesh;
      CGAL::scan_OFF(stream, input_mesh, true /* verbose */);
      if(!stream || !input_mesh.is_valid() || input_mesh.empty())
      {
        prompt_message("Unable to read file");
        return FALSE;
      }

      // Convert Polyhedron vertices to point set.
      // Compute vertices' normals from connectivity.
      Polyhedron::Vertex_const_iterator v;
      for (v = input_mesh.vertices_begin(); v != input_mesh.vertices_end(); v++)
      {
        const Point& p = v->point();
        Vector n = compute_vertex_normal<Polyhedron::Vertex,Kernel>(*v);
        m_points.push_back(Point_with_normal(p,n));
      }
    }
    else // Read OFF file as a point cloud
    {
      std::ifstream stream(lpszPathName);
      if( ! stream || 
          ! CGAL::read_off_point_set(stream,
                                     std::back_inserter(m_points)) )
      {
        prompt_message("Unable to read file");
        return FALSE;
      }
    }
  }
  // if .pwn or .xyz extension
  else if(extension.CompareNoCase(".pwn") == 0 ||
          extension.CompareNoCase(".xyz") == 0)
  {
    std::ifstream stream(lpszPathName);
    if( ! stream || 
        ! CGAL::read_xyz_point_set(stream,
                                   std::back_inserter(m_points)) )
    {
      prompt_message("Unable to read file");
      return FALSE;
    }
  }
  // if Gyroviz .pwc extension
  else if (extension.CompareNoCase(".pwc") == 0)
  {
    std::deque<Point> cameras; // temporary container of cameras to read
    if( ! read_pwc_point_set(lpszPathName,
                             std::back_inserter(m_points),
                             std::back_inserter(cameras)) )
    {
      prompt_message("Unable to read file");
      return FALSE;
    }
  }
  // if Gyroviz .g23 extension
  else if (extension.CompareNoCase(".g23") == 0)
  {
    std::string movie_file_name;
    std::map<int,Point> cameras; // (indexed) cameras
    if( ! read_g23_point_set<Point>(lpszPathName,
                                    std::back_inserter(m_points),
                                    &cameras,
                                    &movie_file_name) )
    {
      prompt_message("Unable to read file");
      return FALSE;
    }
  }
  // if unknown extension
  else
  {
    prompt_message("File format not supported");
    return FALSE;
  }

  // Save original normals for visual comparison
  for (int i=0; i<m_points.size(); i++)
    m_points[i].original_normal() = m_points[i].normal();

  m_points.invalidate_bounds();
  m_edit_mode = POINT_SET;

  status_message("Load point set...done (%.2lf s)", task_timer.time());
  update_status();
  UpdateAllViews(NULL);
  return TRUE;
}

// Save input point set as...  callback
void CPoissonDoc::OnFileSaveAs()
{
  // file filters
  CString szFilter;
  szFilter  = "Points With Normals Files (*.pwn)|*.pwn|";
  szFilter += "XYZ Files (*.xyz)|*.xyz|";
  szFilter += "OFF Files (*.off)|*.off|";
  szFilter += "All Files (*.*)|*.*||";

  // create the Save As dialog
  CFileDialog dlgExport(false, "pwn", NULL,
                        OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, AfxGetMainWnd());

  // dialog title
  dlgExport.m_ofn.lpstrTitle = "Save reconstructed surface to file";

  // show the dialog
  if (dlgExport.DoModal() == IDOK)
  {
    // get extension
    CString file = dlgExport.m_ofn.lpstrFile;
    CString extension = dlgExport.m_ofn.lpstrFile;
    extension = extension.Right(4);
    extension.MakeLower();

    // set current path
    CString path = dlgExport.m_ofn.lpstrFile;
    path = path.Left(path.ReverseFind('\\'));
    SetCurrentDirectory(path);

    // Save normals?
    assert(m_points.begin() != m_points.end());
    bool points_have_normals = (m_points.begin()->normal() != CGAL::NULL_VECTOR);
    bool save_normals = points_have_normals || 
                        extension.CompareNoCase(".pwn"); // pwn means "point with normal"

    // if .pwn or .xyz extension
    if(extension.CompareNoCase(".pwn") == 0 ||
       extension.CompareNoCase(".xyz") == 0)
    {
      std::ofstream stream(dlgExport.m_ofn.lpstrFile);
      if( ! stream || 
          ! CGAL::write_xyz_point_set(stream,
                                      m_points.begin(), m_points.end(),
                                      save_normals) )
      {
        prompt_message("Unable to save file");
        return;
      }
    }
    // if .off extension
    else if(extension.CompareNoCase(".off") == 0)
    {
      std::ofstream stream(dlgExport.m_ofn.lpstrFile);
      if( ! stream || 
          ! CGAL::write_off_point_set(stream,
                                      m_points.begin(), m_points.end()) )
      {
        prompt_message("Unable to save file");
        return;
      }
    }
    else
    {
      prompt_message("File format not supported");
      return;
    }
  }
}

// Save m_points[] only if it is the form visible on screen
void CPoissonDoc::OnUpdateFileSaveAs(CCmdUI *pCmdUI)
{
  assert(m_points.begin() != m_points.end());
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

// Save reconstructed surface as...  callback
void CPoissonDoc::OnFileSaveSurface()
{
  // file filters
  static char szFilter[] = "OFF Files (*.off)|*.off; *.off|All Files (*.*)|*.*||";

  // create the Save As dialog
  CFileDialog dlgExport(false, "off", NULL,
                        OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, AfxGetMainWnd());

  // dialog title
  dlgExport.m_ofn.lpstrTitle = "Save reconstructed surface to file";

  // show the dialog
  if (dlgExport.DoModal() == IDOK)
  {
    // get extension
    CString file = dlgExport.m_ofn.lpstrFile;
    CString extension = dlgExport.m_ofn.lpstrFile;
    extension = extension.Right(4);
    extension.MakeLower();

    // set current path
    CString path = dlgExport.m_ofn.lpstrFile;
    path = path.Left(path.ReverseFind('\\'));
    SetCurrentDirectory(path);

    // if .off extension
    if(extension.CompareNoCase(".off") == 0)
    {
      std::ofstream out((char *)dlgExport.m_ofn.lpstrFile);
      if( !out )
      {
        prompt_message("Unable to save file");
        return;
      }

      CGAL::output_surface_facets_to_off(out, m_surface_mesher_c2t3);
    }
    else
    {
      prompt_message("File format not supported");
      return;
    }
  }
}

// Enable "Save reconstructed surface as..." if surface is computed
void CPoissonDoc::OnUpdateFileSaveSurface(CCmdUI *pCmdUI)
{
  pCmdUI->Enable((m_edit_mode == POISSON || m_edit_mode == APSS)
                 && m_surface_mesher_dt.number_of_vertices() > 0);
}

// Update the number of vertices and tetrahedra in the status bar
// and write them to cerr.
void CPoissonDoc::update_status()
{
  CStatusBar* pStatus =
    (CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(AFX_IDW_STATUS_BAR);
  ASSERT(pStatus != NULL);

  if (m_edit_mode == POINT_SET || m_edit_mode == APSS)
  {
    CString points;
    points.Format("%d points",m_points.size());
    CString selected_points;
    selected_points.Format("%d selected",m_points.nb_selected_points());

    // write message to cerr
    std::cerr << "=> " << points << " (" << selected_points << "), "
                       << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated, "
                       << "largest free block=" << (CGAL::Peak_memory_sizer().largest_free_block()>>20) << " Mb, "
                       << "#blocks over 100 Mb=" << CGAL::Peak_memory_sizer().count_free_memory_blocks(100*1048576)
                       << std::endl;

    // Update status bar
    pStatus->SetPaneText(1,points);
    pStatus->SetPaneText(2,selected_points);
    pStatus->UpdateWindow();
  }
  else if (m_edit_mode == POISSON)
  {
    assert(m_poisson_function != NULL);

    CString vertices;
    vertices.Format("%d vertices",m_poisson_function->triangulation().number_of_vertices());
    CString tets;
    tets.Format("%d tets",m_poisson_function->triangulation().number_of_cells());

    // write message to cerr
    std::cerr << "=> " << vertices << ", " << tets << ", "
                       << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated, "
                       << "largest free block=" << (CGAL::Peak_memory_sizer().largest_free_block()>>20) << " Mb, "
                       << "#blocks over 100 Mb=" << CGAL::Peak_memory_sizer().count_free_memory_blocks(100*1048576)
                       << std::endl;

    // Update status bar
    pStatus->SetPaneText(1,vertices);
    pStatus->SetPaneText(2,tets);
    pStatus->UpdateWindow();
  }
}

// Write user message in status bar and cerr
void CPoissonDoc::status_message(char* fmt,...)
{
  // format message in 'buffer'
  char buffer[256];
  va_list argptr;
  va_start(argptr,fmt);
  vsprintf(buffer,fmt,argptr);
  va_end(argptr);

  // write message to cerr
  std::cerr << buffer << std::endl;

  // write message in status bar
  CStatusBar* pStatus =
    (CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(AFX_IDW_STATUS_BAR);
  ASSERT(pStatus != NULL);
  pStatus->SetPaneText(0,buffer);
  pStatus->UpdateWindow();
}

// Write user message in message box and cerr
void CPoissonDoc::prompt_message(char* fmt,...)
{
  // format message in 'buffer'
  char buffer[256];
  va_list argptr;
  va_start(argptr,fmt);
  vsprintf(buffer,fmt,argptr);
  va_end(argptr);

  // write message to cerr
  std::cerr << buffer << std::endl;

  // write message in message box
  AfxMessageBox(buffer);
}

// Display Options dialog
void CPoissonDoc::OnEditOptions()
{
  CDialogOptions dlg;
  dlg.m_sm_angle_poisson = m_sm_angle_poisson;
  dlg.m_sm_radius_poisson = m_sm_radius_poisson;
  dlg.m_sm_distance_poisson = m_sm_distance_poisson;
  dlg.m_sm_angle_apss = m_sm_angle_apss;
  dlg.m_sm_radius_apss = m_sm_radius_apss;
  dlg.m_sm_distance_apss = m_sm_distance_apss;
  dlg.m_nb_neighbors_avg_spacing = m_nb_neighbors_avg_spacing;
  dlg.m_nb_neighbors_remove_outliers = m_nb_neighbors_remove_outliers;
  dlg.m_nb_neighbors_smooth_jet_fitting = m_nb_neighbors_smooth_jet_fitting;
  dlg.m_nb_neighbors_pca_normals = m_nb_neighbors_pca_normals;
  dlg.m_nb_neighbors_jet_fitting_normals = m_nb_neighbors_jet_fitting_normals;
  dlg.m_nb_neighbors_mst = m_nb_neighbors_mst;
  dlg.m_nb_neighbors_apss = m_nb_neighbors_apss;
  dlg.m_min_cameras_cone_angle = m_min_cameras_cone_angle;
  dlg.m_threshold_percent_avg_knn_sq_dst = m_threshold_percent_avg_knn_sq_dst;
  dlg.m_clustering_step = m_clustering_step;
  dlg.m_random_simplification_percentage = m_random_simplification_percentage;

  if(dlg.DoModal() == IDOK)
  {
    m_sm_angle_poisson = dlg.m_sm_angle_poisson;
    m_sm_radius_poisson = dlg.m_sm_radius_poisson;
    m_sm_distance_poisson = dlg.m_sm_distance_poisson;
    m_sm_angle_apss = dlg.m_sm_angle_apss;
    m_sm_radius_apss = dlg.m_sm_radius_apss;
    m_sm_distance_apss = dlg.m_sm_distance_apss;
    m_nb_neighbors_avg_spacing = dlg.m_nb_neighbors_avg_spacing;
    m_nb_neighbors_remove_outliers = dlg.m_nb_neighbors_remove_outliers;
    m_nb_neighbors_smooth_jet_fitting = dlg.m_nb_neighbors_smooth_jet_fitting;
    m_nb_neighbors_pca_normals = dlg.m_nb_neighbors_pca_normals;
    m_nb_neighbors_jet_fitting_normals = dlg.m_nb_neighbors_jet_fitting_normals;
    m_nb_neighbors_mst = dlg.m_nb_neighbors_mst;
    m_nb_neighbors_apss = dlg.m_nb_neighbors_apss;
    m_min_cameras_cone_angle = dlg.m_min_cameras_cone_angle;
    m_threshold_percent_avg_knn_sq_dst = dlg.m_threshold_percent_avg_knn_sq_dst;
    m_clustering_step = dlg.m_clustering_step;
    m_random_simplification_percentage = dlg.m_random_simplification_percentage;

    UpdateAllViews(NULL);
    EndWaitCursor();
  }
}

// Check the accuracy of normals direction estimation.
// If original normals are available, compare with them and select normals with large deviation.
// @return true on success.
bool CPoissonDoc::verify_normal_direction()
{
  bool success = true;

  m_points.select(m_points.begin(), m_points.end(), false);

  assert(m_points.begin() != m_points.end());
  bool points_have_original_normals = (m_points.begin()->original_normal() != CGAL::NULL_VECTOR);
  if (points_have_original_normals)
  {
    std::cerr << "Compare with original normals:" << std::endl;

    double min_normal_deviation = DBL_MAX; // deviation / original normal
    double max_normal_deviation = DBL_MIN;
    double avg_normal_deviation = 0;
    int invalid_normals = 0; // #normals with large deviation
    for (Point_set::iterator p = m_points.begin(); p != m_points.end(); p++)
    {
      // Compute normal deviation.
      // Orient each normal like the original one (but keep the "non oriented" flag).
      Vector v1 = p->original_normal(); // input normal
      double norm1 = std::sqrt( v1*v1 );
      assert(norm1 != 0.0);
      Vector v2 = p->normal(); // computed normal
      double norm2 = std::sqrt( v2*v2 );
      assert(norm2 != 0.0);
      double cos_normal_deviation = (v1*v2)/(norm1*norm2);
      if (cos_normal_deviation < 0)
      {
        cos_normal_deviation = -cos_normal_deviation;
        assert( ! p->normal().is_oriented() );
        p->normal() = Normal(-v2, false /*non oriented*/);
      }
      double normal_deviation = std::acos(cos_normal_deviation);
      assert(normal_deviation >= 0 && normal_deviation <= M_PI/2.);

      // statistics about normals deviation
      min_normal_deviation = (std::min)(min_normal_deviation, normal_deviation);
      max_normal_deviation = (std::max)(max_normal_deviation, normal_deviation);
      avg_normal_deviation += normal_deviation;

      // count and select normal if large deviation
      bool valid = (normal_deviation <= M_PI/3.); // valid if deviation <= 60 degrees
      if ( ! valid )
      {
        m_points.select(&*p);
        invalid_normals++;
      }
    }
    avg_normal_deviation /= double(m_points.size());

    std::cerr << "  Min normal deviation=" << min_normal_deviation*180.0/M_PI << " degrees\n";
    std::cerr << "  Max normal deviation=" << max_normal_deviation*180.0/M_PI << " degrees\n";
    std::cerr << "  Avg normal deviation=" << avg_normal_deviation*180.0/M_PI << " degrees\n";
    if (invalid_normals > 0)
    {
      std::cerr << "  Error: " << invalid_normals << " normals have a deviation > 60 degrees\n";
      success = false;
    }
  }

  return success;
}

// Compute normals direction by Principal Component Analysis
void CPoissonDoc::OnAlgorithmsEstimateNormalsByPCA()
{
  BeginWaitCursor();
  CGAL::Timer task_timer; task_timer.start();

  // percentage -> number of neighbors
  int nb_neighbors = int(double(m_points.size()) * m_nb_neighbors_pca_normals / 100.0);
  if (nb_neighbors < 7)
    nb_neighbors = 7;
  if ((unsigned int)nb_neighbors > m_points.size()-1)
    nb_neighbors = m_points.size()-1;

  status_message("Estimate Normals Direction by PCA (k=%.2lf%%=%d)...",
                 m_nb_neighbors_pca_normals, nb_neighbors);

  CGAL::pca_estimate_normals(m_points.begin(), m_points.end(),
                             m_points.normals_begin(),
                             nb_neighbors);

  status_message("Estimate Normals Direction by PCA...done (%.2lf s)", task_timer.time());

  // Check the accuracy of normals direction estimation.
  // If original normals are available, compare with them.
  verify_normal_direction();

  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsEstimateNormalsByPCA(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

// Compute normals direction by Jet Fitting
void CPoissonDoc::OnAlgorithmsEstimateNormalsByJetFitting()
{
  BeginWaitCursor();
  CGAL::Timer task_timer; task_timer.start();

  // percentage -> number of neighbors
  int nb_neighbors = int(double(m_points.size()) * m_nb_neighbors_jet_fitting_normals / 100.0);
  if (nb_neighbors < 7)
    nb_neighbors = 7;
  if ((unsigned int)nb_neighbors > m_points.size()-1)
    nb_neighbors = m_points.size()-1;

  status_message("Estimate Normals Direction by Jet Fitting (k=%.2lf%%=%d)...",
                 m_nb_neighbors_jet_fitting_normals, nb_neighbors);

  CGAL::jet_estimate_normals(m_points.begin(), m_points.end(),
                             m_points.normals_begin(),
                             nb_neighbors);

  status_message("Estimate Normals Direction by Jet Fitting...done (%.2lf s)", task_timer.time());

  // Check the accuracy of normals direction estimation.
  // If original normals are available, compare with them.
  verify_normal_direction();

  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsEstimateNormalByJetFitting(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

// Check the accuracy of normal orientation.
// Count and select non-oriented normals.
// If original normals are available, compare with them and select flipped normals.
bool CPoissonDoc::verify_normal_orientation()
{
  bool success = true;

  m_points.select(m_points.begin(), m_points.end(), false);

  // Count and select non-oriented normals
  int unoriented_normals = 0;
  for (Point_set::iterator p = m_points.begin(); p != m_points.end(); p++)
  {
    // Check unit vector
    Vector v = p->normal();
    double norm = std::sqrt( v*v );
    assert(norm > 0.99 || norm < 1.01);

    // Check orientation
    if ( ! p->normal().is_oriented() )
    {
      m_points.select(&*p);
      unoriented_normals++;
    }
  }
  if (unoriented_normals > 0)
  {
    std::cerr << "Error: " << unoriented_normals << " normals are unoriented\n";
    success = false;
  }

  // If original normals are available, compare with them and select flipped normals
  assert(m_points.begin() != m_points.end());
  bool points_have_original_normals = (m_points.begin()->original_normal() != CGAL::NULL_VECTOR);
  if (points_have_original_normals)
  {
    std::cerr << "Compare with original normals:" << std::endl;

    int flipped_normals = 0; // #normals with wrong orientation
    for (Point_set::iterator p = m_points.begin(); p != m_points.end(); p++)
    {
      Vector v1 = p->original_normal(); // input normal
      double norm1 = std::sqrt( v1*v1 );
      assert(norm1 != 0.0);
      Vector v2 = p->normal(); // computed normal
      double norm2 = std::sqrt( v2*v2 );
      assert(norm2 != 0.0);
      double cos_normal_deviation = (v1*v2)/(norm1*norm2);
      if (cos_normal_deviation < 0 // if flipped
       && p->normal().is_oriented()) // unoriented normals are already reported
      {
        m_points.select(&*p);
        flipped_normals++;
      }
    }

    if (flipped_normals == 0)
      std::cerr << "  ok\n";
    else
      std::cerr << "  Error: " << flipped_normals << " normal(s) are flipped\n";
  }

  return success;
}

// Orient the normals using a minimum spanning tree.
void CPoissonDoc::OnAlgorithmsOrientNormalsWithMST()
{
  BeginWaitCursor();
  status_message("Orient Normals with a Minimum Spanning Tree (k=%d)...", m_nb_neighbors_mst);
  CGAL::Timer task_timer; task_timer.start();

  // Mark all normals as unoriented
  for (Point_set::iterator p = m_points.begin(); p != m_points.end(); p++)
    p->normal().set_orientation(false);

  CGAL::mst_orient_normals(m_points.begin(), m_points.end(),
                           get(boost::vertex_index, m_points),
                           get(CGAL::vertex_point, m_points),
                           get(boost::vertex_normal, m_points),
                           m_nb_neighbors_mst);

  status_message("Orient Normals with a Minimum Spanning Tree...done (%.2lf s)", task_timer.time());

  // Check the accuracy of normal orientation.
  // If original normals are available, compare with them.
  verify_normal_orientation();

  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsOrientNormalsWithMST(CCmdUI *pCmdUI)
{
  assert(m_points.begin() != m_points.end());
  bool points_have_normals = (m_points.begin()->normal() != CGAL::NULL_VECTOR);
  pCmdUI->Enable(m_edit_mode == POINT_SET && points_have_normals);
}

// Specific to Gyroviz: orient the normals w.r.t. the position of cameras
// that reconstructed the points by photogrammetry.
void CPoissonDoc::OnAlgorithmsOrientNormalsWrtCameras()
{
  BeginWaitCursor();
  status_message("Orient Normals wrt Cameras...");
  CGAL::Timer task_timer; task_timer.start();

  normal_orientation_wrt_cameras(m_points.begin(), m_points.end(),
                                 get(CGAL::vertex_point, m_points),
                                 get(boost::vertex_normal, m_points),
                                 get(boost::vertex_cameras, m_points));

  status_message("Orient Normals wrt Cameras...done (%.2lf s)", task_timer.time());

  // Check the accuracy of normal orientation.
  // If original normals are available, compare with them.
  verify_normal_orientation();

  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsOrientNormalsWrtCameras(CCmdUI *pCmdUI)
{
  assert(m_points.begin() != m_points.end());
  bool points_have_normals = (m_points.begin()->normal() != CGAL::NULL_VECTOR);
  bool points_have_cameras = (m_points.begin()->cameras_begin() != m_points.begin()->cameras_end());
  pCmdUI->Enable(m_edit_mode == POINT_SET && points_have_normals && points_have_cameras);
}

// Smooth point set using jet fitting + projection
void CPoissonDoc::OnAlgorithmsSmoothUsingJetFitting()
{
  BeginWaitCursor();
  CGAL::Timer task_timer; task_timer.start();

  // percentage -> number of neighbors
  int nb_neighbors = int(double(m_points.size()) * m_nb_neighbors_smooth_jet_fitting / 100.0);
  if (nb_neighbors < 7)
    nb_neighbors = 7;
  if ((unsigned int)nb_neighbors > m_points.size()-1)
    nb_neighbors = m_points.size()-1;

  status_message("Smooth Point Set (k=%.2lf%%=%d)...",
                 m_nb_neighbors_smooth_jet_fitting, nb_neighbors);

  // Smooth points in m_points[]
  CGAL::jet_smooth_point_set(m_points.begin(), m_points.end(),
                             nb_neighbors);
  m_points.invalidate_bounds();

  status_message("Smooth Point Set...done (%.2lf s)", task_timer.time());
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsSmoothUsingJetFitting(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

// Clean up current mode
void CPoissonDoc::CloseMode()
{
  // Nothing to do if m_edit_mode == POINT_SET

  // If m_edit_mode == POISSON
  delete m_poisson_function; m_poisson_function = NULL;

  // If m_edit_mode == APSS
  delete m_apss_function; m_apss_function = NULL;

  m_edit_mode = NO_EDIT_MODE;
}

// Edit >> Mode >> Point set callback
void CPoissonDoc::OnModePointSet()
{
  // No need to convert Poisson triangulation back to point set (yet)

  // Clean up previous mode
  CloseMode();

  m_edit_mode = POINT_SET;

  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateModePointSet(CCmdUI *pCmdUI)
{
  // Edit >> Mode >> Point set is always enabled
  pCmdUI->SetCheck(m_edit_mode == POINT_SET);
}

// "Edit >> Mode >> Poisson" is an alias to
// "Reconstruction >> Poisson Reconstruction w/ Normalized Divergence".
void CPoissonDoc::OnModePoisson()
{
  OnOneStepPoissonReconstructionWithNormalizedDivergence();
}

void CPoissonDoc::OnUpdateModePoisson(CCmdUI *pCmdUI)
{
  OnUpdateOneStepPoissonReconstructionWithNormalizedDivergence(pCmdUI);
  pCmdUI->SetCheck(m_edit_mode == POISSON);
}

// "Reconstruction >> Poisson Reconstruction w/ Normalized Divergence" callback: 
// - Create Poisson Delaunay Triangulation
// - Delaunay refinement
// - Solve Poisson Equation
// - Surface Meshing
void CPoissonDoc::OnOneStepPoissonReconstructionWithNormalizedDivergence()
{
  BeginWaitCursor();
  status_message("1-step Poisson reconstruction with normalized divergence...");
  CGAL::Timer total_timer; total_timer.start();

  // Clean up previous mode
  CloseMode();

  CGAL::Timer task_timer; task_timer.start();

  //***************************************
  // Compute implicit function
  //***************************************

  status_message("Create Poisson triangulation...");

  // Create implicit function and insert vertices.
  assert(m_poisson_function == NULL);
  m_poisson_function = new Poisson_reconstruction_function(m_points.begin(), m_points.end());

  // Print status
  status_message("Create Poisson triangulation...done (%.2lf s)", task_timer.time());
  task_timer.reset();

  status_message("Compute implicit function...");

  // Compute the Poisson indicator function f()
  // at each vertex of the triangulation.
  if ( ! m_poisson_function->compute_implicit_function() )
  {
    status_message("Error: cannot compute implicit function");
    return;
  }

  // Print status
  status_message("Compute implicit function...done (%.2lf s)", task_timer.time());
  task_timer.reset();

  //***************************************
  // Surface mesh generation
  //***************************************

  status_message("Surface meshing...");

  // Clear previous call
  m_surface_mesher_dt.clear();
  m_surface_mesher_c2t3.clear();
  m_surface.clear();

  // Get point inside the implicit surface
  Point inner_point = m_poisson_function->get_inner_point();
  FT inner_point_value = (*m_poisson_function)(inner_point);
  if(inner_point_value >= 0.0)
  {
    status_message("Error: unable to seed (%lf at inner_point)",inner_point_value);
    return;
  }

  // Get implicit function's radius
  Sphere bounding_sphere = m_poisson_function->bounding_sphere();
  FT size = sqrt(bounding_sphere.squared_radius());

  // defining the implicit surface = implicit function + bounding sphere centered at inner_point
  typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;
  Point sm_sphere_center = inner_point;
  FT    sm_sphere_radius = size + std::sqrt(CGAL::squared_distance(bounding_sphere.center(),inner_point));
  sm_sphere_radius *= 1.01; // <= the Surface Mesher fails if the sphere does not contain the surface
  Surface_3 surface(*m_poisson_function,
                    Sphere(sm_sphere_center,sm_sphere_radius*sm_sphere_radius));

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<STr> criteria(m_sm_angle_poisson,  // Min triangle angle (degrees)
                                                      m_sm_radius_poisson*size,  // Max triangle radius
                                                      m_sm_distance_poisson*size); // Approximation error

  // meshing surface
  CGAL::make_surface_mesh(m_surface_mesher_c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

  // Print status
  status_message("Surface meshing...done (%d output vertices, %.2lf s)",
                 m_surface_mesher_dt.number_of_vertices(), task_timer.time());
  task_timer.reset();
    
  // get output surface
  std::deque<Triangle> triangles;
  CGAL::output_surface_facets_to_triangle_soup(m_surface_mesher_c2t3, std::back_inserter(triangles));
  m_surface.insert(m_surface.end(), triangles.begin(), triangles.end());

  // Record new mode
  m_edit_mode = POISSON;

  status_message("1-step Poisson reconstruction with normalized divergence...done (%.2lf s)", total_timer.time());
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateOneStepPoissonReconstructionWithNormalizedDivergence(CCmdUI *pCmdUI)
{
  assert(m_points.begin() != m_points.end());
  bool points_have_normals = (m_points.begin()->normal() != CGAL::NULL_VECTOR);
  bool normals_are_oriented = m_points.begin()->normal().is_oriented();
  pCmdUI->Enable((m_edit_mode == POINT_SET || m_edit_mode == POISSON)
                 && points_have_normals /* && normals_are_oriented */);
}

// Remove vertices / cameras cone's angle is low
void CPoissonDoc::OnAlgorithmsOutlierRemovalWrtCamerasConeAngle()
{
  BeginWaitCursor();
  status_message("Remove outliers / cameras cone's angle...");
  CGAL::Timer task_timer; task_timer.start();

  // Select points to delete
  m_points.select(m_points.begin(), m_points.end(), false);
  Point_set::iterator first_iterator_to_remove =
    remove_outliers_wrt_camera_cone_angle(
            m_points.begin(), m_points.end(),
            m_min_cameras_cone_angle*M_PI/180.0);
  m_points.select(first_iterator_to_remove, m_points.end(), true);

  status_message("Remove outliers / cameras cone's angle...done (%.2lf s)", task_timer.time());
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsOutlierRemovalWrtCamerasConeAngle(CCmdUI *pCmdUI)
{
  assert(m_points.begin() != m_points.end());
  bool points_have_cameras = (m_points.begin()->cameras_begin() != m_points.begin()->cameras_end());
  pCmdUI->Enable(m_edit_mode == POINT_SET && points_have_cameras);
}

// Remove outliers:
// - compute average squared distance to the K nearest neighbors,
// - remove threshold_percent worst points.
void CPoissonDoc::OnOutlierRemoval()
{
  BeginWaitCursor();
  CGAL::Timer task_timer; task_timer.start();

  // percentage -> number of neighbors
  int nb_neighbors = int(double(m_points.size()) * m_nb_neighbors_remove_outliers / 100.0);
  if (nb_neighbors < 7)
    nb_neighbors = 7;
  if ((unsigned int)nb_neighbors > m_points.size()-1)
    nb_neighbors = m_points.size()-1;

  status_message("Remove outliers wrt average squared distance to k nearest neighbors (remove %.2lf%%, k=%.2lf%%=%d)...",
                 m_threshold_percent_avg_knn_sq_dst, m_nb_neighbors_remove_outliers, nb_neighbors);

  // Select points to delete
  m_points.select(m_points.begin(), m_points.end(), false);
  Point_set::iterator first_iterator_to_remove =
    CGAL::remove_outliers(
            m_points.begin(), m_points.end(),
            nb_neighbors,
            m_threshold_percent_avg_knn_sq_dst);
  m_points.select(first_iterator_to_remove, m_points.end(), true);

  status_message("Remove outliers wrt average squared distance to k nearest neighbors...done (%.2lf s)", task_timer.time());
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateOutlierRemoval(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

void CPoissonDoc::OnAnalysisAverageSpacing()
{
  BeginWaitCursor();
  status_message("Compute average spacing to k nearest neighbors (k=%d)...", m_nb_neighbors_avg_spacing);
  CGAL::Timer task_timer; task_timer.start();

  double value = CGAL::compute_average_spacing(m_points.begin(),
                                               m_points.end(),
                                               m_nb_neighbors_avg_spacing);

  // write message in message box
  prompt_message("Average spacing: %lf", value);

  status_message("Compute average spacing to k nearest neighbors...done: %lf (%.2lf s)", value, task_timer.time());
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAnalysisAverageSpacing(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

// "Reconstruction >> APSS reconstruction" callback
void CPoissonDoc::OnReconstructionApssReconstruction()
{
    BeginWaitCursor();
    status_message("APSS reconstruction...");
    CGAL::Timer task_timer; task_timer.start();

    // Clean up previous mode
    CloseMode();

    // Clear previous call
    m_surface_mesher_dt.clear();
    m_surface_mesher_c2t3.clear();
    m_surface.clear();

    // Create implicit function
    m_apss_function = new APSS_reconstruction_function(m_points.begin(), m_points.end(),
                                                       m_nb_neighbors_apss);

    // Get point inside the implicit surface
    Point inner_point = m_apss_function->get_inner_point();
    FT inner_point_value = (*m_apss_function)(inner_point);
    if(inner_point_value >= 0.0)
    {
      status_message("Error: unable to seed (%lf at inner_point)",inner_point_value);
      return;
    }

    // Get implicit function's radius
    Sphere bounding_sphere = m_apss_function->bounding_sphere();
    FT size = sqrt(bounding_sphere.squared_radius());

    // defining the implicit surface = implicit function + bounding sphere centered at inner_point
    typedef CGAL::Implicit_surface_3<Kernel, APSS_reconstruction_function> Surface_3;
    Point sm_sphere_center = inner_point;
    FT    sm_sphere_radius = size + std::sqrt(CGAL::squared_distance(bounding_sphere.center(),inner_point));
    sm_sphere_radius *= 1.01; // <= the Surface Mesher fails if the sphere does not contain the surface
    Surface_3 surface(*m_apss_function,
                      Sphere(sm_sphere_center,sm_sphere_radius*sm_sphere_radius)); 

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(m_sm_angle_apss,  // Min triangle angle (degrees)
                                                        m_sm_radius_apss*size,  // Max triangle radius
                                                        m_sm_distance_apss*size); // Approximation error

    // meshing surface
    CGAL::make_surface_mesh(m_surface_mesher_c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

    // get output surface
    std::deque<Triangle> triangles;
    CGAL::output_surface_facets_to_triangle_soup(m_surface_mesher_c2t3, std::back_inserter(triangles));
    m_surface.insert(m_surface.end(), triangles.begin(), triangles.end());

    // Record new mode
    m_edit_mode = APSS;

    // Print status
    status_message("APSS reconstruction...done (%d vertices, %.2lf s)",
                   m_surface_mesher_dt.number_of_vertices(), task_timer.time());
    update_status();
    UpdateAllViews(NULL);
    EndWaitCursor();
}

// Enable "Reconstruction >> APSS reconstruction" if normals are computed and oriented.
void CPoissonDoc::OnUpdateReconstructionApssReconstruction(CCmdUI *pCmdUI)
{
  assert(m_points.begin() != m_points.end());
  bool points_have_normals = (m_points.begin()->normal() != CGAL::NULL_VECTOR);
  bool normals_are_oriented = m_points.begin()->normal().is_oriented();
  pCmdUI->Enable((m_edit_mode == POINT_SET || m_edit_mode == APSS)
                 && points_have_normals /* && normals_are_oriented */);
}

// "Edit >> Mode >> APSS" is an alias to "Reconstruction >> APSS reconstruction".
void CPoissonDoc::OnModeAPSS()
{
  OnReconstructionApssReconstruction();
}

// "Edit >> Mode >> APSS" is an alias to "Reconstruction >> APSS reconstruction".
void CPoissonDoc::OnUpdateModeAPSS(CCmdUI *pCmdUI)
{
  OnUpdateReconstructionApssReconstruction(pCmdUI);
  pCmdUI->SetCheck(m_edit_mode == APSS);
}


void CPoissonDoc::OnEditDelete()
{
    BeginWaitCursor();
    status_message("Delete selected points");

    m_points.delete_selection();

    update_status();
    UpdateAllViews(NULL);
    EndWaitCursor();
}

void CPoissonDoc::OnUpdateEditDelete(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

void CPoissonDoc::OnEditResetSelection()
{
    BeginWaitCursor();
    status_message("Reset selection");

    m_points.select(m_points.begin(), m_points.end(), false);

    update_status();
    UpdateAllViews(NULL);
    EndWaitCursor();
}

void CPoissonDoc::OnUpdateEditResetSelection(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

void CPoissonDoc::OnPointCloudSimplificationByClustering()
{
  BeginWaitCursor();
  status_message("Point cloud simplification by clustering...");
  CGAL::Timer task_timer; task_timer.start();

  // Get point set's radius
  Sphere bounding_sphere = m_points.bounding_sphere();
  FT size = sqrt(bounding_sphere.squared_radius());

  // Select points to delete
  m_points.select(m_points.begin(), m_points.end(), false);
  Point_set::iterator first_iterator_to_remove =
    CGAL::merge_simplify_point_set(
            m_points.begin(), m_points.end(),
            m_clustering_step*size);
  m_points.select(first_iterator_to_remove, m_points.end(), true);

  status_message("Point cloud simplification by clustering...done (%.2lf s)", task_timer.time());
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdatePointCloudSimplificationByClustering(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

void CPoissonDoc::OnPointCloudSimplificationRandom()
{
  BeginWaitCursor();
  status_message("Random point cloud simplification...");
  CGAL::Timer task_timer; task_timer.start();

  // Select points to delete
  m_points.select(m_points.begin(), m_points.end(), false);
  Point_set::iterator first_iterator_to_remove =
    CGAL::random_simplify_point_set(
            m_points.begin(), m_points.end(),
            m_random_simplification_percentage);
  m_points.select(first_iterator_to_remove, m_points.end(), true);

  status_message("Random point cloud simplification...done (%.2lf s)", task_timer.time());
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdatePointCloudSimplificationRandom(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

void CPoissonDoc::OnRadialNormalOrientation()
{
  BeginWaitCursor();
  status_message("Radial Normal Orientation...");
  CGAL::Timer task_timer; task_timer.start();

  CGAL::radial_orient_normals(m_points.begin(), m_points.end(),
                              get(boost::vertex_index, m_points),
                              get(CGAL::vertex_point, m_points),
                              get(boost::vertex_normal, m_points));

  status_message("Radial Normal Orientation...done (%.2lf s)", task_timer.time());

  // Check the accuracy of normal orientation.
  // If original normals are available, compare with them.
  verify_normal_orientation();

  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateRadialNormalOrientation(CCmdUI *pCmdUI)
{
  assert(m_points.begin() != m_points.end());
  bool points_have_normals = (m_points.begin()->normal() != CGAL::NULL_VECTOR);
  pCmdUI->Enable(m_edit_mode == POINT_SET && points_have_normals);
}

void CPoissonDoc::OnFlipNormals()
{
  BeginWaitCursor();
  status_message("Flip Normals...");
  CGAL::Timer task_timer; task_timer.start();

  // Flip normals
  for (int i=0; i<m_points.size(); i++)
    m_points[i].normal() = -m_points[i].normal();

  status_message("Flip Normals...done (%.2lf s)", task_timer.time());

  // Check the accuracy of normal orientation.
  // If original normals are available, compare with them.
  verify_normal_orientation();

  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateFlipNormals(CCmdUI *pCmdUI)
{
  assert(m_points.begin() != m_points.end());
  bool points_have_normals = (m_points.begin()->normal() != CGAL::NULL_VECTOR);
  pCmdUI->Enable(m_edit_mode == POINT_SET && points_have_normals);
}
