// PoissonDoc.cpp : implementation of the CPoissonDoc class
//

// This demo
#include "stdafx.h"
#include "Poisson.h"
#include "DialogOptions.h"
#include "PoissonDoc.h"
#include "enriched_polyhedron.h"
#include "surface_reconstruction_read_pwc.h"

// CGAL
#define CGAL_C2T3_USE_POLYHEDRON
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/assertions.h>

// This package
#include <CGAL/IO/surface_reconstruction_read_off_point_cloud.h>
#include <CGAL/IO/surface_reconstruction_write_off_point_cloud.h>
#include <CGAL/IO/surface_reconstruction_read_pnb.h>
#include <CGAL/IO/surface_reconstruction_read_pwn.h>
#include <CGAL/IO/surface_reconstruction_write_pwn.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>
#include <CGAL/IO/surface_reconstruction_write_xyz.h>
#include <CGAL/IO/surface_reconstruction_output.h>
#include <CGAL/remove_outliers_wrt_avg_knn_sq_distance_3.h>
#include <CGAL/estimate_normals_jet_fitting_3.h>
#include <CGAL/smooth_jet_fitting_3.h>
#include <CGAL/estimate_normals_pca_3.h>
#include <CGAL/average_spacing_3.h>

// STL
#include <iostream>
#include <fstream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CPoissonDoc

IMPLEMENT_DYNCREATE(CPoissonDoc, CDocument)

BEGIN_MESSAGE_MAP(CPoissonDoc, CDocument)
    ON_COMMAND(ID_RECONSTRUCTION_DELAUNAYREFINEMENT, OnReconstructionDelaunayRefinement)
    ON_COMMAND(ID_RECONSTRUCTION_POISSON, OnReconstructionPoisson)
    ON_COMMAND(ID_ALGORITHMS_REFINEINSHELL, OnAlgorithmsRefineInShell)
    ON_COMMAND(ID_RECONSTRUCTION_SURFACEMESHING, OnReconstructionSurfaceMeshing)
    ON_COMMAND(ID_EDIT_OPTIONS, OnEditOptions)
    ON_UPDATE_COMMAND_UI(ID_RECONSTRUCTION_POISSON, OnUpdateReconstructionPoisson)
    ON_UPDATE_COMMAND_UI(ID_RECONSTRUCTION_SURFACEMESHING, OnUpdateReconstructionSurfaceMeshing)
    ON_COMMAND(ID_ALGORITHMS_MARCHINGTETCONTOURING, OnAlgorithmsMarchingTetContouring)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_MARCHINGTETCONTOURING, OnUpdateAlgorithmsMarchingTetContouring)
    ON_COMMAND(ID_FILE_SAVE_SURFACE, OnFileSaveSurface)
    ON_UPDATE_COMMAND_UI(ID_FILE_SAVE_SURFACE, OnUpdateFileSaveSurface)
    ON_COMMAND(ID_FILE_SAVE_AS, OnFileSaveAs)
    ON_UPDATE_COMMAND_UI(ID_FILE_SAVE_AS, OnUpdateFileSaveAs)
    ON_COMMAND(ID_ALGORITHMS_EXTRAPOLATENORMALS, OnAlgorithmsExtrapolatenormals)
    ON_COMMAND(ID_ALGORITHMS_POISSONSTATISTICS, OnAlgorithmsPoissonStatistics)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_POISSONSTATISTICS, OnUpdateAlgorithmsPoissonStatistics)
    ON_COMMAND(ID_ALGORITHMS_ESTIMATENORMALSBYPCA, OnAlgorithmsEstimateNormalsByPCA)
    ON_COMMAND(ID_ALGORITHMS_ESTIMATENORMALBYJETFITTING, OnAlgorithmsEstimateNormalsByJetFitting)
    ON_COMMAND(ID_ALGORITHMS_ORIENTNORMALSCAMERAS, OnAlgorithmsOrientNormalsWrtCameras)
    ON_COMMAND(ID_ALGORITHMS_ORIENTNORMALSMST, &CPoissonDoc::OnAlgorithmsOrientNormalsWithMST)
	ON_COMMAND(ID_ALGORITHMS_SMOOTHUSINGJETFITTING, &CPoissonDoc::OnAlgorithmsSmoothUsingJetFitting)
    ON_COMMAND(ID_MODE_POINT_SET, &CPoissonDoc::OnModePointSet)
    ON_UPDATE_COMMAND_UI(ID_MODE_POINT_SET, &CPoissonDoc::OnUpdateModePointSet)
    ON_COMMAND(ID_MODE_POISSON, &CPoissonDoc::OnModePoisson)
    ON_UPDATE_COMMAND_UI(ID_MODE_POISSON, &CPoissonDoc::OnUpdateModePoisson)
    ON_COMMAND(ID_CREATE_POISSON_TRIANGULATION, &CPoissonDoc::OnCreatePoissonTriangulation)
    ON_UPDATE_COMMAND_UI(ID_CREATE_POISSON_TRIANGULATION, &CPoissonDoc::OnUpdateCreatePoissonTriangulation)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_SMOOTHUSINGJETFITTING, &CPoissonDoc::OnUpdateAlgorithmsSmoothUsingJetFitting)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_ESTIMATENORMALSBYPCA, &CPoissonDoc::OnUpdateAlgorithmsEstimateNormalsByPCA)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_ESTIMATENORMALBYJETFITTING, &CPoissonDoc::OnUpdateAlgorithmsEstimateNormalByJetFitting)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_ORIENTNORMALSMST, &CPoissonDoc::OnUpdateAlgorithmsOrientNormalsWithMST)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_ORIENTNORMALSCAMERAS, &CPoissonDoc::OnUpdateAlgorithmsOrientNormalsWrtCameras)
    ON_UPDATE_COMMAND_UI(ID_RECONSTRUCTION_DELAUNAYREFINEMENT, &CPoissonDoc::OnUpdateReconstructionDelaunayRefinement)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_REFINEINSHELL, &CPoissonDoc::OnUpdateAlgorithmsRefineInShell)
    ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_EXTRAPOLATENORMALS, &CPoissonDoc::OnUpdateAlgorithmsExtrapolateNormals)
    ON_COMMAND(ID_PROCESSING_REMOVEOUTLIERS, &CPoissonDoc::OnRemoveOutliers)
    ON_UPDATE_COMMAND_UI(ID_PROCESSING_REMOVEOUTLIERS, &CPoissonDoc::OnUpdateRemoveOutliers)
    ON_COMMAND(ID_ANALYSIS_AVERAGE_SPACING, &CPoissonDoc::OnAnalysisAverageSpacing)
    ON_COMMAND(ID_ONE_STEP_POISSON_RECONSTRUCTION, &CPoissonDoc::OnOneStepPoissonReconstruction)
    ON_UPDATE_COMMAND_UI(ID_ONE_STEP_POISSON_RECONSTRUCTION, &CPoissonDoc::OnUpdateOneStepPoissonReconstruction)
    ON_UPDATE_COMMAND_UI(ID_ANALYSIS_AVERAGE_SPACING, &CPoissonDoc::OnUpdateAnalysisAverageSpacing)
END_MESSAGE_MAP()


// CPoissonDoc construction/destruction

CPoissonDoc::CPoissonDoc()
: m_poisson_function(m_poisson_dt), 
  m_surface_mesher_c2t3(m_surface_mesher_dt)
{
  m_edit_mode = NO_EDIT_MODE; // No points yet
  
  m_triangulation_refined = false; // Need to apply Delaunay refinement
  m_poisson_solved = false; // Need to solve Poisson equation

  // Surface mesher options
  m_sm_angle = 20.0; // LR: 30 is OK
  m_sm_radius = 0.1; // as suggested by LR (was 0.01)
  m_sm_distance = 0.001; // was 0.01

  // Delaunay refinement options
  m_dr_shell_size = 0.01;
  m_dr_sizing = 0.5 * m_dr_shell_size;
  m_dr_max_vertices = (unsigned int)5e6;

  // Surface mesher and marching tet common options
  m_contouring_value = 0.0; // 0 by default

  // Normal estimation options
  m_number_of_neighbours = 7; // by default

  // Outlier removal
  m_outlier_percentage = 10.0;
}

CPoissonDoc::~CPoissonDoc()
{
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
  double init = clock();

  if (!CDocument::OnOpenDocument(lpszPathName))
    return FALSE;

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
      AfxMessageBox("Unable to read file");
      return FALSE;
    }
    bool is_mesh = (header.size_of_facets() > 0);
    header_stream.close();

    // Read OFF file as a mesh and compute normals from connectivity
    if (is_mesh)
    {
      // read file in polyhedron 
      typedef Enriched_polyhedron<Kernel> Polyhedron;
      Polyhedron input_mesh;
      std::ifstream file_stream(lpszPathName);
      CGAL::scan_OFF(file_stream, input_mesh, true /* verbose */); 
      if(!file_stream || !input_mesh.is_valid() || input_mesh.empty())
      {
        AfxMessageBox("Unable to read file");
        return FALSE;
      }

      // Compute normals using mesh connectivity
      input_mesh.compute_normals();

      // Copy points to m_points
      Polyhedron::Vertex_iterator v;
      for(v = input_mesh.vertices_begin();
          v != input_mesh.vertices_end();
          v++)
      {
        const Point& p = v->point();
        const Vector& n = v->normal();
        m_points.push_back(Point_with_normal(p,n));
      }
    }
    else // Read OFF file as a point cloud
    {
      if( ! CGAL::surface_reconstruction_read_off_point_cloud(lpszPathName, 
                                                              std::back_inserter(m_points)) )
      {
        AfxMessageBox("Unable to read file");
        return FALSE;
      }
    }
  }
  // if .pwn extension
  else if(extension.CompareNoCase(".pwn") == 0)
  {
    if( ! CGAL::surface_reconstruction_read_pwn(lpszPathName, 
                                                std::back_inserter(m_points)) )
    {
      AfxMessageBox("Unable to read file");
      return FALSE;
    }
  }
  // if .xyz extension
  else if(extension.CompareNoCase(".xyz") == 0)
  {
    if( ! CGAL::surface_reconstruction_read_xyz(lpszPathName, 
                                                std::back_inserter(m_points)) )
    {
      AfxMessageBox("Unable to read file");
      return FALSE;
    }
    // 
  }
  // if .pnb extension
  else if(extension.CompareNoCase(".pnb") == 0)
  {
    if( ! CGAL::surface_reconstruction_read_pnb(lpszPathName, 
                                                std::back_inserter(m_points)) )
    {
      AfxMessageBox("Unable to read file");
      return FALSE;
    }
  }
  // if Gyroviz .pwc extension
  else if (extension.CompareNoCase(".pwc") == 0)
  {
    std::vector<Point> cameras; // temporary container of cameras to read
    if( ! surface_reconstruction_read_pwc(lpszPathName, 
                                          std::back_inserter(m_points),
                                          std::back_inserter(cameras)) )
    {
      AfxMessageBox("Unable to read file");
      return FALSE;
    }
  }
  // if unknown extension
  else
  {
    AfxMessageBox("File format not supported");
    return FALSE;
  }

  m_points.invalidate_bounding_box();
  m_edit_mode = POINT_SET;

  status_message("Load point set (%lf s)",duration(init));
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

    // if .pwn extension
    if(extension.CompareNoCase(".pwn") == 0)
    {
      if( ! CGAL::surface_reconstruction_write_pwn(dlgExport.m_ofn.lpstrFile,
                                                   m_points.begin(), m_points.end()) )
      {
        AfxMessageBox("Unable to save file");
        return;
      }
    }
    // if .xyz extension
    else if(extension.CompareNoCase(".xyz") == 0)
    {
      if( ! CGAL::surface_reconstruction_write_xyz(dlgExport.m_ofn.lpstrFile,
                                                   m_points.begin(), m_points.end()) )
      {
        AfxMessageBox("Unable to save file");
        return;
      }
    }
    // if .off extension
    else if(extension.CompareNoCase(".off") == 0)
    {
      if( ! CGAL::surface_reconstruction_write_off_point_cloud(dlgExport.m_ofn.lpstrFile,
                                                               m_points.begin(), m_points.end()) )
      {
        AfxMessageBox("Unable to save file");
        return;
      }
    }
    else 
    {
      AfxMessageBox("File format not supported");
      return;
    }
  }
}

// Save m_points[] only if it is the form visible on screen
void CPoissonDoc::OnUpdateFileSaveAs(CCmdUI *pCmdUI)
{
  CGAL_assertion(m_points.begin() != m_points.end());
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
        AfxMessageBox("Unable to save file");
        return;
      }
        
      CGAL::output_surface_facets_to_off(out, m_surface_mesher_c2t3);
    }
    else
    {
      AfxMessageBox("File format not supported");
      return;
    }
  }
}

// Enable "Save reconstructed surface as..." if surface is computed
void CPoissonDoc::OnUpdateFileSaveSurface(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POISSON && m_surface_mesher_dt.number_of_vertices() > 0);
}

// Update the number of vertices and tetrahedra in the status bar
void CPoissonDoc::update_status()
{   
  CWinApp *pApp = AfxGetApp();
  if(pApp->m_pMainWnd != NULL) 
  { 
    CStatusBar* pStatus = 
      (CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(AFX_IDW_STATUS_BAR);
    if (pStatus != NULL) 
    {
      // Update status bar
      if (m_edit_mode == POINT_SET)
      {
        CString points;
        points.Format("%d points",m_points.size());
        pStatus->SetPaneText(1,points);
        pStatus->UpdateWindow(); 
      }
      else if (m_edit_mode == POISSON)
      {
        CString vertices;
        vertices.Format("%d vertices",m_poisson_dt.number_of_vertices());
        CString tets;
        tets.Format("%d tets",m_poisson_dt.number_of_cells());
        pStatus->SetPaneText(1,vertices);
        pStatus->SetPaneText(2,tets);
        pStatus->UpdateWindow(); 
      }
    }
  }
}

// Set user message in status bar
void CPoissonDoc::status_message(char* fmt,...)
{   
  CWinApp *pApp = AfxGetApp();
  if(pApp->m_pMainWnd != NULL) 
  { 
    char buffer[256];
    CStatusBar* pStatus = 
      (CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(
      AFX_IDW_STATUS_BAR);
    
    // fill buffer
    va_list argptr;      
    va_start(argptr,fmt);
    vsprintf(buffer,fmt,argptr);
    va_end(argptr);
    
    if(pStatus != NULL) 
    {
      pStatus->SetPaneText(0,buffer);
      pStatus->UpdateWindow(); 
    }
  }
  return;
}

// Display Options dialog
void CPoissonDoc::OnEditOptions()
{
  CDialogOptions dlg;
  dlg.m_sm_angle = m_sm_angle;
  dlg.m_sm_radius = m_sm_radius;
  dlg.m_sm_distance = m_sm_distance;

  dlg.m_dr_sizing = m_dr_sizing;
  dlg.m_dr_shell_size = m_dr_shell_size;
  dlg.m_dr_max_vertices = m_dr_max_vertices;

  dlg.m_contouring_value = m_contouring_value;

  dlg.m_number_of_neighbours = m_number_of_neighbours;
  dlg.m_outlier_percentage = m_outlier_percentage;

  if(dlg.DoModal() == IDOK)
  {
    m_sm_angle = dlg.m_sm_angle;
    m_sm_radius = dlg.m_sm_radius;
    m_sm_distance = dlg.m_sm_distance;

    m_dr_sizing = dlg.m_dr_sizing;
    m_dr_shell_size = dlg.m_dr_shell_size;
    m_dr_max_vertices = dlg.m_dr_max_vertices;

    m_contouring_value = dlg.m_contouring_value;

    m_number_of_neighbours = dlg.m_number_of_neighbours;

    m_outlier_percentage = dlg.m_outlier_percentage;
    
    update_status();
    UpdateAllViews(NULL);
    EndWaitCursor();
  }
}

// Utility: compute elapsed time
double CPoissonDoc::duration(double time_init)
{
  return (clock() - time_init)/CLOCKS_PER_SEC;
}

// Compute normals direction by Principal Component Analysis
void CPoissonDoc::OnAlgorithmsEstimateNormalsByPCA()
{
  BeginWaitCursor();
  status_message("Estimate Normals Direction...");
  double init = clock();

  CGAL::estimate_normals_pca_3(m_points.begin(), m_points.end(), 
                               m_points.normals_begin(), 
                               m_number_of_neighbours);

  status_message("Estimate Normals Direction...done (%lf s)",duration(init));
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
  status_message("Estimate Normals Direction...");
  double init = clock();

  CGAL::estimate_normals_jet_fitting_3(m_points.begin(), m_points.end(), 
                                       m_points.normals_begin(), 
                                       m_number_of_neighbours);

  status_message("Estimate Normals Direction...done (%lf s)",duration(init));
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsEstimateNormalByJetFitting(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

// Orient the normals using a minimum spanning tree.
void CPoissonDoc::OnAlgorithmsOrientNormalsWithMST()
{
  BeginWaitCursor();
  status_message("Orient Normals with MST...");
  double init = clock();

  CGAL::orient_normals_minimum_spanning_tree_3(m_points.begin(), m_points.end(),
                                               get(boost::vertex_index, m_points),
                                               get(CGAL::vertex_point, m_points),
                                               get(boost::vertex_normal, m_points),
                                               m_number_of_neighbours);

  status_message("Orient Normals with MST...done (%lf s)",duration(init));
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsOrientNormalsWithMST(CCmdUI *pCmdUI)
{
  CGAL_assertion(m_points.begin() != m_points.end());
  bool points_have_normals = (m_points.begin()->normal().get_vector() != CGAL::NULL_VECTOR);
  pCmdUI->Enable(m_edit_mode == POINT_SET && points_have_normals);
}

// Specific to Gyroviz: orient the normals w.r.t. the position of cameras
// that reconstructed the points by photogrammetry.
void CPoissonDoc::OnAlgorithmsOrientNormalsWrtCameras()
{
  BeginWaitCursor();
  status_message("Orient Normals wrt Cameras...");
  double init = clock();

  CGAL::orient_normals_wrt_cameras_3(m_points.begin(), m_points.end(),
                                     get(CGAL::vertex_point, m_points),
                                     get(boost::vertex_normal, m_points),
                                     get(boost::vertex_cameras, m_points));

  status_message("Orient Normals wrt Cameras...done (%lf s)",duration(init));
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsOrientNormalsWrtCameras(CCmdUI *pCmdUI)
{
  CGAL_assertion(m_points.begin() != m_points.end());
  bool points_have_normals = (m_points.begin()->normal().get_vector() != CGAL::NULL_VECTOR);
  bool points_have_cameras = (m_points.begin()->cameras_begin() != m_points.begin()->cameras_end());
  pCmdUI->Enable(m_edit_mode == POINT_SET && points_have_normals && points_have_cameras);
}

// Reconstruction >> Poisson >> Create Poisson Triangulation callback:
// Create m_poisson_dt from m_points[] and enter in Poisson edit mode.
// m_points is not editable in this mode.
void CPoissonDoc::OnCreatePoissonTriangulation()
{
  BeginWaitCursor();
  status_message("Create Poisson Triangulation...");
  double init = clock();

  // Copy points to m_poisson_dt
  m_poisson_function.clear();
  m_poisson_function.insert(m_points.begin(), m_points.end());
  
  m_edit_mode = POISSON;
  m_triangulation_refined = false; // Need to apply Delaunay refinement
  m_poisson_solved = false; // Need to solve Poisson equation

  status_message("Create Poisson Triangulation...done (%lf s)",duration(init));
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateCreatePoissonTriangulation(CCmdUI *pCmdUI)
{
  CGAL_assertion(m_points.begin() != m_points.end());
  bool points_have_normals = (m_points.begin()->normal().get_vector() != CGAL::NULL_VECTOR);
  bool normals_are_oriented = m_points.begin()->normal().is_normal_oriented();
  pCmdUI->Enable(m_edit_mode == POINT_SET && points_have_normals && normals_are_oriented);

  pCmdUI->SetCheck(m_edit_mode == POISSON);
}

// Uniform Delaunay refinement
void CPoissonDoc::OnReconstructionDelaunayRefinement()
{
  BeginWaitCursor();

  status_message("Delaunay refinement...");
  const double quality = 2.5;
  const unsigned int max_vertices = (unsigned int)1e7; // max 10M vertices
  const double enlarge_ratio = 1.5;
  double init = clock();
  unsigned int nb_vertices_added = m_poisson_function.delaunay_refinement(quality,max_vertices,enlarge_ratio,50000);
  m_triangulation_refined = true;

  status_message("Delaunay refinement...done (%lf s, %d vertices inserted)",duration(init),nb_vertices_added);
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateReconstructionDelaunayRefinement(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POISSON);
}

// Delaunay refinement in a surface's shell
void CPoissonDoc::OnAlgorithmsRefineInShell()
{
  BeginWaitCursor();

  status_message("Delaunay refinement...");
  const double quality = 2.5;
  const unsigned int max_vertices = (unsigned int)1e7; // max 10M vertices
  const double enlarge_ratio = 1.5;
  double init = clock();
  unsigned int nb_vertices_added = m_poisson_function.delaunay_refinement_shell(m_dr_shell_size,m_dr_sizing,m_dr_max_vertices);
  m_triangulation_refined = true;

  status_message("Delaunay refinement...done (%lf s, %d vertices inserted)",duration(init),nb_vertices_added);
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsRefineInShell(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POISSON);
}

// Extrapolate the normals field:
// compute null normals by averaging neighbour normals.
void CPoissonDoc::OnAlgorithmsExtrapolatenormals()
{
  BeginWaitCursor();
    m_poisson_function.extrapolate_normals();
  EndWaitCursor();
  UpdateAllViews(NULL);
}

void CPoissonDoc::OnUpdateAlgorithmsExtrapolateNormals(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POISSON && m_triangulation_refined);
}

// Solve Poisson equation callback
void CPoissonDoc::OnReconstructionPoisson()
{
  BeginWaitCursor();
  status_message("Solve Poisson equation...");
  double init = clock();

  // Solve Poisson equation such that:
  // - m_poisson_function() = 0 on the input points,
  // - m_poisson_function() < 0 inside the surface.
  double duration_assembly, duration_factorization, duration_solve;
  m_poisson_solved = m_poisson_function.solve_poisson(&duration_assembly,
                                                      &duration_factorization,
                                                      &duration_solve);
  m_poisson_function.set_contouring_value(m_poisson_function.median_value_at_input_vertices());
  m_contouring_value = 0.0;

  double total_duration = duration(init);
  if (!m_poisson_solved)
      status_message("Poisson reconstruction...solver failed");
  else
      status_message("Solve Poisson equation...done (%lf s)", total_duration);
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

// Enable "Solve Poisson equation" if Delaunay refinement is applied
void CPoissonDoc::OnUpdateReconstructionPoisson(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POISSON && m_triangulation_refined);
}

// Surface Meshing callback
void CPoissonDoc::OnReconstructionSurfaceMeshing()
{
    typedef CGAL::Implicit_surface_3<Kernel, Poisson_implicit_function&> Surface_3;

    // Clear previous call
    m_surface_mesher_dt.clear();
    m_surface_mesher_c2t3.clear();
    
    // Apply contouring value defined in Options dialog
    m_poisson_function.set_contouring_value(m_contouring_value);

    // Get inner point
    Point sink = m_poisson_function.sink();
    FT f_sink = m_poisson_function(sink);
    if(f_sink >= 0.0)
    {
      status_message("Unable to seed (%lf at sink)",f_sink);
      return;
    }

    BeginWaitCursor();
    double init = clock();

    // Get implicit surface's size
    Sphere bounding_sphere = m_poisson_function.bounding_sphere();
    FT size = sqrt(bounding_sphere.squared_radius());

    // defining the surface
    Surface_3 surface(m_poisson_function,                    
                      Sphere(sink,4*size*size)); // bounding sphere

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(m_sm_angle,  // lower bound of facets angles (degrees)
                                                        m_sm_radius*size,  // upper bound of Delaunay balls radii
                                                        m_sm_distance*size); // upper bound of distance to surface

    // meshing surface
    status_message("Surface meshing...");
    make_surface_mesh(m_surface_mesher_c2t3, surface, criteria, CGAL::Non_manifold_tag());

    status_message("Surface meshing...done (%d vertices, %lf s)",
                   m_surface_mesher_dt.number_of_vertices(),duration(init));

    // get output surface
    std::list<Triangle> triangles;
    CGAL::output_surface_facets<C2t3,Triangle>(triangles,m_surface_mesher_c2t3);
    m_poisson_dt.set_surface(triangles);
    
    // Reset contouring value 
    m_poisson_function.set_contouring_value(0.0);

    UpdateAllViews(NULL);
    EndWaitCursor();
}

// Enable "Surface Meshing" if Poisson equation is solved
void CPoissonDoc::OnUpdateReconstructionSurfaceMeshing(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POISSON && m_poisson_solved);
}

// Marching Tet Contouring callback
void CPoissonDoc::OnAlgorithmsMarchingTetContouring()
{
  int nb = m_poisson_dt.marching_tet(m_contouring_value);
  status_message("Marching tet contouring...done (%d triangles)",nb);
  UpdateAllViews(NULL);
}

// Enable "Marching Tet Contouring" if Poisson equation is solved
void CPoissonDoc::OnUpdateAlgorithmsMarchingTetContouring(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POISSON && m_poisson_solved);
}

void CPoissonDoc::OnAlgorithmsPoissonStatistics()
{
  BeginWaitCursor();
    char buffer[1000];
    sprintf(buffer, "Poisson implicit function:\n- Median value at input vertices = %lf\n- Average value at input vertices = %lf\n- Min value at input vertices = %lf\n- Max value at input vertices = %lf\n- Median value at convex hull = %lf\n- Average value at convex hull = %lf\n- Min value = %lf", 
                    m_poisson_function.median_value_at_input_vertices(), 
                    m_poisson_function.average_value_at_input_vertices(), 
                    m_poisson_function.min_value_at_input_vertices(), 
                    m_poisson_function.max_value_at_input_vertices(), 
                    m_poisson_function.median_value_at_convex_hull(), 
                    m_poisson_function.average_value_at_convex_hull(), 
                    m_poisson_function(m_poisson_function.sink()));
    AfxMessageBox(buffer);
  EndWaitCursor();
}

// Enable "Poisson Statistics" if Poisson equation is solved
void CPoissonDoc::OnUpdateAlgorithmsPoissonStatistics(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POISSON && m_poisson_solved);
}

// Smooth point set using jet fitting + projection
void CPoissonDoc::OnAlgorithmsSmoothUsingJetFitting()
{
  BeginWaitCursor();
  status_message("Smooth Point Set...");
  double init = clock();

  // Smooth points in m_points[]
  std::vector<Point> output; // intermediate buffer
  CGAL::smooth_jet_fitting_3(m_points.begin(),m_points.end(),std::back_inserter(output),m_number_of_neighbours);
  for (int i=0; i < m_points.size(); i++)
      m_points[i].set_position(output[i]);
  m_points.invalidate_bounding_box();

  status_message("Smooth Point Set...done (%lf s)",duration(init));
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAlgorithmsSmoothUsingJetFitting(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

// Edit >> Mode >> Point set callback
void CPoissonDoc::OnModePointSet()
{
  // No need to convert Poisson triangulation back to point set (yet)
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

// Edit >> Mode >> Poisson is an alias to 
// Reconstruction >> Poisson >> Create Poisson Triangulation.
void CPoissonDoc::OnModePoisson()
{
  OnCreatePoissonTriangulation();
}

void CPoissonDoc::OnUpdateModePoisson(CCmdUI *pCmdUI)
{
  OnUpdateCreatePoissonTriangulation(pCmdUI);
  pCmdUI->SetCheck(m_edit_mode == POISSON);
}

// Reconstruction >> Poisson reconstruction callback: call in 1 step:
// - Create Poisson Delaunay Triangulation
// - Delaunay refinement
// - Solve Poisson Equation
// - Surface Meshing
void CPoissonDoc::OnOneStepPoissonReconstruction()
{
  OnCreatePoissonTriangulation();
  OnReconstructionDelaunayRefinement();
  OnReconstructionPoisson();
  OnReconstructionSurfaceMeshing();
}

void CPoissonDoc::OnUpdateOneStepPoissonReconstruction(CCmdUI *pCmdUI)
{
  OnUpdateCreatePoissonTriangulation(pCmdUI);
}

void CPoissonDoc::OnRemoveOutliers()
{
  BeginWaitCursor();
  status_message("Remove outliers (%3.1lf%%)...",m_outlier_percentage);
  double init = clock();

  // todo: use mutating version when ready
  Point_set output;
  CGAL::remove_outliers_wrt_avg_knn_sq_distance_3(
          m_points.begin(), m_points.end(),
          std::back_inserter(output),
          m_number_of_neighbours,
          m_outlier_percentage);

  m_points.clear();
  std::copy(output.begin(),output.end(),std::back_inserter(m_points));

  status_message("Remove outliers...done (%lf s)",duration(init));
  update_status();
  UpdateAllViews(NULL);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateRemoveOutliers(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}

void CPoissonDoc::OnAnalysisAverageSpacing()
{
  BeginWaitCursor();
    double value = CGAL::average_spacing_3(m_points.begin(),
                                           m_points.end(),
                                           m_number_of_neighbours,
                                           Kernel());
    //status_message("Average spacing: %lf...", value);
    char buffer[100];
    sprintf(buffer, "Average spacing: %lf", value);
    AfxMessageBox(buffer);
  EndWaitCursor();
}

void CPoissonDoc::OnUpdateAnalysisAverageSpacing(CCmdUI *pCmdUI)
{
  pCmdUI->Enable(m_edit_mode == POINT_SET);
}
