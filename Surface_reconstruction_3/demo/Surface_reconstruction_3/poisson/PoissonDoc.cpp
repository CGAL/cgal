// PoissonDoc.cpp : implementation of the CPoissonDoc class
//

// This demo
#include "stdafx.h"
#include "Poisson.h"
#include "DialogOptions.h"
#include "PoissonDoc.h"
#include "enriched_polyhedron.h"

// CGAL
#define CGAL_C2T3_USE_POLYHEDRON
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_scanner_OFF.h>

// This package
#include <CGAL/surface_reconstruction_output.h>

// STL
#include <iostream>
#include <fstream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CPoissonDoc

IMPLEMENT_DYNCREATE(CPoissonDoc, CDocument)

BEGIN_MESSAGE_MAP(CPoissonDoc, CDocument)
	ON_COMMAND(ID_RECONSTRUCTION_DELAUNAYREFINEMENT, OnReconstructionDelaunayrefinement)
	ON_COMMAND(ID_RECONSTRUCTION_POISSON, OnReconstructionPoisson)
	ON_COMMAND(ID_ALGORITHMS_REFINEINSHELL, OnAlgorithmsRefineinshell)
	ON_COMMAND(ID_RECONSTRUCTION_SURFACEMESHING, OnReconstructionSurfacemeshing)
	ON_COMMAND(ID_EDIT_OPTIONS, OnEditOptions)
  ON_UPDATE_COMMAND_UI(ID_RECONSTRUCTION_POISSON, OnUpdateReconstructionPoisson)
  ON_UPDATE_COMMAND_UI(ID_RECONSTRUCTION_SURFACEMESHING, OnUpdateReconstructionSurfacemeshing)
  ON_COMMAND(ID_ALGORITHMS_MARCHINGTETCONTOURING, OnAlgorithmsMarchingtetcontouring)
  ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_MARCHINGTETCONTOURING, OnUpdateAlgorithmsMarchingtetcontouring)
  ON_COMMAND(ID_FILE_SAVE_SURFACE, OnFileSaveSurface)
  ON_UPDATE_COMMAND_UI(ID_FILE_SAVE_SURFACE, OnUpdateFileSaveSurface)
  ON_COMMAND(ID_FILE_SAVE_AS, OnFileSaveAs)
  ON_UPDATE_COMMAND_UI(ID_FILE_SAVE_AS, OnUpdateFileSaveAs)
	ON_COMMAND(ID_ALGORITHMS_EXTRAPOLATENORMALS, OnAlgorithmsExtrapolatenormals)
	ON_COMMAND(ID_ALGORITHMS_POISSONSTATISTICS, OnAlgorithmsPoissonStatistics)
  ON_UPDATE_COMMAND_UI(ID_ALGORITHMS_POISSONSTATISTICS, OnUpdateAlgorithmsPoissonstatistics)
  ON_COMMAND(ID_ALGORITHMS_ESTIMATENORMALSBYPCA, OnAlgorithmsEstimateNormalsByPCA)
  ON_COMMAND(ID_ALGORITHMS_ESTIMATENORMALBYJETFITTING, OnAlgorithmsEstimateNormalsByJetFitting)
  ON_COMMAND(ID_ALGORITHMS_ORIENTNORMALSCAMERAS, OnAlgorithmsOrientNormalsWrtCameras)
  ON_COMMAND(ID_ALGORITHMS_ORIENTNORMALSMST, &CPoissonDoc::OnAlgorithmsOrientNormalsWithMST)
END_MESSAGE_MAP()


// CPoissonDoc construction/destruction

CPoissonDoc::CPoissonDoc()
: m_poisson_function(m_poisson_dt), 
  m_surface_mesher_c2t3(m_surface_mesher_dt)
{
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
      typedef Enriched_polyhedron<Kernel,Enriched_items> Polyhedron;
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

      // Copy points to m_poisson_dt
      typedef Dt3::Point_with_normal Point_with_normal;
		  std::vector<Point_with_normal> pwns;
		  Polyhedron::Vertex_iterator v;
		  for(v = input_mesh.vertices_begin();
				  v != input_mesh.vertices_end();
				  v++)
		  {
			  const Point& p = v->point();
			  const Vector& n = v->normal();
			  pwns.push_back(Point_with_normal(p,n));
		  }
      m_poisson_dt.insert(pwns.begin(), pwns.end(), Dt3::INPUT);
    }
    else // Read OFF file as a point cloud
    {
		  if(!m_poisson_dt.read_off_point_cloud(lpszPathName))
		  {
			  AfxMessageBox("Unable to read file");
			  return FALSE;
		  }
    }
	}
  // if .pwn extension
	else if(extension.CompareNoCase(".pwn") == 0)
	{
		if(!m_poisson_dt.read_pwn((char *)lpszPathName))
		{
			AfxMessageBox("Unable to read file");
			return FALSE;
		}
	}
  // if .xyz extension
	else if(extension.CompareNoCase(".xyz") == 0)
	{
		if(!m_poisson_dt.read_xyz(lpszPathName))
		{
			AfxMessageBox("Unable to read file");
			return FALSE;
		}
	}
  // if .pnb extension
	else if(extension.CompareNoCase(".pnb") == 0)
	{
		if(!m_poisson_dt.read_pnb(lpszPathName))
		{
			AfxMessageBox("Unable to read file");
			return FALSE;
		}
	}
	// if Gyroviz .pwc extension
	else if (extension.CompareNoCase(".pwc") == 0)
	{
		if(!m_poisson_dt.read_pwc(lpszPathName))
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

  status_message("Delaunay triangulation (%lf s)",duration(init));
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
		  if(!m_poisson_dt.save_pwn(dlgExport.m_ofn.lpstrFile))
	    {
		    AfxMessageBox("Unable to save file");
		    return;
	    }
	  }
	  // if .xyz extension
    else if(extension.CompareNoCase(".xyz") == 0)
	  {
		  if(!m_poisson_dt.save_xyz(dlgExport.m_ofn.lpstrFile))
	    {
		    AfxMessageBox("Unable to save file");
		    return;
	    }
	  }
	  // if .off extension
    else if(extension.CompareNoCase(".off") == 0)
	  {
		  if(!m_poisson_dt.save_off_point_cloud(dlgExport.m_ofn.lpstrFile))
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

// Disable "Save input point set as..." if Delaunay refinement was applied
// to avoid jeopardizing the input file.
void CPoissonDoc::OnUpdateFileSaveAs(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(!m_triangulation_refined);
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
	pCmdUI->Enable(m_surface_mesher_dt.number_of_vertices() > 0);
}

// Update the number of vertices and tetrahedra in the status bar
void CPoissonDoc::update_status()
{   
	CWinApp *pApp = AfxGetApp();
	if(pApp->m_pMainWnd != NULL) 
	{ 
		CStatusBar* pStatus = 
			(CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(
			AFX_IDW_STATUS_BAR);
		
		if(pStatus != NULL) 
		{
			CString vertices;
			vertices.Format("%d vertices",m_poisson_dt.number_of_vertices());

			CString tets;
			tets.Format("%d tets",m_poisson_dt.number_of_cells());

			// Update status bar
			pStatus->SetPaneText(1,vertices);
			pStatus->SetPaneText(2,tets);
			pStatus->UpdateWindow(); 
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
	}
}

// Utility: compute elapsed time
double CPoissonDoc::duration(const double time_init)
{
  return (clock() - time_init)/CLOCKS_PER_SEC;
}

// Compute normals direction by Principal Component Analysis
void CPoissonDoc::OnAlgorithmsEstimateNormalsByPCA()
{
	BeginWaitCursor();
	status_message("Estimate Normals Direction...");
	double init = clock();

	m_poisson_function.estimate_normals_pca(m_number_of_neighbours);

	status_message("Estimate Normals Direction...done (%lf s)",duration(init));
  update_status();
	UpdateAllViews(NULL);
	EndWaitCursor();
}

// Compute normals direction by Jet Fitting
void CPoissonDoc::OnAlgorithmsEstimateNormalsByJetFitting()
{
	BeginWaitCursor();
	status_message("Estimate Normals Direction...");
	double init = clock();

	m_poisson_function.estimate_normals_jet_fitting(m_number_of_neighbours);

	status_message("Estimate Normals Direction...done (%lf s)",duration(init));
  update_status();
	UpdateAllViews(NULL);
	EndWaitCursor();
}

/// Orient the normals using a minimum spanning tree.
void CPoissonDoc::OnAlgorithmsOrientNormalsWithMST()
{
	BeginWaitCursor();
	status_message("Orient Normals with MST...");
	double init = clock();

	m_poisson_dt.orient_normals_minimum_spanning_tree(m_number_of_neighbours);

	status_message("Orient Normals with MST...done (%lf s)",duration(init));
  update_status();
	UpdateAllViews(NULL);
	EndWaitCursor();
}

/// Specific to Gyroviz: orient the normals w.r.t. the position of cameras
/// that reconstructed the points by photogrammetry.
void CPoissonDoc::OnAlgorithmsOrientNormalsWrtCameras()
{
	BeginWaitCursor();
	status_message("Orient Normals wrt Cameras...");
	double init = clock();

	m_poisson_dt.orient_normals_wrt_cameras();

	status_message("Orient Normals wrt Cameras...done (%lf s)",duration(init));
  update_status();
	UpdateAllViews(NULL);
	EndWaitCursor();
}

// Uniform Delaunay refinement
void CPoissonDoc::OnReconstructionDelaunayrefinement()
{
	BeginWaitCursor();

	status_message("Delaunay refinement...");
	const double quality = 2.5;
	const unsigned int max_vertices = (unsigned int)1e7; // max 10M vertices
	const double enlarge_ratio = 1.5;
	double init = clock();
	unsigned int nb_vertices_added = m_poisson_function.delaunay_refinement(quality,max_vertices,enlarge_ratio,50000);
	status_message("Delaunay refinement...done (%lf s, %d vertices inserted)",duration(init),nb_vertices_added);
	m_triangulation_refined = true;

  update_status();
	UpdateAllViews(NULL);
	EndWaitCursor();
}

// Delaunay refinement in a surface's shell
void CPoissonDoc::OnAlgorithmsRefineinshell()
{
	BeginWaitCursor();

	status_message("Delaunay refinement...");
	const double quality = 2.5;
	const unsigned int max_vertices = (unsigned int)1e7; // max 10M vertices
	const double enlarge_ratio = 1.5;
	double init = clock();
	unsigned int nb_vertices_added = m_poisson_function.delaunay_refinement_shell(m_dr_shell_size,m_dr_sizing,m_dr_max_vertices);
	status_message("Delaunay refinement...done (%lf s, %d vertices inserted)",duration(init),nb_vertices_added);
	m_triangulation_refined = true;

  update_status();
	UpdateAllViews(NULL);
	EndWaitCursor();
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
	pCmdUI->Enable(m_triangulation_refined);
}

// Surface Meshing callback
void CPoissonDoc::OnReconstructionSurfacemeshing()
{
    typedef CGAL::Implicit_surface_3<Kernel, Poisson_implicit_function&> Surface_3;

    // Clear previous call
    m_surface_mesher_dt.clear();
    m_surface_mesher_c2t3.clear();
    
    // Apply contouring value defined in Options dialog and reset it
		m_poisson_function.set_contouring_value(m_contouring_value);
		m_contouring_value = 0.0;

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

	  UpdateAllViews(NULL);
		EndWaitCursor();
}

// Enable "Surface Meshing" if Poisson equation is solved
void CPoissonDoc::OnUpdateReconstructionSurfacemeshing(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_poisson_solved);
}

// Marching Tet Contouring callback
void CPoissonDoc::OnAlgorithmsMarchingtetcontouring()
{
		int nb = m_poisson_dt.marching_tet(m_contouring_value);
		status_message("Marching tet contouring...done (%d triangles)",nb);
	  UpdateAllViews(NULL);
}

// Enable "Marching Tet Contouring" if Poisson equation is solved
void CPoissonDoc::OnUpdateAlgorithmsMarchingtetcontouring(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_poisson_solved);
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
void CPoissonDoc::OnUpdateAlgorithmsPoissonstatistics(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_poisson_solved);
}
