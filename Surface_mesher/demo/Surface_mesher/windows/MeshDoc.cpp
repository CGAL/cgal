// MeshDoc.cpp : implementation of the CMeshDoc class
//

#include "stdafx.h"
#include "Mesh.h"
#include "MeshDoc.h"
#include "MainFrm.h"

// STL stuff
#include <iostream>
#include <fstream>
#include "DialogOptions.h"
#include ".\meshdoc.h"

#ifdef _DEBUG
	#define new DEBUG_NEW
#endif

// CMeshDoc

IMPLEMENT_DYNCREATE(CMeshDoc, CDocument)

BEGIN_MESSAGE_MAP(CMeshDoc, CDocument)
	ON_COMMAND(ID_EDIT_OPTIONS, OnEditOptions)
	ON_COMMAND(ID_MESHING_RUN, OnMeshingRun)
END_MESSAGE_MAP()


// CMeshDoc construction/destruction
CMeshDoc::CMeshDoc()
{
	m_criterion_uniform_size = 0.1;
	m_max_nb_vertices = 5000;
	m_refresh_each = 50;
	m_init_nb_vertices = 10;
}

CMeshDoc::~CMeshDoc()
{
}


// CMeshDoc serialization

void CMeshDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{}
	else
	{}
}


// CMeshDoc diagnostics

#ifdef _DEBUG
void CMeshDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CMeshDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CMeshDoc commands

// open document file
BOOL CMeshDoc::OnOpenDocument(LPCTSTR lpszPathName)
{
	if(!CDocument::OnOpenDocument(lpszPathName))
		return FALSE;

	// get extension
	CString file = lpszPathName;
	CString extension = lpszPathName;
	extension = extension.Right(4);
	extension.MakeLower();

	// set current path
	// path "c:\path\file.wrl" -> c:\path
	CString path = lpszPathName;
	path = path.Left(path.ReverseFind('\\'));
	SetCurrentDirectory(path);

	// off extension
	if(extension == ".off")
	{
		// read from stream
		std::ifstream stream(lpszPathName);
		if(!stream)
		{
			AfxMessageBox("Unable to open file");
			return false;
		}
		stream >> m_mesh;

		// fill oracle with polyhedron
		m_oracle.init(stream);
		stream.close();
	}
	else
		{
			AfxMessageBox("Unknown extension");
			return false;
		}

	return TRUE;
}

// save file
BOOL CMeshDoc::OnSaveDocument(LPCTSTR lpszPathName)
{
	// Extension-based checking
	CString file = lpszPathName;

	// Extension
	CString extension = lpszPathName;
	extension = extension.Right(4);
	extension.MakeLower();

	// Path "c:\path\file.wrl" -> c:\path
	CString path = lpszPathName;
	path = path.Left(path.ReverseFind('\\'));

	// Current path
	SetCurrentDirectory(path);

	TRACE("\nOpening document\n");
	TRACE("File      : %s\n",lpszPathName);
	TRACE("Path      : %s\n",path);
	TRACE("Extension : %s\n",extension);

	// save off file
	if(extension == ".off")
	{
		// OFF extension
		TRACE("attempt to save OFF file %s...",lpszPathName);
		std::ofstream stream(lpszPathName);
		if(!stream)
		{
			TRACE("failed\n");
			return false;
		}

		// save the mesh
		TRACE("ok\nsave mesh...");
		stream << m_mesh;
		TRACE("..ok\n");
	}
	return TRUE;
}



//*******************************************
// User message in status bar
//*******************************************
void CMeshDoc::StatusMessage(char* fmt,...)
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

void CMeshDoc::OnEditOptions()
{
	CDialogOptions dlg;
	dlg.m_criterion_uniform_size = m_criterion_uniform_size;
	dlg.m_max_nb_vertices = m_max_nb_vertices;
	dlg.m_refresh_each = m_refresh_each;
	dlg.m_init_nb_vertices = m_init_nb_vertices;
	if(dlg.DoModal())
	{
		m_criterion_uniform_size = dlg.m_criterion_uniform_size;
		m_max_nb_vertices = dlg.m_max_nb_vertices;
		m_refresh_each = dlg.m_refresh_each;
		m_init_nb_vertices = dlg.m_init_nb_vertices;
	}
}

void CMeshDoc::OnMeshingRun()
{
  BeginWaitCursor();

	// this is fairly buggy
	// m_oracle.random_points(initial_point_sample,30);
  // random_points is ad hoc for the implicit oracle only

	// cleanup Delaunay triangulation
	m_del.clear();

  // Initial point sample
	std::list<Point> points;
	generate_initial_point_sample(points,m_init_nb_vertices);
  m_del.insert(points.begin(),points.end());
	StatusMessage("Initial number of points: %d",m_del.number_of_vertices());

  // Meshing criteria
  CGAL::Surface_mesher::Curvature_size_criterion<Del> c_s_crit (10000);
                                         // bound on Hausdorff distance
                                         // does not play any role if
                                         // bigger than the square of
                                         // the Uniform_size_criterion

  CGAL::Surface_mesher::Uniform_size_criterion<Del> u_s_crit (m_criterion_uniform_size);
                           // bound on radii of surface Delaunay balls

  CGAL::Surface_mesher::Aspect_ratio_criterion<Del> a_r_crit (30);
                          // lower bound on minimum angle in degrees

  std::vector<Criterion*> crit_vect;
  crit_vect.push_back (&a_r_crit);
  crit_vect.push_back (&u_s_crit);
  crit_vect.push_back (&c_s_crit);

  Criteria criteria(crit_vect);


  // 2D-complex in 3D-Delaunay triangulation
  C2t3 Co2(m_del);

  // Surface meshing
	StatusMessage("Refine...");
  Surface_mesher mesher(m_del, Co2, m_oracle, criteria);

	//mesher.refine_mesh(true);

	unsigned int index = m_init_nb_vertices;
	mesher.init();
	while(!mesher.is_algorithm_done() && index < m_max_nb_vertices)
	{
		int nb_facets_to_refine;
		mesher.make_one_step(nb_facets_to_refine);
		StatusMessage("%d vertices (%d facets to refine)",index++,nb_facets_to_refine);

		if(index % m_refresh_each == 0)
			refresh();
	}

	StatusMessage("Refine...(%d points)",m_del.number_of_vertices());

	// output
	char *filename = "d:\\out.off";
  std::ofstream stream(filename);
  output_surface_facets_to_off(stream,m_del);
  stream.close();

  // refresh client
	UpdateAllViews(NULL);
  EndWaitCursor();
}

void CMeshDoc::refresh(void)
{
  CMainFrame *pFrame = (CMainFrame *)AfxGetApp()->m_pMainWnd;
  CMDIChildWnd *pChild = (CMDIChildWnd *)pFrame->GetActiveFrame();
  CView *pView = (CView *)pChild->GetActiveView();
  pView->RedrawWindow(NULL,NULL,RDW_INVALIDATE | RDW_UPDATENOW);
}


void CMeshDoc::generate_initial_point_sample(std::list<Point>& points,
																						 const unsigned int nb)
{
	srand(0);
  for(unsigned int i = 0;i<nb;i++)
    points.push_back(m_mesh.random_point());
}
