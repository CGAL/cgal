// MeshDoc.cpp : implementation of the CMeshDoc class
//

#include "stdafx.h"
#include "Mesh.h"
#include "MeshDoc.h"

// STL stuff
#include <iostream>
#include <fstream>

// CGAL (off) + obj IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include "Lib/parser_obj.h"

// subdivision
#include "Lib/sqrt3.h"
#include "Lib/quad-triangle.h"
#include "Lib/SurfLab/Polyhedron_subdivision.h"
#include ".\meshdoc.h"

#ifdef _DEBUG
	#define new DEBUG_NEW
#endif


// CMeshDoc

IMPLEMENT_DYNCREATE(CMeshDoc, CDocument)

BEGIN_MESSAGE_MAP(CMeshDoc, CDocument)
	ON_COMMAND(ID_SUBDIVISION_SQRT3, OnSubdivisionSqrt3)
	ON_COMMAND(ID_SUBDIVISION_SQRT3_TWICE, OnSubdivisionSqrt3Twice)
	ON_COMMAND(ID_SUBDIVISION_QUAD, OnSubdivisionQuad)
	ON_UPDATE_COMMAND_UI(ID_SUBDIVISION_QUAD, OnUpdateSubdivisionQuad)
	ON_UPDATE_COMMAND_UI(ID_SUBDIVISION_SQRT3, OnUpdateSubdivisionSqrt3)
	ON_UPDATE_COMMAND_UI(ID_SUBDIVISION_SQRT3_TWICE, OnUpdateSubdivisionSqrt3Twice)
	ON_COMMAND(ID_SUBDIVISION_DOO, OnSubdivisionDoo)
	ON_UPDATE_COMMAND_UI(ID_SUBDIVISION_DOO, OnUpdateSubdivisionDoo)
	ON_COMMAND(ID_SUBDIVISION_LOOP, OnSubdivisionLoop)
	ON_UPDATE_COMMAND_UI(ID_SUBDIVISION_LOOP, OnUpdateSubdivisionLoop)
	ON_COMMAND(ID_SUBDIVISION_CATMULL32868, OnSubdivisionCatmull32868)
	ON_UPDATE_COMMAND_UI(ID_SUBDIVISION_CATMULL32868, OnUpdateSubdivisionCatmull32868)
END_MESSAGE_MAP()


// CMeshDoc construction/destruction
CMeshDoc::CMeshDoc()
{
	m_pMesh = NULL;
}

CMeshDoc::~CMeshDoc()
{
	delete m_pMesh;
	ResetMeshProperties();
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
	
	// allocate a new mesh
	m_pMesh = new Enriched_polyhedron<Enriched_kernel,Enriched_items>;
	CGAL_assertion(m_pMesh != NULL);

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
		stream >> *m_pMesh;
	}
	else
		if(extension == ".obj")
		{
			Parser_obj<Enriched_kernel,Enriched_items> parser;
			parser.read(lpszPathName,m_pMesh);
		}
		else
		{
			AfxMessageBox("Unknown extension");
			return false;
		}

  // update mesh properties in the status bar 
	// and compute normals
	m_pMesh->compute_type();
	UpdateMeshProperties(true,true);
	m_pMesh->compute_normals();

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
		stream << *m_pMesh;
		TRACE("..ok\n");
	}

	// save obj file 
	if(extension == ".obj")
	{
		StatusMessage("Saving OBJ file...");
		m_pMesh->write_obj((char *)lpszPathName);
		StatusMessage("Saving OBJ file...done");
	}
	return TRUE;
}

// sqrt3 subdivision
void CMeshDoc::OnSubdivisionSqrt3()
{
	BeginWaitCursor();
	StatusMessage("Subdivide...");
	CSubdivider_sqrt3<Enriched_polyhedron<Enriched_kernel,Enriched_items>,Enriched_kernel> subdivider;
	subdivider.subdivide(*m_pMesh,1); // one iteration
	m_pMesh->compute_normals();
	UpdateMeshProperties();
	UpdateAllViews(NULL);
	EndWaitCursor();
	StatusMessage("Subdivide...done");
}
void CMeshDoc::OnUpdateSubdivisionSqrt3(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pMesh != NULL && m_pMesh->is_pure_triangle());
}

// sqrt3 subdivision applied twice 
// (boundary edges are trisected)
void CMeshDoc::OnSubdivisionSqrt3Twice()
{
	BeginWaitCursor();
	StatusMessage("Subdivide...");
	CSubdivider_sqrt3<Enriched_polyhedron<Enriched_kernel,Enriched_items>,Enriched_kernel> subdivider;
	subdivider.subdivide(*m_pMesh,2); // two iterations
	m_pMesh->compute_normals();
	UpdateMeshProperties();
	UpdateAllViews(NULL);
	EndWaitCursor();
	StatusMessage("Subdivide...done");
}
void CMeshDoc::OnUpdateSubdivisionSqrt3Twice(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pMesh != NULL && m_pMesh->is_pure_triangle());
}

// quad/triangle subdivision
void CMeshDoc::OnSubdivisionQuad()
{
	double start = clock();

	// subdivision engine
	CSubdivider_quad_triangle<Enriched_polyhedron<Enriched_kernel,Enriched_items>,Enriched_kernel> subdivider;

	BeginWaitCursor();
	StatusMessage("Quad/triangle subdivision...");

	// alloc a new mesh
  Enriched_polyhedron<Enriched_kernel,Enriched_items> *pNewMesh = 
		new Enriched_polyhedron<Enriched_kernel,Enriched_items>;

	// subdivide once
	subdivider.subdivide(*m_pMesh,*pNewMesh,true);

	// copy bounding box (approximate, but fast)
	pNewMesh->copy_bounding_box(m_pMesh);

	// delete previous mesh
	delete m_pMesh;

	// set new mesh
	m_pMesh = pNewMesh;

	// compute normals
	m_pMesh->compute_normals();
	UpdateMeshProperties();

	EndWaitCursor();
	float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
	StatusMessage("Quad/triangle subdivision...done (%g s)",duration);
	UpdateAllViews(NULL);
}
void CMeshDoc::OnUpdateSubdivisionQuad(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pMesh != NULL);
}

// Doo-Sabin
void CMeshDoc::OnSubdivisionDoo()
{
	double start = clock();

	BeginWaitCursor();
	StatusMessage("Doo-Sabin subdivision...");

  Polyhedron_subdivision<Enriched_polyhedron<Enriched_kernel,Enriched_items> >::DooSabin_subdivision(*m_pMesh,1);

	// compute normals
	m_pMesh->compute_normals();
	UpdateMeshProperties();

	EndWaitCursor();
	float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
	StatusMessage("Doo-Sabin subdivision...done (%g s)",duration);
	UpdateAllViews(NULL);
}
void CMeshDoc::OnUpdateSubdivisionDoo(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pMesh != NULL);
}

// Loop subdivision
void CMeshDoc::OnSubdivisionLoop()
{
	double start = clock();

	BeginWaitCursor();
	StatusMessage("Loop subdivision...");

  Polyhedron_subdivision<Enriched_polyhedron<Enriched_kernel,Enriched_items> >::Loop_subdivision(*m_pMesh,1);

	// compute normals
	m_pMesh->compute_normals();
	UpdateMeshProperties();

	EndWaitCursor();
	float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
	StatusMessage("Loop subdivision...done (%g s)",duration);
	UpdateAllViews(NULL);
}
void CMeshDoc::OnUpdateSubdivisionLoop(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pMesh->is_pure_triangle());
}

// Catmull-Clark
void CMeshDoc::OnSubdivisionCatmull32868()
{
	double start = clock();

	BeginWaitCursor();
	StatusMessage("Catmull-Clark subdivision...");

  Polyhedron_subdivision<Enriched_polyhedron<Enriched_kernel,Enriched_items> >::CatmullClark_subdivision(*m_pMesh,1);

	// compute normals
	m_pMesh->compute_normals();
	UpdateMeshProperties();

	EndWaitCursor();
	float duration = (float)((clock()-start)/CLOCKS_PER_SEC);
	StatusMessage("Catmull-Clark...done (%g s)",duration);
	UpdateAllViews(NULL);
}
void CMeshDoc::OnUpdateSubdivisionCatmull32868(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pMesh != NULL);
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

//*******************************************
// update mesh properties in status bar
//*******************************************
void CMeshDoc::UpdateMeshProperties(bool update_component,
																		bool update_boundary)
{   
	CWinApp *pApp = AfxGetApp();
	if(pApp->m_pMainWnd != NULL) 
	{ 
		CStatusBar* pStatus = 
			(CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(
			AFX_IDW_STATUS_BAR);
		
		if(pStatus != NULL) 
		{
			static unsigned int c = 0;
			if(update_component)
			  c = m_pMesh->nb_components();
			static unsigned int b = 0;
			if(update_boundary)
				b = m_pMesh->nb_boundaries();
			unsigned int v = m_pMesh->size_of_vertices();
			unsigned int e = m_pMesh->size_of_halfedges()/2;
			unsigned int f = m_pMesh->size_of_facets();
			unsigned int g = m_pMesh->genus(c,v,f,e,b);

			// components
			CString components;
			components.Format("%d component%c",c,(c>1) ? 's' : ' ');

			// vertices
			CString vertices;
			vertices.Format("%d vertices",v);

			// facets
			CString facets;
			facets.Format("%d facets",f);

			// edges
			CString edges;
			edges.Format("%d edges",e);

			// boundaries
			CString boundaries;
			if(b == 0)
				boundaries.Format("no boundary");
			else
				if(b == 1)
					boundaries.Format("1 boundary");
				else
					boundaries.Format("%d boundaries",b);

			// genus
			CString genus;
			genus.Format("genus %d",g);
			pStatus->SetPaneText(1,components);
			pStatus->SetPaneText(2,vertices);
			pStatus->SetPaneText(3,facets);
			pStatus->SetPaneText(4,edges);
			pStatus->SetPaneText(5,boundaries);
			pStatus->SetPaneText(6,genus);
			pStatus->UpdateWindow(); 
		}
  }
}

//*******************************************
// reset mesh properties in status bar
//*******************************************
void CMeshDoc::ResetMeshProperties()
{   
	if(AfxGetApp()->m_pMainWnd != NULL) 
	{ 
		CStatusBar* pStatus = 
			(CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(
			AFX_IDW_STATUS_BAR);
		if(pStatus != NULL) 
		{
			pStatus->SetPaneText(1,CString(""));
			pStatus->SetPaneText(2,CString(""));
			pStatus->SetPaneText(3,CString(""));
			pStatus->SetPaneText(4,CString(""));
			pStatus->SetPaneText(5,CString(""));
			pStatus->SetPaneText(6,CString(""));
			pStatus->UpdateWindow(); 
		}
  }
}



