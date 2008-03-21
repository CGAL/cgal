// GyrovizDoc_3.cpp : implementation of the CGyrovizDoc_3 class
//

// This demo
#include "stdafx.h"
#include "GyrovizDoc_3.h"
#include "Gyroviz.h"
#include "DialogOptions_3.h"

// STL
#include <iostream>
#include <fstream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CGyrovizDoc_3

IMPLEMENT_DYNCREATE(CGyrovizDoc_3, CDocument)

BEGIN_MESSAGE_MAP(CGyrovizDoc_3, CDocument)
	ON_COMMAND(ID_EDIT_OPTIONS, OnEditOptions)
  ON_COMMAND(ID_FILE_SAVE_SURFACE, OnFileSaveSurface)
  ON_UPDATE_COMMAND_UI(ID_FILE_SAVE_SURFACE, OnUpdateFileSaveSurface)
  ON_COMMAND(ID_FILE_SAVE_AS, OnFileSaveAs)
  ON_UPDATE_COMMAND_UI(ID_FILE_SAVE_AS, OnUpdateFileSaveAs)
END_MESSAGE_MAP()


// CGyrovizDoc_3 construction/destruction

CGyrovizDoc_3::CGyrovizDoc_3()
{
	// options
	m_point_size = 2; // OpenGL point size
}

CGyrovizDoc_3::~CGyrovizDoc_3()
{
}

// CGyrovizDoc_3 diagnostics

#ifdef _DEBUG
void CGyrovizDoc_3::AssertValid() const
{
	CDocument::AssertValid();
}

void CGyrovizDoc_3::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// File >> Open implementation
BOOL CGyrovizDoc_3::OnOpenDocument(LPCTSTR lpszPathName)
{
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

  // if .pnt extension
	if(extension.CompareNoCase(".pnt") == 0)
	{
		double init = clock();
		if(!m_gyroviz_dt.read_pnt((char *)lpszPathName))
		{
			AfxMessageBox("Unable to open file");
			return FALSE;
		}
		status_message("Delaunay triangulation (%lf s)",duration(init));
	}
	
	// if .pwc extension
	else if (extension.CompareNoCase(".pwc") == 0)
	{
	double init = clock();
		if(!m_gyroviz_dt3.read_pwc((char *)lpszPathName))
		{
			AfxMessageBox("Unable to open file");
			return FALSE;
		}
		status_message("3D Delaunay triangulation (%lf s)",duration(init));
	}
	
  else
	{
		AfxMessageBox("File format not supported");
		return FALSE;
	}

  update_status();
	UpdateAllViews(NULL);
  return TRUE;
}

// Save input point set as...  callback
void CGyrovizDoc_3::OnFileSaveAs()
{
}

// Disable "Save input point set as..." if Delaunay refinement was applied
// to avoid jeopardizing the input file.
void CGyrovizDoc_3::OnUpdateFileSaveAs(CCmdUI *pCmdUI)
{
/*
	pCmdUI->Enable(!m_gyroviz_solved);
*/
}

// Save reconstructed surface as...  callback
void CGyrovizDoc_3::OnFileSaveSurface()
{
}

// Enable "Save reconstructed surface as..." if surface is computed
void CGyrovizDoc_3::OnUpdateFileSaveSurface(CCmdUI *pCmdUI)
{
/*
	pCmdUI->Enable(m_surface_mesher_dt.number_of_vertices() > 0);
*/
}

// Update the number of vertices and faces in the status bar
void CGyrovizDoc_3::update_status()
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
			vertices.Format("%d vertices",m_gyroviz_dt.number_of_vertices());

			CString faces;
			faces.Format("%d faces",m_gyroviz_dt.number_of_faces());

			// Update status bar
			pStatus->SetPaneText(1,vertices);
			pStatus->SetPaneText(2,faces);
			pStatus->UpdateWindow(); 
		}
  }
}

// Set user message in status bar
void CGyrovizDoc_3::status_message(char* fmt,...)
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
void CGyrovizDoc_3::OnEditOptions()
{
	CDialogOptions_3 dlg;
	dlg.m_point_size = m_point_size;

	if(dlg.DoModal() == IDOK)
	{
		m_point_size = dlg.m_point_size;
	}
}

double CGyrovizDoc_3::duration(const double time_init)
{
  return (clock() - time_init)/CLOCKS_PER_SEC;
}
