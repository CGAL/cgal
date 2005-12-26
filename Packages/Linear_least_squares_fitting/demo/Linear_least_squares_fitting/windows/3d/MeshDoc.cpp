// MeshDoc.cpp : implementation of the CMeshDoc class
//

#include "stdafx.h"
#include "Mesh.h"
#include "MeshDoc.h"

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
  ON_COMMAND(ID_FIT_PLANE, OnFitPlane)
END_MESSAGE_MAP()


// CMeshDoc construction/destruction
CMeshDoc::CMeshDoc()
{
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

    // add mesh points to point set
    Mesh::Point_iterator it;
    for(it = m_mesh.points_begin();
        it != m_mesh.points_end();
        it++)
      m_points.push_back(*it);
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
	if(dlg.DoModal())
	{
		// update variable options here
	}
}

void CMeshDoc::OnFitPlane()
{
  linear_least_squares_fitting_3(m_points.begin(),
                                 m_points.end(),
                                 m_fitting_plane,
                                 m_centroid);
}
