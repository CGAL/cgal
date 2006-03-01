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
  ON_COMMAND(ID_FIT_FITPOINTSET, OnFitFitpointset)
  ON_COMMAND(ID_FIT_TRIANGLESET32899, OnFitTriangleset32899)
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
    for(Mesh::Point_iterator it = m_mesh.points_begin();
        it != m_mesh.points_end();
        it++)
      m_points.push_back(*it);

    // add mesh triangles to triangle set
    for(Mesh::Facet_iterator f = m_mesh.facets_begin();
        f != m_mesh.facets_end();
        f++)
    {
      const Point& a = f->halfedge()->vertex()->point();
      const Point& b = f->halfedge()->next()->vertex()->point();
      const Point& c = f->halfedge()->next()->next()->vertex()->point();
      m_triangles.push_back(Triangle(a,b,c));
    }
	}
	else
		{
			AfxMessageBox("Unknown extension");
			return false;
		}

  OnFitFitpointset();


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

void CMeshDoc::gl_draw_fitting_primitives()
{
  ::glLineWidth(2.0f);
  ::glColor3ub(255,0,0);

  Vector b1 = m_fitting_plane.base1();
  Vector b2 = m_fitting_plane.base2();
  b1 = b1 / std::sqrt(b1 * b1);
  b2 = b2 / std::sqrt(b2 * b2);
  Point c = m_centroid;
  Point p1 = c - b1 - b2;
  Point p2 = c + b1 - b2;
  Point p3 = c + b1 + b2;
  Point p4 = c - b1 + b2;

  ::glBegin(GL_LINE_LOOP);
    ::glVertex3d(p1.x(),p1.y(),p1.z());
    ::glVertex3d(p2.x(),p2.y(),p2.z());
    ::glVertex3d(p3.x(),p3.y(),p3.z());
    ::glVertex3d(p4.x(),p4.y(),p4.z());
  ::glEnd();

  ::glColor3ub(0,255,0);

  Point p5 = m_fitting_line.point(1);
  Point p6 = m_fitting_line.point(-1);
  ::glBegin(GL_LINES);
    ::glVertex3d(p5.x(),p5.y(),p5.z());
    ::glVertex3d(p6.x(),p6.y(),p6.z());
  ::glEnd();
}




void CMeshDoc::OnFitFitpointset()
{
  linear_least_squares_fitting_3(m_points.begin(),
                                 m_points.end(),
                                 m_fitting_line);
  linear_least_squares_fitting_3(m_points.begin(),
                                 m_points.end(),
                                 m_fitting_line,
                                 m_centroid);

  linear_least_squares_fitting_3(m_points.begin(),
                                 m_points.end(),
                                 m_fitting_plane);
  linear_least_squares_fitting_3(m_points.begin(),
                                 m_points.end(),
                                 m_fitting_plane,
                                 m_centroid);
  UpdateAllViews(NULL);
}

void CMeshDoc::OnFitTriangleset32899()
{
  linear_least_squares_fitting_3(m_triangles.begin(),
                                 m_triangles.end(),
                                 m_fitting_line);
  linear_least_squares_fitting_3(m_triangles.begin(),
                                 m_triangles.end(),
                                 m_fitting_line,
                                 m_centroid);

  linear_least_squares_fitting_3(m_triangles.begin(),
                                 m_triangles.end(),
                                 m_fitting_plane);
  linear_least_squares_fitting_3(m_triangles.begin(),
                                 m_triangles.end(),
                                 m_fitting_plane,
                                 m_centroid);
  UpdateAllViews(NULL);
}
