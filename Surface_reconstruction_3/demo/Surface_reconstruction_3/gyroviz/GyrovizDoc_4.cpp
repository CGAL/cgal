// GyrovizDoc_4.cpp : implementation of the CGyrovizDoc_4 class
//

#include "stdafx.h"
#include "GyrovizDoc_4.h"
#include "MainFrm.h"
#include "DialogOptions_4.h"

using namespace cimg_library;

#include "Gyroviz_segmentation2.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CGyrovizDoc_4

IMPLEMENT_DYNCREATE(CGyrovizDoc_4, CDocument)

BEGIN_MESSAGE_MAP(CGyrovizDoc_4, CDocument)
END_MESSAGE_MAP()


// CGyrovizDoc_4 construction/destruction

CGyrovizDoc_4::CGyrovizDoc_4()
{
}

CGyrovizDoc_4::~CGyrovizDoc_4()
{
}

BOOL CGyrovizDoc_4::OnNewDocument()
{
  if (!CDocument::OnNewDocument())
    return FALSE;
  // TODO: add reinitialization code here
  // (SDI documents will reuse this document)
  return TRUE;
}

// CGyrovizDoc_4 serialization

void CGyrovizDoc_4::Serialize(CArchive& ar)
{
  if (ar.IsStoring())
  {
    // TODO: add storing code here
  }
  else
  {
    // TODO: add loading code here
  }
}


// CGyrovizDoc_4 diagnostics

#ifdef _DEBUG
void CGyrovizDoc_4::AssertValid() const
{
  CDocument::AssertValid();
}

void CGyrovizDoc_4::Dump(CDumpContext& dc) const
{
  CDocument::Dump(dc);
}
#endif //_DEBUG

// CGyrovizDoc_4 commands


// Update the number of vertices and faces in the status bar
void CGyrovizDoc_4::update_status()
{   
	CWinApp *pApp = AfxGetApp();
	if(pApp->m_pMainWnd != NULL) 
	{ 
		CStatusBar* pStatus = 
			(CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(
			AFX_IDW_STATUS_BAR);
		
		if(pStatus != NULL) 
		{
			//CString vertices;
			//vertices.Format("%d vertices",m_gyroviz_dt.number_of_vertices());

			//CString faces;
			//faces.Format("%d faces",m_gyroviz_dt.number_of_faces());

			//// Update status bar
			//pStatus->SetPaneText(1,vertices);
			//pStatus->SetPaneText(2,faces);
			//pStatus->UpdateWindow(); 
		}
  }
}

// User message in status bar
void CGyrovizDoc_4::status_message(char* fmt,...)
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


BOOL CGyrovizDoc_4::OnSaveDocument(LPCTSTR lpszPathName)
{
  // save pslg to a file
  //return m_pslg.save((char *)lpszPathName);
  return CDocument::OnSaveDocument(lpszPathName);
}


BOOL CGyrovizDoc_4::OnOpenDocument(LPCTSTR lpszPathName)
{
  //
  // Read bitmap image
  //
  
  if (!CDocument::OnOpenDocument(lpszPathName))
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

  if(extension == ".bmp")
  {// read bmp format image
    m_cimg_interm_image = CImg<unsigned char>((char *)lpszPathName);
    m_original_image = cimg_image_multiplexer_char(m_cimg_interm_image);
    
    m_cimg_gray_image = to_grayscale(m_cimg_interm_image);
    m_grayscaled_image = cimg_image_multiplexer_char(m_cimg_gray_image);

    m_cimg_filt_image = gauss3(m_cimg_gray_image);
    m_filtered_image = cimg_image_multiplexer_char(m_cimg_filt_image);

    m_cimg_seg_image = grad_freiChen(m_cimg_filt_image,15); // Edge detection on filtered image def : 30
    m_segmented_image = cimg_image_multiplexer_char(m_cimg_seg_image);
  }
  else
	{
		AfxMessageBox("File format not supported");
		return FALSE;
	}

  //
  // Read corresponding Voodoo 2D feature file (.pnt)
  //
  
  // file filters
  static char szFilter[] = "Voodoo 2D Feature Files (*.pnt)|*.pnt|All Files (*.*)|*.*||";

  // create the Open dialog
  CFileDialog dlgOpen(true, "pnt", NULL, 
                      OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST, szFilter, AfxGetMainWnd());
                                
  // dialog title
  dlgOpen.m_ofn.lpstrTitle = "Select corresponding Voodoo 2D feature file";

  // show the dialog
  if (dlgOpen.DoModal() == IDOK)
  {
	  // get extension
	  CString file = dlgOpen.m_ofn.lpstrFile;
	  CString extension = dlgOpen.m_ofn.lpstrFile;
	  extension = extension.Right(4);
	  extension.MakeLower();
  	
	  // set current path
	  CString path = dlgOpen.m_ofn.lpstrFile;
	  path = path.Left(path.ReverseFind('\\'));
	  SetCurrentDirectory(path);

    // if .pnt extension
	  if(extension.CompareNoCase(".pnt") == 0)
	  {
		  double init = clock();
		  if(!m_gyroviz_dt.read_pnt((char *)dlgOpen.m_ofn.lpstrFile))
		  {
			  AfxMessageBox("Unable to open file");
			  return FALSE;
		  }
		  m_gyroviz_dt.nw_add_constraints(m_cimg_seg_image, 5);
		  status_message("Constrained Delaunay triangulation (%lf s)",duration(init));
	  
	  }
    else
	  {
		  AfxMessageBox("File format not supported");
		  return FALSE;
	  }
	}

  update_status();
  UpdateAllViews(NULL);
  return TRUE;
}

double CGyrovizDoc_4::duration(const double time_init)
{
  return (clock() - time_init)/CLOCKS_PER_SEC;
}

