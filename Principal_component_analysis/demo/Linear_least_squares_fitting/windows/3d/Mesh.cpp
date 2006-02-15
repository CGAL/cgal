// Mesh.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "Mesh.h"
#include "MainFrm.h"

#include "ChildFrm.h"
#include "MeshDoc.h"
#include "MeshView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CMeshApp

BEGIN_MESSAGE_MAP(CMeshApp, CWinApp)
	// Standard file based document commands
	ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
	ON_COMMAND(ID_FILE_OPEN, CWinApp::OnFileOpen)
	// Standard print setup command
	ON_COMMAND(ID_FILE_PRINT_SETUP, CWinApp::OnFilePrintSetup)
	ON_COMMAND(ID_CAMERA_COPYVIEWPOINT, OnCameraCopyviewpoint)
	ON_COMMAND(ID_CAMERA_PASTEVIEWPOINT, OnCameraPasteviewpoint)
	ON_UPDATE_COMMAND_UI(ID_CAMERA_COPYVIEWPOINT, OnUpdateCameraCopyviewpoint)
	ON_UPDATE_COMMAND_UI(ID_CAMERA_PASTEVIEWPOINT, OnUpdateCameraPasteviewpoint)
END_MESSAGE_MAP()


// CMeshApp construction

CMeshApp::CMeshApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}


// The one and only CMeshApp object

CMeshApp theApp;

// CMeshApp initialization

BOOL CMeshApp::InitInstance()
{
	// InitCommonControls() is required on Windows XP if an application
	// manifest specifies use of ComCtl32.dll version 6 or later to enable
	// visual styles.  Otherwise, any window creation will fail.
	InitCommonControls();

	CWinApp::InitInstance();

	// Standard initialization
	// If you are not using these features and wish to reduce the size
	// of your final executable, you should remove from the following
	// the specific initialization routines you do not need
	// Change the registry key under which our settings are stored
	// TODO: You should modify this string to be something appropriate
	// such as the name of your company or organization
	SetRegistryKey(_T("Local AppWizard-Generated Applications"));
	LoadStdProfileSettings(16);  // Load standard INI file options (including MRU)
	// Register the application's document templates.  Document templates
	//  serve as the connection between documents, frame windows and views
	CMultiDocTemplate* pDocTemplate;
	pDocTemplate = new CMultiDocTemplate(IDR_MeshTYPE,
		RUNTIME_CLASS(CMeshDoc),
		RUNTIME_CLASS(CChildFrame), // custom MDI child frame
		RUNTIME_CLASS(CMeshView));
	if (!pDocTemplate)
		return FALSE;
	AddDocTemplate(pDocTemplate);
	// create main MDI Frame window
	CMainFrame* pMainFrame = new CMainFrame;
	if (!pMainFrame || !pMainFrame->LoadFrame(IDR_MAINFRAME))
		return FALSE;
	m_pMainWnd = pMainFrame;

	// call DragAcceptFiles only if there's a suffix
	//  In an MDI app, this should occur immediately after setting m_pMainWnd
	// Enable drag/drop open
	m_pMainWnd->DragAcceptFiles();
	// Enable DDE Execute open
	EnableShellOpen();
	RegisterShellFileTypes(TRUE);
	// Parse command line for standard shell commands, DDE, file open
	CCommandLineInfo cmdInfo;
	cmdInfo.m_nShellCommand = CCommandLineInfo::FileNothing;
	ParseCommandLine(cmdInfo);
	
	// Dispatch commands specified on the command line.  Will return FALSE if
	// app was launched with /RegServer, /Register, /Unregserver or /Unregister.
	if (!ProcessShellCommand(cmdInfo))
		return FALSE;
	// The main window has been initialized, so show and update it
	pMainFrame->ShowWindow(SW_SHOW);
	pMainFrame->UpdateWindow();
	return TRUE;
}

// CMeshApp message handlers

//************************************
// copy viewpoint
//************************************
void CMeshApp::OnCameraCopyviewpoint()
{
	CMainFrame *pFrame = (CMainFrame *)AfxGetApp()->m_pMainWnd;
	CMDIChildWnd *pChild = (CMDIChildWnd *)pFrame->GetActiveFrame();
	CMeshView *pView = (CMeshView *)pChild->GetActiveView();
	if(pView == NULL)
		return;
	m_Arcball.Set(pView->arcball());
}
void CMeshApp::OnUpdateCameraCopyviewpoint(CCmdUI *pCmdUI)
{
}

//************************************
// paste viewpoint
//************************************
void CMeshApp::OnCameraPasteviewpoint()
{
	CMainFrame *pFrame = (CMainFrame *)AfxGetApp()->m_pMainWnd;
	CMDIChildWnd *pChild = (CMDIChildWnd *)pFrame->GetActiveFrame();
	CMeshView *pView = (CMeshView *)pChild->GetActiveView();
	if(pView == NULL)
		return;
	pView->arcball()->Set(m_Arcball);
	pView->GetDocument()->UpdateAllViews(NULL);
}
void CMeshApp::OnUpdateCameraPasteviewpoint(CCmdUI *pCmdUI)
{
}
