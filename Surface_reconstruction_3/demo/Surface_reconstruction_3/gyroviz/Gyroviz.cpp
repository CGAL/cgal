// Gyroviz.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "Gyroviz.h"
#include "MainFrm.h"

#include "ChildFrm.h"
#include "GyrovizDoc_4.h"
#include "GyrovizView_4.h"
#include "GyrovizDoc_3.h"
#include "GyrovizView_3.h"
#include "GyrovizDoc_2.h"
#include "GyrovizView_2.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CGyrovizApp

BEGIN_MESSAGE_MAP(CGyrovizApp, CWinApp)
	ON_COMMAND(ID_APP_ABOUT, OnAppAbout)
	// Standard file based document commands
	ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
	ON_COMMAND(ID_FILE_OPEN, CWinApp::OnFileOpen)
	// Standard print setup command
	ON_COMMAND(ID_FILE_PRINT_SETUP, CWinApp::OnFilePrintSetup)
END_MESSAGE_MAP()


// CGyrovizApp construction

CGyrovizApp::CGyrovizApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}


// The one and only CGyrovizApp object

CGyrovizApp theApp;

// CGyrovizApp initialization

BOOL CGyrovizApp::InitInstance()
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
	CMultiDocTemplate* pDocTemplate3;
	pDocTemplate3 = new CMultiDocTemplate(IDR_GyrovizTYPE_3,
		RUNTIME_CLASS(CGyrovizDoc_3),
		RUNTIME_CLASS(CChildFrame), // custom MDI child frame
		RUNTIME_CLASS(CGyrovizView_3));
	if (!pDocTemplate3)
		return FALSE;
	AddDocTemplate(pDocTemplate3);
	//	
	CMultiDocTemplate* pDocTemplate2;
	pDocTemplate2 = new CMultiDocTemplate(IDR_GyrovizTYPE_2,
		RUNTIME_CLASS(CGyrovizDoc_2),
		RUNTIME_CLASS(CChildFrame), // custom MDI child frame
		RUNTIME_CLASS(CGyrovizView_2));
	if (!pDocTemplate2)
		return FALSE;
	AddDocTemplate(pDocTemplate2);
	//
	CMultiDocTemplate* pDocTemplate4;
	pDocTemplate4 = new CMultiDocTemplate(IDR_GyrovizTYPE_4,
		RUNTIME_CLASS(CGyrovizDoc_4),
		RUNTIME_CLASS(CChildFrame), // custom MDI child frame
		RUNTIME_CLASS(CGyrovizView_4));
	if (!pDocTemplate4)
		return FALSE;
	AddDocTemplate(pDocTemplate4);
	
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
	pMainFrame->ShowWindow(m_nCmdShow);
	pMainFrame->UpdateWindow();
	return TRUE;
}



// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()

// App command to run the dialog
void CGyrovizApp::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
}


// CGyrovizApp message handlers

