#pragma once

#include "resource.h"


// Boîte de dialogue CDialogOptions_4

class CDialogOptions_4 : public CDialog
{
	DECLARE_DYNAMIC(CDialogOptions_4)

public:
	CDialogOptions_4(CWnd* pParent = NULL);   // constructeur standard
	virtual ~CDialogOptions_4();

// Données de boîte de dialogue
	enum { IDD = IDD_DIALOG_OPTIONS_4 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // Prise en charge DDX/DDV

	DECLARE_MESSAGE_MAP()

public:
	double m_point_size; // OpenGL point size
};
