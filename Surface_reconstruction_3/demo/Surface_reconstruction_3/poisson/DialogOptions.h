#pragma once


// Boîte de dialogue CDialogOptions

class CDialogOptions : public CDialog
{
	DECLARE_DYNAMIC(CDialogOptions)

public:
	CDialogOptions(CWnd* pParent = NULL);   // constructeur standard
	virtual ~CDialogOptions();

// Données de boîte de dialogue
	enum { IDD = IDD_DIALOG_OPTIONS };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // Prise en charge DDX/DDV

	DECLARE_MESSAGE_MAP()

public:
	double m_sm_distance;
	double m_sm_radius;
	double m_sm_angle;
	unsigned int m_dr_max_vertices;
	double m_dr_shell_size;
	double m_dr_sizing;
	double m_contouring_value;
	unsigned int m_number_of_neighbours;
};
