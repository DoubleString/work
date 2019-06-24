/*
 * RtkBamboo.h
 *
 *  Created on: 2019/4/30
 *      Author: jtao
 */

#ifndef INCLUDE_RTK_NEW_RTKBAMBOO_H_
#define INCLUDE_RTK_NEW_RTKBAMBOO_H_

#include "../Algo/SrifFilter.h"
#include "Line_amb.h"
using namespace std;
namespace bamboo{
#define ATOM_CUT 13
class RtkLogInfo{
public:
	RtkLogInfo(){
		m_reset();
	}
	void m_reset(){
		int isys;
		memset(Raw,0,sizeof(Raw));
		memset(useRaw,0,sizeof(useRaw));
		memset(ifab,0,sizeof(ifab));
		memset(pifab,0,sizeof(pifab));
		memset(status,0,sizeof(status));
		memset(elev,0,sizeof(elev));
		memset(ion_wl,0,sizeof(ion_wl));

		memset(threRaw,0,sizeof(threRaw));
		memset(threIF,0,sizeof(threIF));

		memset(ion,0,sizeof(ion));
		memset(trp,0,sizeof(trp));
		memset(resiIF,0,sizeof(resiIF));
		memset(resiRaw,0,sizeof(resiRaw));

		memset(ion_new,0,sizeof(ion_new));
		memset(trp_new,0,sizeof(trp_new));
		memset(resiIF_new,0,sizeof(resiIF_new));
		memset(resiRaw_new,0,sizeof(resiRaw_new));
		for(isys = 0;isys < MAXSYS;isys++)
			refsat[isys] = -1;
		memset(idif,0,sizeof(idif));
	}
	int Raw[MAXSAT][MAXFREQ],useRaw[MAXSAT][MAXFREQ],ifab[MAXSAT],pifab[MAXSAT],status[MAXSAT],refsat[MAXSYS],idif[MAXSAT];
	double elev[MAXSAT],ion_wl[MAXSAT];

	double threRaw[MAXSAT][MAXFREQ],threIF[MAXSAT];
	double ion[MAXSAT],trp[MAXSAT],resiIF[MAXSAT],resiRaw[MAXSAT][MAXFREQ];
	double ion_new[MAXSAT],trp_new[MAXSAT],resiIF_new[MAXSAT],resiRaw_new[MAXSAT][MAXFREQ];
};
class RtkBamboo : public SrifFilter{
public:
	RtkBamboo(vector<Station>& stalist,double dintv);
	~RtkBamboo();
	virtual void v_onFilterInit();
	virtual void v_onFilterEstimate();
	virtual void v_onAddObsMatrix();
	virtual void v_onPrepare();
	virtual void v_onFinish();
	virtual void v_onAddAmbt();
	virtual void v_onTimeUpdate();
	virtual void v_onOutput();
	Line_amb m_atom;
	Line_net m_line;
protected:
	int m_fixamb_raw(Line_net& line,Line_amb& m_atom);
	int m_fixamb_if(Line_net& line,Line_amb& m_atom);
	void m_addMWCon(double* mat,int ndim,Line_net& line,Line_amb& m_atom);
	void m_addMWCon_(double* mat,int ndim,Line_net& line,Line_amb& m_atom);
	void m_smoothMW(Line_amb& m_atom);
	void m_smoothMWIon(Line_net& line,Line_amb& m_atom);
	void m_holdamb(Line_net& line,Line_amb& m_atom,double* xa);
	bool m_validAmb(bool,Line_net& line,Line_amb& atom,double* xa,int isat,int ifreq,int abfx,int& oAbfx,int threshold,double& res);
	bool m_validMW(Line_amb& atom,int isat,int refsat);
	void m_calAtom(Line_net& line,Line_amb& m_atom);
	void m_calAtom_pam(Line_net& line,Line_amb& atom);
	void m_getAtomStatus(Line_net& line,int psat,int refsat,double& ddion,double& ddtrp);
	void m_calAtom_bak(Line_net& line,Line_amb& m_atom);
	void m_calAtomWL(Line_net& line,Line_amb& atom);
	int m_validIFAmb(bool ledit,Line_net& line,Line_amb& amb,double* xf,int isat,int& oN1,int& oN2,int& subON1,int& subON2,double& res);
	void m_genRMatrix(Line_net& line);
	void m_fillAmat(Line_net& line);
	bool m_getRange(Line_net& iline,int isat,int ifreq,double& range1,double& range0,double& var1,double& var0);
	bool m_getPhase(Line_net& iline,int isat,int ifreq,double& phase1,double& phase0,double& var1,double& var0);
	void m_pkRefsat(Line_net& line,Line_amb& atom);
	void m_onDiaGenSen(int iobs, int pos, double wgt);
	void m_onDia();
	void m_onDiaUpInfs(int nb,int* iset);
	void m_detectOmc(Line_net& line,int isys,int flag[MAXSAT]);
	int m_getRefByLinearCmb(Line_net& line,Line_amb& atom,int* idx, double* elevn, int nx);
	time_t m_lastAct;
	RtkLogInfo loginfo;
	double m_bslength,m_tt;
	int m_nfix,m_total;
};
}
#endif /* INCLUDE_RTK_NEW_RTKBAMBOO_H_ */
