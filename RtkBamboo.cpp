/*
 * Rtkbamboo.cpp
 *
 *  Created on: 2019/4/30
 *      Author: jtao
 */

/// create by jtao on 2019/4/30,support the RAW & IF & SF mode
/// ambiguity fixed support only raw now,support hold,instantaneous,continuous mode
/// check with MW observation
/// support GLONASS rtk here
/// support GLONASS ambiguities fixed here
#include "../../include/Rtk/RtkBamboo.h"
#include "../../include/Model/Model.h"
#include "../../include/Com/CMat.h"
#include "../../include/Com/Com.h"
#include "../../include/Tropo/Tropo_vmf1.h"
#include "../../include/Bamboo/Bamboo.h"
#include "../../include/Ambfix/Ambfix_rtk.h"
#include "../../include/Model/RtkModelItrs.h"
using namespace std;
using namespace bamboo;
//// version change: maxdel done , iteration fix part done , atom calculation,holding function done,ddion smoother
//// find whether nfix > 4,or think the fixed ambiguities is not right
//// do not use holding function to hold ambiguities
RtkBamboo::RtkBamboo(vector<Station>& stalist,double dintv):SrifFilter(stalist,dintv){
	int i;
	double dx[3];
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	vector<Station>::iterator sItr;
	if(m_Sta.size() != 2){
		LOGPRINT("station number = %d,rtk mode should input two stations!",m_Sta.size());
		exit(1);
	}
	for(sItr = this->m_Sta.begin();sItr != this->m_Sta.end();++sItr){
		if((*sItr).gModel != NULL)
			delete (*sItr).gModel;
		(*sItr).gModel = (Model*)new RtkModelItrs(this,&(*sItr));
	}
	m_line.isit = 0;
	m_line.refsit = 1;
	m_atom.isit = 0;
	m_atom.refsit = 1;
	for(i = 0; i < 3;i++){
		dx[i] = m_Sta[0].x[i] - m_Sta[1].x[i];
	}
	m_bslength = sqrt(dot(3,dx,dx));
	m_tt = dintv;
	m_lastAct = 0;
	m_nfix = m_total = 0;
}
RtkBamboo::~RtkBamboo(){
	vector<Station>::iterator sItr;
	for(sItr = this->m_Sta.begin();sItr != this->m_Sta.end();++sItr){
		if((*sItr).gModel != NULL)
			delete (*sItr).gModel;
	}
}
void RtkBamboo::v_onFilterInit(){
	// Get The Filter Configure
	int i,amnum,ifreq,ip,isit,nsit,pmnum,isys,jsys,rClkMind,isat,nsatest,nsitest,nion,nf,nxyz,ntrop;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	const char *xyz = "XYZ";
	char buf[256];
	double minut,fact;
	vector<Station>::iterator sItr;
	vector<Satellite>::iterator pSatItr;
	SrifPamt pM;
	nf = 1;
	nion = 0;
	nsit = this->m_Sta.size();
	nsatest = 0;
	nsitest = 0;
	for(isit = 0,nsitest = 0,nxyz = 0,ntrop = 0;isit < nsit;isit++){
		if(m_Sta[isit].skd[0] != 'F' && m_Sta[isit].skd[0] != 'R')
			nxyz = nxyz + 1;       //
		if(!strstr(dly->ztdmod,"NONE"))
			ntrop = ntrop + 1;     //
	}
	if(dly->promode == PRO_RTK)
		nsit = nsit - 1;
	if(strstr(dly->ionmod,"ION_EST") && !strstr(dly->cobs,"IF"))
		nion = nsit;  //
	if(strstr(dly->cobs,"RAW"))
		nf = 2;                   //
	LOGSTAT(true,"nsatclk(%d),nsitclk(%d),nxyz(%d),ntrop(%d),nion(%d),nf(%d)",nsatest,nsitest,nxyz,ntrop,nion,nf);
	/****************************** ALLOCATE MEMROY **************************************/
	this->nmtxX = nsatest + ntrop + nion * dly->nprn + nsitest * dly->nsys + nxyz * 3 + 3 * (dly->nsys - 1) + dly->nprn * nsit * nf * 2 + 5; // satclk ztd recclk xyz amb * 2
	this->nmtxY = nsatest + ntrop + nion * dly->nprn + nsitest * dly->nsys + nxyz * 3 + 3 * (dly->nsys - 1) + dly->nprn * nsit * nf * 2 +    // amb * 2 for INSTANTANEOUS mode
			dly->nprn * nsit * 2 * nf + 2;
	this->infs = (double*) calloc(this->nmtxX * this->nmtxY,sizeof(double));
	if(this->infs == NULL){
		LOGPRINT("row = %d,ncol = %d,can't allocate memory for Infs Filter!",nmtxY,nmtxX);
		LOGE("row = %d,ncol = %d,can't allocate memory for Infs Filter!",nmtxY,nmtxX);
		exit(1);
	}
	this->iptx = (int*) calloc(this->nmtxX,sizeof(int));
	if(this->iptx == NULL){
		LOGPRINT("nrow = %d,can't allocate memory for iptx Filter!",nmtxX);
		LOGE("nrow = %d,can't allocate memory for iptx Filter!",nmtxX);
		exit(1);
	}
	for(i = 0;i < this->nmtxX;i++){
		this->iptx[i] = this->nmtxY * i;
	}
	this->nptk = dly->nprn * 2 * nf * nsit;
	this->Hk = (double*)calloc(this->nmtxX * this->nptk,sizeof(double));
	if(this->Hk == NULL){
		LOGPRINT("row = %d,ncol = %d,can't allocate memory for Hk Filter!",nptk,nmtxX);
		LOGE("row = %d,ncol = %d,can't allocate memory for Hk Filter!",nptk,nmtxX);
		exit(1);
	}
	this->iphk = (int*)calloc(this->nmtxX,sizeof(int));
	if(this->iphk == NULL){
		LOGPRINT("nrow = %d,can't allocate memory for iphk Filter!",nmtxX);
		LOGE("nrow = %d,can't allocate memory for iphk Filter!",nmtxX);
		exit(1);
	}
	for(i = 0; i < this->nmtxX;i++){
		this->iphk[i] = this->nptk * i;
	}
	LOGSTAT(true,"infs(%d,%d),Hk(%d,%d)",nmtxY,nmtxX,nptk,nmtxX);
	// Allocate Memroy for other Variables
	this->weig = (double*) calloc(dly->nprn * nf * 2 * nsit,sizeof(double));
	this->resi = (double*) calloc(dly->nprn * nf * 2 * nsit,sizeof(double));
	this->ipob = (int(*)[4])calloc(dly->nprn * nf * 2 * 4 * nsit,sizeof(int));
	this->pkresi = (double*)calloc(dly->nprn * nf * 2 * nsit,sizeof(double));
	amnum = nsit * MAXSAT * 2; // amb * 2 for INSTANTANEOUS mode
	pmnum = nsatest + ntrop + nion * dly->nprn + nsitest * dly->nsys + nxyz * 3 + 3 * (dly->nsys - 1); // satclk ztd ion clk xyz
	LOGSTAT(true,"amnum(%d),pmnum(%d)",amnum,pmnum);
	this->amtag = (int(*)[MAXFREQ])calloc(amnum * MAXFREQ,sizeof(int)); // for find over delay convenience
	this->pmtag = (int*)calloc(pmnum,sizeof(int));
	this->ptrAmtag = (SrifPamt* (*)[MAXFREQ])calloc(amnum * MAXFREQ,sizeof(SrifPamt*));
	for(i = 0;i < amnum;i++){
		for(ifreq = 0;ifreq < MAXFREQ;ifreq++){
			this->amtag[i][ifreq] = -1;
			this->ptrAmtag[i][ifreq] = NULL;
		}
	}
	/***************** Initial Filter ************************/
	ip = 0;
	this->nc = this->np = this->nps = this->ns = this->imtx = 0;
	// ztd
	if(!strstr(dly->ztdmod,"NONE")){
		lestZtd = true;
		for(sItr = m_Sta.begin(),isit = 0;sItr != this->m_Sta.end();++isit,++sItr){
			pM.iobs = 0;
			pM.iattach = -1;
			pM.ind = ip;
			sprintf(pM.pname, "%s_%03d", dly->ztdmod, isit);
			if (strstr(dly->ztdmod, ":")) {
				sscanf(strstr(dly->ztdmod, ":") + 1, "%lf", &minut);
				fact = sqrt(minut / 60.0);
			} else {
				fact = 1.0;
			}
			pM.rw = 1.0 / ((*sItr).qztd * fact);
			pM.zw = 0.0;
			pM.map = 1.0;

			pM.xini = 0.0;
			pM.xcor = 0.0;
			pM.xest = 0.0;
			pM.xsig = 0.0;

			pM.psat = -1;
			pM.isit = isit;

			this->pmtag[pM.ind] = ip;
			this->infs[this->iptx[ip] + ip] = 1.0 / ((*sItr).dztd0);
			this->m_Pm.push_back(pM);
			ip++;
		}
	}
	// ion
	if(strstr(dly->ionmod,"ION_EST")){
		for(sItr = m_Sta.begin(),isit = 0;sItr != this->m_Sta.end();++isit,++sItr){
			if((*sItr).skd[0] == 'R')
				continue;
			for (isat = 0; isat < dly->nprn; isat++) {
				isys = index_string(SYS, dly->cprn[isat][0]);
				pM.iobs = 0;
				pM.ind = ip;
				pM.iattach = -1;
				sprintf(pM.pname, "ION_%03d_%03d", isit, isat);
				pM.map = 1.0;
				pM.rw = 1.0 / (*sItr).qion;   // should be changed after
				pM.zw = 0.0;
				pM.xini = 1e-6;
				pM.xcor = 0.0;
				pM.xest = 0.0;
				pM.xsig = 0.0;
				pM.isit = isit;
				pM.psat = isat;
				this->pmtag[pM.ind] = ip;
				this->infs[this->iptx[ip] + ip] = 1.0
						/ ((*sItr).dion0 * m_bslength / 1E4);
				this->m_Pm.push_back(pM);
				ip++;
			}
		}
	}
	/// station coordinate
	for (sItr = this->m_Sta.begin(), isit = 0; sItr != this->m_Sta.end();
			++isit, ++sItr) {
		if ((*sItr).skd[0] == 'K' || (*sItr).skd[0] == 'D') {
			lestPos = true;
			for (i = 0; i < 3; i++) {
				pM.iobs = 0;
				pM.ind = ip;
				pM.iattach = -1;
				sprintf(pM.pname, "STAP_%03d_%c", isit, xyz[i]);
				pM.rw = 1.0 / ((*sItr).qx[i] * dly->scal_qx * sqrt(intv));
				pM.zw = 0.0;
				pM.map = 0.0;

				pM.xini = (*sItr).x[i];
				pM.xcor = 0.0;
				pM.xest = 0.0;
				pM.xsig = 0.0;

				pM.psat = -1;
				pM.isit = isit;

				this->pmtag[pM.ind] = ip;
				this->infs[this->iptx[ip] + ip] = 1.0
						/ ((*sItr).dx0[i] * dly->scal_dx);
				this->m_Pm.push_back(pM);
				ip++;
			}
		}
	}
	/*********************** Static Position Parameters *****************************/
	ns = 0;
	for (sItr = this->m_Sta.begin(), isit = 0; sItr != this->m_Sta.end();
			++sItr, ++isit) {
		if ((*sItr).skd[0] == 'S') {
			lestPos = true;
			for (i = 0; i < 3; i++) {
				pM.iobs = 0;
				pM.ind = ip + ns;
				pM.iattach = -1;
				sprintf(pM.pname, "STAP_%03d_%c", isit, xyz[i]);

				pM.rw = 1.0 / ((*sItr).qx[i] * dly->scal_qx);
				pM.zw = 0.0;
				pM.map = 1.0;

				pM.psat = -1;
				pM.isit = isit;

				pM.xini = (*sItr).x[i];
				pM.xcor = 0.0;
				pM.xest = 0.0;
				pM.xsig = 0.0;

				this->pmtag[pM.ind] = ip + ns;
				this->infs[this->iptx[ip + ns] + ip + ns] = 1.0
						/ ((*sItr).dx0[i] * dly->scal_dx);
				this->m_Pm.push_back(pM);
				ns++;
			}
		}
	}
	/// static pco parameter
//	double pest[] = {-0.001  ,     -0.021,       -0.015};
//	double pcov[] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}; //{0.10,0.10,0.10,0.0001,0.0001,0.0001};
	double pest[] = {0.00,0.00,0.00};
	double pcov[] = {0.0201,0.0201,0.0501,0.0001,0.0001,0.0001}; //{0.10,0.10,0.10,0.0001,0.0001,0.0001};
	for(sItr = this->m_Sta.begin(),isit = 0;sItr != this->m_Sta.end() && dly->lpcoest;++sItr,++isit){
		if((*sItr).skd[0] == 'R')
			continue;
		for(jsys = 0;jsys < dly->nsys;jsys++){
			isys = index_string(SYS,dly->system[jsys]);
			if(isys == dly->iref)
				continue;
			for(i = 0; i < 3;i++){
				pM.iobs = 0;
				pM.ind = ip + ns;
				pM.iattach = -1;
				sprintf(pM.pname,"PCO_%04d_%c_%c",isit,SYS[isys],xyz[i]);

				pM.rw = 1 / pcov[3 + i];
				pM.zw = 0.0;
				pM.map = 1.0;

				pM.psat = -1;
				pM.isit = isit;

				pM.xini = 0.0;
				pM.xcor = 0.0;
				pM.xest = pest[i];
				pM.xsig = 0.0;


				this->pmtag[pM.ind] = ip + ns;
				this->infs[this->iptx[ip + ns] + ip + ns] = 1.0 / pcov[i];
				this->m_Pm.push_back(pM);
				ns++;
			}
		}
	}
	this->np = ip;
	this->ns = ns;
	this->nc = nc;
	this->nps = this->np + this->ns;
	this->imtx = this->nps + this->nc;
	this->isInit = true;
	for(i = 0; i < nps;i++){
		infs[iptx[imtx] + i] = m_Pm[i].xest * infs[iptx[i] + i];
	}
	LOGSTAT(true,"index                  param-name   isat isit      std      process-std      value ");
	for(ip = 0;ip < imtx;ip++){
		LOGSTAT(true,"%03d %30s %4d %4d %12.4lf %12.4lf %12.4lf",ip,m_Pm[ip].pname,m_Pm[ip].psat,m_Pm[ip].isit,
				1.0 / infs[iptx[ip] + ip],1.0 / m_Pm[ip].rw,infs[iptx[imtx] + ip]);
	}
	LOGSTAT(true,"np(%d),ns(%d),nps(%d),nc(%d),imtx(%d)",np,ns,nps,nc,imtx);
}
void RtkBamboo::v_onPrepare(){
	time_t now_t;
	int isat,isit,isys,ifreq;
	vector<Satellite>::iterator satItr;
	vector<Station>::iterator staItr;
	list<SrifPamt*>::iterator ambItr;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	for (satItr = dly->SAT.begin(), isat = 0; satItr != dly->SAT.end();
			++satItr, ++isat) {
		if ((*satItr).isabled == false) {
			LOGE("cprn = %s,satellite is unavailable",(*satItr).cprn);
			for (staItr = m_Sta.begin(); staItr != m_Sta.end(); ++staItr) {
				memset((*staItr).ob.obs[isat], 0, sizeof(double) * 2 * MAXFREQ);
			}
		}
	}
	memset(m_sclock, 0, sizeof(double) * MAXSAT);
	for (isat = 0, satItr = dly->SAT.begin();
			satItr != dly->SAT.end() && clkAdapter != NULL; ++satItr, ++isat) {
		if (!(clkAdapter->v_lRead(isat))) {
			if (!clkAdapter->v_readClk((*satItr).cprn, this->mjd, this->sod,
					m_sclock + isat)){
				LOGPRINT("cprn = %s,no ephemeris for satellite clock",(*satItr).cprn);
			}
			clkAdapter->v_sRead(isat, 1);
		}
	}
	if (!strncmp(dly->ambmod, "INSTANTANEOUS", strlen(dly->ambmod))) {
		for(isit = 0;isit < m_Sta.size();isit++){
			for (isat = 0; isat < dly->nprn; isat++) {
				isys = index_string(SYS, dly->cprn[isat][0]);
				for (ifreq = 0; ifreq < dly->nfq[isys]; ifreq++) {
					if (m_Sta[isit].ob.obs[isat][ifreq] != 0) {
						m_Sta[isit].ob.flag[isat][ifreq] = 1;
					}
				}
			}
		}
	}
	/// add for rtk here
	m_line.m_reset();
	now_t = mjd2time(mjd,sod);
	m_tt = m_lastAct == 0 ? intv :now_t - m_lastAct;
	m_lastAct = now_t;
}
void RtkBamboo::v_onFinish(){
	int isat,isit,isys,ifreq,i;
	vector<Station>::iterator sItr;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	for(isat = 0; isat < dly->nprn && clkAdapter != NULL;isat++){
		if(this->clkAdapter->v_lRead(isat)){
			this->clkAdapter->v_sRead(isat,0);
		}
	}
	for (isat = 0; isat < dly->nprn; isat++) {
		isys = index_string(SYS, dly->cprn[isat][0]);
		if (isys == -1)
			continue;
		for (ifreq = 0; ifreq < dly->nfq[isys]; ifreq++) {
			if (m_Sta[0].ob.omc[isat][ifreq] == 0
					|| m_Sta[1].ob.omc[isat][ifreq] == 0)
				continue;
			this->ptrAmtag[isat][ifreq]->ptime[1] = mjd + sod / 86400.0;
		}
	}
	/******************* Reset Ambiguity Parameters Including Reference Satellite***************/
	for (isat = 0; isat < dly->nprn; isat++) {
		isys = index_string(SYS, dly->cprn[isat][0]);
		if (isys == -1)
			continue;
		for (ifreq = 0; ifreq < dly->nfq[isys]; ifreq++) {
			if (m_Sta[0].ob.omc[isat][ifreq] == 0
					|| m_Sta[1].ob.omc[isat][ifreq] == 0)
				continue;
			if (m_Sta[0].ob.flag[isat][ifreq] != 0)
				m_Sta[0].ob.flag[isat][ifreq] = 0;
			if (m_Sta[1].ob.flag[isat][ifreq] != 0)
				m_Sta[1].ob.flag[isat][ifreq] = 0;
		}
	}
}
void RtkBamboo::v_onAddAmbt(){
	int nbias = 0, isat, ifreq, isys, i, j;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	SrifPamt* aM;
	list<SrifPamt*>::iterator aItr;
	if (strstr(dly->uobs, "PHASE") || strstr(dly->uobs, "GRAPHIC")) {
		for (isat = 0; isat < dly->nprn; isat++) {
			isys = index_string(SYS, dly->cprn[isat][0]);
			for (ifreq = 0; ifreq < dly->nfq[isys]; ifreq++) {
				if (m_Sta[0].ob.omc[isat][ifreq] == 0
						|| m_Sta[1].ob.omc[isat][ifreq] == 0)
					continue;
				if (m_Sta[0].ob.flag[isat][ifreq] != 0
						|| m_Sta[1].ob.flag[isat][ifreq] != 0) {
					aM = new SrifPamt;
					sprintf(aM->pname, "AMBL_%02d", ifreq);
					sprintf(aM->staname, "%s-%s", m_Sta[0].name.c_str(),
							m_Sta[1].name.c_str());
					aM->xini = 0;
					aM->xcor = 0;
					aM->xest = 0;
					aM->xsig = 0;

					aM->map = 1.0;
					aM->rw = 1.0 / 0.001;
					aM->zw = 0.0;

					aM->ifreq = ifreq;
					aM->psat = isat;

					aM->isit = 0;
					aM->elev = 0.0;
					aM->swPtr = NULL;

					aM->ipob = -1;
					aM->iobs = 0;
					aM->ptime[0] = aM->ptime[1] =  this->mjd + this->sod / 86400.0;

					aM->ifab = AMB_INIT;
					aM->xrwl = 0;
					aM->weig = 0;
					aM->xswl = 0;

					this->m_Am.push_back(aM);
					this->ptrAmtag[isat][ifreq] = aM;
					this->amtag[isat][ifreq] = this->imtx + nbias;
					LOGI(LEVEL_ALGO,"mjd = %d,sod = %lf,sta = %s,cprn = %s,ifreq = %d,adding ambiguity",mjd,sod,m_Sta[0].name.c_str(),dly->cprn[isat].c_str(),ifreq);
					nbias++;
				}

			}
		}
	}
	/*move the X to the corresponding column*/
	if (nbias != 0) {
		for (i = 0; i < imtx; i++) {
			infs[iptx[imtx + nbias] + i] = infs[iptx[imtx] + i];
		}
		// Reset the memory
		for (i = imtx; i < imtx + nbias; i++) {
			for (j = 0; j < imtx + nbias + 1; j++) {
				infs[iptx[i] + j] = 0.0;
				infs[iptx[j] + i] = 0.0;
			}
		}
		for (i = 0; i < nbias; i++) {
			infs[iptx[imtx + i] + imtx + i] = 1e-4;    // set the diag
		}
		imtx = imtx + nbias;
		nc = nc + nbias;
	}
	if (imtx + 1 > nmtxX) {
		cout << "***ERROR(onAddAmbt):Matrix nmtxX is too small! nmtxX = "
				<< nmtxX << "imtx + 1 = " << imtx + 1 << endl;
		exit(1);
	}
}
void RtkBamboo::v_onTimeUpdate(){
	int ind, indext, i, ip, ipext, j, isat, isit, isys, ifreq, lmvup = false;
	double beta, sum, minut,rw,xsig,fact;
	int infnp = np + np + 1;
	vector<Station>::iterator sItr;
	vector<SrifPamt> oldPam;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	if(infs == NULL)  return;// in case of none implement,NULL pointer will make it crash
	double* infs_u = (double*) calloc(infnp * (np + imtx + 1), sizeof(double));
	for (i = 0; i < imtx + 1; i++) {
		if (i < np) {
			memcpy(&infs_u[infnp * i], &infs[iptx[i]], sizeof(double) * np);
		} else {
			memcpy(&infs_u[infnp * (i + np)], &infs[iptx[i]],
					sizeof(double) * np);
		}
	}
	for (ind = 0; ind < np; ++ind) {
		lmvup = false;
		ip = this->pmtag[ind];
		rw = m_Pm[ip].rw;
		oldPam.push_back(m_Pm[ip]);
		if (!strncmp(m_Pm[ip].pname, "STAP", 4)) {
			m_Pm[ip].zw = -rw * m_Pm[ip].map * m_Pm[ip].xcor;
			m_Pm[ip].xini = m_Pm[ip].xini + m_Pm[ip].xcor;
		}
		if (!strncmp(m_Pm[ip].pname, "SATCLK", 6)) {
			m_Pm[ip].zw = -rw * m_Pm[ip].map * m_Pm[ip].xcor;
			m_Pm[ip].xini = m_Pm[ip].xini + m_Pm[ip].xcor;
		}
		if (!strncmp(m_Pm[ip].pname, "RECCLK", 6)) {
			m_Pm[ip].zw = 0.0;
		}
		if (strstr(m_Pm[ip].pname, "ZTD")) {
			if (mjd * 86400.0 + sod - m_Pm[ip].ptime[1] * 86400.0 < -1E-2) {
				lmvup = true;
				// Update Xini
			} else {
				if (!strstr(m_Pm[ip].pname, ":")) {
					rw = rw * 1.0 / sqrt(m_tt); /// update every epoch
					m_Pm[ip].zw = -rw * m_Pm[ip].map * m_Pm[ip].xcor;
					m_Pm[ip].xini = m_Pm[ip].xini + m_Pm[ip].xcor;

					m_Pm[ip].ptime[0] = this->mjd + this->sod / 86400.0;
					m_Pm[ip].ptime[1] = this->mjd + this->sod / 86400.0;
				} else {
					m_Pm[ip].zw = -rw * m_Pm[ip].map * m_Pm[ip].xcor; // update every timeval
					m_Pm[ip].xini = m_Pm[ip].xini + m_Pm[ip].xcor;

					i = index_string(m_Pm[ip].pname, ':');
					minut = atof(m_Pm[ip].pname + i + 1);
					m_Pm[ip].ptime[0] = this->mjd + this->sod / 86400.0;
					m_Pm[ip].ptime[1] = this->mjd + this->sod / 86400.0
							+ minut / 1440.0;
				}
			}
		}
		if(strstr(m_Pm[ip].pname,"ION")){
			m_Pm[ip].zw = 0.0;
			i = index_string(m_Pm[ip].pname,'_');
			isit = atoi(m_Pm[ip].pname + i + 1);
			j = index_string(m_Pm[ip].pname + i + 1,'_');
			isat = atoi(m_Pm[ip].pname + i + 1 + j + 1);
			fact = 1.0;
			if(m_Sta[isit].ob.elev[isat] > 0.0)
				fact = cos(m_Sta[isit].ob.elev[isat]);
			rw = rw * 1.0 / (m_bslength / 1E4 * fact * sqrt(m_tt));
			m_Pm[ip].zw = -rw * m_Pm[ip].map * m_Pm[ip].xcor;
			m_Pm[ip].xini = m_Pm[ip].xini + m_Pm[ip].xcor;
			if((mjd - m_Pm[ip].ptime[1]) * 86400.0 + sod > 1200.0)
				m_Pm[ip].xini = 0.0;
		}
		m_Pm[ip].iobs = 0;
		if (lmvup) {
			infs_u[infnp * ind + np + ind] = -1E6;
			infs_u[infnp * (np + ind) + np + ind] = 1E6;
			infs_u[infnp * (imtx + np) + np + ind] = 0.0;
			continue;
		}
		infs_u[infnp * ind + np + ind] = -rw * m_Pm[ip].map;
		infs_u[infnp * (np + ind) + np + ind] = rw;
		infs_u[infnp * (imtx + np) + np + ind] = m_Pm[ip].zw;
	}
	CMat::CMat_Householder(infs_u, infnp, np + imtx + 1, np + np, np + imtx + 1,
			false);
	for (i = 0; i < imtx + 1; i++) {
		memcpy(&infs[iptx[i]], &infs_u[infnp * (np + i) + np],
				sizeof(double) * np);
	}
	/////////////////////////////// update the smSig here /////////////////////////
	if(nobs != 0)
		this->smAdt.m_smooth(esig,this->smSig,xsig);
	free(infs_u);
}
bool RtkBamboo::m_getRange(Line_net& iline,int isat,int ifreq,double& range1,double& range0,double& var1,double& var0){
	int isit,refsit;
	var1 = var0 = range1 = range0 = 0.0;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	if(!strstr(dly->uobs,"CODE"))
		return true;
	isit = iline.isit;
	refsit = iline.refsit;
	if(m_Sta[isit].ob.omc[isat][MAXFREQ + ifreq] == 0 || m_Sta[refsit].ob.omc[isat][MAXFREQ + ifreq] == 0)
		return false;
	if(strstr(dly->cobs,"IF")){
		range1 = m_Sta[isit].ob.omc[isat][MAXFREQ] * dly->SAT[isat].fac[0] -
				m_Sta[isit].ob.omc[isat][MAXFREQ + 1] * dly->SAT[isat].fac[1];
		range0 = m_Sta[refsit].ob.omc[isat][MAXFREQ] * dly->SAT[isat].fac[0] -
				m_Sta[refsit].ob.omc[isat][MAXFREQ + 1] * dly->SAT[isat].fac[1];
		var1 = pow(dly->SAT[isat].fac[0],2) * m_Sta[isit].ob.var[isat][MAXFREQ] +
				pow(dly->SAT[isat].fac[1],2) * m_Sta[isit].ob.var[isat][MAXFREQ + 1];
		var0 = pow(dly->SAT[isat].fac[0],2) * m_Sta[refsit].ob.var[isat][MAXFREQ] +
				pow(dly->SAT[isat].fac[1],2) * m_Sta[refsit].ob.var[isat][MAXFREQ + 1];
	}else if(strstr(dly->cobs,"SF")){
		range1 = m_Sta[isit].ob.omc[isat][MAXFREQ];
		range0 = m_Sta[refsit].ob.omc[isat][MAXFREQ];

		var1 = m_Sta[isit].ob.var[isat][MAXFREQ];
		var0 = m_Sta[refsit].ob.var[isat][MAXFREQ];
	}else if(strstr(dly->cobs,"RAW")){
		range1 = m_Sta[isit].ob.omc[isat][MAXFREQ + ifreq];
		range0 = m_Sta[refsit].ob.omc[isat][MAXFREQ + ifreq];

		var1 = m_Sta[isit].ob.var[isat][MAXFREQ + ifreq];
		var0 = m_Sta[refsit].ob.var[isat][MAXFREQ + ifreq];
	}
	return true;
}
bool RtkBamboo::m_getPhase(Line_net& iline,int isat,int ifreq,double& phase1,double& phase0,double& var1,double& var0){
	int isit,refsit;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	phase1 = phase0 = 0.0;
	if(!strstr(dly->uobs,"PHASE") && !strstr(dly->uobs,"GRAPHIC"))
		return true;
	isit = iline.isit;
	refsit = iline.refsit;
	if(m_Sta[isit].ob.omc[isat][ifreq] == 0 || m_Sta[refsit].ob.omc[isat][ifreq] == 0)
		return false;
	if(strstr(dly->cobs,"IF")){
		phase1 = m_Sta[isit].ob.omc[isat][0] * dly->SAT[isat].fac[0] -
				m_Sta[isit].ob.omc[isat][1] * dly->SAT[isat].fac[1];
		phase0 = m_Sta[refsit].ob.omc[isat][0] * dly->SAT[isat].fac[0] -
				m_Sta[refsit].ob.omc[isat][1] * dly->SAT[isat].fac[1];
		var1 = pow(dly->SAT[isat].fac[0],2) * m_Sta[isit].ob.var[isat][0] +
				pow(dly->SAT[isat].fac[1],2) * m_Sta[isit].ob.var[isat][1];

		var0 = pow(dly->SAT[isat].fac[0],2) * m_Sta[refsit].ob.var[isat][0] +
				pow(dly->SAT[isat].fac[1],2) * m_Sta[refsit].ob.var[isat][1];

	}else if(strstr(dly->cobs,"SF")){
		phase1 = m_Sta[isit].ob.omc[isat][0];
		phase0 = m_Sta[refsit].ob.omc[isat][0];

		var1 = m_Sta[isit].ob.var[isat][0];
		var0 = m_Sta[refsit].ob.var[isat][0];
	}else if(strstr(dly->cobs,"RAW")){
		phase1 = m_Sta[isit].ob.omc[isat][ifreq];
		phase0 = m_Sta[refsit].ob.omc[isat][ifreq];

		var1 = m_Sta[isit].ob.var[isat][ifreq];
		var0 = m_Sta[refsit].ob.var[isat][ifreq];
	}
	return true;
}

void RtkBamboo::m_genRMatrix(Line_net& line){
	int j,jsys,isys,ifreq,isat,isit,refsit,nref,i,jsat;
	double range1_ref,range0_ref,phase1_ref,phase0_ref;
	double range1,range0,phase1,phase0,var_1,var_0;
	double var[MAXSAT*2*MAXFREQ] = {0};
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	isit = line.isit;
	refsit = line.refsit;
	line.nobs = 0;
	for(jsys = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		if(line.refsat[isys] == -1)
			continue;
		for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
			if (!m_getRange(line, line.refsat[isys], ifreq,range1_ref,
					range0_ref,var_1,var_0)
					|| !m_getPhase(line, line.refsat[isys],ifreq,
							phase1_ref, phase0_ref,var_1,var_0))
				continue;
			for (isat = 0; isat < dly->nprn; isat++) {
				if (dly->cprn[isat][0] != SYS[isys]
						|| isat == line.refsat[isys])
					continue;
				if (strstr(dly->uobs, "PHASE")
						|| strstr(dly->uobs, "GRAPHIC")) {
					if (m_getPhase(line, isat, ifreq,phase1, phase0,var_1,var_0)){
						var[line.nobs] = var_1 + var_0;
						line.ipt[line.nobs][0] = 0;  // phase
						line.ipt[line.nobs][1] = isat; //
						line.ipt[line.nobs][2] = ifreq;
						line.ipt[line.nobs][3] = isit;
						line.nobs = line.nobs + 1;
					}
				}
				if (strstr(dly->uobs, "CODE")) {
					if (m_getRange(line, isat, ifreq,range1, range0,var_1,var_0)){
						var[line.nobs] = var_1 + var_0;
						line.ipt[line.nobs][0] = 1; // range
						line.ipt[line.nobs][1] = isat; //
						line.ipt[line.nobs][2] = ifreq;
						line.ipt[line.nobs][3] = isit;
						line.nobs = line.nobs + 1;
					}
				}
			}
		}
	}
	/// TODO:
	for(jsys = 0,nref = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		if(line.refsat[isys] == -1)
			continue;
		for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
			if(m_getPhase(line, line.refsat[isys],ifreq, phase1_ref, phase0_ref,var_1,var_0)){
				var[line.nobs + nref] = var_1 + var_0;

				line.ipt[line.nobs + nref][0] = 0;
				line.ipt[line.nobs + nref][1] = line.refsat[isys];
				line.ipt[line.nobs + nref][2] = ifreq;
				nref = nref + 1;
 			}
			if (m_getRange(line, line.refsat[isys], ifreq,range1_ref, range0_ref,var_1,var_0)){
				var[line.nobs + nref] = var_1 + var_0;

				line.ipt[line.nobs + nref][0] = 1;
				line.ipt[line.nobs + nref][1] = line.refsat[isys];
				line.ipt[line.nobs + nref][2] = ifreq;
				nref = nref + 1;
			}
		}
	}
	if(line.nobs > 0){
		double *B = (double*)calloc((line.nobs + nref) * (line.nobs + nref),sizeof(double));
		double *P = (double*)calloc((line.nobs + nref) * (line.nobs + nref),sizeof(double));
		double* tmp = (double*)calloc(line.nobs * (line.nobs + nref),sizeof(double));
		line.R = (double*)calloc(line.nobs * line.nobs,sizeof(double));
		for(i = 0; i < line.nobs + nref;i++){
			P[(line.nobs + nref) * i + i] = var[i];
		}
		for(i = 0; i < line.nobs + nref; i++){
			isat = line.ipt[i][1];
			B[(line.nobs + nref) * i + i] = 1.0; ///
			if(i >= line.nobs){
				//
				for(j = 0; j < line.nobs + nref;j++){
					jsat = line.ipt[j][1];
					///
					if(line.ipt[j][0] == line.ipt[i][0] && line.ipt[j][2] == line.ipt[i][2] && dly->cprn[isat][0] == dly->cprn[jsat][0]){
						B[(line.nobs + nref) * i + j] = -1; // i col,j rowl
					}
				}
			}
		}
		CMat::CMat_Matmul("NN",line.nobs,line.nobs + nref,line.nobs + nref,1.0,B,line.nobs + nref,P,line.nobs + nref,0.0,tmp,line.nobs); //B*P
		CMat::CMat_Matmul("NT",line.nobs,line.nobs,line.nobs + nref,1.0,tmp,line.nobs,B,line.nobs + nref,0.0,line.R,line.nobs);
		CMat::CMat_Inverse(line.R,line.nobs,line.nobs); //
		///
		CMat::CMat_Cholesky('U',line.R,line.nobs,line.nobs); // cholesky
		///
		for(i = 0; i < line.nobs - 1;i++){
			memset(line.R + line.nobs * i + i + 1,0,sizeof(double) * (line.nobs - i - 1));
		}
		free(B);
		free(P);
		free(tmp);
	}
}
void RtkBamboo::m_fillAmat(Line_net& line){
	/// TODO:
	int isat,isit,refsit,isys,ref,ifreq,ip,curnp,ix;
	char ztdtag_0[256],ztdtag_1[256],pco_tag[256];
	char iontag_1_ref[256],iontag_0_ref[256],iontag_1[256],iontag_0[256];
	char stap_1[256],stap_0[256];
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	isit = line.isit;
	refsit = line.refsit;
	//
	sprintf(ztdtag_1,"%s_%03d",dly->ztdmod,line.isit); ///
	sprintf(ztdtag_0,"%s_%03d",dly->ztdmod,line.refsit);  ///
	///
	sprintf(stap_1,"STAP_%03d",isit);
	sprintf(stap_0,"STAP_%03d",refsit);
	for(isat = 0;isat < dly->nprn;isat++){
		///
		if(m_Sta[isit].ob.omc[isat][MAXFREQ] == 0 || m_Sta[refsit].ob.omc[isat][MAXFREQ] == 0)
			continue;
		isys = index_string(SYS,dly->cprn[isat][0]);
		ref = line.refsat[isys];
		if(ref == -1 || isat == ref)
			continue;
		///
		sprintf(iontag_1_ref,"ION_%03d_%03d",isit,ref);
		sprintf(iontag_1,"ION_%03d_%03d",isit,isat);

		sprintf(pco_tag,"PCO_%04d_%c",isit,SYS[isys]);
		curnp = 0;
		for (ip = 0; ip < m_Pm.size(); ip++) {
			if (strstr(m_Pm[ip].pname, stap_1)) {
				///
				for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
					for(ix = 0;ix < 3;ix++){
						line.d_coef[isat][curnp + ix][ifreq] = m_Sta[isit].dloudx[isat][ix] - m_Sta[isit].dloudx[ref][ix];
						line.d_coef[isat][curnp + ix][MAXFREQ + ifreq] = m_Sta[isit].dloudx[isat][ix] - m_Sta[isit].dloudx[ref][ix];

						line.i_coef[isat][curnp + ix][ifreq] = m_Pm[ip].ind + ix;
						line.i_coef[isat][curnp + ix][MAXFREQ + ifreq] = m_Pm[ip].ind + ix;
					}
				}
				curnp = curnp + 3;
				ip = ip + 2;
			} else if(strstr(m_Pm[ip].pname,pco_tag)){
				for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
					for(ix = 0;ix < 3;ix++){
						line.d_coef[isat][curnp + ix][ifreq] = m_Sta[isit].dloudx[isat][ix] - m_Sta[isit].dloudx[ref][ix];
						line.d_coef[isat][curnp + ix][MAXFREQ + ifreq] = m_Sta[isit].dloudx[isat][ix] - m_Sta[isit].dloudx[ref][ix];

						line.i_coef[isat][curnp + ix][ifreq] = m_Pm[ip].ind + ix;
						line.i_coef[isat][curnp + ix][MAXFREQ + ifreq] = m_Pm[ip].ind + ix;
					}
				}
				curnp = curnp + 3;
				ip = ip + 2;
			}else if (strstr(m_Pm[ip].pname, stap_0)) {
				///
				for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
					for(ix = 0;ix < 3;ix++){
						line.d_coef[isat][curnp + ix][ifreq] = -1 * (m_Sta[refsit].dloudx[isat][ix] - m_Sta[refsit].dloudx[ref][ix]);
						line.d_coef[isat][curnp + ix][MAXFREQ + ifreq] = -1 * (m_Sta[refsit].dloudx[isat][ix] - m_Sta[refsit].dloudx[ref][ix]);

						line.i_coef[isat][curnp + ix][ifreq] = m_Pm[ip].ind + ix;
						line.i_coef[isat][curnp + ix][MAXFREQ + ifreq] = m_Pm[ip].ind + ix;
					}
				}
				curnp = curnp + 3;
				ip = ip + 2;
			} else if (strstr(m_Pm[ip].pname, iontag_1)) {
				///
				for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
					line.d_coef[isat][curnp][ifreq] = -1.0
							* m_Sta[isit].ob.ionpart[isat]
							* pow(dly->SAT[isat].freq[0], 2)
							/ pow(dly->SAT[isat].freq[ifreq], 2);
					line.d_coef[isat][curnp][MAXFREQ + ifreq] = 1.0
							* m_Sta[isit].ob.ionpart[isat]
							* pow(dly->SAT[isat].freq[0], 2)
							/ pow(dly->SAT[isat].freq[ifreq], 2);

					line.i_coef[isat][curnp][ifreq] = m_Pm[ip].ind;
					line.i_coef[isat][curnp][MAXFREQ + ifreq] = m_Pm[ip].ind;
				}
				curnp = curnp + 1;
			} else if (strstr(m_Pm[ip].pname, iontag_1_ref)) {
				///
				for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
					line.d_coef[isat][curnp][ifreq] = -1.0 * -1.0
							* m_Sta[isit].ob.ionpart[ref]
							* pow(dly->SAT[ref].freq[0], 2)
							/ pow(dly->SAT[ref].freq[ifreq], 2);
					line.d_coef[isat][curnp][MAXFREQ + ifreq] = -1.0 * 1.0
							* m_Sta[isit].ob.ionpart[ref]
							* pow(dly->SAT[ref].freq[0], 2)
							/ pow(dly->SAT[ref].freq[ifreq], 2);
					line.i_coef[isat][curnp][ifreq] = m_Pm[ip].ind;
					line.i_coef[isat][curnp][MAXFREQ + ifreq] = m_Pm[ip].ind;
				}
				curnp = curnp + 1;
			}  else if (strstr(m_Pm[ip].pname, ztdtag_1)) {
				for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
					line.d_coef[isat][curnp][ifreq] = m_Sta[isit].ob.zmap[isat] - m_Sta[isit].ob.zmap[ref];
					line.d_coef[isat][curnp][MAXFREQ + ifreq] = m_Sta[isit].ob.zmap[isat] - m_Sta[isit].ob.zmap[ref];

					line.i_coef[isat][curnp][ifreq] = m_Pm[ip].ind;
					line.i_coef[isat][curnp][ifreq] = m_Pm[ip].ind;
				}
				curnp = curnp + 1;
			} else if (strstr(m_Pm[ip].pname, ztdtag_0)) {
				for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
					line.d_coef[isat][curnp][ifreq] = -1 * (m_Sta[refsit].ob.zmap[isat] - m_Sta[refsit].ob.zmap[ref]);
					line.d_coef[isat][curnp][MAXFREQ + ifreq] = -1 * (m_Sta[refsit].ob.zmap[isat] - m_Sta[refsit].ob.zmap[ref]);

					line.i_coef[isat][curnp][ifreq] = m_Pm[ip].ind;
					line.i_coef[isat][curnp][ifreq] = m_Pm[ip].ind;
				}
				curnp = curnp + 1;
			}
		}
		///
		//
		for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
			line.d_coef[isat][curnp][ifreq] = strstr(dly->cobs,"IF") ? 1.0 : dly->SAT[isat].lamda[ifreq];
			line.d_coef[isat][curnp][MAXFREQ + ifreq] = 0.0;
			line.i_coef[isat][curnp][ifreq] = amtag[isit * MAXSAT + isat][ifreq];
			line.i_coef[isat][curnp][MAXFREQ + ifreq] = amtag[isit * MAXSAT + isat][ifreq];
		}
		curnp = curnp + 1;
		//
		for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
			line.d_coef[isat][curnp][ifreq] = strstr(dly->cobs,"IF") ? -1.0 : -dly->SAT[isat].lamda[ifreq]; /// glonass will correct the bias here
											//strstr(dly->cobs,"IF") ? -1.0 : -dly->SAT[ref].lamda[ifreq];
			line.d_coef[isat][curnp][MAXFREQ + ifreq] = 0.0;
			line.i_coef[isat][curnp][ifreq] = amtag[isit * MAXSAT + ref][ifreq];
			line.i_coef[isat][curnp][MAXFREQ + ifreq] = amtag[isit * MAXSAT + ref][ifreq];
		}
		curnp = curnp + 1;
		line.n_np[isat] = curnp;
	}
}
void RtkBamboo::v_onFilterEstimate(){
	int i, ip, lDebug = false,lfix,mjd0,mjd1;
	double sod0,sod1;
	char tm0[256],tm1[256];
	vector<SrifPamt>::iterator pItr;
	list<SrifPamt*>::iterator aItr;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	if(infs == NULL) return;
	CMat::CMat_Householder(infs, nmtxY, imtx, imtx + nobs, imtx + 1, false);
	cout << "imtx: " << imtx << " nobs: " << nobs << " " << imtx + nobs << " - " << nobs << " np: " << np << endl;
	this->esig = 0.0;
	for (i = 0; i < nobs; i++) {
		resi[i] = infs[iptx[imtx] + imtx + i];
		esig += resi[i] * resi[i];
		infs[iptx[imtx] + imtx + i] = 0.0;
	}
	if (nobs != 0) {
		esig = sqrt(esig / nobs);
		this->m_onGetResi();
	}
	// Solve Linear System Equation
	///////////////////////////////////////////////////////////
	double* est = (double*) calloc(imtx, sizeof(double));
	memcpy(est, &infs[iptx[imtx]], sizeof(double) * imtx);
	CMat::CMat_SolveLinear(imtx, infs, nmtxY, est, imtx);
	f_x.clear();
	//LOGSTAT(true,"sod = %lf",sod);
	for (i = 0, pItr = m_Pm.begin(), aItr = m_Am.begin(); i < imtx; i++) {
		if (i < nps) {
			ip = this->pmtag[i];
			(*pItr).xsig = 0.0;
			(*pItr).xcor = est[i];
			(*pItr).xest = (*pItr).xini + (*pItr).xcor;
			f_x.push_back((*pItr).xest);
			//LOGSTAT(true,"%03d %04d %16.3lf %16.3lf %16.3lf",i,(*pItr).iobs,(*pItr).xini,(*pItr).xcor,(*pItr).xest);
			++pItr;
		} else {
			(*aItr)->xsig = 0.0;
			(*aItr)->xcor = est[i];
			(*aItr)->xest = (*aItr)->xini + (*aItr)->xcor;
			++aItr;
		}
	}
	free(est);
	m_getDopValue(); // compute dop here
	if(dly->liar){
		lfix = strstr(dly->cobs,"IF") ? m_fixamb_if(m_line,m_atom):m_fixamb_raw(m_line,m_atom);
	}
	m_calAtom_pam(m_line,m_atom);
}
void RtkBamboo::v_onAddObsMatrix(){
	int refsit,i,isat, refsat,ifreq, jp, i_sitobs, ind, isys, isit, jnd, jsys,nsumsit,nuseob;
	double phase, range, range1, phase1, range1_ref, phase1_ref, range0, phase0,
			range0_ref, phase0_ref,var1,var0,lambda,lambda_ref;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	// STEP1,smooth the MW ION here
	m_smoothMW(m_atom); // smooth mw observation
//	m_smoothMWIon(m_line,m_atom);
	// begin to sort the
	nsumsit = nobs = 0;
	isit = m_line.isit;
	refsit = m_line.refsit;
	m_pkRefsat(m_line,m_atom);
	m_fillAmat(m_line);  //
	m_genRMatrix(m_line); //
	///
	for (i = 0; i < imtx + 1; i++) {
		memset(infs + iptx[i] + imtx + nobs, 0, sizeof(double) * m_line.nobs);
		memset(Hk + iphk[i] + nobs, 0, sizeof(double) * m_line.nobs);
	}
	range1_ref = range0_ref = phase1_ref = phase0_ref = 0.0;
	for (jsys = 0, i_sitobs = 0; jsys < dly->nsys; jsys++) {
		isys = index_string(SYS, dly->system[jsys]);
		refsat = m_line.refsat[isys];
		if (refsat == -1)
			continue;
		for (ifreq = 0; ifreq < dly->nfq[isys]; ifreq++) {
			if (!m_getRange(m_line, refsat, ifreq, range1_ref, range0_ref, var1,
					var0)
					|| !m_getPhase(m_line, refsat, ifreq, phase1_ref, phase0_ref,
							var1, var0))
				continue;
			// TODO:
			lambda_ref = strstr(dly->cobs,"IF") ? 1.0 : dly->SAT[refsat].lamda[ifreq];
			if (m_Sta[isit].ob.flag[refsat][ifreq] || m_Sta[refsit].ob.flag[refsat][ifreq])
				ptrAmtag[isit * MAXSAT + refsat][ifreq]->xest = ptrAmtag[isit * MAXSAT + refsat][ifreq]->xini = (phase1_ref
						- range1_ref - phase0_ref + range0_ref) / lambda_ref;
			m_line.omc[refsat][ifreq] = phase1_ref - phase0_ref;
			m_line.omc[refsat][MAXFREQ + ifreq] = range1_ref - range0_ref;
			for (isat = 0; isat < dly->nprn; isat++) {
				if (dly->cprn[isat][0] != SYS[isys]
						|| isat == m_line.refsat[isys])
					continue;
				lambda = strstr(dly->cobs,"IF") ? 1.0 : dly->SAT[isat].lamda[ifreq];
				nuseob = 0;
				if (strstr(dly->uobs, "PHASE") || strstr(dly->uobs, "GRAPHIC")) {
					if (!m_getPhase(m_line, isat, ifreq, phase1, phase0, var1,
							var0)) //
						continue;
					ipob[nobs][0] = isat;
					ipob[nobs][1] = 1;
					ipob[nobs][2] = ifreq;
					ipob[nobs][3] = isit;
					weig[nobs] = 1.0 / sqrt(var1 + var0);
					phase = phase1 - phase0 - phase1_ref + phase0_ref;
//					if(SYS[isys] == 'R'){
//						//// will remove the reference constant value here
//						phase = phase - (lambda - lambda_ref) * ptrAmtag[isit * MAXSAT+refsat][ifreq]->xest;
//					}
					m_line.omc[isat][ifreq] = phase1 - phase0;
					nuseob++;
					nobs++;
					i_sitobs++;
				}
				if (strstr(dly->uobs, "CODE")) {
					if (!m_getRange(m_line, isat, ifreq, range1, range0, var1,
							var0))
						continue;
					ipob[nobs][0] = isat;
					ipob[nobs][1] = 2;
					ipob[nobs][2] = ifreq;
					ipob[nobs][3] = isit;
					weig[nobs] = 1.0 / sqrt(var1 + var0);
					range = range1 - range0 - range1_ref + range0_ref;
					m_line.omc[isat][MAXFREQ + ifreq] = range1 - range0;
					nuseob++;
					nobs++;
					i_sitobs++;
				}
				//TODO:
				if (m_Sta[isit].ob.flag[isat][ifreq] || m_Sta[refsit].ob.flag[isat][ifreq])
					ptrAmtag[isit * MAXSAT + isat][ifreq]->xini = (phase1 - range1 - phase0 + range0) / lambda;
				for (jp = 0; jp < m_line.n_np[isat]; jp++) {
					jnd = m_line.i_coef[isat][jp][ifreq];
					if ((strstr(dly->uobs, "PHASE") || strstr(dly->uobs, "GRAPHIC"))
							&& strstr(dly->uobs, "CODE")) {
						/// ADD the HK here
						Hk[iphk[jnd] + nobs - 2] = m_line.d_coef[isat][jp][ifreq];
						Hk[iphk[jnd] + nobs - 1] = m_line.d_coef[isat][jp][MAXFREQ
								+ ifreq];
						Hk[iphk[imtx] + nobs - 2] = phase
								- (ptrAmtag[isit * MAXSAT + isat][ifreq]->xini * lambda
								- ptrAmtag[isit * MAXSAT + refsat][ifreq]->xini * lambda_ref);
						Hk[iphk[imtx] + nobs - 1] = range;
						if(jnd >= nps)
							continue;
						m_Pm[pmtag[jnd]].iobs += 2;
					} else if (strstr(dly->uobs, "PHASE")
							|| strstr(dly->uobs, "GRAPHIC")) {
						Hk[iphk[jnd] + nobs - 1] = m_line.d_coef[isat][jp][ifreq];
						Hk[iphk[imtx] + nobs - 1] = phase
								- (ptrAmtag[isit * MAXSAT + isat][ifreq]->xini * lambda
								- ptrAmtag[isit * MAXSAT + refsat][ifreq]->xini * lambda_ref);
						if(jnd >= nps)
							continue;
						m_Pm[pmtag[jnd]].iobs += 1;
					} else {
						Hk[iphk[jnd] + nobs - 1] = m_line.d_coef[isat][jp][MAXFREQ
								+ ifreq];
						Hk[iphk[imtx] + nobs - 1] = range;
						if(jnd >= nps)
							continue;
						m_Pm[pmtag[jnd]].iobs += 1;
					}
				}
				for (ind = 0; ind < m_line.nobs; ind++) {
					for (jp = 0; jp < m_line.n_np[isat]; jp++) {
						jnd =m_line.i_coef[isat][jp][ifreq];
						for (i = 0; i < nuseob; i++) {
							infs[iptx[jnd] + imtx + nsumsit + ind] += Hk[iphk[jnd]
									+ nobs - nuseob + i]
									* m_line.R[(i_sitobs - nuseob + i)
											* m_line.nobs + ind];
						}
					}
					for (i = 0; i < nuseob; i++) {
						infs[iptx[imtx] + imtx + nsumsit + ind] += Hk[iphk[imtx]
								+ nobs - nuseob + i]
								* m_line.R[(i_sitobs - nuseob + i) * m_line.nobs
										+ ind];
					}
				}
			}
		}
		if (strstr(dly->uobs, "PHASE") || strstr(dly->uobs, "GRAPHIC")) {
			for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
				for(isat = 0;isat < dly->nprn;isat++){
					if(dly->cprn[isat][0] != SYS[isys] || m_Sta[isit].ob.omc[isat][ifreq] == 0 ||
							m_Sta[refsit].ob.omc[isat][ifreq] == 0)
						continue;
					ptrAmtag[isit * MAXSAT + isat][ifreq]->iobs++;
					ptrAmtag[isit * MAXSAT + isat][ifreq]->elev =  ptrAmtag[isit * MAXSAT + isat][ifreq]->elev +
							m_Sta[isit].ob.elev[isat];
				}
			}
		}
	}
	nsumsit = nsumsit + i_sitobs;
	//CMat::CMat_PrintMatrix(Hk,nptk,nobs,imtx+1,"obs");
	//LOGMAT(infs,nmtxY,nobs + imtx ,imtx+1,9,3,true,"infs");
}
int RtkBamboo::m_fixamb_raw(Line_net& line,Line_amb& atom){
	///TODO:FIX THE AMBIGUITY HERE,USING LAMBDA METHOD
	int nfix_sys[MAXSYS] = {0};
	int ifreq,isat,iamb,jsys,isys,refsat,i,isit,refsit,ind,refind[MAXFREQ],nb,ntot,rind,nfix,lcont,insit,ldebug = false,stat,outamb,ndmw;
	double nlratio;
	nb=m_Am.size();
	isit = line.isit;
	refsit = line.refsit;
	line.xa = (double*)calloc(imtx,sizeof(double));
	line.xf = (double*)calloc(imtx,sizeof(double));
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	ndmw = imtx + MAXSAT;
	double* xi = (double*)calloc(imtx,sizeof(double)); /// initial x-states
	double* x_ = (double*)calloc(imtx,sizeof(double)); /// iteration x-states
	double* xa = (double*)calloc(imtx,sizeof(double)); /// resume fixed x-states
	double* xt = (double*)calloc(imtx,sizeof(double)); /// temporary structure
	double* P = (double*)calloc(imtx * imtx,sizeof(double)); /// variance structure
	double* tmp = (double*)calloc(imtx * imtx,sizeof(double)); /// temporary structure
	double* D = (double*)calloc(imtx * imtx,sizeof(double));   /// D-transformation matrix
	double* MW_sav = (double*)calloc(ndmw * (imtx + 1),sizeof(double));
	SrifAmbd* AB = (SrifAmbd*)calloc(nb,sizeof(SrifAmbd));     /// fix structure
	list<SrifPamt*>::iterator ambItr;
	atom.m_init(); /// initial value
	for(i = 0; i < imtx + 1;i++)
		memcpy(MW_sav + ndmw * i,infs + iptx[i],sizeof(double) * imtx);
	m_addMWCon(MW_sav,ndmw,m_line,m_atom); // add MW observation to constrain ionosphere
	//STEP 1,generate variance
	CMat::CMat_Matmul("TN",imtx,imtx,imtx,1.0,MW_sav,ndmw,MW_sav,ndmw,0.0,P,imtx);
	CMat::CMat_Inverse(P,imtx,imtx);
	/// initial x_float values
	for(i = 0; i < nps;i++)
		xi[i] = m_Pm[i].xest;
	for(ambItr = m_Am.begin(),i=0;ambItr != m_Am.end();++ambItr,++i){
		isit = (*ambItr)->isit;
		isat = (*ambItr)->psat;
		ifreq = (*ambItr)->ifreq;
		AB[i].abst = -1; // -2: will be deleted -1:no need to be fixed 0: to be fixed 1: fixed
		if(amtag[isit * MAXSAT + isat][ifreq] - nps != i)
			AB[i].abst = -2;  // will be deleted,should not use for ambiguities
		AB[i].pab=*ambItr;
		xi[i + nps] = (*ambItr)->xest;
	}
	if (ldebug) {
		LOGSTAT(true,"initial states:");
		LOGSTAT(true,"sod = %lf",sod);
		for (i = 0; i < imtx; i++) {
			if (i < nps) {
				LOGSTAT(true, "X(%02d) = %13.3lf ", i, x_[i]);
			} else {
				LOGSTAT(true,
						"X(%02d) = %13.3lf AB(cprn,abfr,abfx,psat,ifreq,abst,elev) = %s %9.3lf %9.3lf %02d %02d %5.1lf", i,
						xi[i], dly->cprn[AB[i - nps].pab->psat].c_str(),AB[i - nps].abfr, AB[i - nps].abfx, AB[i - nps].pab->ifreq,AB[i - nps].abst,m_Sta[refsit].ob.elev[AB[i - nps].pab->psat] * RAD2DEG);
			}
		}
	}
	memcpy(x_,xi,sizeof(double) * imtx); // copy the initial value into x_
	/// STEP 2,begin to fix the ambiguities for each system
	nlratio = 0.0;
	for(jsys = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		nfix_sys[isys] = 0;
		LOGSTAT(true,"begin to fix %c ambiguities",SYS[isys]);
		// generate D transformation matrix
		memset(D,0,sizeof(double) * imtx * imtx);
		for(i = 0; i < imtx;i++){
			D[i * imtx + i] = 1.0;
		}
		refsat = line.refsat[isys];
		// since the AB will be re-ordered,iteration for each system to find the index of reference
		for(ifreq = 0;ifreq < dly->nfq[isys] && refsat != -1;ifreq++){
			for(i = 0; i < nb;i++){
				if(AB[i].pab->ifreq == ifreq && AB[i].pab->psat == refsat && AB[i].abst != -2) // in case of the reference is replaced by new ambiguities
					refind[ifreq] = i;
			}
		}
		for(i = 0;i < nb && refsat != -1;i++){
			isat = AB[i].pab->psat;
			ifreq = AB[i].pab->ifreq;
			if(dly->cprn[isat][0] != SYS[isys] || m_Sta[isit].ob.omc[isat][0] == 0.0 ||
					m_Sta[refsit].ob.omc[isat][0] == 0.0 || refsat == isat)
				continue;
			if(AB[i].abst == -2)  // the ambiguity is replaced,will not fix this ambiguity
				continue;
			if(m_Sta[refsit].ob.elev[isat] < dly->cutoff || m_Sta[refsit].ob.elev[refsat] < dly->cutoff)
				continue;
			if(m_Sta[isit].ob.elev[isat] < dly->cutoff || m_Sta[isit].ob.elev[refsat] < dly->cutoff)
				continue;
			if((m_atom.ptime[0][isat][1] - m_atom.ptime[0][isat][0]) * 86400.0 < dly->rawseccommon)
				continue;
			if((m_atom.ptime[1][isat][1] - m_atom.ptime[1][isat][0]) * 86400.0 < dly->rawseccommon)
				continue;
			AB[i].elev = m_Sta[refsit].ob.elev[isat];
			AB[i].abst = 0;
			AB[i].abfr = x_[nps + i] - x_[refind[ifreq] + nps]; /// only for output here
			D[(refind[ifreq] + nps)* imtx + i + nps] = -1.0;
		}
		// transformed float values
		CMat::CMat_Matmul("NN",imtx,1,imtx,1.0,D,imtx,x_,imtx,0.0,xt,imtx);
		memcpy(x_,xt,sizeof(double) * imtx);
		// generate the variance
		CMat::CMat_Matmul("NN",imtx,imtx,imtx,1.0,D,imtx,P,imtx,0.0,tmp,imtx);
		CMat::CMat_Matmul("NT",imtx,imtx,imtx,1.0,tmp,imtx,D,imtx,0.0,P,imtx);
		// STEP 3,using partial lambda to fix ambiguities
		lcont = true;
		while(lcont){
			// begin to sort the P,change the order here,return the number of ambiguities without fixed ones
			ntot = Ambfix_rtk::s_plcfloat(AB,P,imtx,x_,imtx,nps);
			// begin to using lambda to fix the ambiguities
			lcont = Ambfix_rtk::s_resolveLambda(AB,P,imtx,x_,ntot,nps,dly->nlmaxdel[isys],nlratio) && (strstr(dly->fixmode,"ITERATION") || strstr(dly->fixmode,"iteration"));// iteration to fix the ambiguities
		}
		// change the unfixed ambiguities into -1
		for(i = 0; i < nb;i++){
			if(AB[i].abst == 0) AB[i].abst = -1;
		}
		// restore the ambiguities
		for (i = 0; i < nb; i++) {
			if (AB[i].abst == 1 && dly->cprn[AB[i].pab->psat][0] == SYS[isys]) { // fixed
				ifreq = AB[i].pab->ifreq;
				isat = AB[i].pab->psat;
				isys = index_string(SYS, dly->cprn[isat][0]);
				refsat = line.refsat[isys];
				rind = amtag[isit * MAXSAT + refsat][ifreq];
				ind = amtag[isit * MAXSAT + isat][ifreq];
				x_[nps + i] = xi[rind] + AB[i].abfx;
				nfix_sys[isys] = ifreq == 0 ? nfix_sys[isys] + 1 : nfix_sys[isys];
			}else if(dly->cprn[AB[i].pab->psat][0] == SYS[isys]){
				// when the slip is detected,the m_Am will reserve two ambiguities for this satellite
				// in this case,the x_ may be incorrect for the replaced ambiguity
				ifreq = AB[i].pab->ifreq;
				isat = AB[i].pab->psat;
				isys = index_string(SYS, dly->cprn[isat][0]);
				ind = amtag[isit * MAXSAT + isat][ifreq];
				x_[nps + i] = xi[ind];               // initial value
			}
		}
		nfix_sys[isys] = nfix_sys[isys] > 0 ? nfix_sys[isys] + 1 : nfix_sys[isys];
		// another system
		if (ldebug) {
			LOGSTAT(true,"sod = %lf ratio = %lf",sod,m_fixAdapter.ratio);
			for (i = 0; i < imtx; i++) {
				if (i < nps) {
					LOGSTAT(true, "X(%02d) = %13.3lf ", i, x_[i]);
				} else {
					LOGSTAT(true,
							"X(%02d) = %13.3lf AB(cprn,abfr,abfx,psat,ifreq,abst,elev) = %s %9.3lf %9.3lf %02d %02d %5.1lf", i,
							x_[i], dly->cprn[AB[i - nps].pab->psat].c_str(),AB[i - nps].abfr, AB[i - nps].abfx, AB[i - nps].pab->ifreq,AB[i - nps].abst,m_Sta[refsit].ob.elev[AB[i - nps].pab->psat] * RAD2DEG);
				}
			}
		}
		LOGSTAT(true,"after fix %c ambiguities",SYS[isys]);
	}
	// resume to fixed states here,will continue for lower fixed satellites
	memcpy(xa,x_,sizeof(double) * nps);
	Line_amb amb;
	bool lfix = false;
	for(jsys = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		for(i = 0; i < nb;i++){
			ifreq = AB[i].pab->ifreq;
			isat = AB[i].pab->psat;
			insit = AB[i].pab->isit;
			if(dly->cprn[isat][0] != SYS[isys])
				continue;
			ind = amtag[insit * MAXSAT + isat][ifreq];
			xa[ind] = x_[nps + i];
			if(m_Sta[isit].ob.flag[isat][ifreq])
				atom.m_resetAmb(isit,isat,ifreq); /// reset ambiguities
			if(m_Sta[refsit].ob.flag[isat][ifreq])
				atom.m_resetAmb(refsit,isat,ifreq); /// reset ambiguities
			if(AB[i].abst == 1 && nfix_sys[isys] >= dly->nlminsav){
				// should look deeper
				amb.ifab[isat][ifreq] = AMB_FIX;
				amb.ifab[line.refsat[isys]][ifreq] = AMB_FIX;

				amb.Raw[isat][ifreq] = NINT(AB[i].abfx);
				amb.Raw[line.refsat[isys]][ifreq] = 0;
				amb.refsat[isys][ifreq] = line.refsat[isys];
				lfix = true;
			}
		}
	}
	this->m_fixAdapter.nfix = sum(nfix_sys,MAXSYS);
	this->m_fixAdapter.ratio = nlratio;
	memcpy(line.xf,xi,sizeof(double) * imtx);
	memcpy(line.xa,xa,sizeof(double) * imtx);
	atom.m_setAmbiguity(amb);
	// debug here
	if(lfix)
		m_nfix = m_nfix + 1;
	m_total = m_total + 1;
	LOGPRINT("%lf\n",1.0 * m_nfix / m_total * 100.0);
	LOGI(LEVEL_ALGO,"mjd = %d,sod = %lf,fixrate = %lf",mjd,sod,1.0 * m_nfix / m_total * 100.0);
	if (ldebug) {
		LOGMAT(xi,imtx,imtx,1,13,3,true,"float-states");
		LOGMAT(xa,imtx,imtx,1,13,3,true,"fix-states");
	}
	// decide whether to hold,detect resi and wl,using hold number to hold ambiguities
	m_holdamb(line,atom,xa);
	// update the xest here
	for(i = 0;i < nps;i++){
		m_Pm[i].xest = x_[i];
		if(lfix){
			m_Pm[i].xcor = m_Pm[i].xest - m_Pm[i].xini; // for time update
		}
	}
	free(P);
	free(D);
	free(AB);
	free(x_);
	free(xi);
	free(xt);
	free(xa);
	free(tmp);
	free(MW_sav);
}
void RtkBamboo::m_smoothMW(Line_amb& atom){
	int i,isit,isat;
	double value,wwl,rwl;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	/// begin to valid ambiguities using MW combination
	for(i = 0; i < 2;i++){
		if(i == 0) isit = atom.isit;
		else if(i == 1) isit = atom.refsit;
		for(isat = 0;isat < dly->nprn;isat++){
			if(m_Sta[isit].ob.omc[isat][MAXFREQ] == 0.0)
				continue;
			if(m_Sta[isit].ob.flag[isat][0] || m_Sta[isit].ob.flag[isat][1]){
				atom.abin[i][isat] = atom.xrwl[i][isat] = atom.weig[i][isat] = atom.xswl[i][isat] = atom.ptime[i][isat][0] = atom.ptime[i][isat][1] = 0.0;
				atom.abwl[isat] = atom.abnl[isat] = 999;
			}
			value = m_Sta[isit].ob.obs[isat][0] - m_Sta[isit].ob.obs[isat][1]
					- (dly->SAT[isat].g * m_Sta[isit].ob.obs[isat][MAXFREQ]
							+ m_Sta[isit].ob.obs[isat][MAXFREQ + 1])
							/ (1.0 + dly->SAT[isat].g)
							/ dly->SAT[isat].lamdw;
			wwl = 1.0;
			if(m_Sta[isit].ob.elev[isat] * RAD2DEG <= 30.0)
				wwl = wwl * 2.0 * sin(m_Sta[isit].ob.elev[isat]);
			if(atom.ptime[i][isat][0] == 0.0){
				atom.abin[i][isat] = NINT(value);
				atom.ptime[i][isat][0] = mjd + sod / 86400.0;
			}
			rwl = value - atom.abin[i][isat];
			atom.xrwl[i][isat] = atom.xrwl[i][isat] + wwl * rwl;
			atom.weig[i][isat] = atom.weig[i][isat] + wwl;
			atom.xswl[i][isat] = atom.xswl[i][isat] + wwl * rwl * rwl;
			atom.ptime[i][isat][1] = mjd + sod / 86400.0;
		}
	}
}
void RtkBamboo::m_holdamb(Line_net& line,Line_amb& atom,double* xa){
	/// hold ambiguities here,using hold count,resi,wl-test to decide whether to hold this ambiguities
	int i,isat,isit,refsit,ifreq,isys,jsys,ind,refind,m_nobs;
	list<SrifPamt*>::iterator ambItr;
	const double WEIG_HOLD=1.0 / sqrt(0.001);
	isit = line.isit;
	refsit = line.refsit;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	for(i = 0,ambItr = m_Am.begin(); i < imtx;i++){
		if(i < nps){
			xa[i] = xa[i] - m_Pm[i].xini;
		}else{
			xa[i] = xa[i] - (*ambItr)->xini;
			++ambItr;
		}
	}
	m_nobs = 0;
	for(jsys = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		for(isat = 0;isat < dly->nprn;isat++){
			if(dly->cprn[isat][0] != SYS[isys] || isat == line.refsat[isys]) continue;
			if(m_Sta[refsit].ob.elev[isat] > dly->holdcutoff){
				if(atom.ifab[isat][0] < AMB_FIX || atom.ihold[isat][0] * intv < dly->holdtimelim)
					continue;
				if(atom.ifab[isat][0] < AMB_FIX || atom.ihold[isat][0] * intv < dly->holdtimelim)
					continue;
				if(m_validMW(atom,isat,line.refsat[isys])){
					atom.pifab[isat][0] = AMB_HOLD;
					atom.pifab[isat][1] = AMB_HOLD;
					if(strstr(dly->ambmod,"HOLD")){
						ind = amtag[isit * MAXSAT + isat][0];
						refind = amtag[isit * MAXSAT + line.refsat[isys]][0];

						infs[iptx[ind] + imtx + m_nobs] = 1.0 * WEIG_HOLD;
						infs[iptx[refind] + imtx + m_nobs] = -1.0 * WEIG_HOLD;
						infs[iptx[imtx] + imtx + m_nobs] = (xa[ind] - xa[refind]) * WEIG_HOLD;
						m_nobs = m_nobs + 1;

						ind = amtag[isit * MAXSAT + isat][1];
						refind = amtag[isit * MAXSAT + line.refsat[isys]][1];

						infs[iptx[ind] + imtx + m_nobs] = 1.0 * WEIG_HOLD;
						infs[iptx[refind] + imtx + m_nobs] = -1.0 * WEIG_HOLD;
						infs[iptx[imtx] + imtx + m_nobs] = (xa[ind] - xa[refind]) * WEIG_HOLD;
						m_nobs = m_nobs + 1;
					}
				}
			}
		}
	}
	if(m_nobs > 0)
		CMat::CMat_Householder(infs, nmtxY, imtx, imtx + m_nobs, imtx + 1, false);
}
void RtkBamboo::m_calAtomWL(Line_net& line,Line_amb& atom){
}
void RtkBamboo::m_calAtom(Line_net& line,Line_amb& atom){
	/// calculate the atomsphere using fixed ambiguities
	RtkDDIon ion,trp;
	double maxiondif = 0.08;
	int i,isys,jsys,refsat,ifreq,isit,refsit,isat,oN1,oN2,ledit,lfixed,subON1,subON2,status,nfix;
	char ztdtag_1[256],ztdtag_0[256];
	double obs0,obs1,obsref0,obsref1,dion,N1,N2,obsif0,obsif1,obsifref0,obsifref1,ddtrp,ddion,ddion_wl,ddion1,ddion2;
	double ztdcor_1,ztdcor_0,resd;
	loginfo.m_reset();
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	isit = atom.isit;
	refsit = atom.refsit;
	Line_net cline;
	cline.isit = line.isit;
	cline.refsit = line.refsit;
	atom.m_inquireRef(this,cline);
	m_fillAmat(cline);
	m_genRMatrix(cline);
	memcpy(cline.omc,line.omc,sizeof(line.omc));
	for(jsys = 0,nfix = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		if(-1 == (refsat = cline.refsat[isys]))
			continue;
		// no observation for reference satellite
		if(m_Sta[isit].ob.omc[refsat][0] == 0 || m_Sta[refsit].ob.omc[refsat][0] == 0)
			continue;
		obsref0 = m_Sta[refsit].ob.obs[refsat][0] * dly->SAT[refsat].lamda[0] - m_Sta[refsit].ob.obs[refsat][1] * dly->SAT[refsat].lamda[1];
		obsref1 = m_Sta[isit].ob.obs[refsat][0] * dly->SAT[refsat].lamda[0] - m_Sta[isit].ob.obs[refsat][1] * dly->SAT[refsat].lamda[1];

		obsifref0 = m_Sta[refsit].ob.obs[refsat][0] * dly->SAT[refsat].lamda[0] * dly->SAT[refsat].fac[0] - m_Sta[refsit].ob.obs[refsat][1] * dly->SAT[refsat].lamda[1] * dly->SAT[refsat].fac[1];
		obsifref1 = m_Sta[isit].ob.obs[refsat][0] * dly->SAT[refsat].lamda[0] * dly->SAT[refsat].fac[0] - m_Sta[isit].ob.obs[refsat][1] * dly->SAT[refsat].lamda[1] * dly->SAT[refsat].fac[1];

		for(isat = 0;isat < dly->nprn;isat++){
			if(dly->cprn[isat][0] != SYS[isys])
				continue;
			// no observation
			if(m_Sta[isit].ob.omc[isat][0] == 0 || m_Sta[refsit].ob.omc[isat][0] == 0)
				continue;
			if(m_Sta[isit].ob.elev[isat] < m_Sta[isit].cutoff || m_Sta[refsit].ob.elev[isat] < m_Sta[refsit].cutoff)
				continue;
			/// 10 degrees for cut,only output 10 degrees
			if(m_Sta[isit].ob.elev[isat] < ATOM_CUT * DEG2RAD || m_Sta[refsit].ob.elev[isat] < ATOM_CUT * DEG2RAD)
				continue;
			// not fixed in the memory
			if((atom.pifab[isat][0] < AMB_FIX || atom.pifab[isat][1] < AMB_FIX))
				continue;
			lfixed = true;
			N1 = atom.Raw[isat][0] - atom.Raw[refsat][0]; //  single difference
			N2 = atom.Raw[isat][1] - atom.Raw[refsat][1]; //  single difference
			obs0 = m_Sta[refsit].ob.obs[isat][0] * dly->SAT[isat].lamda[0] - m_Sta[refsit].ob.obs[isat][1] * dly->SAT[isat].lamda[1];
			obs1 = m_Sta[isit].ob.obs[isat][0] * dly->SAT[isat].lamda[0] - m_Sta[isit].ob.obs[isat][1] * dly->SAT[isat].lamda[1];

			obsif0 = m_Sta[refsit].ob.obs[isat][0] * dly->SAT[isat].lamda[0] * dly->SAT[isat].fac[0] - m_Sta[refsit].ob.obs[isat][1] * dly->SAT[isat].lamda[1] * dly->SAT[isat].fac[1];
			obsif1 = m_Sta[isit].ob.obs[isat][0] * dly->SAT[isat].lamda[0] * dly->SAT[isat].fac[0] - m_Sta[isit].ob.obs[isat][1] * dly->SAT[isat].lamda[1] * dly->SAT[isat].fac[1];
			// logging part
			loginfo.Raw[isat][0] = N1;
			loginfo.Raw[isat][1] = N2;
			// current not fixed or elevation is smaller than 25 degrees,using the valid function to find whether the ambiguity is ok
			if(atom.ifab[isat][0] < AMB_FIX || atom.ifab[isat][1] < AMB_FIX){
				oN1 = NINT(N1);
				oN2 = NINT(N2);
				/// not hold and elev is smaller than 20 degrees,will continue
				if((atom.pifab[isat][0] < AMB_HOLD || atom.pifab[isat][0] < AMB_HOLD) && m_Sta[refsit].ob.elev[isat] < 20.0 * DEG2RAD)
					continue;
				/// hold or elevation is more than 20 degrees,will search here
				/// no wl ion for compare,will continue
				if(!atom.smion.m_inquireIon(mjd + sod / 86400.0,isat,refsat,ddion_wl))
					continue;
				ddion_wl = ddion_wl * dly->SAT[isat].freq[1] / dly->SAT[isat].freq[0];
				/// valid IF residual failed,will continue if the residual is big
				if(!(status = m_validIFAmb(true,cline,atom,line.xa,isat,oN1,oN2,subON1,subON2,resd))){ // using float solution instead
					lfixed = false;
				}else if(status == 2){
					/// dd-ztd is not distinguished enough,so will search the ambiguities
					ddion1 = dly->SAT[isat].fac[1] * (obs1 - obs0 - obsref1 + obsref0 - (dly->SAT[isat].lamda[0] * oN1 - dly->SAT[isat].lamda[1] * oN2));
					ddion2 = dly->SAT[isat].fac[1] * (obs1 - obs0 - obsref1 + obsref0 - (dly->SAT[isat].lamda[0] * subON1 - dly->SAT[isat].lamda[1] * subON2));
					/// compare the two diff,and decide which is the best
					if(fabs(ddion_wl - ddion1) < fabs(ddion_wl - ddion2) && fabs(ddion_wl - ddion1) > maxiondif){
						lfixed = false;
					}else if(fabs(ddion_wl - ddion1) > fabs(ddion_wl - ddion2)){
						if(fabs(ddion_wl - ddion2) > maxiondif)
							lfixed = false;
						else{
							oN1 = subON1;
							oN2 = subON2;
						}
					}
				}else{
					/// dd-ztd is distinguished enough,verify the ion with wl-ion
					ddion1 = dly->SAT[isat].fac[1] * (obs1 - obs0 - obsref1 + obsref0 - (dly->SAT[isat].lamda[0] * oN1 - dly->SAT[isat].lamda[1] * oN2));
					if(fabs(ddion_wl - ddion1) > maxiondif)
						lfixed = false;
				}
				N1 = oN1;
				N2 = oN2;
			}else{
				oN1 = NINT(N1);
				oN2 = NINT(N2);
				/// check check check,if the elevation is lower than 20,and verified failed,just continue,and will not search
				if(m_Sta[isit].ob.elev[isat] < 20.0 * DEG2RAD || m_Sta[refsit].ob.elev[isat] < 20.0 * DEG2RAD){
					if(!m_validIFAmb(false,cline,atom,line.xa,isat,oN1,oN2,subON1,subON2,resd)){// using float solution instead,if the residual is not distinguished enough,make it not fix
						continue;
					}else{
						/// check the ion
						if(!atom.smion.m_inquireIon(mjd + sod / 86400.0,isat,refsat,ddion_wl))
							continue;
						ddion_wl = ddion_wl * dly->SAT[isat].freq[1] / dly->SAT[isat].freq[0];
						/// dd-ion
						ddion = dly->SAT[isat].fac[1] * (obs1 - obs0 - obsref1 + obsref0 - (dly->SAT[isat].lamda[0] * oN1 - dly->SAT[isat].lamda[1] * oN2));
						if(fabs(ddion - ddion_wl) > maxiondif){
							/// distinguish quite enough
							continue;
						}
					}
				}
			}
			if(lfixed){
				/// dd-ion
				ddion = dly->SAT[isat].fac[1] * (obs1 - obs0 - obsref1 + obsref0 - (dly->SAT[isat].lamda[0] * N1 - dly->SAT[isat].lamda[1] * N2));
				/// dd-trop
				ddtrp = obsif1 - obsif0 - obsifref1 + obsifref0 - N1 * dly->SAT[isat].lamda[0] * dly->SAT[isat].fac[0] + N2 * dly->SAT[isat].lamda[1] * dly->SAT[isat].fac[1] - (
									m_Sta[isit].ob.rleng[isat] - m_Sta[isit].ob.rleng[refsat] - m_Sta[refsit].ob.rleng[isat] + m_Sta[refsit].ob.rleng[refsat]);
			}else{
				// search failed or verified failed
				if(m_Sta[refsit].ob.elev[isat] < 24.0 * DEG2RAD)
					continue;
				if (!atom.smion.m_inquireIon(mjd + sod / 86400.0, isat, refsat,
						ddion_wl))
					continue;
				ddion = ddion_wl * dly->SAT[isat].freq[1]
						/ dly->SAT[isat].freq[0];

				sprintf(ztdtag_1, "%s_%03d", dly->ztdmod, isit);
				sprintf(ztdtag_0, "%s_%03d", dly->ztdmod, refsit);
				for (i = 0; i < nps; i++) {
					if (strstr(m_Pm[i].pname, ztdtag_0)) {
						ztdcor_0 = line.xa[m_Pm[i].ind];
					}
					if (strstr(m_Pm[i].pname, ztdtag_1)) {
						ztdcor_1 = line.xa[m_Pm[i].ind];
					}
				}
				ddtrp = m_Sta[isit].ob.trpModel[isat] + m_Sta[isit].ob.zmap[isat] * ztdcor_1
						- m_Sta[refsit].ob.trpModel[isat] - m_Sta[refsit].ob.zmap[isat] * ztdcor_0
						- (m_Sta[isit].ob.trpModel[refsat]
								+ m_Sta[isit].ob.zmap[refsat] * ztdcor_1
								- m_Sta[refsit].ob.trpModel[refsat]
								- m_Sta[refsit].ob.zmap[refsat] * ztdcor_0);
			}
			nfix = nfix + 1;
			/// update structure
			ion.elev[isat] = MIN(m_Sta[isit].ob.elev[isat],m_Sta[refsit].ob.elev[isat]);
			ion.ddion[isat] = ddion;
			ion.ifix[isat] = atom.pifab[isat][0];
			ion.refsat[isys] = refsat;
			ion.ptime = mjd + sod / 86400.0;

			trp.elev[isat] = MIN(m_Sta[isit].ob.elev[isat],m_Sta[refsit].ob.elev[isat]);
			trp.ddion[isat] = ddtrp;
			trp.ifix[isat] = atom.pifab[isat][0];
			trp.refsat[isys] = refsat;
			trp.ptime = mjd + sod / 86400.0;

			/// update logging part
			loginfo.ifab[isat] = atom.ifab[isat][0];
			loginfo.pifab[isat] = atom.pifab[isat][0];

			loginfo.elev[isat] = m_Sta[isit].ob.elev[isat];
			loginfo.ion[isat] = ion.ddion[isat];
			loginfo.trp[isat] = trp.ddion[isat];
			loginfo.refsat[isys] = refsat;
			if (lfixed) {
				loginfo.useRaw[isat][0] = N1;
				loginfo.useRaw[isat][1] = N2;
				loginfo.status[isat] = 1; // valid success
			} else {
				loginfo.useRaw[isat][0] = -1;
				loginfo.useRaw[isat][1] = -1;
				loginfo.status[isat] = 2; // valid failed,using wl-ion
			}
			if (loginfo.useRaw[isat][0] != loginfo.Raw[isat][0]
					|| loginfo.useRaw[isat][1] != loginfo.Raw[isat][1]) {
				loginfo.idif[isat] = 1;
			}

			if (loginfo.status[isat] != 2) {
				m_validIFAmb(false, cline, atom, line.xa, isat,
						loginfo.useRaw[isat][0], loginfo.useRaw[isat][1],
						subON1, subON2, loginfo.resiIF_new[isat]);
				m_validAmb(false, cline, atom, line.xa, isat, 0,
						loginfo.useRaw[isat][0], oN1, 2.8,
						loginfo.resiRaw_new[isat][0]);
				m_validAmb(false, cline, atom, line.xa, isat, 1,
						loginfo.useRaw[isat][1], oN2, 2.8,
						loginfo.resiRaw_new[isat][1]);
			}
			if (loginfo.idif[isat]) {
				m_validIFAmb(false, cline, atom, line.xa, isat,
						loginfo.Raw[isat][0], loginfo.Raw[isat][1], subON1,
						subON2, loginfo.resiIF[isat]);
				m_validAmb(false, cline, atom, line.xa, isat, 0,
						loginfo.Raw[isat][0], oN1, 2.8,
						loginfo.resiRaw[isat][0]);
				m_validAmb(false, cline, atom, line.xa, isat, 1,
						loginfo.Raw[isat][1], oN2, 2.8,
						loginfo.resiRaw[isat][1]);

				loginfo.ion_new[isat] = dly->SAT[isat].fac[1]
						* (obs1 - obs0 - obsref1 + obsref0
								- (dly->SAT[isat].lamda[0]
										* loginfo.Raw[isat][0]
										- dly->SAT[isat].lamda[1]
												* loginfo.Raw[isat][1]));
				loginfo.trp_new[isat] = obsif1 - obsif0 - obsifref1 + obsifref0
						- loginfo.Raw[isat][0] * dly->SAT[isat].lamda[0]
								* dly->SAT[isat].fac[0]
						+ loginfo.Raw[isat][1] * dly->SAT[isat].lamda[1]
								* dly->SAT[isat].fac[1]
						- (m_Sta[isit].ob.rleng[isat]
								- m_Sta[isit].ob.rleng[refsat]
								- m_Sta[refsit].ob.rleng[isat]
								+ m_Sta[refsit].ob.rleng[refsat]);
			}
		}
		for (isat = 0; isat < dly->nprn; isat++) {
			if (dly->cprn[isat][0] != SYS[isys])
				continue;
			// no observation
			if (m_Sta[isit].ob.omc[isat][0] == 0
					|| m_Sta[refsit].ob.omc[isat][0] == 0)
				continue;
			if (m_Sta[isit].ob.elev[isat] < m_Sta[isit].cutoff
					|| m_Sta[refsit].ob.elev[isat] < m_Sta[refsit].cutoff)
				continue;
			///////////////////////////// update logging part /////////////////////////
			if (atom.smion.m_inquireIon(mjd + sod / 86400.0, isat, refsat,
					ddion_wl)) {
				loginfo.ion_wl[isat] = ddion_wl * dly->SAT[isat].freq[1]
						/ dly->SAT[isat].freq[0];
			}
		}
	}
	if(nfix > 0){
		/// ion
		atom.sml1ion.m_ObsUpdate(mjd + sod / 86400.0,ion.ddion,ion.elev,ion.ifix,ion.refsat);
		for(isat = 0;isat < dly->nprn;isat++){
			isys = index_string(SYS,dly->cprn[isat][0]);
			if(atom.sml1ion.m_inquireIon(mjd + sod / 86400.0,isat,ion.refsat[isys],ddion)){
				/// if smooth success,use the smooth value instead
				ion.ddion[isat] = ddion;
			}
		}
		/// ztd
		atom.sml1ztd.m_ObsUpdate(mjd + sod / 86400.0,trp.ddion,trp.elev,trp.ifix,trp.refsat);
		for (isat = 0; isat < dly->nprn; isat++) {
			isys = index_string(SYS, dly->cprn[isat][0]);
			if (atom.sml1ztd.m_inquireIon(mjd + sod / 86400.0, isat,
					trp.refsat[isys], ddtrp)) {
				/// if smooth success,use the smooth value instead
				trp.ddion[isat] = ddtrp;
			}
		}
	}
	atom.ddblion.m_setData(ion);
	atom.ddbltrp.m_setData(trp);
}

void RtkBamboo::v_onOutput(){
	int isys,isit,refsit,jsys,nfix[MAXSYS],nfix_atom[MAXSYS],ntotal[MAXSYS],isat,iobs,ldebug = false;
	double pkresi0,pkresi1;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	// output the current status here
	m_output();
	isit = m_line.isit;
	refsit = m_line.refsit;
	const char* s_logAmb = "amb_atom.sav";
	if(!Bamboo::s_getInstance()->logger.m_lexist(s_logAmb)){
		Bamboo::s_getInstance()->logger.m_openLog(s_logAmb,86400);
	}
	for (jsys = 0; jsys < dly->nsys; jsys++) {
		isys = index_string(SYS, dly->system[jsys]);
		nfix_atom[isys] = nfix[isys] = ntotal[isys] = 0;
		for (isat = 0; isat < dly->nprn; isat++) {
			if (dly->cprn[isat][0] != SYS[isys])
				continue;
			if (m_Sta[isit].ob.obs[isat][0] == 0
					|| m_Sta[refsit].ob.obs[isat][0] == 0)
				continue;
			if (m_Sta[isit].ob.elev[isat] < m_Sta[isit].cutoff
					|| m_Sta[refsit].ob.elev[isat]
							< m_Sta[refsit].cutoff)
				continue;

			ntotal[isys]++;
			if (m_atom.ifab[isat][0] >= AMB_FIX && m_atom.ifab[isat][1] >= AMB_FIX)
				nfix[isys]++;
			if(loginfo.status[isat] != 0)
				nfix_atom[isys]++;
		}
	}
	Bamboo::s_getInstance()->logger.m_wtMsg("@%s %05d %7.1lf %s - %s ",
			s_logAmb, mjd, sod, m_Sta[isit].name.c_str(),
			m_Sta[refsit].name.c_str());
	for (jsys = 0; jsys < dly->nsys; jsys++) {
		isys = index_string(SYS, dly->system[jsys]);
		Bamboo::s_getInstance()->logger.m_wtMsg("@%s %02d/%02d/%02d ", s_logAmb,
				nfix[isys],nfix_atom[isys],ntotal[isys]);
	}
	Bamboo::s_getInstance()->logger.m_wtMsg("@%s %9.3lf ",s_logAmb,esig);

	Bamboo::s_getInstance()->logger.m_wtMsg("@%s \n", s_logAmb);
	for (isat = 0; isat < dly->nprn; isat++) {
		isys = index_string(SYS, dly->cprn[isat][0]);
		if (m_Sta[isit].ob.obs[isat][0] == 0
				|| m_Sta[refsit].ob.obs[isat][0] == 0)
			continue;
		if (m_Sta[isit].ob.elev[isat] < m_Sta[isit].cutoff
				|| m_Sta[refsit].ob.elev[isat] < m_Sta[refsit].cutoff)
			continue;
		for(iobs = 0,pkresi0 = pkresi1 = 0;iobs < nobs;iobs++){
			if(ipob[iobs][0] == isat && ipob[iobs][1] == 1 && ipob[iobs][3] == isit){
				if(ipob[iobs][2] == 0)
					pkresi0 = pkresi[iobs];
				if(ipob[iobs][2] == 1)
					pkresi1 = pkresi[iobs];
			}
		}
		Bamboo::s_getInstance()->logger.m_wtMsg(
				"@%s %s-%s %s-%s %02d %02d %02d %02d %02d %6.3lf %6.3lf %6.3lf %05d %5.1lf | %05d %05d %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %-10d %-10d",
				s_logAmb, m_Sta[isit].name.c_str(),
				m_Sta[refsit].name.c_str(), dly->cprn[isat].c_str(),
				loginfo.refsat[isys] != -1 ? dly->cprn[loginfo.refsat[isys]].c_str(): " * ",
				m_atom.ifab[isat][0],loginfo.pifab[isat],loginfo.status[isat],m_Sta[isit].ob.flag[isat][0] || m_Sta[refsit].ob.flag[isat][0],
				m_Sta[isit].ob.flag[isat][1] || m_Sta[refsit].ob.flag[isat][1],
				pkresi0,pkresi1,loginfo.ion_wl[isat],m_atom.abwl[isat],m_Sta[isit].ob.elev[isat] * RAD2DEG,
				loginfo.useRaw[isat][0],loginfo.useRaw[isat][1],
				loginfo.ion[isat],loginfo.trp[isat],loginfo.resiRaw_new[isat][0],loginfo.resiRaw_new[isat][1],
				loginfo.resiIF_new[isat],NINT((m_atom.ptime[0][isat][1] - m_atom.ptime[0][isat][0]) * 86400.0 / intv),
				NINT((m_atom.ptime[1][isat][1] - m_atom.ptime[1][isat][0]) * 86400.0 / intv));
		if(loginfo.idif[isat]){
			Bamboo::s_getInstance()->logger.m_wtMsg(
					"@%s | %05d %05d %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf ",s_logAmb,
					loginfo.Raw[isat][0],loginfo.Raw[isat][1],loginfo.ion_new[isat],loginfo.trp_new[isat],
					loginfo.resiRaw[isat][0],loginfo.resiRaw[isat][1],
					loginfo.resiIF[isat]);
		}
		Bamboo::s_getInstance()->logger.m_wtMsg("@%s \n",s_logAmb);
	}
	if(ldebug){
		const char* s_int_ion = "rtk_ion.txt";
		const char* s_int_trp = "rtk_trp.txt";
		if(!Bamboo::s_getInstance()->logger.m_lexist(s_int_trp))
			Bamboo::s_getInstance()->logger.m_openLog(s_int_trp);
		if(!Bamboo::s_getInstance()->logger.m_lexist(s_int_ion))
			Bamboo::s_getInstance()->logger.m_openLog(s_int_ion);

		for(isat = 0; isat < dly->nprn;isat++){
			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %6.3lf",s_int_ion,loginfo.ion[isat]);
			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %6.3lf",s_int_trp,loginfo.ion_wl[isat]);
		}
		Bamboo::s_getInstance()->logger.m_wtMsg("@%s \n",s_int_ion);
		Bamboo::s_getInstance()->logger.m_wtMsg("@%s \n",s_int_trp);
	}
}
/// it is not right since the reference is not right
bool RtkBamboo::m_validAmb(bool ledit,Line_net& line,Line_amb& atom,double* xa,int isat,int ifreq,int abfx,int& oabfx,int threshold,double& res){
	int i,jp,jnd,isit,refsit,refsat,isys,iobs,ip,it_abfx,idx;
	double omc_sum,lambda,axini,raxini,dres,res_tmp,phase,varc,minres;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	int i_i[] = {-1,1};
	bool lfixed;
	res = 0.0;
	///begin to compute residual using fixed states
	isys = index_string(SYS,dly->cprn[isat][0]);
	isit = line.isit;
	refsit = line.refsit;
	refsat = line.refsat[isys];
	omc_sum = 0.0;
	lambda = strstr(dly->cobs,"IF") ? 1.0 : dly->SAT[refsat].lamda[ifreq];
	for(jp = 0;jp < line.n_np[isat];jp++){
		jnd = line.i_coef[isat][jp][ifreq];
		ip = pmtag[jnd];
		if(jnd < nps){
			// parameters like xyz,ztd,ion
			omc_sum += line.d_coef[isat][jp][ifreq] * (xa[jnd] - m_Pm[ip].xini);
		}
	}
	// ambiguities parameters
	oabfx = abfx;
	lfixed = true;
	for(iobs = 0;iobs < line.nobs;iobs++){
		if(line.ipt[iobs][0] == 0 && line.ipt[iobs][1] == isat && line.ipt[iobs][2] == ifreq && line.ipt[iobs][3] == isit){
			varc = 1.0 / line.R[iobs * line.nobs + iobs];
			phase = line.omc[isat][ifreq] - line.omc[refsat][ifreq];
			minres = dres = fabs(phase - omc_sum - lambda * oabfx);
			/// compute the minimum ddres
			idx = 0;
			if(ledit){
				for(i = 0; i < sizeof(i_i) / sizeof(int);i++){
					res_tmp = fabs(phase - omc_sum - lambda * (oabfx + i_i[i]));
					if(res_tmp < minres){
						idx = i_i[i];
						minres = res_tmp;
					}
				}
			}
			lfixed = false;
			if(minres < varc * threshold){
				lfixed = true;
				oabfx = oabfx + idx;
			}
			break;
		}
	}
	res = dres;
	return lfixed;
}
int RtkBamboo::m_validIFAmb(bool ledit,Line_net& line,Line_amb& amb,double* xf,int isat,int& iN1,int& iN2,int& subON1,int& subON2,double& res){
	int lfixed;
	int i,j,isys,isit,refsit,refsat,jp,ifreq,jnd,ip,iobs,itN1,itN2,indx,nsav,idx[9][2],imin,imin2;
	double omc_sum[MAXFREQ] = {0.0},phase0,phase1,ddres;
	double ddres_sav[9],minres,std,var0,var1,minres2;
	subON1 = subON2 = 999;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	int i_i[] = {-1,0,1};
	isys = index_string(SYS,dly->cprn[isat][0]);
	isit = line.isit;
	refsit = line.refsit;
	refsat = line.refsat[isys];
	if(isat == refsat)  // do not deal with the reference satellite
		return true;
	for(jp = 0;jp < line.n_np[isat];jp++){
		for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
			jnd = line.i_coef[isat][jp][ifreq];
			ip = pmtag[jnd];
			if(jnd < nps){
				omc_sum[ifreq] += line.d_coef[isat][jp][ifreq] * (xf[jnd] - m_Pm[ip].xini);
			}
		}
	}
	/// begin to form IF observation
	lfixed = true;
	for(iobs = 0;iobs < line.nobs;iobs++){
		if(line.ipt[iobs][0] == 0 && line.ipt[iobs][1] == isat && line.ipt[iobs][2] == 0 && line.ipt[iobs][3] == isit){
			phase0 = line.omc[isat][0] - line.omc[refsat][0];
			var0 = 1.0 / line.R[iobs * line.nobs + iobs];
		}
		if(line.ipt[iobs][0] == 0 && line.ipt[iobs][1] == isat && line.ipt[iobs][2] == 1 && line.ipt[iobs][3] == isit){
			phase1 = line.omc[isat][1] - line.omc[refsat][1];
			var1 = 1.0 / line.R[iobs * line.nobs + iobs];
		}
	}
	nsav = 0;
	for(i = 0;i < sizeof(i_i) / sizeof(int);i++){
		itN1 = iN1 + i_i[i];
		for(j = 0;j < sizeof(i_i) / sizeof(int);j++){
			itN2 = iN2 + i_i[j];
			idx[nsav][0] = i;
			idx[nsav][1] = j;
			ddres_sav[nsav++] = fabs(phase0 * dly->SAT[isat].fac[0] - phase1 * dly->SAT[isat].fac[1] -
					dly->SAT[isat].lamda[0] * dly->SAT[isat].fac[0] * itN1 + dly->SAT[isat].lamda[1] * dly->SAT[isat].fac[1] * itN2);
		}
	}
	imin = -1;
	imin2 = -1;
	minres = 9e9;
	minres2 = 9e9;
	for(i = 0; i < 9;i++){
		if(ddres_sav[i] < minres){
			// update memory first
			imin2 = imin;
			minres2 = minres;
			// update
			minres = ddres_sav[i];
			imin = i;
		}else if(ddres_sav[i] < minres2){
			minres2 = ddres_sav[i];
			imin2 = i;
		}
	}
	std = sqrt(pow(dly->SAT[isat].fac[0],2) * pow(var0,2) + pow(dly->SAT[isat].fac[1],2) * pow(var1,2));
	lfixed = false;
	if(minres < 2.8 * std){
		lfixed = true;
		if(ledit){
			if(fabs(minres2 - minres) > 0.05){
				/// at least 5 cm difference,(0.19 * 2.54 - 0.24 * 1.54) = 0.11
				iN1 = iN1 + i_i[idx[imin][0]];
				iN2 = iN2 + i_i[idx[imin][1]];
			}else{
				subON1 = iN1 + i_i[idx[imin2][0]];
				subON2 = iN2 + i_i[idx[imin2][1]];


				iN1 = iN1 + i_i[idx[imin][0]];
				iN2 = iN2 + i_i[idx[imin][1]];
				lfixed = 2;
			}
		}else if(fabs(minres2 - minres) < 0.05 || i_i[idx[imin][0]] != 0 || i_i[idx[imin][1]] != 0){
			/// is not distinguish enough
			lfixed = false;
		}
	}
	res = ddres_sav[4];
	return lfixed;
}
void RtkBamboo::m_smoothMWIon(Line_net& line,Line_amb& atom){
	int isat,refsat,refsit,isit,lset[MAXSAT],nfix,jsys,isys,i,ifreq;
	double fac_wl[2],ddion;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	double ztdcor_0,ztdcor_1,lmw0,lmw1;
	double abwl0,abwl1,ddabwl,abwl0_ref,abwl1_ref,DDamb,atoms[MAXSAT],elevs[MAXSAT];
	double omc0[MAXFREQ],omc1[MAXFREQ],refomc0[MAXFREQ],refomc1[MAXFREQ],ioncor[MAXSAT][MAXFREQ];
	char ztdtag_0[256],ztdtag_1[256];
	isit = atom.isit;
	refsit = atom.refsit;

	ztdcor_1 = m_Sta[isit].ztdcor;
	ztdcor_0 = m_Sta[refsit].ztdcor;
	for(isat = 0;isat < dly->nprn;isat++){
		if(m_Sta[isit].ob.omc[isat][0] == 0 || m_Sta[refsit].ob.omc[isat][0] == 0)
			continue;
		isys = index_string(SYS,dly->cprn[isat][0]);
		for(ifreq = 0;ifreq < dly->nfq[isys];ifreq++){
			ioncor[isat][ifreq] = m_Sta[isit].ob.ioncor[isat] * m_Sta[isit].ob.ionpart[isat] * pow(dly->SAT[isat].freq[0], 2)
					/ pow(dly->SAT[isat].freq[ifreq], 2);
		}
	}
	memset(atoms,0,sizeof(atoms));
	memset(elevs,0,sizeof(elevs));
	for(jsys = 0,nfix = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		if(SYS[isys] == 'R')
			continue;
		if(-1 == (refsat = line.refsat[isys]))
			continue;
		if(m_Sta[isit].ob.omc[refsat][0] == 0 || m_Sta[refsit].ob.omc[refsat][0] == 0)
			continue;
		refomc0[0] = m_Sta[refsit].ob.omc[refsat][0] - ztdcor_0 * m_Sta[refsit].ob.zmap[refsat];
		refomc0[1] = m_Sta[refsit].ob.omc[refsat][1] - ztdcor_0 * m_Sta[refsit].ob.zmap[refsat];
		refomc1[0] = m_Sta[isit].ob.omc[refsat][0] - ztdcor_1 * m_Sta[isit].ob.zmap[refsat] - ioncor[refsat][0];
		refomc1[1] = m_Sta[isit].ob.omc[refsat][1] - ztdcor_1 * m_Sta[isit].ob.zmap[refsat] - ioncor[refsat][1];
		for(isat = 0;isat < dly->nprn;isat++){
			if(dly->cprn[isat][0] != SYS[isys])
				continue;
			lset[isat] = 0;
			if((atom.ptime[0][isat][1] - atom.ptime[0][isat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[1][isat][1] - atom.ptime[1][isat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[0][refsat][1] - atom.ptime[0][refsat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[1][refsat][1] - atom.ptime[1][refsat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if(m_Sta[isit].ob.omc[isat][MAXFREQ] == 0.0 || m_Sta[refsit].ob.omc[isat][MAXFREQ] == 0.0)
				continue;
			abwl0 = atom.xrwl[0][isat] / atom.weig[0][isat] + atom.abin[0][isat];
			abwl1 = atom.xrwl[1][isat] / atom.weig[1][isat] + atom.abin[1][isat];

			abwl0_ref = atom.xrwl[0][refsat] / atom.weig[0][refsat] + atom.abin[0][refsat];
			abwl1_ref = atom.xrwl[1][refsat] / atom.weig[1][refsat] + atom.abin[1][refsat];
			ddabwl = abwl0 - abwl1 - abwl0_ref + abwl1_ref;
			if(fabs(ddabwl - NINT(ddabwl)) > dly->wlmaxdev)
				continue;
			DDamb = NINT(ddabwl);
			omc0[0] = m_Sta[refsit].ob.omc[isat][0] - ztdcor_0 * m_Sta[refsit].ob.zmap[isat];
			omc0[1] = m_Sta[refsit].ob.omc[isat][1] - ztdcor_0 * m_Sta[refsit].ob.zmap[isat];
			omc1[0] = m_Sta[isit].ob.omc[isat][0] - ztdcor_1 * m_Sta[isit].ob.zmap[isat] - ioncor[isat][0];
			omc1[1] = m_Sta[isit].ob.omc[isat][1] - ztdcor_1 * m_Sta[isit].ob.zmap[isat] - ioncor[isat][1];

			fac_wl[0] = dly->SAT[isat].freq[0] / (dly->SAT[isat].freq[0] - dly->SAT[isat].freq[1]);
			fac_wl[1] = dly->SAT[isat].freq[1] / (dly->SAT[isat].freq[0] - dly->SAT[isat].freq[1]);


			lmw0=fac_wl[0] * (omc0[0] - refomc0[0]) - fac_wl[1] * (omc0[1] - refomc0[1]);
			lmw1=fac_wl[0] * (omc1[0] - refomc1[0]) - fac_wl[1] * (omc1[1] - refomc1[1]);

			elevs[isat] = m_Sta[refsit].ob.elev[isat];
			atoms[isat] = lmw1 - lmw0 - DDamb * dly->SAT[isat].lamdw; //
			lset[isat] = true;
			if(isat != refsat)
				nfix = nfix + 1;
		}
	}
	if(nfix > 0){
		atom.smion.m_ObsUpdate(mjd + sod / 86400.0,atoms,elevs,lset,line.refsat);
		/// detect the wl-ion,especially the bds ion,which maybe incorrect
		for(isat = 0;isat < dly->nprn;isat++){
			isys = index_string(SYS,dly->cprn[isat][0]);
			if(atom.smion.m_inquireIon(mjd + sod / 86400.0,isat,line.refsat[isys],ddion)){
				ddion = ddion * dly->SAT[isat].freq[1] / dly->SAT[isat].freq[0];
				if(fabs(ddion) > 0.5){
					/// which is maybe incorrect so remove this observation
					memset(m_Sta[isit].ob.omc[isat],0,sizeof(double) * 2 * MAXFREQ);
					memset(m_Sta[refsit].ob.omc[isat],0,sizeof(double) * 2 * MAXFREQ);
				}
			}
		}
	}
	//// for debug,output ion and its smooth value
//	const char* ion_raw = "ion_raw";
//	const char* ion_smooth = "ion_smooth";
//	const char* ion_elev = "ion_elev";
//	const char* ion_l1 = "ion_l1";
//	RtkDDIon ion_dd = atom.ddblion.m_inquireIon(line.refsat,mjd + sod / 86400.0);
//	if(!Bamboo::s_getInstance()->logger.m_lexist(ion_raw)){
//		Bamboo::s_getInstance()->logger.m_openLog(ion_raw);
//	}
//	if(!Bamboo::s_getInstance()->logger.m_lexist(ion_smooth)){
//		Bamboo::s_getInstance()->logger.m_openLog(ion_smooth);
//	}
//	if(!Bamboo::s_getInstance()->logger.m_lexist(ion_elev)){
//		Bamboo::s_getInstance()->logger.m_openLog(ion_elev);
//	}
//	if(!Bamboo::s_getInstance()->logger.m_lexist(ion_l1)){
//		Bamboo::s_getInstance()->logger.m_openLog(ion_l1);
//	}
//	for(isat = 0;isat < dly->nprn;isat++){
//		isys = index_string(SYS,dly->cprn[isat][0]);
//		if(lset[isat]){
//			double ion_sm = 0.0;
//			atom.smion.m_inquireIon(mjd + sod / 86400.0,isat,line.refsat[isys],ion_sm);
//			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %7.3lf ",ion_raw,atoms[isat] * dly->SAT[isat].freq[1] / dly->SAT[isat].freq[0]);
//			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %7.3lf ",ion_elev,elevs[isat] * RAD2DEG);
//			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %7.3lf ",ion_smooth,ion_sm * dly->SAT[isat].freq[1] / dly->SAT[isat].freq[0]);
//		}else{
//			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %7.3lf ",ion_raw,0.0);
//			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %7.3lf ",ion_elev,0.0);
//			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %7.3lf ",ion_smooth,0.0);
//		}
//		if(ion_dd.ifix[isat]){
//			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %7.3lf ",ion_l1,ion_dd.ddion[isat]);
//		}else{
//			Bamboo::s_getInstance()->logger.m_wtMsg("@%s %7.3lf ",ion_l1,0.0);
//		}
//
//
//	}
//	Bamboo::s_getInstance()->logger.m_wtMsg("@%s \n",ion_raw);
//	Bamboo::s_getInstance()->logger.m_wtMsg("@%s \n",ion_elev);
//	Bamboo::s_getInstance()->logger.m_wtMsg("@%s \n",ion_smooth);
//	Bamboo::s_getInstance()->logger.m_wtMsg("@%s \n",ion_l1);
}
bool RtkBamboo::m_validMW(Line_amb& atom,int isat,int refsat){
	int isit = atom.isit,refsit = atom.refsit,lsuccess,isys;
	int N1 = atom.Raw[isat][0];
	int N2 = atom.Raw[isat][1];
	double abwl0,abwl1,ddabwl,abwl0_ref,abwl1_ref;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	if((atom.ptime[0][isat][1] - atom.ptime[0][isat][0]) * 86400.0 < dly->minseccommon)
		return false;
	if((atom.ptime[1][isat][1] - atom.ptime[1][isat][0]) * 86400.0 < dly->minseccommon)
		return false;

	if((atom.ptime[0][refsat][1] - atom.ptime[0][refsat][0]) * 86400.0 < dly->minseccommon)
		return false;
	if((atom.ptime[1][refsat][1] - atom.ptime[1][refsat][0]) * 86400.0 < dly->minseccommon)
		return false;
	/// smaller than resolution cutoff,return false,make it no hold any more
	if(m_Sta[isit].ob.elev[isat] < dly->cutoff || m_Sta[refsit].ob.elev[isat] < dly->cutoff)
		return false;
	abwl0 = atom.xrwl[0][isat] / atom.weig[0][isat] + atom.abin[0][isat];
	abwl1 = atom.xrwl[1][isat] / atom.weig[1][isat] + atom.abin[1][isat];

	abwl0_ref = atom.xrwl[0][refsat] / atom.weig[0][refsat] + atom.abin[0][refsat];
	abwl1_ref = atom.xrwl[1][refsat] / atom.weig[1][refsat] + atom.abin[1][refsat];
	ddabwl = abwl0 - abwl1 - abwl0_ref + abwl1_ref;
	if(fabs(ddabwl - NINT(ddabwl)) > dly->wlmaxdev)
		return false;
	lsuccess = fabs(N1 - N2 - NINT(ddabwl)) < MAXWND;
	return lsuccess;
}
int RtkBamboo::m_fixamb_if(Line_net& line,Line_amb& m_atom){
}
int RtkBamboo::m_getRefByLinearCmb(Line_net& line,Line_amb& atom,int* idx, double* elevn, int nx) {
	int isat, isys, ix, nfix, ref_in, maxnfix, imax,isit,refsit;
	double mw0, mw1, ddmw, maxelev,abwl0_ref,abwl1_ref,abwl0,abwl1,ddabwl;
	vector<double> elevs_it;
	vector<int> idx_it;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	isit = line.isit;
	refsit = line.refsit;
	for (ix = 0, imax = -1, maxnfix = 0; ix < nx && !strstr(dly->cobs, "SF");
			ix++) {
		ref_in = idx[ix];
		if((atom.ptime[0][ref_in][1] - atom.ptime[0][ref_in][0]) * 86400.0 < dly->minseccommon)
			continue;
		if((atom.ptime[1][ref_in][1] - atom.ptime[1][ref_in][0]) * 86400.0 < dly->minseccommon)
			continue;

		abwl0_ref = atom.xrwl[0][ref_in] / atom.weig[0][ref_in] + atom.abin[0][ref_in];
		abwl1_ref = atom.xrwl[1][ref_in] / atom.weig[1][ref_in] + atom.abin[1][ref_in];

		isys = index_string(SYS, dly->cprn[ref_in][0]);
		for (isat = 0, nfix = 0; isat < dly->nprn; isat++) {
			if (SYS[isys] != dly->cprn[isat][0])
				continue;
			if (m_Sta[isit].ob.omc[isat][0] == 0.0
					|| m_Sta[refsit].ob.omc[isat][0] == 0.0)
				continue;
			if (isat == ref_in)
				continue;
			if((atom.ptime[0][isat][1] - atom.ptime[0][isat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[1][isat][1] - atom.ptime[1][isat][0]) * 86400.0 < dly->minseccommon)
				continue;

			abwl0 = atom.xrwl[0][isat] / atom.weig[0][isat] + atom.abin[0][isat];
			abwl1 = atom.xrwl[1][isat] / atom.weig[1][isat] + atom.abin[1][isat];

			ddabwl = abwl0 - abwl1 - abwl0_ref + abwl1_ref;

			if (fabs(ddmw - NINT(ddmw)) < dly->wlmaxdev) {
				nfix = nfix + 1;
			}
		}
		if (nfix > maxnfix) {
			maxnfix = nfix;
			imax = ref_in;
			idx_it.clear();
			elevs_it.clear();
		}else if(nfix == maxnfix && nfix != 0){
			idx_it.push_back(ref_in);
			elevs_it.push_back(elevn[ix]);
		}
	}
	if (imax == -1) {
		for (ix = 0, maxelev = 0.0; ix < nx; ix++) {
			if (elevn[ix] > maxelev) {
				maxelev = elevn[ix];
				imax = idx[ix];
			}
		}
	}else{
		for(ix = 0,maxelev = 0.0;ix < elevs_it.size();ix++){
			if (elevs_it[ix] > maxelev) {
				maxelev = elevs_it[ix];
				imax = idx_it[ix];
			}
		}
	}
	return imax;
}
void RtkBamboo::m_detectOmc(Line_net& line,int isys,int* flag){
	int nx,ix,isat,itg[MAXSAT],flg[MAXSAT],k,ip,isit,refsit;
	double rx[MAXSAT],wgt[MAXSAT],mean,rms,sig;
	isit = line.isit;
	refsit = line.refsit;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	memset(flag,0,sizeof(int) * MAXSAT);
	for(isat = 0,nx = 0;isat < dly->nprn;isat++){
		if(index_string(SYS,dly->cprn[isat][0]) != isys)
			continue;
		if(m_Sta[isit].ob.omc[isat][0] == 0 || m_Sta[refsit].ob.omc[isat][0] == 0)
			continue;
		if(m_Sta[isit].ob.elev[isat] < 30.0 * DEG2RAD || m_Sta[refsit].ob.elev[isat] < 30.0 * DEG2RAD)
			continue;
		if(!strstr(dly->cobs,"SF")){
			rx[nx] = dly->SAT[isat].fac[0] * m_Sta[isit].ob.omc[isat][0] - dly->SAT[isat].fac[1] * m_Sta[isit].ob.omc[isat][1] - dly->SAT[isat].fac[0] * m_Sta[refsit].ob.omc[isat][0] +
					dly->SAT[isat].fac[1] * m_Sta[refsit].ob.omc[isat][1];
		}else{
			rx[nx] = m_Sta[isit].ob.omc[isat][0] - m_Sta[refsit].ob.omc[isat][0];
		}
		if(!strstr(dly->cobs,"RAW")){
			rx[nx] = rx[nx] - ptrAmtag[isit * MAXSAT + isat][0]->xest;
		}else{
			rx[nx] = rx[nx] - (dly->SAT[isat].lamda[0] * dly->SAT[isat].fac[0] * ptrAmtag[isit * MAXSAT + isat][0]->xest - dly->SAT[isat].lamda[1] * dly->SAT[isat].fac[1] * ptrAmtag[isit * MAXSAT + isat][1]->xest);
		}

		wgt[nx] = 1.0;
		itg[nx] = isat;
		flg[nx] = 0;
		nx = nx + 1;
	}
	Qc::Qc_getWgtMean(true, rx, flg, wgt, nx, &k, &mean, &rms, &sig, 0.10);
	for(ix = 0;ix < nx;ix++){
		if(flg[ix] != 0){
			flag[itg[ix]] = 1;
		}
	}
}
void RtkBamboo::m_pkRefsat(Line_net& line,Line_amb& atom){
	char csys;
	int isat,jsys,isys,flag[MAXSAT],nx,idx[MAXSAT],nsat,nsys,isit,refsit;
	double elevn[MAXSAT],cut = 30.0,maxelev = 0.0;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	isit = line.isit;
	refsit = line.refsit;
	for(jsys = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		memset(flag,0,sizeof(flag));
		m_detectOmc(line,isys,flag);
		for(isat = 0,nx = 0;isat < dly->nprn;isat++){
			if(dly->cprn[isat][0] != SYS[isys])
				continue;
			if(m_Sta[isit].ob.omc[isat][0] == 0.0 || m_Sta[refsit].ob.omc[isat][0] == 0.0)
				continue;
			if(m_Sta[isit].ob.elev[isat] * RAD2DEG >= cut && m_Sta[refsit].ob.elev[isat] * RAD2DEG >= cut && flag[isat] == 0){
				elevn[nx] = MIN(m_Sta[isit].ob.elev[isat],m_Sta[refsit].ob.elev[isat]);
				idx[nx++] = isat;
			}
		}
		if(SYS[isys] != 'R')
			line.refsat[isys] = m_getRefByLinearCmb(line,atom,idx,elevn,nx);
		if(line.refsat[isys] == -1 || SYS[isys] == 'R'){
			maxelev = 0.0;
			for(isat = 0,nx = 0;isat < dly->nprn;isat++){
				if(dly->cprn[isat][0] != SYS[isys])
					continue;
				if(m_Sta[isit].ob.omc[isat][0] == 0.0 || m_Sta[refsit].ob.omc[isat][0] == 0.0)
					continue;
				if(m_Sta[isit].ob.elev[isat] >  maxelev){
					line.refsat[isys] = isat;
					maxelev = m_Sta[isit].ob.elev[isat];
				}
			}
		}
	}
	for(isat = 0,nsat = 0,nsys = 0,csys = ' ';isat < dly->nprn;isat++){
		if(m_Sta[isit].ob.omc[isat][0] != 0 && m_Sta[refsit].ob.omc[isat][0] != 0){
			nsat = nsat + 1;
			if(csys != dly->cprn[isat][0]){
				csys = dly->cprn[isat][0];
				nsys = nsys + 1;
			}
		}
	}
	if(nsat < 3 + nsys){
		for(isys = 0;isys < MAXSYS;isys++)
			line.refsat[isys] = -1;
	}
}
void RtkBamboo::m_addMWCon(double* mat,int ndim,Line_net& line,Line_amb& atom){
	/// decide whether got wl observation here
	int jsys,isys,isit,refsit,nfix,refsat,isat,ifreq,nobs_wl = 0,ip,ind,i;
	isit = atom.isit;
	refsit = atom.refsit;
	vector<SrifPamt>::iterator pItr;
	list<SrifPamt*>::iterator aItr;
	double fac_wl[2],lmw0,lmw1,ddomc,weig,var1,var0,ambcon = 1e3;
	double abwl0,abwl1,ddabwl,abwl0_ref,abwl1_ref,DDamb,atoms[MAXSAT],elevs[MAXSAT];
	double omc0[MAXFREQ],omc1[MAXFREQ],refomc0[MAXFREQ],refomc1[MAXFREQ],ioncor[MAXSAT][MAXFREQ];
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	for(jsys = 0,nfix = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		if(SYS[isys] == 'R')
			continue;
		if(-1 == (refsat = line.refsat[isys]))
			continue;
		if(m_Sta[isit].ob.omc[refsat][0] == 0 || m_Sta[refsit].ob.omc[refsat][0] == 0)
			continue;
		refomc0[0] = m_Sta[refsit].ob.omc[refsat][0];
		refomc0[1] = m_Sta[refsit].ob.omc[refsat][1];
		refomc1[0] = m_Sta[isit].ob.omc[refsat][0];
		refomc1[1] = m_Sta[isit].ob.omc[refsat][1];
		for(isat = 0;isat < dly->nprn;isat++){
			if(dly->cprn[isat][0] != SYS[isys] || isat == refsat)
				continue;
			if((atom.ptime[0][isat][1] - atom.ptime[0][isat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[1][isat][1] - atom.ptime[1][isat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[0][refsat][1] - atom.ptime[0][refsat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[1][refsat][1] - atom.ptime[1][refsat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if(m_Sta[isit].ob.omc[isat][MAXFREQ] == 0.0 || m_Sta[refsit].ob.omc[isat][MAXFREQ] == 0.0)
				continue;
			if(m_Sta[refsit].ob.elev[isat] < ATOM_CUT * DEG2RAD || m_Sta[isit].ob.elev[isat] <  ATOM_CUT * DEG2RAD)
				continue;
			abwl0 = atom.xrwl[0][isat] / atom.weig[0][isat] + atom.abin[0][isat];
			abwl1 = atom.xrwl[1][isat] / atom.weig[1][isat] + atom.abin[1][isat];

			abwl0_ref = atom.xrwl[0][refsat] / atom.weig[0][refsat] + atom.abin[0][refsat];
			abwl1_ref = atom.xrwl[1][refsat] / atom.weig[1][refsat] + atom.abin[1][refsat];
			ddabwl = abwl0 - abwl1 - abwl0_ref + abwl1_ref;
			if(fabs(ddabwl - NINT(ddabwl)) > dly->wlmaxdev)
				continue;
			atom.ifab[refsat][0] = AMB_WLFIX;
			atom.abwl[refsat] = 0;

			atom.ifab[isat][0] = AMB_WLFIX;
			atom.abwl[isat] = DDamb = NINT(ddabwl);
			omc0[0] = m_Sta[refsit].ob.omc[isat][0];
			omc0[1] = m_Sta[refsit].ob.omc[isat][1];
			omc1[0] = m_Sta[isit].ob.omc[isat][0];
			omc1[1] = m_Sta[isit].ob.omc[isat][1];

			fac_wl[0] = dly->SAT[isat].freq[0] / (dly->SAT[isat].freq[0] - dly->SAT[isat].freq[1]);
			fac_wl[1] = dly->SAT[isat].freq[1] / (dly->SAT[isat].freq[0] - dly->SAT[isat].freq[1]);

			lmw0=fac_wl[0] * (omc0[0] - refomc0[0]) - fac_wl[1] * (omc0[1] - refomc0[1]);
			lmw1=fac_wl[0] * (omc1[0] - refomc1[0]) - fac_wl[1] * (omc1[1] - refomc1[1]);

			ddomc = lmw1 - lmw0 - DDamb * dly->SAT[isat].lamdw; //

			var1 = pow(fac_wl[0],2) * (m_Sta[isit].ob.var[isat][0] + m_Sta[isit].ob.var[refsat][0]) +
					pow(fac_wl[1],2) * (m_Sta[isit].ob.var[isat][1] + m_Sta[isit].ob.var[refsat][1]);

			var0 = pow(fac_wl[0],2) * (m_Sta[refsit].ob.var[isat][0] + m_Sta[refsit].ob.var[refsat][0]) +
					pow(fac_wl[1],2) * (m_Sta[refsit].ob.var[isat][1] + m_Sta[refsit].ob.var[refsat][1]);
			weig = 1.0 / sqrt(var1 + var0);

			for (ip = 0; ip < line.n_np[isat]; ip++) {
				ind = line.i_coef[isat][ip][0];
				if(ind < nps){
					if(!strstr(m_Pm[pmtag[ind]].pname,"ION")){
						mat[ndim * ind + imtx + nobs_wl] = line.d_coef[isat][ip][0] * weig;
					}
					else{
						mat[ndim * ind + imtx + nobs_wl] = -dly->SAT[isat].freq[0] / dly->SAT[isat].freq[1] * line.d_coef[isat][ip][0] * weig;
					}
				}
			}
			mat[ndim * imtx + imtx + nobs_wl] = ddomc * weig;
			nobs_wl = nobs_wl + 1;
			/// add ambiguities contraint
			mat[ndim * amtag[isit * MAXSAT + isat][0] + imtx + nobs_wl] = 1.0 * ambcon;
			mat[ndim * amtag[isit * MAXSAT + isat][1] + imtx + nobs_wl] = -1.0 * ambcon;

			mat[ndim * amtag[isit * MAXSAT + refsat][0] + imtx + nobs_wl] = -1.0 * ambcon;
			mat[ndim * amtag[isit * MAXSAT + refsat][1] + imtx + nobs_wl] = 1.0 * ambcon;

			mat[ndim * imtx + imtx + nobs_wl] = (DDamb - ptrAmtag[isit * MAXSAT + isat][0]->xini + ptrAmtag[isit * MAXSAT + isat][1]->xini +
					ptrAmtag[isit * MAXSAT + refsat][0]->xini - ptrAmtag[isit * MAXSAT + refsat][1]->xini) * ambcon;
			nobs_wl = nobs_wl + 1;
		}
	}
	if(nobs_wl > 0){
		CMat::CMat_Householder(mat, ndim, imtx, imtx + nobs_wl, imtx + 1, false);
		double* est = (double*) calloc(imtx, sizeof(double));
		memcpy(est, mat + ndim * imtx, sizeof(double) * imtx);
		CMat::CMat_SolveLinear(imtx, mat, ndim, est, imtx);
		//LOGSTAT(true,"sod = %lf",sod);
		for (i = 0, pItr = m_Pm.begin(), aItr = m_Am.begin(); i < imtx; i++) {
			if (i < nps) {
				ip = this->pmtag[i];
				(*pItr).xsig = 0.0;
				(*pItr).xest = (*pItr).xini + est[i];
				//LOGSTAT(true,"%03d %04d %16.3lf %16.3lf %16.3lf",i,(*pItr).iobs,(*pItr).xini,(*pItr).xcor,(*pItr).xest);
				++pItr;
			} else {
				(*aItr)->xsig = 0.0;
				(*aItr)->xest = (*aItr)->xini + est[i];
				//LOGSTAT(true,"%03d %04d %16.3lf %16.3lf %16.3lf",i,(*aItr)->iobs,(*aItr)->xini,(*aItr)->xcor,(*aItr)->xest);
				++aItr;
			}
		}
		free(est);
	}
}
void RtkBamboo::m_addMWCon_(double* mat,int ndim,Line_net& line,Line_amb& atom){
	/// 1) using MW observation equation to find whether the residual is ok
	int jsys,isys,isit,refsit,nfix,refsat,isat,ifreq,nobs_wl = 0,ip,ind,i,ipt[MAXSAT],ldebug = true;
	double wlresi[MAXSAT],wlresi_all[MAXSAT] = {0};
	isit = atom.isit;
	refsit = atom.refsit;
	vector<SrifPamt>::iterator pItr;
	list<SrifPamt*>::iterator aItr;
	double fac_wl[2],lmw0,lmw1,ddomc,weig,var1,var0,ambcon = 1e3;
	double abwl0,abwl1,ddabwl,abwl0_ref,abwl1_ref,DDamb,atoms[MAXSAT],elevs[MAXSAT];
	double omc0[MAXFREQ],omc1[MAXFREQ],refomc0[MAXFREQ],refomc1[MAXFREQ],ioncor[MAXSAT][MAXFREQ];
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	for(jsys = 0,nfix = 0,nobs_wl = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		if(SYS[isys] == 'R')
			continue;
		if(-1 == (refsat = line.refsat[isys]))
			continue;
		if(m_Sta[isit].ob.omc[refsat][0] == 0 || m_Sta[refsit].ob.omc[refsat][0] == 0)
			continue;
		refomc0[0] = m_Sta[refsit].ob.omc[refsat][0];
		refomc0[1] = m_Sta[refsit].ob.omc[refsat][1];
		refomc1[0] = m_Sta[isit].ob.omc[refsat][0];
		refomc1[1] = m_Sta[isit].ob.omc[refsat][1];
		for(isat = 0;isat < dly->nprn;isat++){
			if(dly->cprn[isat][0] != SYS[isys] || isat == refsat)
				continue;
			if((atom.ptime[0][isat][1] - atom.ptime[0][isat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[1][isat][1] - atom.ptime[1][isat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[0][refsat][1] - atom.ptime[0][refsat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if((atom.ptime[1][refsat][1] - atom.ptime[1][refsat][0]) * 86400.0 < dly->minseccommon)
				continue;
			if(m_Sta[isit].ob.omc[isat][MAXFREQ] == 0.0 || m_Sta[refsit].ob.omc[isat][MAXFREQ] == 0.0)
				continue;
			if(m_Sta[refsit].ob.elev[isat] < ATOM_CUT * DEG2RAD || m_Sta[isit].ob.elev[isat] <  ATOM_CUT * DEG2RAD)
				continue;
			abwl0 = atom.xrwl[0][isat] / atom.weig[0][isat] + atom.abin[0][isat];
			abwl1 = atom.xrwl[1][isat] / atom.weig[1][isat] + atom.abin[1][isat];

			abwl0_ref = atom.xrwl[0][refsat] / atom.weig[0][refsat] + atom.abin[0][refsat];
			abwl1_ref = atom.xrwl[1][refsat] / atom.weig[1][refsat] + atom.abin[1][refsat];
			ddabwl = abwl0 - abwl1 - abwl0_ref + abwl1_ref;
			if(fabs(ddabwl - NINT(ddabwl)) > dly->wlmaxdev)
				continue;
			ipt[nobs_wl] = isat;
			atom.ifab[isat][0] = AMB_WLFIX;
			atom.abwl[isat] = DDamb = NINT(ddabwl);
			DDamb = atom.abwl[isat];
			/// add ambiguities contraint
			mat[ndim * amtag[isit * MAXSAT + isat][0] + imtx + nobs_wl] = 1.0 * ambcon;
			mat[ndim * amtag[isit * MAXSAT + isat][1] + imtx + nobs_wl] = -1.0 * ambcon;

			mat[ndim * amtag[isit * MAXSAT + refsat][0] + imtx + nobs_wl] = -1.0 * ambcon;
			mat[ndim * amtag[isit * MAXSAT + refsat][1] + imtx + nobs_wl] = 1.0 * ambcon;

			mat[ndim * imtx + imtx + nobs_wl] = (DDamb - ptrAmtag[isit * MAXSAT + isat][0]->xini + ptrAmtag[isit * MAXSAT + isat][1]->xini +
					ptrAmtag[isit * MAXSAT + refsat][0]->xini - ptrAmtag[isit * MAXSAT + refsat][1]->xini) * ambcon;
			nobs_wl = nobs_wl + 1;
		}
	}
	if(nobs_wl == 0)
		return;
	double* est = (double*) calloc(imtx, sizeof(double));
	/// compute the residual here
	CMat::CMat_Householder(mat, ndim, imtx, imtx + nobs_wl, imtx + 1, false);
	memcpy(est, mat + ndim * imtx, sizeof(double) * imtx);
	CMat::CMat_SolveLinear(imtx, mat, ndim, est, imtx);
	for (i = 0, pItr = m_Pm.begin(), aItr = m_Am.begin(); i < imtx; i++) {
		if (i < nps) {
			ip = this->pmtag[i];
			(*pItr).xsig = 0.0;
			(*pItr).xest = (*pItr).xini + est[i];
			//LOGSTAT(true,"%03d %04d %16.3lf %16.3lf %16.3lf",i,(*pItr).iobs,(*pItr).xini,(*pItr).xcor,(*pItr).xest);
			++pItr;
		} else {
			(*aItr)->xsig = 0.0;
			(*aItr)->xest = (*aItr)->xini + est[i];
			//LOGSTAT(true,"%03d %04d %16.3lf %16.3lf %16.3lf",i,(*aItr)->iobs,(*aItr)->xini,(*aItr)->xcor,(*aItr)->xest);
			++aItr;
		}
	}
	free(est);
}
void RtkBamboo::m_getAtomStatus(Line_net& line,int psat,int ref,double& ddion,double& ddtrp){
	char ztdtag_1[1024],ztdtag_0[1024],iontag_1_ref[1024],iontag_1[1024];
	int isit,refsit,ip;
	double ztdest1,ztdest0,ion1,ion1_ref;
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	isit = line.isit;
	refsit = line.refsit;
	sprintf(ztdtag_1,"%s_%03d",dly->ztdmod,isit); ///
	sprintf(ztdtag_0,"%s_%03d",dly->ztdmod,refsit);  ///
	sprintf(iontag_1_ref,"ION_%03d_%03d",isit,ref);
	sprintf(iontag_1,"ION_%03d_%03d",isit,psat);
	for(ip = 0; ip < nps;ip++){
		if(strstr(m_Pm[ip].pname,ztdtag_1)){
			ztdest1 = m_Pm[ip].xest;
		}
		if(strstr(m_Pm[ip].pname,ztdtag_0)){
			ztdest0 = m_Pm[ip].xest;
		}
		if(strstr(m_Pm[ip].pname,iontag_1)){
			ion1 = m_Pm[ip].xest;
		}
		if(strstr(m_Pm[ip].pname,iontag_1_ref)){
			ion1_ref = m_Pm[ip].xest;
		}
	}
	ddion = ion1 * m_Sta[isit].ob.ionpart[psat] - ion1_ref * m_Sta[isit].ob.ionpart[ref];
	ddtrp = (m_Sta[isit].ob.trpModel[psat] + ztdest1 * m_Sta[isit].ob.zmap[psat]) - (m_Sta[isit].ob.trpModel[ref] + ztdest1 * m_Sta[isit].ob.zmap[ref]) -
			(m_Sta[refsit].ob.trpModel[psat] + ztdest0 * m_Sta[refsit].ob.zmap[psat]) + (m_Sta[refsit].ob.trpModel[ref] + ztdest0 * m_Sta[refsit].ob.zmap[ref]);
}
void RtkBamboo::m_calAtom_pam(Line_net& line,Line_amb& atom){
	/// use the ion parameters in the filter to service,find out whether it is ok
	/// calculate the atomsphere using fixed ambiguities
	RtkDDIon ion,trp;
	double maxiondif = 0.08;
	int i,isys,jsys,refsat,ifreq,isit,refsit,isat,oN1,oN2,ledit,lfixed,subON1,subON2,status,nfix,ip;
	char ztdtag_1[256],ztdtag_0[256];
	double obs0,obs1,obsref0,obsref1,dion,N1,N2,obsif0,obsif1,obsifref0,obsifref1,ddtrp,ddion,ddion_wl,ddion1,ddion2;
	double ztdcor_1,ztdcor_0,resd;
	loginfo.m_reset();
	Deploy* dly = &Bamboo::s_getInstance()->dly;
	isit = atom.isit;
	refsit = atom.refsit;
	Line_net cline;
	cline.isit = line.isit;
	cline.refsit = line.refsit;
	atom.m_inquireRef(this,cline);
	m_fillAmat(cline);
	m_genRMatrix(cline);
	memcpy(cline.omc,line.omc,sizeof(line.omc));
	for(jsys = 0,nfix = 0;jsys < dly->nsys;jsys++){
		isys = index_string(SYS,dly->system[jsys]);
		if(-1 == (refsat = cline.refsat[isys]))
			continue;
		// no observation for reference satellite
		if(m_Sta[isit].ob.omc[refsat][0] == 0 || m_Sta[refsit].ob.omc[refsat][0] == 0)
			continue;
		obsref0 = m_Sta[refsit].ob.obs[refsat][0] * dly->SAT[refsat].lamda[0] - m_Sta[refsit].ob.obs[refsat][1] * dly->SAT[refsat].lamda[1];
		obsref1 = m_Sta[isit].ob.obs[refsat][0] * dly->SAT[refsat].lamda[0] - m_Sta[isit].ob.obs[refsat][1] * dly->SAT[refsat].lamda[1];

		obsifref0 = m_Sta[refsit].ob.obs[refsat][0] * dly->SAT[refsat].lamda[0] * dly->SAT[refsat].fac[0] - m_Sta[refsit].ob.obs[refsat][1] * dly->SAT[refsat].lamda[1] * dly->SAT[refsat].fac[1];
		obsifref1 = m_Sta[isit].ob.obs[refsat][0] * dly->SAT[refsat].lamda[0] * dly->SAT[refsat].fac[0] - m_Sta[isit].ob.obs[refsat][1] * dly->SAT[refsat].lamda[1] * dly->SAT[refsat].fac[1];

		for(isat = 0;isat < dly->nprn;isat++){
			if(dly->cprn[isat][0] != SYS[isys])
				continue;
			// no observation
			if(m_Sta[isit].ob.omc[isat][0] == 0 || m_Sta[refsit].ob.omc[isat][0] == 0)
				continue;
			if(m_Sta[isit].ob.elev[isat] < m_Sta[isit].cutoff || m_Sta[refsit].ob.elev[isat] < m_Sta[refsit].cutoff)
				continue;
			/// 10 degrees for cut,only output 10 degrees
			if(m_Sta[isit].ob.elev[isat] < ATOM_CUT * DEG2RAD || m_Sta[refsit].ob.elev[isat] < ATOM_CUT * DEG2RAD)
				continue;
			// not fixed in the memory,if the wl is fixed,can also use the ion information
			if(atom.ifab[isat][0] < AMB_WLFIX && SYS[isys] != 'R')
				continue;
			// not fixed in the memory for glonass and the common time is smaller than the min-common time,will continue for this case
			if(atom.ifab[isat][0] < AMB_WLFIX && ((m_atom.ptime[0][isat][1] - m_atom.ptime[0][isat][0]) * 86400.0) < dly->minseccommon)
				continue;
			lfixed = false;
			if(atom.ifab[isat][0] >= AMB_FIX && atom.ifab[isat][1] >= AMB_FIX){
				lfixed = true;
				N1 = atom.Raw[isat][0] - atom.Raw[refsat][0]; //  single difference
				N2 = atom.Raw[isat][1] - atom.Raw[refsat][1]; //  single difference

				loginfo.Raw[isat][0] = N1;
				loginfo.Raw[isat][1] = N2;
			}
			nfix = nfix + 1;
			m_getAtomStatus(line,isat,refsat,ddion,ddtrp);
			/// update structure
			ion.elev[isat] = MIN(m_Sta[isit].ob.elev[isat],m_Sta[refsit].ob.elev[isat]);
			ion.ddion[isat] = ddion;
			ion.ifix[isat] = AMB_FIX;
			ion.refsat[isys] = refsat;
			ion.ptime = mjd + sod / 86400.0;

			trp.elev[isat] = MIN(m_Sta[isit].ob.elev[isat],m_Sta[refsit].ob.elev[isat]);
			trp.ddion[isat] = ddtrp;
			trp.ifix[isat] = AMB_FIX;
			trp.refsat[isys] = refsat;
			trp.ptime = mjd + sod / 86400.0;

			/// update logging part
			loginfo.ifab[isat] = atom.ifab[isat][0];
			loginfo.pifab[isat] = atom.pifab[isat][0];

			loginfo.elev[isat] = m_Sta[isit].ob.elev[isat];
			loginfo.ion[isat] = ion.ddion[isat];
			loginfo.trp[isat] = trp.ddion[isat];
			loginfo.refsat[isys] = refsat;
			if (lfixed) {
				loginfo.useRaw[isat][0] = N1;
				loginfo.useRaw[isat][1] = N2;
				loginfo.status[isat] = 1; // valid success
			} else {
				loginfo.useRaw[isat][0] = -1;
				loginfo.useRaw[isat][1] = -1;
				loginfo.status[isat] = 2; // valid failed,using wl-ion
			}
			if (lfixed && (loginfo.useRaw[isat][0] != loginfo.Raw[isat][0]
					|| loginfo.useRaw[isat][1] != loginfo.Raw[isat][1])) {
				loginfo.idif[isat] = 1;
			}
			if (loginfo.status[isat] != 2) {
				m_validIFAmb(false, cline, atom, line.xa, isat,
						loginfo.useRaw[isat][0], loginfo.useRaw[isat][1],
						subON1, subON2, loginfo.resiIF_new[isat]);
				m_validAmb(false, cline, atom, line.xa, isat, 0,
						loginfo.useRaw[isat][0], oN1, 2.8,
						loginfo.resiRaw_new[isat][0]);
				m_validAmb(false, cline, atom, line.xa, isat, 1,
						loginfo.useRaw[isat][1], oN2, 2.8,
						loginfo.resiRaw_new[isat][1]);
			}
			if (loginfo.idif[isat]) {
				m_validIFAmb(false, cline, atom, line.xa, isat,
						loginfo.Raw[isat][0], loginfo.Raw[isat][1], subON1,
						subON2, loginfo.resiIF[isat]);
				m_validAmb(false, cline, atom, line.xa, isat, 0,
						loginfo.Raw[isat][0], oN1, 2.8,
						loginfo.resiRaw[isat][0]);
				m_validAmb(false, cline, atom, line.xa, isat, 1,
						loginfo.Raw[isat][1], oN2, 2.8,
						loginfo.resiRaw[isat][1]);

				loginfo.ion_new[isat] = dly->SAT[isat].fac[1]
						* (obs1 - obs0 - obsref1 + obsref0
								- (dly->SAT[isat].lamda[0]
										* loginfo.Raw[isat][0]
										- dly->SAT[isat].lamda[1]
												* loginfo.Raw[isat][1]));
				loginfo.trp_new[isat] = obsif1 - obsif0 - obsifref1 + obsifref0
						- loginfo.Raw[isat][0] * dly->SAT[isat].lamda[0]
								* dly->SAT[isat].fac[0]
						+ loginfo.Raw[isat][1] * dly->SAT[isat].lamda[1]
								* dly->SAT[isat].fac[1]
						- (m_Sta[isit].ob.rleng[isat]
								- m_Sta[isit].ob.rleng[refsat]
								- m_Sta[refsit].ob.rleng[isat]
								+ m_Sta[refsit].ob.rleng[refsat]);
			}
		}
		for (isat = 0; isat < dly->nprn; isat++) {
			if (dly->cprn[isat][0] != SYS[isys] || SYS[isys] == 'R')
				continue;
			// no observation
			if (m_Sta[isit].ob.omc[isat][0] == 0
					|| m_Sta[refsit].ob.omc[isat][0] == 0)
				continue;
			if (m_Sta[isit].ob.elev[isat] < m_Sta[isit].cutoff
					|| m_Sta[refsit].ob.elev[isat] < m_Sta[refsit].cutoff)
				continue;
			///////////////////////////// update logging part /////////////////////////
			if (atom.smion.m_inquireIon(mjd + sod / 86400.0, isat, refsat,
					ddion_wl)) {
				loginfo.ion_wl[isat] = ddion_wl * dly->SAT[isat].freq[1]
						/ dly->SAT[isat].freq[0];
			}
		}
	}
	atom.ddblion.m_setData(ion);
	atom.ddbltrp.m_setData(trp);
}
/// smooth the ion
