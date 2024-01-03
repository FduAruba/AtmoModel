#include "ztdmodel.h"

void ZtdModel::setBasicOption(IN ProOption& opt, IN int res)
{
	_useres = res;
	_fittype = opt._troptype;
	switch (_fittype)
	{
	case MODEL_OFC:		{ _ncoeff = 10; break; }
	case MODEL_OFC_MSL: { _ncoeff =  6; break; }
	case MODEL_MOFC:    { _ncoeff =  7; break; }
	default:            { _ncoeff =  0; break; }
	}
}

void ZtdModel::settime(IN const Gtime t)
{
	this->_tnow = t;
	if (timediff(t, this->_ztdModCur._time) > 360.0) {
		this->_ztdModCur._valid = false;
	}
}

bool ZtdModel::ztdModSaas(IN double* pos, IN double* azel, IN double humi, OUT double* ztd)
{
	const double temp0 = 15.0;	/* temparature at sea level */
	double hgt = 0, pres = 0, temp = 0, e = 0, z = 0, zhd = 0, zwd = 0;

	if (pos[2] < -100.0 || 1.0 / 2.2557E-5 <= pos[2] || azel[1] <= 0) {
		return false;
	}
	hgt = pos[2] < 0.0 ? 0.0 : pos[2];

	/* standard atmosphere */
	pres = 1013.25 * pow(1.0 - 2.2557E-5 * hgt, 5.2568);
	temp = temp0 - 6.5E-3 * hgt + 273.16;
	e = 6.108 * humi * exp((17.15 * temp - 4684.0) / (temp - 38.45));
	/* saastamoninen model */
	z = PI / 2.0 - azel[1];
	zhd = 0.0022768 * pres / (1.0 - 0.00266 * cos(2.0 * pos[0]) - 0.00028 * hgt / 1E3) / cos(z);
	zwd = 0.002277 * (1255.0 / temp + 0.05) * e / cos(z);

	ztd[0] = zhd; ztd[1] = zwd;
	return true;
}

void ZtdModel::ztd2msl(IN int bmean, IO ZtdInfos& ztds)
{
	Dvec staztd;
	double azel[2] = { 0.0, PI / 2.0 };
	double iztd[2] = { 0.0 }, pos[3] = { 0 };
	double vmea = 0, vsig = 0;
	double ztd0 = 0, zhd0 = 0, scal = 0;
	bool loop = true;

	/* 1.计算区域zhd0均值zhdm */
	if (bmean) {
		for (auto& iSta : ztds) {
			memset(iztd, 0, sizeof(iztd));
			memcpy(pos, iSta.second._blh, sizeof(pos));
			pos[2] = 0.0;
			if (ztdModSaas(pos, azel, REL_HUMI, iztd)) {
				staztd.push_back(iztd[0]);
			}
		}

		while (loop) {
			loop = false;
			calcMeanStd(staztd, vmea, vsig);
			for (auto it = staztd.begin(); it != staztd.end();) {
				if (fabs(*it - vmea) > 2.0 * vsig) {
					it = staztd.erase(it);
					loop = true;
					continue;
				}
				++it;
			}
		}
		zhd0 = vmea;
	}
	
	/* 2.计算各站ztd0、zwd0 */
	for (auto& iSta : ztds) {
		scal = exp(-1.3137e-4 * iSta.second._blh[2]);
		ztd0 = (iSta.second._zhd + iSta.second._zwd) / scal;

		if (!bmean) {
			memset(iztd, 0, sizeof(iztd));
			memcpy(pos, iSta.second._blh, sizeof(pos));
			pos[2] = 0.0;
			if (ztdModSaas(pos, azel, REL_HUMI, iztd)) {
				zhd0 = iztd[0];
			}
		}

		iSta.second._zhd0 = zhd0;
		iSta.second._zwd0 = ztd0 - zhd0;
		//printf("%s %10.6f %8.4f\n", iSta.first.c_str(), zhd0, iSta.second._zwd0);
	}

	return;
}

VecXd ZtdModel::buildRowA(IN int nx, IN int np, IN double* blh)
{
	VecXd xrow = VecXd::Zero(nx);
	int maxorder = (int)(log(nx) / log(2));

	switch (np)
	{
	case 2: {
		for (int iord = 0, id = 0; iord <= maxorder && id < nx; ++iord) {
			polynomial2(blh[0], blh[1], iord, nx, id, xrow);
		}
		break;
	}
	case 3: {
		for (int iord = 0, id = 0; iord <= maxorder && id < nx; ++iord) {
			polynomial3(blh[0], blh[1], blh[2], iord, nx, id, xrow);
		}
		break;
	}
	default:
		break;
	}
	//cout << xrow.transpose() << endl;

	return xrow;
}

int ZtdModel::bulidMatrix(IN int nx, IN ZtdInfos& ztds, IN double* c, IN unordered_set<string> ext,
	OUT VecXd* y, OUT VecXd* P, OUT MatXd* A, OUT map<int, string>& idx)
{
	int ny = 0;
	int np = nx <= 7 ? 2 : 3;
	double blh[3] = { 0.0 };
	VecXd arow;
	
	idx.clear();

	for (auto& iSta : ztds) {
		string sta = iSta.first;
		if (ext.find(sta) != ext.end()) {
			continue;
		}
		
		idx.emplace(ny, sta);
		blh[0] = iSta.second._blh[0] - c[0];
		blh[0] *= R2D;
		blh[1] = iSta.second._blh[1] - c[1];
		blh[1] *= R2D;
		blh[2] = iSta.second._blh[2] * 1.0E-3;

		arow = buildRowA(nx, np, blh);
		A->row(ny) = arow.transpose();
		(*y)[ny] = iSta.second._zwd0;
		(*P)[ny] = 1.0E-3 / iSta.second._std_zwd;
		//(*P)[ny] = 1;

		ny++;
	}
	y->conservativeResize(ny);
	A->conservativeResize(ny, nx);
	P->conservativeResize(ny);
	//cout << (*A) << endl;
	return ny;
}

double ZtdModel::estres(int nx, int ny, VecXd* x, VecXd* y, MatXd* A, VecXd* v)
{
	double sig = 0;
	(*v) = (*y) - (*A) * (*x);
	sig = (*v).dot(*v);
	sig = sqrt(sig / (ny - nx));

	return sig;
}

bool ZtdModel::ofcModel(IN Gtime t, IN ZtdInfos& ztds, IN GridInfo& grids, IN int bmean)
{
	int nx = _ncoeff, ny = (int)ztds.size();
	double c[3] = { 0.0 }, sig = 0;
	double maxres = 0;
	int imax = -1, iter = 0;
	map<int, string> idx;
	unordered_set<string> ext;
	bool stat = true;

	if (ny <= nx) { return false; }
	c[0] = grids._center[0] * D2R;
	c[1] = grids._center[1] * D2R;

	/* eigen 初始化 */
	VecXd*  x = new VecXd(nx);		x->setZero();	// 待估系数
	VecXd* dx = new VecXd(nx);	   dx->setZero();	// 迭代状态
	VecXd*  y = new VecXd(ny);		y->setZero();	// 观测量
	VecXd*  v = new VecXd(ny);		v->setZero();	// 观测残差
	VecXd*  P = new VecXd(ny);		P->setZero();	// 观测协方差
	MatXd*  A = new MatXd(ny, nx);	A->setZero();	// 观测矩阵
	MatXd*  Q = new MatXd(nx, nx);  Q->setZero();	// 系数协方差

	/* 循环迭代估计系数 */
	for (iter = 0; iter < 10; iter++) {
		/* 矩阵构建 */
		ny = bulidMatrix(nx, ztds, c, ext, y, P, A, idx);
		if (ny <= nx) {
			stat = false; 
			break;
		}

		/* 先验残差估计 */
		sig = estres(nx, ny, x, y, A, v);

		/* LSQ估计 (*/
		(*dx) = LSQ(*A, *v, *P, Q);
		(*x) += (*dx);

		/* 后验残差剔除 */
		sig = estres(nx, ny, x, y, A, v);
		maxres = v->cwiseAbs().maxCoeff(&imax);
		
		if (fabs(maxres) > 2.5 * sig) {
			ext.emplace(idx[imax]);
			continue;
		}
		if (sqrt((*dx).dot(*dx)) < 1.0E-4) {
			break;
		}
	}

	if (!stat || iter >= 10) {
		printf("out iter %6.4f!\n", sqrt((*dx).dot(*dx)));
		x->resize(0); dx->resize(0); y->resize(0); v->resize(0); P->resize(0);
		A->resize(0, 0); Q->resize(0, 0);
		delete x; delete dx; delete y; delete v; delete P; delete A; delete Q;
		return false;
	}

	/* 记录建模系数 */
	_ztdModCur._time = t;
	_ztdModCur._qi = sig;
	_ztdModCur._nsta = ny;
	_ztdModCur._zhd = bmean ? ztds.begin()->second._zhd0 : 0.0;
	for (int i = 0; i < _ncoeff; i++) {
		_ztdModCur._coeff[i] = (*x)(i);
		_ztdModCur._coeff_rms[i] = sqrt((*Q)(i, i));
	}
	_ztdModCur._valid = true;

	x->resize(0); dx->resize(0); y->resize(0); v->resize(0); P->resize(0);
	A->resize(0, 0); Q->resize(0, 0);
	delete x; delete dx; delete y; delete v; delete P; delete A; delete Q;
	return true;
}