/*
 * TransitModel.cxx
 *
 *  Created on: 31 Jan 2014
 *      Author: futyand
 */

#include <iostream>
#include <cmath>

#include "TransitModel.hxx"

using namespace std;

TransitModel::TransitModel(double radiusRatio, bool doLimbDarkening, pair<double,double> limbDarkeningCoeffs) {
	m_doLimbDarkening = doLimbDarkening;
	m_radiusRatio=radiusRatio;
	m_limbDarkeningCoeffs = limbDarkeningCoeffs;
}

//////////////////////////////////////////////////////////////////////////////////
/// This function implements the formulae of Mandel and Agol given in
/// http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/ for the case of quadratic limb
/// darkening.  The code is a conversion to C++ of C code written by Laura Kreidberg, available
/// from Eric Agol's website: http://www.astro.washington.edu/users/agol/transit.html
double TransitModel::fluxFactor(double starRadiusFraction) const {

	double p = m_radiusRatio;
	double z = starRadiusFraction;
	double u1 = m_limbDarkeningCoeffs.first;
	double u2 = m_limbDarkeningCoeffs.second;

	double lambdad, etad, lambdae, pi, x1, x2, x3, omega;
	double kap0 = 0.0, kap1 = 0.0, q, Kk, Ek, Pk, n;
	double muo1, mu0;

	if(fabs(p - 0.5) < 1.0e-3) p = 0.5;

	omega=1.0-u1/3.0-u2/6.0;
	pi=M_PI;

	//	#pragma omp parallel for private(z, x1,x2,x3,n,q,Kk,Ek,Pk,kap0,kap1)

	x1 = pow((p - z), 2.0);
	x2 = pow((p + z), 2.0);
	x3 = p*p - z*z;

	//source is unocculted:
	if(z > 1.0 + p)
	{
		//cout << "zone 1\n" << endl;
		lambdad = 0.0;
		etad = 0.0;
		lambdae = 0.0;
		muo1 = 1.0-((1.0-u1-2.0*u2)*lambdae+(u1+2.0*u2)*lambdad+u2*etad)/omega;
		mu0 = 1.0-lambdae;
		return m_doLimbDarkening ? muo1 : mu0;
	}
	//source is completely occulted:
	if(p > 1.0 && z < p - 1.0)
	{
		//cout << "zone 2\n" << endl;
		lambdad = 0.0;
		etad = 0.5;		//error in fortran code corrected here, following Eastman's python code
		lambdae = 1.0;
		muo1 = 1.0-((1.0-u1-2.0*u2)*lambdae+(u1+2.0*u2)*(lambdad + 2.0/3.0)+u2*etad)/omega;
		mu0 = 1.0-lambdae;
		return m_doLimbDarkening ? muo1 : mu0;
	}
	//source is partly occulted and occulting object crosses the limb:
	if(z > fabs(1.0 - p) && z < 1.0 + p)
	{
		//cout << "zone 3\n" << endl;
		kap1 = acos(min((1.0-p*p+z*z)/2.0/z,1.0));
		kap0 = acos(min((p*p+z*z-1.0)/2.0/p/z,1.0));
		lambdae = p*p*kap0+kap1;
		lambdae = (lambdae - 0.50*sqrt(max(4.0*z*z-pow((1.0+z*z-p*p), 2.0),0.0)))/pi;
	}
	//occulting object transits the source but doesn't completely cover it:
	if(z < 1.0 - p)
	{
		//cout << "zone 4\n" << endl;
		lambdae = p*p;
	}
	//edge of the occulting star lies at the origin
	if(fabs(z-p) < 1.0e-4*(z+p))
	{
		//cout << "zone 5\n" << endl;
		if(z > 0.5)
		{
			//cout << "zone 5.1\n" << endl;
			q = 0.5/p;
			Kk = ellk(q);
			Ek = ellec(q);
			lambdad = 1.0/3.0+16.0*p/9.0/pi*(2.0*p*p-1.0)*Ek-
					(32.0*pow(p, 4.0)-20.0*p*p+3.0)/9.0/pi/p*Kk;
			etad = 1.0/2.0/pi*(kap1+p*p*(p*p+2.0*z*z)*kap0-
					(1.0+5.0*p*p+z*z)/4.0*sqrt((1.0-x1)*(x2-1.0)));
			if(p == 0.5)
			{
				lambdad = 1.0/3.0-4.0/pi/9.0;
				etad = 3.0/32.0;
			}
			muo1 = 1.0-((1.0-u1-2.0*u2)*lambdae+(u1+2.0*u2)*lambdad+u2*etad)/omega;
			mu0 = 1.0-lambdae;
			return m_doLimbDarkening ? muo1 : mu0;
		}
		else
		{
			//cout << "zone 5.2\n" << endl;
			q = 2.0*p;
			Kk = ellk(q);
			Ek = ellec(q);
			lambdad = 1.0/3.0+2.0/9.0/pi*(4.0*(2.0*p*p-1.0)*Ek + (1.0-4.0*p*p)*Kk);
			etad = p*p/2.0*(p*p+2.0*z*z);
			muo1 = 1.0-((1.0-u1-2.0*u2)*lambdae+(u1+2.0*u2)*lambdad+u2*etad)/omega;
			mu0 = 1.0-lambdae;
			return m_doLimbDarkening ? muo1 : mu0;
		}
	}
	//occulting star partly occults the source and crosses the limb:
	if((z > 0.5 + fabs(p -0.5) && z < 1.0 + p) || (p > 0.5 && z > fabs(1.0-p)*1.0001
			&& z < p))
	{
		//cout << "zone 6\n" << endl;
		q = sqrt((1.0-pow((p-z), 2.0))/4.0/z/p);
		Kk = ellk(q);
		Ek = ellec(q);
		n = 1.0/x1-1.0;
		Pk = Kk-n/3.0*rj(0.0,1.0-q*q,1.0,1.0+n);
		lambdad = 1.0/9.0/pi/sqrt(p*z)*(((1.0-x2)*(2.0*x2+
				x1-3.0)-3.0*x3*(x2-2.0))*Kk+4.0*p*z*(z*z+
						7.0*p*p-4.0)*Ek-3.0*x3/x1*Pk);
		if(z < p) lambdad += 2.0/3.0;
		etad = 1.0/2.0/pi*(kap1+p*p*(p*p+2.0*z*z)*kap0-
				(1.0+5.0*p*p+z*z)/4.0*sqrt((1.0-x1)*(x2-1.0)));
		muo1 = 1.0-((1.0-u1-2.0*u2)*lambdae+(u1+2.0*u2)*lambdad+u2*etad)/omega;
		mu0 = 1.0-lambdae;
		return m_doLimbDarkening ? muo1 : mu0;
	}
	//occulting star transits the source:
	if(p < 1.0  && z < (1.0 - p)*1.0001)
	{
		//cout << "zone 7\n" << endl;
		q = sqrt((x2-x1)/(1.0-x1));
		Kk = ellk(q);
		Ek = ellec(q);
		n = x2/x1-1.0;
		Pk = Kk-n/3.0*rj(0.0,1.0-q*q,1.0,1.0+n);
		lambdad = 2.0/9.0/pi/sqrt(1.0-x1)*((1.0-5.0*z*z+p*p+
				x3*x3)*Kk+(1.0-x1)*(z*z+7.0*p*p-4.0)*Ek-3.0*x3/x1*Pk);
		if(z < p) lambdad += 2.0/3.0;
		if(fabs(p+z-1.0) < 1.0e-4)
		{
			lambdad = 2.0/3.0/pi*acos(1.0-2.0*p)-4.0/9.0/pi*
					sqrt(p*(1.0-p))*(3.0+2.0*p-8.0*p*p);
		}
		etad = p*p/2.0*(p*p+2.0*z*z);
	}
	muo1 = 1.0-((1.0-u1-2.0*u2)*lambdae+(u1+2.0*u2)*lambdad+u2*etad)/omega;
	mu0 = 1.0-lambdae;

	return m_doLimbDarkening ? muo1 : mu0;

}

double TransitModel::rc(double x, double y) const {
	double rc, ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2, C3,C4;
	ERRTOL=0.04; TINY=1.69e-38; SQRTNY=1.3e-19; BIG=3.0e37;
	TNBG=TINY*BIG; COMP1=2.236/SQRTNY; COMP2=TNBG*TNBG/25.0;
	THIRD=1.0/3.0; C1=0.3; C2=1.0/7.0; C3=0.375; C4=9.0/22.0;

	double alamb,ave,s,w,xt,yt;
	if(x < 0.0 || y == 0.0 || (x+fabs(y)) < TINY || (x+fabs(y)) > BIG || (y < -COMP1 && x > 0 && x < COMP2)){
		cout << "Invalid argument(s) in rc\n" << endl;
		return 0;
	}
	if(y > 0.0)
	{
		xt=x;
		yt=y;
		w=1.0;
	}
	else
	{
		xt=x-y;
		yt=-y;
		w=sqrt(x)/sqrt(xt);
	}
	s = ERRTOL*10.0;
	while(fabs(s) > ERRTOL)
	{
		alamb = 2.0*sqrt(xt)*sqrt(yt)+yt;
		xt = 0.25*(xt+alamb);
		yt  =0.25*(yt+alamb);
		ave = THIRD*(xt+yt+yt);
		s = (yt-ave)/ave;
	}
	rc = w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
	return rc;
}

double TransitModel::rj(double x, double y, double z, double p) const {
	double rj, ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8, tempmax;

	ERRTOL=0.05; TINY=2.5e-13; BIG=9.0e11; C1=3.0/14.0;
	C2=1.0/3.0; C3=3.0/22.0; C4=3.0/26.0; C5=.750*C3;
     	C6=1.50*C4; C7=.50*C2; C8=C3+C3;

	double  a = 0.0,alamb,alpha,ave,b = 0.0,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,
     		fac,pt,rcx = 0.0,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;

	if(x < 0.0 || y < 0.0 || z < 0.0 || (x+y) < TINY || (x+z) < TINY || (y+z) < TINY || fabs(p) < TINY
		|| x > BIG || y > BIG || z > BIG || fabs(p) > BIG)
	{
		cout << "Invalid argument(s) in rj\n" << endl;
		return 0;
	}
	sum=0.0;
	fac=1.0;
	if(p > 0.0)
	{
		xt=x;
		yt=y;
		zt=z;
		pt=p;
	}
	else
	{
		xt = min(x, y);
		xt = min(xt,z);
		zt = max(x, y);
		zt = max(zt, z);
		yt = x+y+z-xt-zt;
		a = 1.0/(yt-p);
		b = a*(zt-yt)*(yt-xt);
		pt = yt+b;
		rho = xt*zt/yt;
		tau = p*pt/yt;
		rcx = rc(rho,tau);
	}
	tempmax = ERRTOL*10.0;
	while(tempmax > ERRTOL)
	{
		sqrtx = sqrt(xt);
		sqrty = sqrt(yt);
		sqrtz = sqrt(zt);
		alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		alpha = pow((pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz), 2.0);
		beta = pt*(pt+alamb)*(pt + alamb);
		sum = sum+fac*rc(alpha,beta);
		fac = 0.25*fac;
		xt = 0.25*(xt+alamb);
		yt = 0.25*(yt+alamb);
		zt = 0.250*(zt+alamb);
		pt = 0.25*(pt+alamb);
		ave = 0.2*(xt+yt+zt+pt+pt);
		delx = (ave-xt)/ave;
		dely = (ave-yt)/ave;
		delz = (ave-zt)/ave;
		delp = (ave-pt)/ave;
		tempmax = max(fabs(delx), fabs(dely));
		tempmax = max(tempmax, fabs(delz));
		tempmax = max(tempmax, fabs(delp));
	}
	ea = delx*(dely+delz)+dely*delz;
	eb = delx*dely*delz;
	ec = delp*delp;
	ed = ea-3.0*ec;
	ee = eb+2.0*delp*(ea-ec);
	rj = 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4)) +
		delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
	if(p < 0.0) rj=a*(b*rj+3.0*(rcx-rf(xt,yt,zt)));
	return rj;
}

double TransitModel::ellec(double k) const {
	double m1,a1,a2,a3,a4,b1,b2,b3,b4,ee1,ee2,ellec;
	// Computes polynomial approximation for the complete elliptic
	// integral of the second kind (Hasting's approximation):
	m1=1.0-k*k;
	a1=0.44325141463;
	a2=0.06260601220;
	a3=0.04757383546;
	a4=0.01736506451;
	b1=0.24998368310;
	b2=0.09200180037;
	b3=0.04069697526;
	b4=0.00526449639;
	ee1=1.0+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
	ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.0/m1);
	ellec=ee1+ee2;
	return ellec;
}

double TransitModel::ellk(double k) const {
	double a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ellk, ek1,ek2,m1;
	// Computes polynomial approximation for the complete elliptic
	// integral of the first kind (Hasting's approximation):
	m1=1.0-k*k;
	a0=1.38629436112;
	a1=0.09666344259;
	a2=0.03590092383;
	a3=0.03742563713;
	a4=0.01451196212;
	b0=0.5;
	b1=0.12498593597;
	b2=0.06880248576;
	b3=0.03328355346;
	b4=0.00441787012;
	ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
	ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1);
	ellk=ek1-ek2;
	return ellk;
}

double TransitModel::rf(double x, double y, double z) const {
	double rf, ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4, tempmax;

	ERRTOL=0.08; TINY=1.5e-38; BIG=3.0e37; THIRD=1.0/3.0;
	C1=1.0/24.0; C2=0.1; C3=3.0/44.0; C4=1.0/14.0;

	double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;

	if(min(x,y) < 0.0 || z < 0.0 || min(x+y, x+z) < TINY || y+z < TINY || max(x,y) > BIG || z > BIG)
	{
		cout << "Invalid argument(s) in rf\n" << endl;
		return 0;
	}
	xt=x;
	yt=y;
	zt=z;
	tempmax = ERRTOL*10.0;
	while(tempmax > ERRTOL)
	{
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=THIRD*(xt+yt+zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
		tempmax = max(fabs(delx), fabs(dely));
		tempmax = max(tempmax, fabs(delz));
	}
	e2=delx*dely-delz*delz;
	e3=delx*dely*delz;
	rf=(1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
	return rf;
}

double TransitModel::max(double a, double b) const {
	if(a > b) return a;
	else return b;
}

double TransitModel::min(double a, double b) const {
	if(a < b) return a;
	else return b;
}
