#include "Lsqr.h"
/*
   Code is based on scipy.sparse library
   https://github.com/scipy/scipy/blob/v1.5.4/scipy/sparse/linalg/isolve/lsqr.py
*/


int sign(double val){
	return (val > 0) - (val < 0);
}


void Lsqr::_sym_ortho(double a, double b) {
	if (b == 0) {
		c = sign(a);
		s = 0;
		r = abs(a);
	}
	else if (a == 0) {
		c = 0;
		s = sign(b);
		r = abs(b);
	}
	else if (abs(b) > abs(a)) {
		double tau = a / b;
		s = sign(b) / sqrt((1 + tau * tau));
		c = s * tau;
		r = b / s;
	}
	else {
		double tau = b / a;
		c = sign(a) / sqrt(1 + tau * tau);
		s = c * tau;
		r = a / c;
	}

}

Array& Lsqr::solve(Applyble& A, Array& b, Array& x) {

	std::cout << "solve\n";
	int* shape = A.get_shape();
	int m = shape[0], n = shape[1];
	Array u(m), v(n), w(n), dk(n);
	u = 0;
	v = 0;
	w = 0;
	dk = 0;
	if (iter_lim == 0)
		iter_lim = 2 * n;
	int itn = 0, istop = 0;
	double ctol = 0, anorm = 0, acond = 0, dampsq = damp * damp, ddnorm = 0, res1, res2 = 0,
		xnorm, xxnorm = 0, z = 0, cs = -1, sn2 = 0, beta, alfa, rhobar, phibar, rnorm,
		r1norm, r2norm, arnorm, sn, rho, theta, bnorm, rhobar1, cs1, cs2 = -1, sn1,
		psi, phi, tau, t1, t2, delta, gambar, rhs, zbar, gamma, r1sq, test1, test2,
		test3, rtol;
	if (conlim > 0)
		ctol = 1 / conlim;
	bnorm = b.norm();
	A.apply(u, x);
	u.iuminus();
	u += b;
	beta = u.norm();

	if (beta > 0) {
		u *= 1 / beta;
		A.applyT(v, u);
		alfa = v.norm();
	}
	
	else {
		v = x;
		alfa = 0;
	}

	if (alfa > 0) {
		v *= 1 / alfa;
		w = v;
	}
	rhobar = alfa;
	phibar = beta;
	rnorm = beta;
	r1norm = rnorm;
	r2norm = rnorm;

	arnorm = alfa * beta;
    
	if (arnorm == 0) {
		return x;
	}
	while (itn < iter_lim) {
		itn++;
		u *= -alfa;

		A.apply(u, v);

		beta = u.norm();

		if (beta > 0) {
			u *= 1 / beta;
			anorm = sqrt(anorm * anorm + alfa * alfa + beta * beta + damp * damp);
			v *= -beta;
			A.applyT(v, u);
			alfa = v.norm();
			if (alfa > 0)
				v *= 1 / alfa;
		}

		rhobar1 = sqrt(rhobar * rhobar + damp * damp);
		cs1 = rhobar / rhobar1;
		sn1 = damp / rhobar1;
		psi = sn1 * phibar;
		phibar = cs1 * phibar;

		// Use a plane rotation to eliminate the subdiagonal element (beta)
		// of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
		_sym_ortho(rhobar1, beta);

		cs = c;
		sn = s;
		rho = r;
		theta = sn * alfa;
		rhobar = -cs * alfa;
		phi = cs * phibar;
		phibar = sn * phibar;
		tau = sn * phi;

		// Update x and w.
		t1 = phi / rho;
		t2 = -theta / rho;
		dk = w;
		dk *= 1 / rho;

		x.iplus_weighted(w, t1);
		w *= t2;
		w += v;
		ddnorm = ddnorm + dk.norm2();

		delta = sn2 * rho;
		gambar = -cs2 * rho;
		rhs = phi - delta * z;
		zbar = rhs / gambar;
		xnorm = sqrt(xxnorm + zbar * zbar);
		gamma = sqrt(gambar * gambar + theta * theta);
		cs2 = gambar / gamma;
		sn2 = theta / gamma;
		z = rhs / gamma;
		xxnorm = xxnorm + z * z;

		// Test for convergence.
		// First, estimate the condition of the matrix  Abar,
		// and the norms of  rbar  and  Abar'rbar.
		acond = anorm * sqrt(ddnorm);
		res1 = phibar * phibar;
		res2 = res2 + psi * psi;
		rnorm = sqrt(res1 + res2);
		arnorm = alfa * abs(tau);

		/*# Distinguish between
		#    r1norm = ||b - Ax|| and
		#    r2norm = rnorm in current code
		#           = sqrt(r1norm^2 + damp^2*||x||^2).
		#    Estimate r1norm from
		#    r1norm = sqrt(r2norm^2 - damp^2*||x||^2).
		# Although there is cancellation, it might be accurate enough.*/

		r1sq = rnorm * rnorm - dampsq * xxnorm;
		r1norm = sqrt(abs(r1sq));
		if (r1sq < 0)
			r1norm = -r1norm;
		r2norm = rnorm;

		//# Now use these norms to estimate certain other quantities,
		//# some of which will be small near a solution.
		test1 = rnorm / bnorm;
		test2 = arnorm / (anorm * rnorm + eps);
		test3 = 1 / (acond + eps);
		t1 = test1 / (1 + anorm * xnorm / bnorm);
		rtol = btol + atol * anorm * xnorm / bnorm;


		if (itn >= iter_lim)
			istop = 7;
		if (1 + test3 <= 1)
			istop = 6;
		if (1 + test2 <= 1)
			istop = 5;
		if (1 + t1 <= 1)
			istop = 4;

		//# Allow for tolerances set by the user.
		if (test3 <= ctol)
			istop = 3;
		if (test2 <= atol)
			istop = 2;
		if (test1 <= rtol)
			istop = 1;
		if (istop != 0)
			break;
	}
	return x;
}


