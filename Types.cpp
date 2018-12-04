#include "Types.h"

using namespace std;

void Barrel::calcHi1() {
	hi1 = 1.0 - (Consts::k - 1.0) * p_mid / Consts::f * (1.0 - Consts::alpha_k * Consts::delta) / Consts::delta;
}

void Barrel::calcForTest(double _B)
{
	B = _B;
	Delta = 500;
	pm = 240e6;
	p_mid = 0.5 * pm;
	double eta_K = 0.5;

	hi1 = 1.0 - (Consts::k - 1.0) * p_mid / Consts::f * (1.0 - Consts::alpha_k * Consts::delta) / Consts::delta;

	psi0 = (1.0 / Delta - 1.0 / Consts::delta) / (Consts::f / p0 - (1.0 - Consts::alpha_k * Consts::delta) / Consts::delta);
	sigma0 = sqrt(1.0 + 4.0 * Consts::lambda / Consts::kapa * psi0);
	z0 = 2.0 * psi0 / (Consts::kapa * (1.0 + sigma0));
	K1 = Consts::kapa * sigma0;

	B1 = 0.5 * (Consts::k - 1) * B - Consts::kapa * Consts::lambda * hi1;
	gamma1 = psi0 * B1 / pow(K1, 2.0);
	alpha1 = 2 * Consts::kapa * Consts::lambda / B1;
	beta1 = 0.5 * hi1 + sqrt(gamma1 + 0.25 * pow(hi1, 2.0));
	beta2 = 0.5 * hi1 - sqrt(gamma1 + 0.25 * pow(hi1, 2.0));

	beta_m = (Consts::k * hi1 - 1) / (2.0 * Consts::k + alpha1);
	beta_m_1 = beta_m / beta1;
	beta_m_2 = beta_m / beta2;

	fi_m = pow(1.0 - beta_m_2, (1.0 + alpha1 * beta2) / (beta1 - beta2) - 1.0) /
		pow(1.0 - beta_m_1, (1.0 + alpha1 * beta1) / (beta1 - beta2) + 1.0);

	pi_m = (1 - beta_m_1) * (1 - beta_m_2) * pow(fi_m, -1.0 / (Consts::k - 1.0));

	calcLambdaK();
	calcForEtaK(eta_K);
}

void Barrel::calcB(double a, double b)
{
	p_mid = 0.5 * pm;
	hi1 = 1.0 - (Consts::k - 1) * p_mid / Consts::f * (1 - Consts::alpha_k * Consts::delta) / Consts::delta;

	psi0 = (1.0 / Delta - 1.0 / Consts::delta) / (Consts::f / p0 - (1.0 - Consts::alpha_k * Consts::delta) / Consts::delta);
	sigma0 = sqrt(1 + 4.0 * Consts::lambda / Consts::kapa * psi0);
	z0 = 2.0 * psi0 / (Consts::kapa * (1.0 + sigma0));
	K1 = Consts::kapa * sigma0;

	halfSearchB(a, b);
}

void Barrel::halfSearchB(double a, double b)
{
	double pi_m_star = pm / p0;

	while (true)
	{
		B = 0.5 * (a + b);

		B1 = 0.5 * (Consts::k - 1) * B - Consts::kapa * Consts::lambda * hi1;
		gamma1 = psi0 * B1 / pow(K1, 2.0);
		alpha1 = 2 * Consts::kapa * Consts::lambda / B1;
		beta1 = 0.5 * hi1 + sqrt(gamma1 + 0.25 * pow(hi1, 2.0));
		beta2 = 0.5 * hi1 - sqrt(gamma1 + 0.25 * pow(hi1, 2.0));

		beta_m = (Consts::k * hi1 - 1) / (2.0 * Consts::k + alpha1);
		beta_m_1 = beta_m / beta1;
		beta_m_2 = beta_m / beta2;

		fi_m = pow(1.0 - beta_m_2, (1.0 + alpha1 * beta2) / (beta1 - beta2) - 1.0) /
			pow(1.0 - beta_m_1, (1.0 + alpha1 * beta1) / (beta1 - beta2) + 1.0);

		pi_m = (1 - beta_m_1) * (1 - beta_m_2) * pow(fi_m, -1.0 / (Consts::k - 1.0));

		if (fabs(pi_m - pi_m_star) < eps) break;					// Условие выхода из цикла

		if (pi_m > pi_m_star) a = B;
		else b = B;
	}
}

void Barrel::calcLambdaK()
{
	beta_k = B1 / K1 * (1.0 - z0);
	beta_k1 = beta_k / beta1;
	beta_k2 = beta_k / beta2;
	fi_k = pow(1.0 - beta_k2, (1 + alpha1 * beta2) / (beta1 - beta2) - 1.0) /
		pow(1.0 - beta_k1, (1 + alpha1 * beta1) / (beta1 - beta2) + 1.0);

	Lambda_K = Consts::f * Delta * psi0 / p0 * (pow(fi_k, 1.0 / (Consts::k - 1.0)) - 1.0) -
		(1.0 - Consts::alpha_k * Consts::delta) / Consts::delta * psi0 * Delta * beta_k / gamma1 *
		(1.0 + 0.5 * alpha1 * beta_k);
	r_K = 0.5 * (Consts::k - 1.0) * B * pow(1.0 - z0, 2.0);
}

void Barrel::calcForEtaK(double eta_K)
{
	Lambda_D = Lambda_K / eta_K;
	r_D = hi1 - (hi1 - r_K) * pow((1.0 - Consts::alpha_k * Delta + Lambda_K) /
		(1.0 - Consts::alpha_k * Delta + Lambda_D), Consts::k - 1.0);

	omega = q * K / (2 * Consts::f / (Consts::k - 1.0) * r_D / pow(vd, 2.0) - 1.0 / 3.0);
	omega_q = omega / q;
	W0 = omega / Delta;
	double S = 0.25 * pi * pow(d, 2.0) * ns;
	L0 = W0 / S;
	LD = Lambda_D * L0;
	W = LD * S;
	fi = K + 1.0 / 3.0 * omega / q;
	Ik = sqrt(Consts::f * omega * fi * q * B) / S;
}

double Barrel::calcZSluh()
{
	Z_Sluh = Consts::C * sqrt(1.0 + Lambda_D) / (pow(omega_q, 1.5) * pow((LD + L0 / hi + 1.5 * d) / d, 4.0));
	return Z_Sluh;
}

double Barrel::calcZSluh(double v1_vd) {
	Z_Sluh = Consts::C * sqrt(1.0 + 1.0 / Lambda_D) / (pow(omega_q, 1.5) * pow((LD + L0 / hi + 1.5 * d) / d, 4.0) * v1_vd);
	return Z_Sluh;
}

Matrix& Matrix::T()
{
	double **a = new double*[y];
	for (int i = 0; i < y; i++)
		a[i] = new double[x];

	for (int i = 0; i < x; i++)
		for (int j = 0; j < y; j++)
			a[j][i] = data[i][j];

	for (int i = 0; i < x; i++)
		delete[] data[i];
	delete[] data;

	int buf = x;
	x = y;

	y = buf;
	data = a;

	return *this;
}

Matrix& Matrix::zeros()
{
	for (int i = 0; i < x; i++)
		for (int j = 0; j < y; j++)
			data[i][j] = 0;

	return *this;
}

void Criterion::calcCriterion(const CResult &res, const CriterionParams &cp)
{
	double alpha_n[4];
	alpha_n[0] = -cp.k[0] * cp.alpha[0] * ksi_l_d(res.L / cp.d, cp.l_d_max, cp.a);
	alpha_n[1] = -cp.k[1] * cp.alpha[1];
	alpha_n[2] = -cp.k[2] * cp.alpha[2];
	alpha_n[3] = -cp.k[3] * cp.alpha[3] * ksi_p(res.p_max, cp.pm_star, cp.b);

	pwd_name = cp.pwd_name;
	Delta = res.Delta;
	w_q = res.w_q;
	Z = pow((res.L / cp.d) / cp.l_d_ref, alpha_n[0]) *
			pow((res.W0 / pow(cp.d, 3.0)) / cp.W0_d_ref, alpha_n[1]) *
			pow(res.w_q / cp.w_q_ref, alpha_n[2]) *
			pow(res.p_max / cp.pm_star, alpha_n[3]);
}

double Criterion::ksi_l_d(double l_d, double l_d_max, double a)
{
	if (l_d > l_d_max)
		return a * l_d / l_d_max;
	return 1.0;
}

double Criterion::ksi_p(double pm, double pm_star, double b)
{
	if (pm > pm_star)
		return b * pm / pm_star;
	return 1.0;
}
