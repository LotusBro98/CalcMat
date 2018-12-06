#include <math.h>
#include <stdio.h>

#define f(x) 1 / cos(x) / cos(x)
#define Pervoobr(x) tan(x)

//ti - корни многочленов Лежандра
//ci - веса квадратуры Гаусса-Лежандра на (-1;1) по трем узлам

double t1 = -0.774596669241483;
double c1 = 0.555555555555555;

double t2 = 0.000000000000000;
double c2 = 0.888888888888889;

double t3 = 0.774596669241483;
double c3 = 0.555555555555555;

double quadraturaGaussa(int n, double a, double b)
{
	double S = 0;
	for (int i = 0; i < n; i++)
	{
		double xa = a + i * (b - a) / n;
		double xb = a + (i + 1) * (b - a) / n;

		double x1 = 0.5 * (xa + xb) + 0.5 * (xb - xa) * t1;
		double x2 = 0.5 * (xa + xb) + 0.5 * (xb - xa) * t2;
		double x3 = 0.5 * (xa + xb) + 0.5 * (xb - xa) * t3;

		double deltaS = c1 * f(x1) + c2 * f(x2) + c3 * f(x3);

		S += 0.5 * (xb - xa) * deltaS;
	}

	return S;
}

double quadraturaSimpsona(int n, double a, double b)
{
	double S = 0;
	for (int i = 0; i < n; i++)
	{
		double xa = a + i * (b - a) / n;
		double xb = a + (i + 1) * (b - a) / n;
		double xm = 0.5 * (xa + xb);

		S += (xb - xa) / 6.0 * (f(xa) + 4 * f(xm) + f(xb));
	}

	return S;
}

int main()
{
	int nValues[4] = {10, 20, 40, 80};
	double a = 0;
	double b = 1;

	for (int i = 0; i < 4; i++)
	{
		int n = nValues[i];

		double Int1 = Pervoobr(b) - Pervoobr(a);
		double Int2 = quadraturaSimpsona(n, a, b);
		double Int3 = quadraturaGaussa(n, a, b);

		printf("n: %d\nInt1: %.15lf\nInt2: %.15lf\nInt3: %.15lf\n|Int2 - Int1|: %lg\n|Int3 - Int1|: %lg\n\n",
				n,
				Int1,
				Int2,
				Int3,
				fabs(Int2 - Int1),
				fabs(Int3 - Int1));
	}

	return 0;
}
