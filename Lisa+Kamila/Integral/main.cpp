#include <iostream>
#include <cmath>
#include <iomanip>

double f(double x)
{
	return exp(x);
}

double F(double x)
{
	return exp(x);
}

double t[3] = {
	-sqrt(3.0 / 5.0),
	0,
	sqrt(3.0 / 5.0)
};

double w[3] = {
	2.0 * (1 + 3 * t[1] * t[2]) / (t[0] - t[1]) / (t[0] - t[2]) / 3.0,
	2.0 * (1 + 3 * t[0] * t[2]) / (t[1] - t[0]) / (t[1] - t[2]) / 3.0,
	2.0 * (1 + 3 * t[0] * t[1]) / (t[2] - t[0]) / (t[2] - t[1]) / 3.0 
};

double gauss(double from, double to, int n)
{
	double Int = 0;
	for (int i = 0; i < n; i++)
	{
		double xa = from + i * (to - from) / n;
		double xb = from + (i + 1) * (to - from) / n;

		double x[3];
		for (int j = 0; j < 3; j++)
			x[j] = 0.5 * (xa + xb) + 0.5 * (xb - xa) * t[j];

		double sum = 0;
		for (int j = 0; j < 3; j++)
			sum += w[j] * f(x[j]);

		Int += 0.5 * (xb - xa) * sum;
	}

	return Int;
}

double simpson(double from, double to, int n)
{
	double Int = 0;
	for (int i = 0; i < n; i++)
	{
		double xa = from + i * (to - from) / n;
		double xb = from + (i + 1) * (to - from) / n;
		double xm = 0.5 * (xa + xb);

		Int += (xb - xa) / 6.0 * (f(xa) + 4 * f(xm) + f(xb));
	}

	return Int;
}

int main()
{
	int nlist[4] = {10, 20, 40, 80};
	double a = 0;
	double b = 1;

	for (int i = 0; i < 4; i++)
	{
		int n = nlist[i];

		double Int1 = F(b) - F(a);
		double Int2 = simpson(a, b, n);
		double Int3 = gauss(a, b, n);

		std::cout << std::setprecision(16) 
			<< "n             = " << n << std::endl
			<< "Int1          = " << Int1 << std::endl
			<< "Int2          = " << Int2 << std::endl
			<< "Int3          = " << Int3 << std::endl
			<< "|Int2 - Int1| = " << fabs(Int2 - Int1) << std::endl
			<< "|Int3 - Int1| = " << fabs(Int3 - Int1) << std::endl
			<< std::endl;
	}

	return 0;
}
