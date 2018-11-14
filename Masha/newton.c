#include <stdio.h>
#include <math.h>

void mulAdd(double * a, double * b, double c, int n)
{
	for (int i = 0; i < n; i++)
		a[i] += c * b[i];
}

void scale(double * a, int n, double c)
{
	for (int i = 0; i < n; i++)
		a[i] *= c;
}

double dotVec(double * a, double * b, int n)
{
	double c = 0;
	for (int i = 0; i < n; i++)
		c += a[i] * b[i];
	return c;
}

double len(double * a, int n)
{
	return sqrt(dotVec(a, a, n));
}

void printVec(double * a, int n)
{
	for (int i = 0; i < n; i++)
		printf("%lf\n", a[i]);
}

double * line(double * A, int m, int i)
{
	return A + m * i;
}

double * aij(double * A, int m, int i, int j)
{
	return A + m * i + j;
}

void printMat(double * A, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			printf("%6.2lf ", *aij(A, m, i, j));
		printf("\n");
	}
	printf("\n");
}

double * getCol(double * A, int n, int m, int j, double * buf)
{
	for (int i = 0; i < n; i++)
		buf[i] = *aij(A, m, i, j);
	return buf;
}

double * setCol(double * A, int n, int m, int j, double * buf)
{
	for (int i = 0; i < n; i++)
		*aij(A, m, i, j) = buf[i];
	return buf;
}

void hstack(double * A, double * B, double * buf, int n, int m1, int m2)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m1; j++)
			*aij(buf, m1 + m2, i, j) = *aij(A, m1, i, j);

		for (int j = 0; j < m2; j++)
			*aij(buf, m1 + m2, i, j + m1) = *aij(B, m2, i, j);
	}
}

void setVec(double * dst, double * src, int n)
{
	for (int i = 0; i < n; i++)
		dst[i] = src[i];
}

//---------------------------Gauss-------------------------------

int findNonzeroInCol(double * A, int n, int m, int j, int i)
{
	for (; i < n; i++)
		if (*aij(A, m, i, j) != 0)
			return i;
	return -1;
}

void swapLines(double * A, int m, int i1, int i2)
{
	double t;
	for (int j = 0; j < m; j++)
	{
		t = *aij(A, m, i1, j);
		*aij(A, m, i1, j) = *aij(A, m, i2, j); 
		*aij(A, m, i2, j) = t; 
	}
}

void Gauss(double * A, int n, int m)
{
	int i = 0;
	int j = 0;
	while (j < m && i < n)
	{
		if (*aij(A, m, i, j) == 0)
		{
			int nz = findNonzeroInCol(A, n, m, j, i);
			if (nz == -1)
			{
				j++;
				continue;
			}
			swapLines(A, m, i, nz);
		}

		scale(line(A, m, i), m, 1 / *aij(A, m, i, j));
		for (int i1 = 0; i1 < n; i1++)
			if (i1 != i)
				mulAdd(line(A, m, i1), line(A, m, i), -*aij(A, m, i1, j), m);

		i++;
		j++;
	}
}

double * solveGauss(double * A, double * b, double * x, int n)
{
	double Ab[n * (n + 1)];
	hstack(A, b, Ab, n, n, 1);	

	Gauss(Ab, n, n + 1);

	for (int i = 0; i < n; i++)
		if (*aij(Ab, n + 1, i, i) != 1)
			return NULL;

	getCol(Ab, n, n + 1, n, x);	

	return x;
}
// -------------------------Newton------------------------------

double * Newton(void (*setJ)(double * J, double * x), void (*setF)(double * F, double * x),
		int n, double * x, double * x0, double eps)
{
	double J[n * n];
	double F[n];
	double u[n];
	
	double err;

	setVec(x, x0, n);	

	do
	{
		setF(F, x);
		setJ(J, x);

		solveGauss(J, F, u, n);
		mulAdd(x, u, -1, n);

		err = len(F, n);
		printf("%lf\n", err);
	}
	while (err > eps);

	return x;
}

//---------------------------Data-------------------------------

void setJ(double * J, double * x)
{
	J[0] = 4 - x[3];
	J[1] = -1;
	J[2] = 1;
	J[3] = -x[0];

	J[4] = 1;
	J[5] = -2;
	J[6] = 3 - x[3];
	J[7] = -x[2];

	J[8] = -1;
	J[9] = 3 - x[3];
	J[10] = -2;
	J[11] = -x[1];

	J[12] = 2 * x[0];
	J[13] = 2 * x[1];
	J[14] = 2 * x[2];
	J[15] = 0;	
}

void setF(double * F, double * x)
{
	F[0] = 4 * x[0] - x[1] + x[2] - x[0] * x[3];
	F[1] = x[0] - 2 * x[1] + 3 * x[2] - x[2] * x[3];
	F[2] = -x[0] + 3 * x[1] - 2 * x[2] - x[1] * x[3];
	F[3] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1;
}

//-------------------------main---------------------------------

int main()
{
	double x[4] = {0.1, 0.1, 0.1, 0.1};

	Newton(setJ, setF, 4, x, x, 1e-4);

	printf("---\n");
	printVec(x, 4);
	printf("---\n");
}
