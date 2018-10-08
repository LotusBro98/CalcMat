#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

class Matrix
{
	public:

    explicit Matrix(int n, int m)
	{
		this->_n = n;
        this->_m = m;
		this->_a = new double[n * m];
	}

/*
    explicit Matrix(std::istream &stream)
    {
        int n;
        stream >> n;

        this->_n = n;
        this->_m = n;
        this->_a = new double[n * n];

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                stream >> a(i, j);
            }
        }
    }
*/

/*
	Matrix(double (*f)(int i, int j, void* arg), int n, int m, void * arg) : Matrix(n, m)
	{
		for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
				a(i, j) = f(i, j, arg);
	}
*/

    friend std::ostream& operator<< (std::ostream & os, Matrix * A)
    {
        for (int i = 0; i < A->_n; ++i) {
            for (int j = 0; j < A->_m; ++j) {
                os << std::setw(3) << A->a(i, j) << ' ';
            }
            os << std::endl;
        }
        os << std::endl;
        return os;
    }

    Matrix * operator | (Matrix & b)
    {
        int n = std::max(_n, b._n);
        int m = _m + b._m;

        Matrix * compound = new Matrix(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                if (i < _n && j < _m)
                    compound->a(i, j) = a(i, j);
                else if (i < b._n && j >= _m)
                    compound->a(i, j) = b.a(i, j - _m);
                else
                    compound->a(i, j) = 0;
            }
        }

        return compound;
    }

    void addLine(int i0, int i1, double c = 1)
    {
        for (int j = 0; j < _m; ++j) {
            a(i0, j) += a(i1, j) * c;
        }
    }

    void subtractLine(int i0, int i1, double c = 1)
    {
        for (int j = 0; j < _m; ++j) {
            a(i0, j) -= a(i1, j) * c;
        }
    }

    void mulLine(int i, double c)
    {
        for (int j = 0; j < _m; ++j) {
            a(i, j) *= c;
        }
    }

    void divLine(int i, double c)
    {
        for (int j = 0; j < _m; ++j) {
            a(i, j) /= c;
        }
    }

    void swapLines(int i1, int i2)
    {
        for (int j = 0; j < _m; ++j) {
            double t = a(i1, j);
            a(i1, j) = a(i2, j);
            a(i2, j) = t;
        }
    }

    void gauss()
    {
        for (int j = 0; j < _m; ++j)
        {
            int in0;
            for (in0 = j; in0 < _n; in0++)
                if (a(in0, j) != 0)
                    break;
            if (in0 == _n)
                continue;

            if (in0 != j)
                swapLines(in0, j);

            divLine(j, a(j, j));

            for (int i = 0; i < _n; ++i) {
                if (i == j || a(i, j) == 0)
                    continue;

                subtractLine(i, j, a(i, j));
            }
        }
    }

    Matrix * subMatrix(int i0, int j0, int n, int m)
    {
        Matrix * sub = new Matrix(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                sub->a(i, j) = a(i0 + i, j0 + j);
            }
        }

        return sub;
    }

    friend Matrix * operator * (Matrix & A, Matrix & B)
    {
        if (A._m != B._n)
            return nullptr;

        Matrix * C = new Matrix(A._n, B._m);
        for (int i = 0; i < C->_n; ++i) {
            for (int j = 0; j < C->_m; ++j) {
                C->a(i, j) = 0;
                for (int k = 0; k < A._m; ++k) {
                    C->a(i, j) += A.a(i, k) * B.a(k, j);
                }
            }
        }

        return C;
    }

    friend Matrix * operator - (Matrix & A, Matrix & B)
    {
        if (A._n != B._n || A._m != B._m)
            return nullptr;

        Matrix * C = new Matrix(A._n, A._m);
        for (int i = 0; i < A._n; ++i) {
            for (int j = 0; j < A._m; ++j) {
                C->a(i, j) = A.a(i, j) - B.a(i, j);
            }
        }

        return C;
    }

    friend Matrix * operator + (Matrix & A, Matrix & B)
    {
        if (A._n != B._n || A._m != B._m)
            return nullptr;

        Matrix * C = new Matrix(A._n, A._m);
        for (int i = 0; i < A._n; ++i) {
            for (int j = 0; j < A._m; ++j) {
                C->a(i, j) = A.a(i, j) + B.a(i, j);
            }
        }

        return C;
    }


    double vectorNormSqrt()
    {
        double len2 = 0;
        if (_n == 1)
            for (int j = 0; j < _m; ++j)
                len2 += a(0, j) * a(0, j);
        else if (_m == 1)
            for (int i = 0; i < _n; ++i)
                len2 += a(i, 0) * a(i, 0);
        else
            return -1;

        return std::sqrt(len2);
    }


    Matrix * copy()
    {
        Matrix * A = new Matrix(_n, _m);
        A->_a = new double[_n * _m];

        std::copy(_a, _a + _n * _m,  A->_a);

        return A;
    }

    ~Matrix()
    {
        delete(_a);
    }

	int n()
    {
        return _n;
    }

    int m() const {
        return _m;
    }

    inline double & a(int i, int j)
    {
        return _a[_m * i + j];
    }

	private:
	int _n;
    int _m;
	double * _a;
};

Matrix * projectOn1DGrid(double (*f)(double x), double * X, int np)
{
	Matrix * vf = new Matrix(np, 1);
	for (int i = 0; i < np; i++)
		vf->a(i, 0) = f(X[i]);
	
	return vf;
}

double scalarProduct(Matrix * va, Matrix * vb)
{
	if (va->n() != vb->n())
		throw;

	double sum = 0;
	for (int i = 0; i < va->n(); i++)
		sum += va->a(i, 0) * vb->a(i, 0);

	return sum;
}

Matrix * approx(double * X, double * F, double (**bf)(double x), int np, int nf)
{
	Matrix ** bfproj = new Matrix*[nf];
	for (int i = 0; i < nf; i++)
		bfproj[i] = projectOn1DGrid(bf[i], X, np);

	Matrix * fproj = new Matrix(np, 1);
	for (int i = 0; i < np; i++)
		fproj->a(i, 0) = F[i];

	Matrix * D = new Matrix(nf, nf);
	for (int i = 0; i < nf; i++)
		for (int j = 0; j < nf; j++)
			D->a(i, j) = scalarProduct(bfproj[i], bfproj[j]);

	Matrix * bff = new Matrix(nf, 1);
	for (int i = 0; i < nf; i++) 
		bff->a(i, 0) = scalarProduct(fproj, bfproj[i]);

	Matrix * Dexp = *D | *bff;

	Dexp->gauss();

	Matrix * u = Dexp->subMatrix(0, Dexp->m() - 1, Dexp->n(), 1);
	return u;
}

double f0(double x)
{
	return 1;
}

double f1(double x)
{
	return x;
}

double f2(double x)
{
	if (x < 20) return 1;
	if (x >= 28) return 0;
	return 1 - (x - 20)/8;
}

double f3(double x)
{
	if ((x < 20) || (x >= 39)) return 0;
	if (x < 28) return (x - 20) / 8;
	return 1 - (x - 28) / 11;
}

double f4(double x)
{
	if ((x < 28) || (x >= 45)) return 0;
	if (x < 39) return (x - 28)/11;
	return 1 - (x - 39)/6;
}

double f5(double x)
{
	if (x < 39) return 0;
	if (x >= 45) return 1;
	return (x - 39)/6;
}


int main(int argc, char * argv[])
{
	int np = 26;
	int nf = 4; //2

	double dataset[np] = {431, 409, 429, 422, 530, 505, 459, 499, 526, 563, 587, 595, 647, 669, 746, 760, 778, 828, 846, 836, 916, 956, 1014, 1076, 1134, 1024};
	double points[np] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45};

	double (*bf[nf])(double x) = {f2, f3, f4, f5};//f0,f1

	Matrix * u = approx(points, dataset, bf, np, nf);

	std::cout << u << "\n"; 

	return 0;
}


// ------------------------------------------------------------
