#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string.h>

class Matrix
{
	public:

    explicit Matrix(int n, int m)
	{
		this->_n = n;
        this->_m = m;
		this->_a = new double[n * m];
	}

    explicit Matrix(std::istream &stream)
    {
        int n;
		int m;
        stream >> n >> m;

        this->_n = n;
        this->_m = m;
        this->_a = new double[n * m];

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                stream >> a(i, j);
            }
        }
    }

	Matrix(double (*f)(int i, int j), int n, int m) : Matrix(n, m)
	{
		for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
				a(i, j) = f(i, j);
	}

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

		std::cout << "SES " << n << " " << m << "\n";

        Matrix * compound = new Matrix(n, m);
		std::cout << "SES\n";
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

	Matrix * operator += (Matrix & A)
	{
		if (A._n != _n || A._m != _m)
            return nullptr;

        for (int i = 0; i < _n; ++i)
            for (int j = 0; j < _m; ++j)
                a(i, j) += A.a(i, j);

		return this;
	}

	Matrix * operator -= (Matrix & A)
	{
		if (A._n != _n || A._m != _m)
            return nullptr;

        for (int i = 0; i < _n; ++i)
            for (int j = 0; j < _m; ++j)
                a(i, j) -= A.a(i, j);

		return this;
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

        std::copy(_a, _a + _n * _m,  A->_a);

        return A;
    }

	double mulLines(int i, Matrix * v, int j = 0)
	{
		if (v->_n != _m)
			throw;

		double sum = 0;
		for (int k = 0; k < _m; k++)
			sum += a(i, k) * v->a(k, j);
		
		return sum;
	}

	void setMul(Matrix * A, Matrix * B)
	{
	   	if (A->_m != B->_n || _n != A->_n || _m != B->_m)
			throw;

        for (int i = 0; i < _n; ++i) {
            for (int j = 0; j < _m; ++j) {
                a(i, j) = 0;
                for (int k = 0; k < A->_m; ++k) {
                    a(i, j) += A->a(i, k) * B->a(k, j);
                }
            }
        }
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

Matrix * solveGauss(Matrix * A, Matrix * b)
{
	Matrix * expanded = *A | *b;
	expanded->gauss();
	Matrix * x = expanded->subMatrix(0, A->m(), A->n(), 1);
	delete expanded;
	return x;
}

Matrix * solveIter(Matrix * A, Matrix * b, double w = 1, double eps = 1e-4)
{
	if (A->n() != b->n())
		throw;
	
	Matrix * B = A->copy();
	Matrix * c = b->copy();
	for (int i = 0; i < A->n(); i++)
	{
		B->divLine(i, A->a(i, i));
		c->a(i, 0) /= A->a(i, i);
	}

	Matrix * x = new Matrix(A->m(), 1);	
	Matrix * r = new Matrix(A->m(), 1);
	do
	{
		for (int i = 0; i < x->n(); i++)
			x->a(i, 0) = x->a(i, 0) - w * ( B->mulLines(i, x) - c->a(i, 0) );

		r->setMul(A, x);
		*r -= *b;
		std::cout << r->vectorNormSqrt() << "\n";
	}
	while (r->vectorNormSqrt() > eps);

	return x; 
}

double matGen(int i, int j)
{
	if (i == j)
		return 10;
	if (i + 1 < j)
		return 0;
	return (1.0 / (i + 1));
}

double stolbGen(int i, int j)
{
	return i + 1;
}


int main()
{
	Matrix * A = new Matrix(matGen, 100, 100);
	Matrix * b = new Matrix(stolbGen, 100, 1);

	Matrix * xGauss = solveGauss(A, b);
	std::cout << xGauss;
	//std::cout << (*(*A * *xGauss) -= *b)->vectorNormSqrt() << "\n";

	std::cout << "---------------\n";

	Matrix * xZeydel = solveIter(A, b);
	std::cout << xZeydel;
}

// ------------------------------------------------------------
