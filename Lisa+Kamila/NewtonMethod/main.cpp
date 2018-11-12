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
        for (int j = 0; j < std::min(_m, _n); ++j)
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

	Matrix * operator /= (double c)
	{
        for (int i = 0; i < _n; ++i)
            for (int j = 0; j < _m; ++j)
                a(i, j) /= c;

		return this;
	}

	Matrix * operator *= (double c)
	{
        for (int i = 0; i < _n; ++i)
            for (int j = 0; j < _m; ++j)
                a(i, j) *= c;

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

	void set(Matrix * A)
	{
		if (A->_n != _n || A->_m != _m)
			throw;

		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				a(i, j) = A->a(i, j);
	}

	Matrix * normalize()
	{
		return *this /= vectorNormSqrt();
	}

	Matrix * newE(int n = 0)
	{
		if (n == 0)
			n = _n;

		return new Matrix([](int i, int j) -> double {return i == j ? 1 : 0;}, n, n); 
	}

	Matrix * inverse()
	{
		auto E = newE();
		auto * AA_1 = *this | *E;
		AA_1->gauss();

		auto A_1 = AA_1->subMatrix(0, _n, _n, _n);

		delete AA_1;
		delete E;

		return A_1;
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

void fillJ(Matrix * x, Matrix * Jx)
{
	Jx->a(0, 0) = 4 - x->a(3, 0);
	Jx->a(0, 1) = -1;
	Jx->a(0, 2) = 1;
	Jx->a(0, 3) = -x->a(0, 0);

	Jx->a(1, 0) = 1;
	Jx->a(1, 1) = -2;
	Jx->a(1, 2) = 3 - x->a(3, 0);
	Jx->a(1, 3) = -x->a(2, 0);

	Jx->a(2, 0) = -1;
	Jx->a(2, 1) = 3 - x->a(3, 0);
	Jx->a(2, 2) = -2;
	Jx->a(2, 3) = -x->a(1, 0);

	Jx->a(3, 0) = 2 * x->a(0, 0);
	Jx->a(3, 1) = 2 * x->a(1, 0);
	Jx->a(3, 2) = 2 * x->a(2, 0);
	Jx->a(3, 3) = 0;
}

void fillF(Matrix * x, Matrix * Fx)
{
	Fx->a(0, 0) = 4*x->a(0, 0)      - x->a(1, 0)      + x->a(2, 0)      - x->a(0, 0) * x->a(3, 0);
	Fx->a(0, 1) = x->a(0, 0)        - 2 * x->a(1, 0)  + 3 * x->a(2, 0)  - x->a(2, 0) * x->a(3, 0);
	Fx->a(0, 2) = -x->a(0, 0)       + 3 * x->a(1, 0)  - 2 * x->a(2, 0)  - x->a(1, 0) * x->a(3, 0);
	Fx->a(0, 3) =  x->a(0, 0)*x->a(0, 0)   + x->a(1, 0)*x->a(1, 0) + x->a(2, 0)*x->a(2, 0) - 1;
}

Matrix * searchNewton(void (*fillJ)(Matrix * x, Matrix * Jx), void (*fillF)(Matrix * x, Matrix * Fx), Matrix * x0, double eps = 1e-4)
{
	Matrix * x = x0->copy();
	Matrix * Jx = new Matrix(4, 4);
	Matrix * Fx = new Matrix(4, 1);

	double df;

	fillF(x, Fx);
	do
	{
		fillJ(x, Jx);

		Matrix * dx = solveGauss(Jx, Fx);
		*x -= *dx;
		delete dx;

		fillF(x, Fx);
		df = Fx->vectorNormSqrt();
		std::cout << df << std::endl;
	}
	while (df > eps);

	delete Jx;
	delete Fx;

	return x;
}



int main()
{
	Matrix * x0 = new Matrix(4, 1);
	x0->a(0, 0) = 0.1;
	x0->a(0, 1) = 0.1;
	x0->a(0, 2) = 0.1;
	x0->a(0, 3) = 0.1;

	Matrix * x = searchNewton(fillJ, fillF, x0, 1e-10);

	std::cout << "-----------------\n" << x;
}

// ------------------------------------------------------------
