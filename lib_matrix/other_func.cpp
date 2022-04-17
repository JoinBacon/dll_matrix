#include "pch.h"
#include "matrix_class.h"

Matrix HadamardProduct(Matrix c1, Matrix c2) {
	if (c1.GetLine() != c2.GetLine() || c1.GetCol() != c2.GetCol()) {
		throw Exception1("Error: the matrices have different sizes.");
	}
	Matrix H(c1.GetLine(), c1.GetCol());
	for (int i = 0; i < c1.GetLine(); ++i) {
		for (int j = 0; j < c1.GetCol(); ++j) {
			H.GetElem(i, j, H) = c1.GetElem(i, j, c1) * c2.GetElem(i, j, c2);
		}
	}
	return H;
}

double Tr(Matrix c1) {
	if (c1.GetLine() != c1.GetCol())
	{
		throw Exception1("Error: the matrix must be square");
	}
	double S = 0;
	for (int i = 0; i < c1.GetLine(); ++i) {
		for (int j = 0; j < c1.GetCol(); ++j) {
			if (i == j) {
				S += c1.GetElem(i, j, c1);
			}
		}
	}
	return S;
}

double Determinate(Matrix c1) {
	if (c1.GetLine() != c1.GetCol()) {
		throw Exception1("Error: to find the determinant, the matrix must be square");
	}
	int swapcount = 0;
	double determinate = 1;
	int n = c1.GetLine();
	for (int i = 0; i < n - 1; ++i) {
		double maxcol = fabs(c1.GetElem(i, i, c1));
		int imax = i;
		for (int j = i + 1; j < n; ++j) {
			if (c1.GetElem(j, i, c1) > maxcol) {
				maxcol = fabs(c1.GetElem(j, i, c1));
				imax = j;
			}
		}
		if (i != imax) {
			for (int k = 0; k < n; ++k) {
				double tmp = c1.GetElem(i, k, c1);
				c1.GetElem(i, k, c1) = c1.GetElem(imax, k, c1);
				c1.GetElem(imax, k, c1) = tmp;
			}
			++swapcount;
		}
		for (int l = i + 1; l < n; ++l) {
			if (c1.GetElem(i, i, c1) != 0) {
				double coef = -c1.GetElem(l, i, c1) / c1.GetElem(i, i, c1);
				for (int m = 0; m < c1.GetLine(); ++m) {
					c1.GetElem(l, m, c1) += c1.GetElem(i, m, c1) * coef;
				}
			}
		}
	}
	for (int i = 0; i < n; ++i) {
		determinate *= c1.GetElem(i, i, c1);
	}
	if (swapcount % 2 == 0 || determinate == 0) {
		return determinate;
	}
	else {
		return -determinate;
	}
}

double Scalar(Matrix c1, Matrix c2) {
	double Sc = 0;
	if ((c1.GetLine() != 1) & (c2.GetCol() != 1)) {
		throw Exception1("Error: dot product is only possible for vectors.");
	}
	if (c1.GetCol() != c2.GetLine()) {
		throw Exception1("Error: with such sizes of vectors, the dot product is impossible.");;
	}
	for (int i = 0; i < c1.GetCol(); ++i) {
		Sc += c1.GetElem(0, i, c1) * c2.GetElem(i, 0, c2);
	}
	return Sc;
}

double EuclideanNorm(Matrix c1) {
	double norm = 0;
	if ((c1.GetLine() != 1) & (c1.GetCol() != 1)) {
		throw Exception1("Error: euclidean norm is only possible for vectors.");
	}
	if (c1.GetLine() == 1) {
		for (int i = 0; i < c1.GetCol(); ++i) {
			norm += c1.GetElem(0, i, c1) * c1.GetElem(0, i, c1);
		}
	}
	else {
		for (int i = 0; i < c1.GetLine(); ++i) {
			norm += c1.GetElem(i, 0, c1) * c1.GetElem(i, 0, c1);
		}
	}
	return sqrt(norm);
}

double MaximumNorm(Matrix c1) {
	double max = 0;
	if ((c1.GetLine() != 1) & (c1.GetCol() != 1)) {
		throw Exception1("Error: maximum norm is only possible for vectors.");
	}
	if (c1.GetLine() == 1) {
		max = fabs(c1.GetElem(0, 0, c1));
		for (int i = 1; i < c1.GetCol(); ++i) {
			if (max < fabs(c1.GetElem(0, i, c1)))
			{
				max = fabs(c1.GetElem(0, i, c1));
			}
		}
	}
	else {
		max = fabs(c1.GetElem(0, 0, c1));
		for (int i = 0; i < c1.GetLine(); ++i) {
			if (max < fabs(c1.GetElem(i, 0, c1)))
			{
				max = fabs(c1.GetElem(i, 0, c1));
			}
		}
	}
	return max;
}

double FrobeniusNorm(Matrix c1) {
	double s = 0;
	for (int i = 0; i < c1.GetLine(); ++i) {
		for (int j = 0; j < c1.GetCol(); ++j) {
			s += c1.GetElem(i, j, c1) * c1.GetElem(i, j, c1);
		}
	}
	return sqrt(s);
}

double Angle(Matrix c1, Matrix c2) {
	double scalar = Scalar(c1, c2);
	double norm1 = EuclideanNorm(c1);
	double norm2 = EuclideanNorm(c2);
	double angle = scalar / (norm1 * norm2);
	return acos(angle);
}

Matrix ReverseM(Matrix A) {
	if (A.GetLine() != A.GetCol())
	{
		throw Exception1("Error: the matrix must be square");
	}
	if (!Determinate(A)) {
		throw Exception1("Error: the determinate of matrix equals 0");
	}
	int n = A.GetLine();
	double N1 = 0, Ninf = 0;
	Matrix A0 = A;
	for (int i = 0; i < n; ++i) {
		double colsum = 0, linesum = 0;
		for (int j = 0; j < n; ++j) {
			linesum += fabs(A0.GetElem(i, j, A0));
			colsum += fabs(A0.GetElem(j, i, A0));
		}
		N1 = max(colsum, N1);
		Ninf = max(linesum, Ninf);
	}
	A0 = Transpose(A0);
	A0 = A0 * (1 / (N1 * Ninf));
	EMatrix E2(n);
	for (int i = 0; i < n; ++i) {
		E2.GetElem(i, i, E2) = 2 * E2.GetElem(i, i, E2);
	}
	Matrix inv = A0;
	double EPS = 0.0001;
	while (fabs(Determinate(A * inv) - 1) >= EPS) {
		Matrix prev = inv;
		inv = A * prev;
		inv = (-1) * inv;
		inv = inv + E2;
		inv = prev * inv;
	}
	return inv;
}


Matrix Transpose(Matrix c1) {
	Matrix c2(c1.GetCol(), c1.GetLine());
	for (int i = 0; i < c1.GetCol(); ++i) {
		for (int j = 0; j < c1.GetLine(); ++j) {
			c2.GetElem(i, j, c2) = c1.GetElem(j, i, c1);
		}
	}
	return c2;
}

int Rank(Matrix c1) {
	const double EPS = 1E-9;
	int n = c1.GetLine();
	int m = c1.GetCol();
	int rank = max(n, m);
	std::vector<char> line_used(n);
	for (int i = 0; i < m; ++i) {
		int j = 0;
		for (j = 0; j < n; ++j) {
			if (!line_used[j] && fabs(c1.GetElem(j, i, c1)) > EPS) break;
		}
		if (j == n) {
			--rank;
		}
		else
		{
			line_used[j] = true;
			for (int p = i + 1; p < m; ++p) {
				c1.GetElem(j, p, c1) /= c1.GetElem(j, i, c1);
			}
			for (int k = 0; k < n; ++k) {
				if (k != j && fabs(c1.GetElem(k, i, c1)) > EPS) {
					for (int p = i + 1; p < m; ++p) {
						c1.GetElem(k, p, c1) -= c1.GetElem(j, p, c1) * c1.GetElem(k, i, c1);
					}
				}
			}
		}
	}
	return rank;
}
