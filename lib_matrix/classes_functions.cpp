#include "pch.h"
#include "matrix_class.h"

Matrix::Matrix(void) :
	line(0), column(0)
{

}

Matrix::Matrix(int m, int n) :
	line(m), column(n)
{
	M.resize(line, std::vector<double>(column));
}

Matrix::Matrix(std::vector<std::vector<double>> vals) {
	for (int i = 0; (unsigned)i < vals.size(); i++) {
		if (vals[0].size() != vals[i].size()) {
			throw Exception1("Error");
		}
	}
	line = vals.size();
	column = vals[0].size();
	M = vals;
}
int Matrix::GetLine() const {
	return line;
}

int Matrix::GetCol() const {
	return column;
}

double& Matrix::GetElem(const int i, const int j, const Matrix& c1) {
	return M[i][j];
}

Matrix operator+(const Matrix& c1, const Matrix& c2) {

	if (c1.line != c2.line || c1.column != c2.column) {
		throw Exception1("Error: the matrices have different sizes.");
	}
	Matrix Sum(c1.line, c1.column);
	for (int i = 0; i < c1.line; ++i) {
		for (int j = 0; j < c1.column; ++j) {
			Sum.M[i][j] = c1.M[i][j] + c2.M[i][j];
		}
	}
	return Sum;
}

Matrix operator-(const Matrix& c1, const Matrix& c2) {
	if (c1.line != c2.line || c1.column != c2.column) {
		throw Exception1("Error: the matrices have different sizes.");
	}
	Matrix Dif(c1.line, c1.column);
	for (int i = 0; i < c1.line; ++i) {
		for (int j = 0; j < c1.column; ++j) {
			Dif.M[i][j] = c1.M[i][j] - c2.M[i][j];
		}
	}
	return Dif;
}

Matrix operator*(const Matrix& c1, const Matrix& c2) {
	if (c1.column != c2.line) {
		throw Exception1("Error: the matrices have different sizes.");
	}
	Matrix Res(c1.line, c2.column);
	for (int i = 0; i < c1.line; ++i) {
		for (int j = 0; j < c2.column; ++j) {
			for (int l = 0; l < c1.column; ++l) {
				Res.M[i][j] += c1.M[i][l] * c2.M[l][j];
			}
		}
	}
	return Res;
}

Matrix operator*(const Matrix& c1, const double value) {
	Matrix Res(c1.line, c1.column);
	for (int i = 0; i < c1.line; ++i) {
		for (int j = 0; j < c1.column; ++j) {
			Res.M[i][j] = c1.M[i][j] * value;
		}
	}
	return Res;
}

Matrix operator*(const double value, const Matrix& c1) {
	Matrix Res(c1.line, c1.column);
	for (int i = 0; i < c1.line; ++i) {
		for (int j = 0; j < c1.column; ++j) {
			Res.M[i][j] = value * c1.M[i][j];
		}
	}
	return Res;
}


std::ostream& operator<<(std::ostream& ostr, const Matrix& c1) {
	for (int i = 0; i < c1.line; ++i)
	{
		for (int j = 0; j < c1.column; ++j) {
			ostr << c1.M[i][j] << "\t";
		}
		ostr << "\n";
	}
	return(ostr);
}
std::ofstream& operator<<(std::ofstream& f, Matrix& c1) {
	for (int i = 0; i < c1.line; ++i)
	{
		for (int j = 0; j < c1.column - 1; ++j)
		{
			f.precision(3); ///!!!
			f << c1.M[i][j] << "\t";
		}
		f << c1.M[i][c1.column - (size_t)1] << "\n";
	}
	return f;
}

std::ifstream& operator>>(std::ifstream& f, Matrix& c1) {
	std::vector<double> values;
	int columns = 0, lines = 0;

	double a;
	while (f >> a)
	{
		values.push_back(a);
		if (f.peek() == '\n' || f.peek() == EOF)
		{
			if (columns != 0 && columns != values.size())
			{
				throw Exception1("Error: the different shapes of lines in matrix");
			}
			++lines;
			columns = values.size();
			c1.M.push_back(values);
			values.clear();
		}
	}
	c1.line = lines;
	c1.column = columns;
	return f;
}





void Matrix::BinaryWrite(Matrix& c1, std::string NameFile) const {
	std::ofstream file(NameFile, std::ios::binary | std::ios::app | std::ios_base::app);
	if (!file)
	{
		throw Exception1("file open error");
	}
	int lines = c1.GetCol();
	int columns = c1.GetCol();
	file.write((char*)&lines, sizeof(int));
	file.write((char*)&columns, sizeof(int));
	for (int i = 0; i < c1.GetLine(); ++i)
	{
		for (int j = 0; j < c1.GetCol(); ++j)
		{
			file.write((char*)&c1.M[i][j], sizeof(double));
		}
	}
	file.close();
}
void Matrix::BinaryRead(Matrix& c1, std::string NameFile) {
	std::ifstream f(NameFile, std::ios::binary | std::ios_base::app);
	f.read((char*)&c1.line, sizeof(int));
	f.read((char*)&c1.column, sizeof(int));
	c1.M.resize(c1.line, std::vector<double>(c1.column));
	if (!f) {
		throw Exception1("file open error");
	}
	else
	{
		for (int i = 0; i < c1.line; ++i)
		{
			for (int j = 0; j < c1.column; ++j)
			{
				f.read((char*)&c1.M[i][j], sizeof(double));
			}
		}
	}
	f.close();
}
EMatrix::EMatrix(int m) :
	Matrix(m, m)
{
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			if (i == j)
				M[i][j] = 1;
			else
				M[i][j] = 0;
		}
	}
}


DiagonalMatrix::DiagonalMatrix(Matrix c1) :
	Matrix(c1.GetLine(), c1.GetLine())
{
	if (c1.GetLine() != c1.GetCol()) {
		throw Exception1("diagonal: matrix must be square");
	}
	for (int i = 0; i < c1.GetLine(); ++i) {
		for (int j = 0; j < c1.GetLine(); ++j) {
			if (i == j) {
				M[i][j] = c1.GetElem(i, j, c1);
			}
			else {
				M[i][j] = 0;
			}
		}
	}

}
DiagonalMatrix::DiagonalMatrix(int m) :
	Matrix(m, m)
{

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			if (i == j) {
				M[i][j] = (double)(rand() % 10000) / 1000;
			}
			else {
				M[i][j] = 0;
			}
		}
	}
}
UpperTriangular::UpperTriangular(Matrix c1) :
	Matrix(c1.GetLine(), c1.GetLine())
{
	if (c1.GetLine() != c1.GetCol()) {
		throw Exception1("upper triangle: matrix must be square");
	}
	for (int i = 0; i < c1.GetLine(); ++i) {
		for (int j = 0; j < c1.GetLine(); ++j) {
			if (i >= j) {
				M[i][j] = c1.GetElem(i, j, c1);
			}
			else {
				M[i][j] = 0;
			}
		}
	}
}

UpperTriangular::UpperTriangular(int m) :
	Matrix(m, m)
{
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			if (i >= j) {
				M[i][j] = (double)(rand() % 10000) / 1000;
			}
			else {
				M[i][j] = 0;
			}
		}
	}
}

LowerTriangular::LowerTriangular(Matrix c1) :
	Matrix(c1.GetLine(), c1.GetLine())
{
	if (c1.GetLine() != c1.GetCol()) {
		throw Exception1("upper triangle: matrix must be square");
	}
	for (int i = 0; i < c1.GetLine(); ++i) {
		for (int j = 0; j < c1.GetLine(); ++j) {
			if (i <= j) {
				M[i][j] = c1.GetElem(i, j, c1);
			}
			else {
				M[i][j] = 0;
			}
		}
	}
}
LowerTriangular::LowerTriangular(int m) :
	Matrix(m, m)
{
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			if (i <= j) {
				M[i][j] = (double)(rand() % 10000) / 1000;
			}
			else {
				M[i][j] = 0;
			}
		}
	}
}
SymmetricalMatrix::SymmetricalMatrix(int m) :
	Matrix(m, m)
{
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			if (j >= i) {
				M[i][j] = (double)(rand() % 10000) / 1000;
			}
			else {
				M[i][j] = M[j][i];
			}
		}
	}
}

//Matrix::Matrix(const std::initializer_list <std::vector <double>> &list):
//	Matrix(2, 2)
//{
//	int k = 0;
//	for (auto& v : list) {
//		for (auto& elem : v) {
//			int m = 0;
//			M[k][m] = elem;
//			++m;
//			
//		}
//		++k;
//	}
//		
//
//}

Matrix RSA::center(Matrix M) {
	int I = M.GetLine();
	int J = M.GetCol();
	Matrix M1 = M;
	for (int j = 0; j < J; ++j) {
		double sum_col = 0;
		for (int i = 0; i < I; ++i) {
			sum_col += M1.GetElem(i, j, M1);
		}
		double mj = sum_col / I;
		for (int i = 0; i < I; ++i) {
			M1.GetElem(i, j, M1) -= mj;
		}
	}
	return M1;
}

Matrix RSA::scaling(Matrix M) {
	int I = M.GetLine();
	int J = M.GetCol();
	Matrix M1 = M;
	for (int j = 0; j < J; ++j) {
		double sum_col = 0;
		for (int i = 0; i < I; ++i) {
			sum_col += M1.GetElem(i, j, M1);
		}
		double mj = sum_col / I;
		double sj = 0;
		for (int i = 0; i < I; ++i) {
			sj += pow(M1.GetElem(i, j, M1) - mj, 2);
		}
		sj = sqrt(sj / (static_cast<unsigned __int64>(I) - 1));
		for (int i = 0; i < I; ++i) {
			M1.GetElem(i, j, M1) = (M1.GetElem(i, j, M1) - mj) / sj;
		}
	}
	return M1;
}

void RSA::NIPALS(Matrix M) {
	double eps = 1e-8;
	int PC = min(M.GetLine(), M.GetCol());
	Matrix D = M;
	D = center(D);
	D = scaling(D);
	Matrix E = D;
	Matrix P(E.GetCol(), PC);
	Matrix T(E.GetLine(), PC);
	for (int h = 1; h <= PC; ++h) {
		Matrix t(E.GetLine(), 1);
		for (int j = 0; j < E.GetLine(); ++j) {
			t.GetElem(j, 0, t) = E.GetElem(j, h - 1, E);
		}
		Matrix d(E.GetLine(), 1);
		Matrix p(E.GetCol(), 1);
		do {
			p = Transpose((Transpose(t) * E) * (1 / Scalar(Transpose(t), t)));
			p = p * (1 / EuclideanNorm(p));
			Matrix t_old = t;
			t = (E * p) * (1 / Scalar(Transpose(p), p));
			d = t_old - t;
		} while (EuclideanNorm(d) > eps);
		E = E - t * Transpose(p);
		for (int i = 0; i < P.GetLine(); ++i) {
			P.GetElem(i, h - 1, P) = p.GetElem(i, 0, p);
		}

		for (int i = 0; i < T.GetLine(); ++i) {
			T.GetElem(i, h - 1, T) = t.GetElem(i, 0, p);
		}
	}
	std::ofstream scores_file("scores.txt");
	std::ofstream loadings_file("loadings.txt");
	std::ofstream remains_file("remains.txt");
	if (scores_file.is_open()) {
		scores_file << T;
	}
	if (loadings_file.is_open()) {
		loadings_file << P;
	}
	if (remains_file.is_open()) {
		remains_file << D - T * Transpose(P);
	}
	scores_file.close();
	loadings_file.close();
	remains_file.close();
}

std::vector<double> RSA::leverage(Matrix T) {
	std::vector<double> la;
	std::vector<double> Hi;
	Matrix G = Transpose(T) * T;
	for (int i = 0; i < G.GetLine(); ++i) {
		la.push_back(G.GetElem(i, i, G));
	}
	for (int i = 0; i < T.GetLine(); ++i) {
		double hi = 0;
		for (int j = 0; j < T.GetCol(); ++j) {
			hi += pow(T.GetElem(i, j, T), 2) / la[j];
		}
		Hi.push_back(hi);
	}
	return Hi;
}

void RSA::deviation(Matrix M, Matrix T, Matrix P) {
	std::vector<double> v;
	double v0 = 0, sum = 0;
	Matrix X = M;
	X = center(X);
	X = scaling(X);
	Matrix E = X - T * Transpose(P);
	for (int i = 0; i < E.GetLine(); ++i) {
		double vi = 0;
		for (int j = 0; j < E.GetCol(); ++j) {
			vi += pow(E.GetElem(i, j, E), 2);
		}
		v.push_back(vi);
	}
	for (auto& vii : v) {
		v0 += vii;
	}
	v0 /= M.GetLine();
	double TRV = v0 / M.GetCol();
	for (int i = 0; i < X.GetLine(); ++i) {
		for (int j = 0; j < X.GetCol(); ++j) {
			sum += pow(X.GetElem(i, j, X), 2);
		}
	}
	double ERV = 1 - ((M.GetLine() * v0) / sum);

}