#pragma once

#define MATHLIBRARY_API  __declspec( dllexport )

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>


class Exception1 : public std::exception {
	std::string text;
public:
	MATHLIBRARY_API  Exception1(std::string str) : text(str) {};
	const char* what() const throw () {
		return text.c_str();
	}
};

class Matrix {
protected:
	int line;
	int column;
	std::vector <std::vector <double>> M;
public:
	MATHLIBRARY_API Matrix(void);
	MATHLIBRARY_API Matrix(int, int);
	MATHLIBRARY_API Matrix(std::vector<std::vector<double>>);
	/*Matrix(const std::initializer_list<std::vector <double>>& list);*/
	MATHLIBRARY_API int GetLine() const;
	MATHLIBRARY_API int GetCol() const;
	MATHLIBRARY_API double& GetElem(const int, const int, const Matrix&);
	MATHLIBRARY_API friend Matrix operator+(const Matrix&, const Matrix&);
	MATHLIBRARY_API friend Matrix operator-(const Matrix&, const Matrix&);
	MATHLIBRARY_API friend Matrix operator*(const Matrix&, const Matrix&);
	MATHLIBRARY_API friend Matrix operator*(const Matrix&, const double);
	MATHLIBRARY_API friend Matrix operator*(const double, const Matrix&);
	MATHLIBRARY_API friend std::ostream& operator<<(std::ostream&, const Matrix&);
	MATHLIBRARY_API friend std::ofstream& operator<<(std::ofstream&, Matrix&);
	MATHLIBRARY_API friend std::ifstream& operator>>(std::ifstream&, Matrix&);
	MATHLIBRARY_API void BinaryWrite(Matrix&, std::string) const;
	MATHLIBRARY_API void BinaryRead(Matrix&, std::string);
};
MATHLIBRARY_API Matrix HadamardProduct(Matrix, Matrix);
MATHLIBRARY_API double Tr(Matrix);
MATHLIBRARY_API double Determinate(Matrix);
MATHLIBRARY_API double Scalar(Matrix, Matrix);
MATHLIBRARY_API double EuclideanNorm(Matrix);
MATHLIBRARY_API double MaximumNorm(Matrix);
MATHLIBRARY_API double FrobeniusNorm(Matrix);
MATHLIBRARY_API int Rank(Matrix);
MATHLIBRARY_API double Angle(Matrix, Matrix);
MATHLIBRARY_API Matrix ReverseM(Matrix);
MATHLIBRARY_API Matrix Transpose(Matrix);

class EMatrix : public Matrix {
public:
	MATHLIBRARY_API EMatrix(int);
};

class DiagonalMatrix : public Matrix {
public:
	MATHLIBRARY_API DiagonalMatrix(Matrix);
	MATHLIBRARY_API DiagonalMatrix(int);
};

class UpperTriangular : public Matrix {
public:
	MATHLIBRARY_API UpperTriangular(Matrix);
	MATHLIBRARY_API UpperTriangular(int);
};

class LowerTriangular : public Matrix {
public:
	MATHLIBRARY_API LowerTriangular(Matrix);
	MATHLIBRARY_API LowerTriangular(int);
};

class SymmetricalMatrix : public Matrix {
public:
	MATHLIBRARY_API  SymmetricalMatrix(int);
};

 class RSA {
public:
	MATHLIBRARY_API Matrix center(Matrix);
	MATHLIBRARY_API Matrix scaling(Matrix);
	MATHLIBRARY_API	void NIPALS(Matrix);
	MATHLIBRARY_API std::vector<double> leverage(Matrix);
	MATHLIBRARY_API void deviation(Matrix, Matrix, Matrix);
};