#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

class Matrix {
private:
    std::vector<std::vector<double>> mat;
    int rows;
    int cols;

public:
    Matrix(int rows, int cols);
    void input();
    void display() const;
    Matrix inverseGauss() const;
    Matrix inverseCramer() const;
    double determinant(const std::vector<std::vector<double>>& m, int n) const;
    std::vector<std::vector<double>> adjoint() const;
    std::vector<std::vector<double>> getCofactor(const std::vector<std::vector<double>>& m, int p, int q, int n) const;
    Matrix multiply(const Matrix& other) const;
    bool isIdentity() const;
    bool isSquare() const { return rows == cols; }
    void roundSmallValues(double threshold = 1e-8);
};

#endif // MATRIX_H