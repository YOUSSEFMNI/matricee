#include "Matrix.h"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

Matrix::Matrix(int rows, int cols) : rows(rows), cols(cols) {
    mat.resize(rows, std::vector<double>(cols));
}

void Matrix::input() {
    std::cout << "Enter elements of " << rows << "x" << cols << " matrix row-wise:" << std::endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cin >> mat[i][j];
        }
    }
}

void Matrix::display() const {
    for (const auto& row : mat) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

Matrix Matrix::inverseGauss() const {
    if (!isSquare()) {
        throw std::logic_error("Non-square matrices cannot be inverted using Gauss.");
    }

    int n = rows;
    std::vector<std::vector<double>> a = mat;
    std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0));

    for (int i = 0; i < n; ++i) {
        inv[i][i] = 1;
    }

    for (int i = 0; i < n; ++i) {
        double diagElement = a[i][i];
        if (diagElement == 0) {
            throw std::logic_error("Matrix is singular and cannot be inverted.");
        }
        for (int j = 0; j < n; ++j) {
            a[i][j] /= diagElement;
            inv[i][j] /= diagElement;
        }
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = a[k][i];
                for (int j = 0; j < n; ++j) {
                    a[k][j] -= factor * a[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }

    Matrix inverseMatrix(n, n);
    inverseMatrix.mat = inv;
    return inverseMatrix;
}

double Matrix::determinant(const std::vector<std::vector<double>>& m, int n) const {
    if (n == 1) {
        return m[0][0];
    }
    if (n == 2) {
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }

    double det = 0;
    for (int x = 0; x < n; x++) {
        std::vector<std::vector<double>> submatrix = getCofactor(m, 0, x, n);
        det += ((x % 2 == 0 ? 1 : -1) * m[0][x] * determinant(submatrix, n - 1));
    }

    return det;
}

std::vector<std::vector<double>> Matrix::adjoint() const {
    if (!isSquare()) {
        throw std::logic_error("Non-square matrices do not have an adjoint.");
    }

    int n = rows;
    if (n == 1) {
        return { {1} };
    }

    std::vector<std::vector<double>> adj(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::vector<std::vector<double>> cofactor = getCofactor(mat, i, j, n);
            adj[j][i] = (determinant(cofactor, n - 1) * ((i + j) % 2 == 0 ? 1 : -1));
        }
    }

    return adj;
}

std::vector<std::vector<double>> Matrix::getCofactor(const std::vector<std::vector<double>>& m, int p, int q, int n) const {
    std::vector<std::vector<double>> temp(n - 1, std::vector<double>(n - 1));
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = m[row][col];
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
    return temp;
}

Matrix Matrix::inverseCramer() const {
    if (!isSquare()) {
        throw std::logic_error("Non-square matrices cannot be inverted using Cramer.");
    }

    int n = rows;
    double det = determinant(mat, n);
    if (det == 0) {
        throw std::logic_error("Matrix is singular and cannot be inverted.");
    }

    std::vector<std::vector<double>> adj = adjoint();
    std::vector<std::vector<double>> inv(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv[i][j] = adj[i][j] / det;
        }
    }

    Matrix inverseMatrix(n, n);
    inverseMatrix.mat = inv;
    return inverseMatrix;
}

Matrix Matrix::multiply(const Matrix& other) const {
    if (cols != other.rows) {
        throw std::logic_error("Matrices dimensions do not match for multiplication.");
    }

    Matrix result(rows, other.cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            result.mat[i][j] = 0;
            for (int k = 0; k < cols; ++k) {
                result.mat[i][j] += mat[i][k] * other.mat[k][j];
            }
        }
    }

    return result;
}

void Matrix::roundSmallValues(double threshold) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (std::abs(mat[i][j]) < threshold) {
                mat[i][j] = 0;
            }
        }
    }
}

bool Matrix::isIdentity() const {
    if (!isSquare()) {
        return false;
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if ((i == j && std::abs(mat[i][j] - 1) >= 1e-8) || (i != j && std::abs(mat[i][j]) >= 1e-8)) {
                return false;
            }
        }
    }

    return true;
}
