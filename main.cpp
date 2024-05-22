#include <iostream>
#include "Matrix.h"

void executeMethod(const Matrix& matrix, const std::string& method) {
    try {
        Matrix inverseMatrix = (method == "gauss") ? matrix.inverseGauss() : matrix.inverseCramer();
        std::cout << "Inverse Matrix using " << ((method == "gauss") ? "Gauss" : "Cramer") << ":" << std::endl;
        inverseMatrix.display();

        // Verify the result
        Matrix product = matrix.multiply(inverseMatrix);
        std::cout << "Product of original matrix and its inverse:" << std::endl;
        product.roundSmallValues();
        product.display();

        if (product.isIdentity()) {
            std::cout << "The inversion is correct. The product is the identity matrix." << std::endl;
        }
        else {
            std::cout << "The inversion is incorrect. The product is not the identity matrix." << std::endl;
        }
    }
    catch (const std::logic_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main() {
    int rows, cols;
    std::cout << "Enter the number of rows of the matrix: ";
    std::cin >> rows;
    std::cout << "Enter the number of columns of the matrix: ";
    std::cin >> cols;

    Matrix matrix(rows, cols);
    matrix.input();
    std::cout << "Original Matrix:" << std::endl;
    matrix.display();

    if (!matrix.isSquare()) {
        std::cerr << "Error: Only square matrices can be inverted." << std::endl;
        return 1;
    }

    std::cout << "Choose inversion method (gauss/cramer): ";
    std::string method;
    std::cin >> method;

    executeMethod(matrix, method);

    return 0;
}