#include <iostream>
#include <Eigen/Dense>
#include <fstream>

// 4 * 2^k divisions
void createMesh(int k) {
    int numVertices = 4 * std::pow(2,k);
    Eigen::MatrixXd vertices = Eigen::MatrixXd::Zero(numVertices,2);
    Eigen::MatrixXi elements = Eigen::MatrixXi::Zero(numVertices,2);

    elements.col(0) = Eigen::VectorXi::LinSpaced(numVertices,0,numVertices-1);
    elements.col(1) = Eigen::VectorXi::LinSpaced(numVertices,1,numVertices);
    elements(numVertices-1,1) = 0;

    elements += Eigen::MatrixXi::Constant(numVertices,2,1);

    Eigen::Vector2d SW(-0.25,-0.25);
    Eigen::Vector2d SE(0.25,-0.25);
    Eigen::Vector2d NE(0.25,0.25);
    Eigen::Vector2d NW(-0.25,0.25);

    for (int i = 0 ; i < numVertices/4 ; ++i ) {
        // Bottom edge
        vertices.row(i) = SW + 4. * ((double)i)/ (double)numVertices  * (SE - SW);
        // Right edge
        vertices.row(numVertices/4 + i) = SE + 4. * ((double)i)/ (double)numVertices  * (NE - SE);
        // Top edge
        vertices.row(2 * numVertices/4 + i) = NE + 4. * ((double)i)/ (double)numVertices  * (NW - NE);
        // Left edge
        vertices.row(3* numVertices/4 + i) = NW + 4. * ((double)i)/ (double)numVertices  * (SW - NW);
    }
    std::ofstream coordinates_out("capacitor" + std::to_string(k) + "_coordinates" + ".dat");
    std::ofstream elements_out("capacitor" + std::to_string(k) + "_elements" + ".dat");

    coordinates_out << vertices << std::endl;
    elements_out << elements << std::endl;

}

int main() {
    for (int i = 0 ; i < 5 ; ++i) {
        createMesh(i);
    }
    return 0;
}