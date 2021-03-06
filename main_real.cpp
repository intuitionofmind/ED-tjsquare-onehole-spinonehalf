#include "head.h"
#include "precondition.hpp"
#include "hamiltonian.hpp"
#include "arssym.h"
// #include "arcomp.h"
// #include "arscomp.h"

int main() {
        std::ofstream file_log("log", std::ios_base::app);
        std::ofstream file_eigvals("eigenvalues.dat", std::ios_base::app | std::ios_base::binary);
        std::ofstream file_eigvecs("eigenvectors.dat", std::ios_base::app | std::ios_base::binary);
        time_t start, end;
        start = time(NULL);
        auto eigValR = new double[numEval];
        auto eigValI = new double[numEval];
        auto eigVec = new double[(numEval+1)*dim];
        auto v1 = new double[dim];
        auto v2 = new double[dim];
        double J = 0.3;
        // PrintHam(J);
        file_log << numSite << std::endl;
        file_log << numSam << std::endl;
        file_log << numEval << std::endl;
        file_log << dim << std::endl;
        file_log << J << std::endl;
        file_log << step << std::endl;
        file_log << sigma << std::endl;

        for (int i = 0; i < numSam; ++i) {
            tjSquareHalf<double> H(dim, J);
            ARSymStdEig<double, tjSquareHalf<double>>
            prob(dim, numEval, &H, &tjSquareHalf<double>::MultVec, "SA");
            int nconv = prob.EigenValVectors(eigVec, eigValR, eigValI);

            for (int j = 0; j < nconv; ++j) {
                file_eigvals.write((char*)(&eigValR[j]), sizeof(double));
                for (int k = 0; k < dim; ++k) {
                    file_eigvecs.write((char*)(&eigVec[j*dim+k]), sizeof(double));
                    double imag = 0.;
                    file_eigvecs.write((char*)(&imag), sizeof(double));
                    }
                }

            // Check the orthogonality of the eigenvectors.
            /* int nPrint = nconv;
            for (int j = 0; j < nPrint; ++j) {
            std::cout << std::setprecision(10) << std::real(eigValR[j]) << std::endl;
            }
            for (int j = 0; j < nPrint; ++j) {
            for (int k = 0; k < nPrint; ++k) {
            for (int l = 0; l < dim; l++) {
            v1[l] = eigVec[j*dim+l];
            v2[l] = eigVec[k*dim+l];
            }
            std::cout << std::abs(H.Dot(v1, v2)) << "  ";
            }
            std::cout << std::endl;
            }
            */
            file_log << (i+1) << std::endl;
            J += step;
            }

        tjSquareHalf<double> H(dim, J);
        ARSymStdEig<double, tjSquareHalf<double>>
        prob(dim, numEval, &H, &tjSquareHalf<double>::MultVec, "SA");
        int nconv = prob.EigenValVectors(eigVec, eigValR, eigValI);

        for (int j = 0; j < nconv; ++j) {
                file_eigvals.write((char*)(&eigValR[j]), sizeof(double));
                for (int k = 0; k < dim; ++k) {
                    file_eigvecs.write((char*)(&eigVec[j*dim+k]), sizeof(double));
                    double imag = 0.;
                    file_eigvecs.write((char*)(&imag), sizeof(double));
                    }
                }

        std::cout << "Now for OTOC." << std::endl;

        int num = 100;
        double timeStep = 0.001;
        int x0 = 1;
        int y0 = 0;
        int x1 = 2;
        int y1 = 1;

        auto v = new arcomplex<double>[dim];
        int s = 0;
        for (int l = 0; l < dim; ++l) {
                arcomplex<double> e (eigValR[s*dim+l], eigValI[s*dim+l])
                v[l] = e;
                }
        // TimeEvolution test.
        tjSquareHalf<arcomplex<double>> A(dim, 1.0);
        auto temp = new arcomplex<double>[dim];
        for (int i = 0; i < num; ++i) {
            int n = i*10;
            A.TimeEvolutionSimple(timeStep, n, v, v1);
            A.TimeEvolutionAG(timeStep, n, v, v2);
            A.VecMinus(v1, v2, temp);
            // std::cout << n << " " << A.Dot(v1, v1) << " " << A.Dot(v2, v2) << " " << A.Dot(temp, temp) << std::endl;
            std::cout << std::setprecision(14) << std::abs(A.Dot(temp, temp)) << std::endl;
            std::cout << std::endl;
            }
/* 
                    for (int i = 0; i < num; ++i) {
                            int n = i*10;
                            double bf = H.TimeSzCommutatorSquare(x0, y0, x1, y1, timeStep, n, v);
                            std::cout << n << " " << bf << std::endl; 
                            }

                    int cut = 100;
                    auto zArray = new double[cut]; // The array for partition function. 
                    for (int i = 0; i < cut; ++i) {
                            double deltaE = std::real(eigValSort[i])-std::real(eigValSort[0]);
                            zArray[i] = std::exp(-1.0*beta*deltaE);
                            }
                    for (int i = 0; i < num; ++i) {
                            int n = i*10;
                            double pf = 0.0; // Partition function. 
                            double tr = 0.0; // Quantity traced over. 
                            for (int j = 0; j < cut; ++j) {
                                double deltaE = std::real(eigValSort[j])-std::real(eigValSort[0]);
                                double z = std::exp(-1.0*beta*deltaE);
                                std::cout << z << std::endl;
                                pf += z;
                                for (int l = 0; l < dim; ++l) {
                                    v[l] = eigVec[order[j]*dim+l];
                                    }
                                double butterfly = H.TimeSzCommutatorSquare(x0, y0, x1, y1, timeStep, n, v);
                                tr += z*butterfly;
                                }
                            tr = tr/pf;
                            std::cout << n << " " << tr << std::endl;
                            }
 */

                    delete [] eigValR;
        delete [] eigValI;
        delete [] eigVec;
        delete [] v1;
        delete [] v2;

        end = time(NULL);

        file_log << "Time: " << (end-start)/60.0 << " min" << std::endl;
        file_eigvals.close();
        file_eigvecs.close();
        file_log.close();
        return 1;
    }

// PrintHam() function is used to print the Hamiltonian matrix explicitly to check whether the code is right or not.
void PrintHam(double JJ) {
        tjSquareHalf<arcomplex<double>> A(dim, JJ);
        arcomplex<double>* v1 = new arcomplex<double>[dim];
        arcomplex<double>* v2 = new arcomplex<double>[dim];
        arcomplex<double>* temp = new arcomplex<double>[dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                A.SetOne(v1, i);
                A.SetOne(v2, j);
                A.MultVec(v2, temp);
                std::cout << A.Dot(v1, temp) << "  ";
                }
            std::cout << std::endl;
            }
        delete [] v1;
        delete [] v2;
        delete [] temp;
        }

