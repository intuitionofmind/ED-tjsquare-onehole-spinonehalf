#include "head.h"
#include "precondition.hpp"
#include "hamiltonian.hpp"
//#include "arssym.h"
#include "arcomp.h"
#include "arscomp.h"

int main() {
/*        
        arcomplex<double> a (1.0, 2.0);  // Test arcomplex class. 
        arcomplex<double> b (3.0, 4.0);
        std::cout << a*b << std::endl;
        std::cout << std::real(a*b) << std::endl;
        std::cout << std::imag(a*b) << std::endl;
        std::cout << std::cos(flux) << std::endl;
        std::cout << std::sin(flux) << std::endl;
        std::cout << std::polar(1.0, flux) << std::endl;
*/
        std::ofstream file_log("log", std::ios_base::app);
        std::ofstream file_eigvals("eigenvalues.dat", std::ios_base::app | std::ios_base::binary);
        std::ofstream file_eigvecs("eigenvectors.dat", std::ios_base::app | std::ios_base::binary);
        std::ofstream file_eigvecs_translated("eigenvectors_translated.dat", std::ios_base::app | std::ios_base::binary);
        time_t start, end;
        start = time(NULL);
        auto eigVal = new arcomplex<double>[numEval];
        auto eigValSort = new arcomplex<double>[numEval];
        auto eigVec = new arcomplex<double>[numEval*dim];
        auto v1 = new arcomplex<double>[dim];
        auto v2 = new arcomplex<double>[dim];
        // double J = 0.3333333333333333;
        double J = 0.3;
        // PrintHam(J);
        file_log << numSite << std::endl;
        file_log << numSam << std::endl;
        file_log << numEval << std::endl;
        file_log << dim << std::endl;
        file_log << J << std::endl;
        file_log << step << std::endl;
        file_log << sigma << std::endl;

//        auto vv1 = new arcomplex<double>[dim];
//        auto vv2 = new arcomplex<double>[dim];
//        auto vv3 = new arcomplex<double>[dim];
//        auto vv4 = new arcomplex<double>[dim];

        for (int i = 0; i < numSam; ++i) {
            std::cout << i << std::endl;
            tjSquareHalf<arcomplex<double>> H(dim, J);
            ARCompStdEig<double, tjSquareHalf<arcomplex<double>>> prob;
            prob.DefineParameters(dim, numEval, &H, &tjSquareHalf<arcomplex<double>>::MultVec, "SR");
            int nconv = prob.EigenValVectors(eigVec, eigVal);
            std::vector<int> order;
            H.SortEval(nconv, eigVal, eigValSort, order);
            for (int j = 0; j < nconv; ++j) {
                double r = std::real(eigValSort[j]);
                file_eigvals.write((char*)(&r), sizeof(double));
                if (0 == i) { std::cout << "eval: " << std::setprecision(14) << (eigValSort[j])*3.0 << std::endl; }
                for (int k = 0; k < dim; ++k) {
                    v1[k] = eigVec[order[j]*dim+k];
                    }
//                if (flagBoun == "PBC") {H.Translation(v1, v2);}
                for (int k = 0; k < dim; ++k) {v2[k] = v1[k];}
                for (int k = 0; k < dim; ++k) {
                    double rr = std::real(v1[k]);
                    double ii = std::imag(v1[k]);
                   if (0 == i && 0 == j) { std::cout << rr << " " << ii << std::endl; }
                    file_eigvecs.write((char*)(&rr), sizeof(double));
                    file_eigvecs.write((char*)(&ii), sizeof(double));
                    rr = std::real(v2[k]);
                    ii = std::imag(v2[k]);
                    file_eigvecs_translated.write((char*)(&rr), sizeof(double));
                    file_eigvecs_translated.write((char*)(&ii), sizeof(double));
                    }

                // Check the Marshall sign rule.
/*                if (0 == j && 4 == i) {
                    auto x0 = v1[0];
                    for (int k = 0; k < dim; ++k) {v1[k] = std::polar(1.0, -1.0*std::arg(x0))*v1[k];}
                    H.Marshall(v1, v2);
                    for (int k = 0; k < dim; ++k) {std::cout <<  "Marshall sign: " << v1[k] << "  " << v2[k] << std::endl;}
                    }
*/

                // Check the correlation function.
/*
                if (0 == j && i == 10) {
                    arcomplex<double> t (0.0, 0.0);
                    for (int k = 0; k < numSite; ++k) {
                        for (int l = 0; l < numSite; ++l) {
                            if (l == k) {t += H.Correlation(v1, k, l);}
                            std::cout << H.Correlation(v1, k, l) << " ";
                            }
                        std::cout << std::endl;
                        }
                    std::cout << "Trace: " << t << std::endl;
                    }
*/

                // Check the commutation relation between T and H.
/*
                if (2 == j) { 
                    std::cout << H.Dot(v1, v1) << " " << H.Dot(v1, v2)  << std::endl;
                    std::cout << H.Dot(v2, v1) << " " << H.Dot(v2, v2)  << std::endl;

                    arcomplex<double>*  v3 = new arcomplex<double>[dim];  
                    arcomplex<double>*  v4 = new arcomplex<double>[dim];
                    H.MultVec(v2, v3); // Compute HT.
                    H.MultVec(v1, v2); // compute TH.
                    H.Translation(v2, v4);
                    for (int k = 0; k < dim; ++k) {std::cout <<  "Commutation check: " << (v3[k]-v4[k]) << std::endl;}
                    delete [] v3;
                    delete [] v4;
                    }
*/

                // Check the translational cyclic. 
/*
                if (0 == j) {
                    arcomplex<double>*  v = new arcomplex<double>[dim];
                    H.Copy(v1, v);
                    for (int k = 0; k < numSite; ++k) {
                        H.Translation(v1, v2);
                        H.Copy(v2, v1);
                        for (int l = 0; l < dim; ++l) {std::cout << v1[l] << " ";}
                        std::cout << std::endl << std::endl;
                        }
                    for (int k = 0; k < dim; ++k) {std::cout << "Translation check: " << (v[k]-v2[k]) << std::endl;}
                    delete [] v;
                    }
*/
/*                if (0 == j) {
                    H.Copy(v1, vv1);
                    H.Copy(v2, vv2);
                    }
                if (1 == j) {
                    H.Copy(v1, vv3);
                    H.Copy(v2, vv4);
                    }*/
                }
//            std::cout << "translation matrix in the gegenerate space:" << std::endl;
//            std::cout << H.Dot(vv1, vv2) << " " << H.Dot(vv1, vv4) << std::endl;
//            std::cout << H.Dot(vv3, vv2) << " " << H.Dot(vv3, vv4) << std::endl;
//            for (int j = 0; j < dim; ++j) {std::cout << vv3[j] << " ";}
//            std::cout << std::endl;
            file_log << (i+1) << std::endl;
            J += 0.1;
            }
        delete [] eigVal;
        delete [] eigValSort;
        delete [] eigVec;
        delete [] v1;
        delete [] v2;

//        delete [] vv1;
//        delete [] vv2;
//        delete [] vv3;
//        delete [] vv4;
 
        end = time(NULL);

        file_log << "Time: " << (end-start)/60.0 << " min" << std::endl;
        file_eigvals.close();
        file_eigvecs.close();
        file_eigvecs_translated.close();
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

