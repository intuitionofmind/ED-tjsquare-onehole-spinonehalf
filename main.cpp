#include "head.h"
#include "precondition.hpp"
#include "hamiltonian.hpp"
// #include "arssym.h"
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
/* 
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

                file_log << (i+1) << std::endl;
                J += 0.1;
                }
            }
 */
        tjSquareHalf<arcomplex<double>> H(dim, J);
        ARCompStdEig<double, tjSquareHalf<arcomplex<double>>> prob;
        prob.DefineParameters(dim, numEval, &H, &tjSquareHalf<arcomplex<double>>::MultVec, "SR");
        int nconv = prob.EigenValVectors(eigVec, eigVal);
        std::vector<int> order;
        H.SortEval(nconv, eigVal, eigValSort, order);

        std::cout << "Now for OTOC." << std::endl;

        int numTime = 1000;
        double timeStep = 0.01;
        int x0 = 1;
        int y0 = 1;
        int x1 = 2;
        int y1 = 1;

        auto v0 = new arcomplex<double>[dim];
        auto v = new arcomplex<double>[dim];
        auto w = new arcomplex<double>[dim];
        for (int i = 0; i < numTime; ++i) {
                int j = 0;
                for (int l = 0; l < dim; ++l) {
                    v0[l] = eigVec[order[j]*dim+l];
                    }
                // w(t)vvw(t)
                H.TimeEvolution(timeStep, i, v0, w);
                H.Sz(x1, y1, w, v);
                H.TimeEvolution(-1.0*timeStep, i, v, w);
                H.Sz(x0, y0, w, v);
                H.Sz(x0, y0, v, w);
                H.TimeEvolution(timeStep, i, w, v);
                H.Sz(x1, y1, v, w);
                H.TimeEvolution(-1.0*timeStep, i, w, v);
                double p0 = std::abs(H.Dot(v0, v));
                // vw(t)w(t)v
                H.Sz(x0, y0, v0, w);
                H.TimeEvolution(timeStep, i, w, v);
                H.Sz(x1, y1, v, w);
                H.Sz(x1, y1, w, v);
                H.TimeEvolution(-1.0*timeStep, i, v, w);
                H.Sz(x0, y0, w, v);
                double p1 = std::abs(H.Dot(v0, v));
                // w(t)vw(t)v
                H.Sz(x0, y0, v0, w);
                H.TimeEvolution(timeStep, i, w, v);
                H.Sz(x1, y1, v, w);
                H.TimeEvolution(-1.0*timeStep, i, w, v);
                H.Sz(x0, y0, v, w);
                H.TimeEvolution(timeStep, i, w, v);
                H.Sz(x1, y1, v, w);
                H.TimeEvolution(-1.0*timeStep, i, w, v);
                double p2 = std::abs(H.Dot(v0, v));
                // vw(t)vw(t)
                H.TimeEvolution(timeStep, i, v0, w);
                H.Sz(x1, y1, w, v);
                H.TimeEvolution(-1.0*timeStep, i, v, w);
                H.Sz(x0, y0, w, v);
                H.TimeEvolution(timeStep, i, v, w);
                H.Sz(x1, y1, w, v);
                H.TimeEvolution(-1.0*timeStep, i, v, w);
                H.Sz(x0, y0, w, v);
                double p3 = std::abs(H.Dot(v0, v));

                double p = p0+p1-p2-p3;
                std::cout << i << " " << p << std::endl;
                }
        delete [] v0;
        delete [] v;
        delete [] w;

        delete [] eigVal;
        delete [] eigValSort;
        delete [] eigVec;
        delete [] v1;
        delete [] v2;

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

