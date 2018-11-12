/*
 * This it the file to define the Hamiltonian operoator and related operations.
 */

#include "arcomp.h"

template<int N>
int RetrieveSpin(std::bitset<N> config, int h) {
        std::bitset<N> spin;
        for (int i = 0; i < N; ++i) {
            if (i < h) { spin[i] = config[i]; }
            else { spin[i] = config[i+1]; }
            }
        return spin.to_ulong();
        }

// note that std::swap() function does not apply to std::bitset
template<int N>
void SwapBit(std::bitset<N>& config, int i, int j) {
        int a = config[i];
        config[i] = config[j];
        config[j] = a;
        }

template<typename T>
class tjSquareHalf {
        private:
        int n;
        T J;
        double xFlux;
        double yFlux;

        double alpha;

        public:
        int Dim();
        T CoupStren();
        double XFlux();
        double YFlux();
        double Alpha();
        T Dot(T* v1, T* v2);
        void VecPlus(T* v1, T* v2, T* w);
        void VecMinus(T* v1, T* v2, T* w);
        void MultNumber(T a, T*v);
        void Copy(T* v1, T* v2);
        void Normalize(T* v);
        void SetOne(T* v, int i);
        void SortEval(int n, T* w1, T* w2, std::vector<int>& order);
        void MultVec(T* v, T* w);
        void TimeEvolutionSimple(double step, int num, T* vInit, T* vFinal);
        void TimeEvolutionAG(double step, int num, T* vInit, T* vFinal);
        void Sz(int x, int y, T* v, T* w);
        double TimeSzCommutatorSquare(int x0, int y0, int x1, int y1, double timeStep, int num, T* vec);
        T Correlation(T* v, int i, int j);
        void Marshall(T* v1, T* v2);
        // constructor
        tjSquareHalf(int d, T j, double fx, double fy, double a) { 
            n = d;
            J = j;
            xFlux = fx;
            yFlux = fy;
            alpha = a;
            }  
        };

template<typename T>
int tjSquareHalf<T>::Dim() { return n; }

template<typename T>
T tjSquareHalf<T>::CoupStren() { return J; }

template<typename T>
double tjSquareHalf<T>::XFlux() { return xFlux; }

template<typename T>
double tjSquareHalf<T>::YFlux() { return yFlux; }

template<typename T>
double tjSquareHalf<T>::Alpha() { return alpha; }

template<typename T>
T tjSquareHalf<T>::Dot(T* v1, T* v2) {
        int len = Dim();
        T res = 0.;
        for (int i = 0; i < len; ++i) { res += std::conj(v1[i])*v2[i]; }
        return res;
        }

template<typename T>
void tjSquareHalf<T>::VecPlus(T* v1, T* v2, T* w) {
        int len = Dim();
        for (int l = 0; l < len; ++l) { w[l] = 0.0; } 
        for (int i = 0; i < len; ++i) { w[i]= v1[i]+v2[i]; }
        }

template<typename T>
void tjSquareHalf<T>::VecMinus(T* v1, T* v2, T* w) {
        int len = Dim();
        for (int l = 0; l < len; ++l) { w[l] = 0.0; } 
        for (int i = 0; i < len; ++i) { w[i]= v1[i]-v2[i]; }
        }

template<typename T>
void tjSquareHalf<T>::MultNumber(T a, T* v) {
        int len = Dim();
        for (int i = 0; i < len; ++i) { v[i] = a*v[i]; }
        }

template<typename T>
void tjSquareHalf<T>::Copy(T* v1, T* v2) {
        int len = Dim();
        for (int i = 0; i < len; ++i) { v2[i] = v1[i]; }
        }

template<typename T>
void tjSquareHalf<T>::Normalize(T* v) {
        T z = Dot(v, v);
        double f = sqrt(std::real(z));
        int len = Dim();
        for (int i = 0; i < len; ++i) { v[i] = v[i]/f; }
        }

template<typename T>
void tjSquareHalf<T>::SetOne(T* v, int i) {
        int len = Dim();
        for (int j = 0; j < len; ++j) { v[j] = 0.0; }
        v[i] = 1.0;
        }

/*
* Note that the ArcomStdEig() in ARPACKPP will not sort the eigenvalues you want while ArsymStdEig() does. SortEval() funtion helps to sort the eigenvalues and its order is stored in the vector "order" such that you can access the i'th smallest eigenvalues according to "order[i]" in the original sequence. 
 */

template<typename T>
void tjSquareHalf<T>::SortEval(int n, T* w1, T* w2, std::vector<int>& order) { 
        std::vector<double> vec;
        for (int i = 0; i < n; ++i) { vec.push_back(std::real(w1[i])); }
        std::sort(vec.begin(), vec.end());
        for (int i = 0; i < n; ++i) { w2[i] = T(vec[i]); }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::vector<int>::iterator it = std::find(order.begin(), order.end(), j); 
                if (std::real(w2[i]) == std::real(w1[j]) && it == order.end()) {
                    order.push_back(j);
                    }
                }
            }
        }

/*
 * One dimensional spin-1/2 one hole doped t-J Hamiltonian. Note that there is a 1/4 shift in the Hamiltonian: H = \sum_{\langle i, j\rangle}(S_{i}S_{j}-1/4)+t*c_{i}^{\dagger}c_{j}.
 * 
 * Required by Arpack++ package handbook:
 * There only requirements make by ARPACK++ are that MultVec musth have two pointers to vectors of type T as paraments and the input vector must precede the output vector.
 */

template<typename T>
void tjSquareHalf<T>::MultVec(T* v, T* w) {
        int len = Dim();
        T J = CoupStren();
        double fx = XFlux();
        double fy = YFlux();
        double alpha = Alpha();

        for (int l = 0; l < len; ++l) { w[l] = 0.0; } 
        for (int l = 0; l < len; ++l) {
            if (0.0 == v[l]) { continue; }
            int s = spinBasis[l % subDim]; // convention is i = h*subDim+s  
            int h = holeBasis[int(l/subDim)];
            std::bitset<numSite> spinConfig(s); // the highest bits for the number of (numSite-numHole) are filled with '0' while do not have real meanings but just for convenience
            std::bitset<numSite> config; // initialize config with all 0s

            // Fuse hole and spin to give rise to a tJ configuration.
            int ii = 0;
            for (int i = 0; i < numSite; ++i) {
                if (h == i) { continue; }
                else {
                    config[i] = spinConfig[ii];
                    ++ii;
                    }
                }

            for (int j = 0; j < numSiteY; ++j) {
                for (int i = 0; i <numSiteX; ++i) {
                    T phase = 1.0; // Sign for fermion hopping and possible flux added on the boudary.
                    int k = j*numSiteX+i;

                    // along x-direction
                    int kx = j*numSiteX+((i+1) % numSiteX);
                    if ("PBC" == flagBounX || kx > k) {
                        // Heisenberg J term
                        if (k != h && kx != h && config[k] != config[kx]) {
                            w[l] -= alpha*0.5*J*v[l];
                            std::bitset<numSite> temp (config);
                            temp.flip(k);
                            temp.flip(kx);
                            int ss = RetrieveSpin<numSite>(temp, h);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            w[h*subDim+std::distance(spinBasis.begin(), iter)] += alpha*0.5*J*v[l];
                            }

                        // Hopping t term.
                        else if (h == k) { // Hole's hopping from k to kx.
                            int hh = kx;
                            int sign = 1; // Sign for sigma t-J model.
                            std::complex<double> phaseX(1.0, 0.0);
                            if (sigma && 0 == config[hh]) { sign = -1; }
                            if (0 == (kx % numSiteX)) { phaseX = std::polar(1.0, fx); } // In case of flux insertion. 
                            std::bitset<numSite> temp (config);
                            SwapBit<numSite>(temp, k, kx);
                            int ss = RetrieveSpin<numSite>(temp, hh);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= alpha*1.0*pow(-1.0, k-kx+1)*sign*phaseX*v[l];
                            }
                        else if (h == kx) {  // Hole's hopping from kx to k. 
                            int hh = k;
                            int sign = 1;
                            std::complex<double> phaseX(1.0, 0.0);
                            if (sigma && 0 == config[hh]) { sign = -1; }
                            if (0 == (kx % numSiteX)) { phaseX = std::polar(1.0, -1.0*fx); }
                            std::bitset<numSite> temp (config);
                            SwapBit<numSite>(temp, k, kx);
                            int ss = RetrieveSpin<numSite>(temp, hh);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= alpha*1.0*pow(-1.0, kx-k+1)*sign*phaseX*v[l];
                            }
                        }

                    // along y-direction
                    int ky = ((j+1) % numSiteY)*numSiteX+i;
                    if ("PBC" == flagBounY || ky > k) {
                        if (k != h && ky != h && config[k] != config[ky]) {
                            w[l] -= 0.5*J*v[l];
                            std::bitset<numSite> temp (config);
                            temp.flip(k);
                            temp.flip(ky);
                            int ss = RetrieveSpin<numSite>(temp, h);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            w[h*subDim+std::distance(spinBasis.begin(), iter)] += 0.5*J*v[l];
                            }

                        else if (h == k) {
                            int hh = ky;
                            int sign = 1;
                            if (sigma && 0 == config[hh]) { sign = -1; }
                            std::complex<double> phaseY(1.0, 0.0);
                            if (0 == (ky % numSiteY)) { phaseY = std::polar(1.0, fy); }
                            std::bitset<numSite> temp (config);
                            SwapBit<numSite>(temp, k, ky);
                            int ss = RetrieveSpin<numSite>(temp, hh);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= 1.0*pow(-1.0, k-ky+1)*sign*phaseY*v[l];
                            }
                        else if (h == ky) {
                            int hh = k;
                            int sign = 1;
                            if (sigma && 0 == config[hh]) { sign = -1; }
                            std::complex<double> phaseY(1.0, 0.0);
                            if (0 == (ky % numSiteY)) { phase = std::polar(1.0, -1.0*fy); }
                            std::bitset<numSite> temp (config);
                            SwapBit<numSite>(temp, k, ky);
                            int ss = RetrieveSpin<numSite>(temp, hh);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= 1.0*pow(-1.0, k-ky+1)*sign*phaseY*v[l];
                            }
                        }
                    }
                }
            }
        }

template<typename T>
void tjSquareHalf<T>::TimeEvolutionSimple(double step, int num, T* vInit, T* vFinal) {
        arcomplex<double> imagUnit (0.0, 1.0);
        auto v = new arcomplex<double>[dim];
        auto w = new arcomplex<double>[dim];
        auto temp = new arcomplex<double>[dim];

        for (int i = 0; i < num; ++i) {
            if (0 == i) {
                Copy(vInit, w);
                Copy(w, v);
                }
            else {
                MultVec(v, temp);
                MultNumber(imagUnit*step, temp);
                VecMinus(v, temp, w);
                Normalize(w);
                Copy(w, v);
                }
            }

        Copy(w, vFinal);
        delete [] v;
        delete [] w;
        delete [] temp;
        }

// Integrate the time evolution operator U=exp(-i*H*t) by num steps with Askar-Goldberg method.
template<typename T>
void tjSquareHalf<T>::TimeEvolutionAG(double step, int num, T* vInit, T* vFinal) {
        arcomplex<double> imagUnit (0.0, 1.0);
        // In askar-Goldberg method, u, v, w denote n-1, n, n+1 states.
        auto u = new arcomplex<double>[dim];
        auto v = new arcomplex<double>[dim];
        auto w = new arcomplex<double>[dim];
        auto ww = new arcomplex<double>[dim];
        auto temp = new arcomplex<double>[dim];
        auto tempTemp = new arcomplex<double>[dim];

        for (int i = 0; i < num; ++i) {
            if (0 == i) {
                Copy(vInit, w);
                Copy(w, v);
                }
            else if (1 == i) {
            // else {
                MultVec(v, temp);
                MultNumber(imagUnit*step, temp);
                VecMinus(v, temp, w);
                Normalize(w);
                // Move a step forward.
                Copy(v, u);
                Copy(w, v);
                }

            else {
                MultVec(v, temp);
                MultVec(temp, tempTemp);
                MultNumber(2.0*imagUnit*step, temp);
                VecMinus(u, temp, ww);
                MultVec(tempTemp, temp);
                MultNumber(0.5*imagUnit*std::pow(step, 3), temp);
                VecPlus(ww, temp, w);
                Copy(v, u);
                Copy(w, v);
                }
 
            }

        Copy(w, vFinal);
        delete [] u;
        delete [] v;
        delete [] w;
        delete [] temp;
        delete [] tempTemp;
        }

template<typename T>
void tjSquareHalf<T>::Sz(int x, int y, T* v, T* w) {
        int r = y*numSiteX+x;
        int len = Dim();
        for (int l = 0; l < len; ++l) { w[l] = 0.0; } 

        for (int l = 0; l < len; ++l) {
            int s = spinBasis[l % subDim];
            int h = holeBasis[int(l/subDim)];
            std::bitset<numSite> spinConfig(s);
            std::bitset<numSite> config; 

            for (int i = 0; i < numSite; ++i) {
                if (i < h) { config[i] = spinConfig[i]; }
                else if (i > h) { config[i] = spinConfig[i-1]; }
                }

            if (h != r && 1 == config[r]) { w[l] = 0.5*v[l]; }
            else if (h != r && 0 == config[r]) { w[l] =-0.5*v[l]; }
            }
        }

template<typename T>
double tjSquareHalf<T>::TimeSzCommutatorSquare(int x0, int y0, int x1, int y1, double timeStep, int num, T* vec) {
        auto v = new arcomplex<double>[dim];
        auto w = new arcomplex<double>[dim];
        int len = Dim();
        arcomplex<double> zero (0.0 ,0.0);
        for (int i = 0; i < len; ++i) {
            v[i] = zero;
            w[i] = zero;
            }
/* 
        // w(t)vvw(t)
        TimeEvolutionSimple(timeStep, num, vec, w);
        Sz(x1, y1, w, v);
        TimeEvolutionSimple(-1.0*timeStep, num, v, w);
        Sz(x0, y0, w, v);
        Sz(x0, y0, v, w);
        TimeEvolutionSimple(timeStep, num, w, v);
        Sz(x1, y1, v, w);
        TimeEvolutionSimple(-1.0*timeStep, num, w, v);
        double p0 = std::abs(Dot(vec, v));
        // vw(t)w(t)v
        Sz(x0, y0, vec, w);
        TimeEvolutionSimple(timeStep, num, w, v);
        Sz(x1, y1, v, w);
        Sz(x1, y1, w, v);
        TimeEvolutionSimple(-1.0*timeStep, num, v, w);
        Sz(x0, y0, w, v);
        double p1 = std::abs(Dot(vec, v));
        // w(t)vw(t)v
        Sz(x0, y0, vec, w);
        TimeEvolutionSimple(timeStep, num, w, v);
        Sz(x1, y1, v, w);
        TimeEvolutionSimple(-1.0*timeStep, num, w, v);
        Sz(x0, y0, v, w);
        TimeEvolutionSimple(timeStep, num, w, v);
        Sz(x1, y1, v, w);
        TimeEvolutionSimple(-1.0*timeStep, num, w, v);
        double p2 = std::abs(Dot(vec, v));
        // vw(t)vw(t)
        TimeEvolutionSimple(timeStep, num, vec, w);
        Sz(x1, y1, w, v);
        TimeEvolutionSimple(-1.0*timeStep, num, v, w);
        Sz(x0, y0, w, v);
        TimeEvolutionSimple(timeStep, num, v, w);
        Sz(x1, y1, w, v);
        TimeEvolutionSimple(-1.0*timeStep, num, v, w);
        Sz(x0, y0, w, v);
        double p3 = std::abs(Dot(vec, v));
 */
        // w(t)vvw(t)
        TimeEvolutionAG(timeStep, num, vec, w);
        Sz(x1, y1, w, v);
        TimeEvolutionAG(-1.0*timeStep, num, v, w);
        Sz(x0, y0, w, v);
        Sz(x0, y0, v, w);
        TimeEvolutionAG(timeStep, num, w, v);
        Sz(x1, y1, v, w);
        TimeEvolutionAG(-1.0*timeStep, num, w, v);
        double p0 = std::abs(Dot(vec, v));
        // vw(t)w(t)v
        Sz(x0, y0, vec, w);
        TimeEvolutionAG(timeStep, num, w, v);
        Sz(x1, y1, v, w);
        Sz(x1, y1, w, v);
        TimeEvolutionAG(-1.0*timeStep, num, v, w);
        Sz(x0, y0, w, v);
        double p1 = std::abs(Dot(vec, v));
        // w(t)vw(t)v
        Sz(x0, y0, vec, w);
        TimeEvolutionAG(timeStep, num, w, v);
        Sz(x1, y1, v, w);
        TimeEvolutionAG(-1.0*timeStep, num, w, v);
        Sz(x0, y0, v, w);
        TimeEvolutionAG(timeStep, num, w, v);
        Sz(x1, y1, v, w);
        TimeEvolutionAG(-1.0*timeStep, num, w, v);
        double p2 = std::abs(Dot(vec, v));
        // vw(t)vw(t)
        TimeEvolutionAG(timeStep, num, vec, w);
        Sz(x1, y1, w, v);
        TimeEvolutionAG(-1.0*timeStep, num, v, w);
        Sz(x0, y0, w, v);
        TimeEvolutionAG(timeStep, num, v, w);
        Sz(x1, y1, w, v);
        TimeEvolutionAG(-1.0*timeStep, num, v, w);
        Sz(x0, y0, w, v);
        double p3 = std::abs(Dot(vec, v));
        double p = -p0-p1+p2+p3;
        // std::cout << p << std::endl;
        delete [] v;
        delete [] w;

        return p;
    }

template<typename T>
T tjSquareHalf<T>::Correlation(T* v, int i, int j) {
        int len = Dim();
        T res = 0.0;
        for (int k = 0; k < len; ++k) {
            int sk = spinBasis[k % subDim];
            int hk = holeBasis[int(k/subDim)];
            for (int l = 0; l < len; ++l) {
                int sl = spinBasis[l % subDim];
                int hl = holeBasis[int(l/subDim)];
                if (j == hk && i == hl && sk == sl) {
                    res += pow(-1, hk+hl+1)*std::conj(v[k])*v[l];
                    }
                }
            }
        return res;
        }

// Change the groud state into the Marshall sign basis.
template<typename T>
void tjSquareHalf<T>::Marshall(T* v1, T* v2) { 
        int len = Dim();
        for (int i = 0; i < len; ++i) {
            int s = spinBasis[i % subDim];
            int h = holeBasis[int(i/subDim)];
            std::bitset<numSite> spinConfig(s);
            int n = 0;
            for (int j = 0; j < numSite; ++j){
                if ( 0 == j % 2) {
                    if (j == h) {++n;} // Only for the electron removed with spin down.
                    else if (j < h && 0 == spinConfig[j]) { ++n; }
                    else if (j > h && 0 == spinConfig[j-1]) { ++n; }
                    }
                }
            v2[i] = pow(-1, n)*pow(-1, h)*v1[i]; // The second pow() only exists for the electron removed with spin down.
            }
        }
