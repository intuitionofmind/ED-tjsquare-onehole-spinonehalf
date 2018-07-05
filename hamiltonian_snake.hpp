/*
 * This it the file to define the Hamiltonian operoator and related operations.
 */

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

        public:
        int Dim();
        T CoupStren();
        T Dot(T* v1, T* v2);
        void Copy(T* v1, T* v2);
        void SetOne(T* v, int i);
        void SortEval(int n, T* w1, T* w2, std::vector<int>& order);
        void MultVec(T* v, T* w);
        T Correlation(T* v, int i, int j);
        void Marshall(T* v1, T* v2);
        // constructor
        tjSquareHalf(int d, T j) { 
            n = d;
            J = j;
            }  
        };

template<typename T>
int tjSquareHalf<T>::Dim() { return n; }

template<typename T>
T tjSquareHalf<T>::CoupStren() { return J; }

template<typename T>
T tjSquareHalf<T>::Dot(T* v1, T* v2) {
        int len = Dim();
        T res = 0.;
        for (int i = 0; i < len; ++i) { res += std::conj(v1[i])*v2[i]; }
        return res;
        }

template<typename T>
void tjSquareHalf<T>::Copy(T* v1, T* v2) {
        int len = Dim();
        for (int i = 0; i < len; ++i) { v2[i] = v1[i]; }
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
        auto snake = new int[numSite];
        bool vertical = 1;
        for (int j = 0; j < numSiteY; ++j) {
            for (int i = 0; i < numSiteX; ++i) {
                    int k = j*numSiteX+i;
                    if (vertical) { snake[k] = k; }
                    else { snake[k] = j*numSiteX+(numSiteX-i-1); }
                }
            vertical = !vertical;
            }
        // for (int i = 0; i < numSite; ++i) { std::cout << snake[i] << " "; }
        // std::cout << std::endl;

        int len = Dim();
        T J = CoupStren();
        for (int l = 0; l < len; ++l) { w[l] = 0.0; } 
        for (int l = 0; l < len; ++l) {
            if (0.0 == v[l]) { continue; }
            int s = spinBasis[l % subDim]; // convention is i = h*subDim+s  
            int h = holeBasis[int(l/subDim)];
            std::bitset<numSite> spinConfig(s); // the highest bits for the number of (numSite-numHole) are filled with '0' while do not have real meanings but just for convenience
            std::bitset<numSite> config; // initialize config with all 0s

            // fuse hole and spin to give rise to a tJ configuration
            for (int i = 0; i < numSite; ++i) {
                if (i < h) { config[i] = spinConfig[i]; }
                else if (i > h) { config[i] = spinConfig[i-1]; }
                }

            for (int j = 0; j < numSiteY; ++j) {
                for (int i = 0; i <numSiteX; ++i) {
                    T phase = 1.0; // sign for fermion hopping and possible flux added on the boudary 
                    int k = snake[j*numSiteX+i];

                    // along x-direction
                    int kx = snake[j*numSiteX+((i+1) % numSiteX)];
                    if ("PBC" == flagBounX || ((i+1) % numSiteX) > i) {
                        // Heisenberg J term
                        if (k != h && kx != h && config[k] != config[kx]) {
                            w[l] -= 0.5*J*v[l];
                            std::bitset<numSite> temp (config);
                            temp.flip(k);
                            temp.flip(kx);
                            int ss = RetrieveSpin<numSite>(temp, h);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            w[h*subDim+std::distance(spinBasis.begin(), iter)] += 0.5*J*v[l];
                            }

                        // hopping t term
                        else if (h == k) {
                            int hh = kx;
                            int sign = 1; // the sign for sigma t-J model  
                            if (sigma && 0 == config[hh]) { sign = -1; }
                            // if (0 == (kx % numSiteX)) { phase = std::polar(1.0, xFlux); } // In case of flux insertion. 
                            std::bitset<numSite> temp (config);
                            SwapBit<numSite>(temp, k, kx);
                            int ss = RetrieveSpin<numSite>(temp, hh);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            // w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= 1.0*pow(-1.0, k-kx+1)*sign*phase*v[l];
                            w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= 1.0*pow(-1.0, k-kx+1)*sign*v[l];
                            }
                        else if (h == kx) {
                            int hh = k;
                            int sign = 1;
                            if (sigma && 0 == config[hh]) { sign = -1; }
                            // if (0 == (kx % numSiteX)) { phase = std::polar(1.0, xFlux); }
                            std::bitset<numSite> temp (config);
                            SwapBit<numSite>(temp, k, kx);
                            int ss = RetrieveSpin<numSite>(temp, hh);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= 1.0*pow(-1.0, k-kx+1)*sign*phase*v[l];
                            }
                        }

                    // along y-direction
                    int ky = snake[((j+1) % numSiteY)*numSiteX+i];
                    if ("PBC" == flagBounY || ((j+1) % numSiteY) > j) {
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
                            // if (0 == (ky % numSiteY)) { phase = std::polar(1.0, yFlux); }
                            std::bitset<numSite> temp (config);
                            SwapBit<numSite>(temp, k, ky);
                            int ss = RetrieveSpin<numSite>(temp, hh);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            // w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= 1.0*pow(-1.0, k-ky+1)*sign*phase*v[l];
                            w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= 1.0*pow(-1.0, k-ky+1)*sign*v[l];
                            }
                        else if (h == ky) {
                            int hh = k;
                            int sign = 1;
                            if (sigma && 0 == config[hh]) { sign = -1; }
                            // if (0 == (ky % numSiteY)) { phase = std::polar(1.0, yFlux); }
                            std::bitset<numSite> temp (config);
                            SwapBit<numSite>(temp, k, ky);
                            int ss = RetrieveSpin<numSite>(temp, hh);
                            std::vector<int>::iterator iter = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                            // w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= 1.0*pow(-1.0, k-ky+1)*sign*phase*v[l];
                            w[hh*subDim+std::distance(spinBasis.begin(), iter)] -= 1.0*pow(-1.0, k-ky+1)*sign*v[l];
                            }
                        }
                    }
                }
            }
            delete [] snake;
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
                    else if (j < h && 0 == spinConfig[j]) {++n;}
                    else if (j > h && 0 == spinConfig[j-1]) {++n;}
                    }
                }
            v2[i] = pow(-1, n)*pow(-1, h)*v1[i]; // The second pow() only exists for the electron removed with spin down.
            }
        }
