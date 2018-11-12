const int numSiteX = 6;
const int numSiteY = 2;
const int numSite = numSiteX*numSiteY;
const int numHole = 1;
const int numEval = 20; // Desired eigenvalues to be found by ARPACK++.
const int numSam = 100;
const double step = 1.0;
// const double xFlux = 0.0*PI;
// const double yFlux = 0.0*PI;
const double zMag = 0.5;
const std::string flagBounX = "PBC";
const std::string flagBounY = "OBC";
const bool sigma = false; // Used to switch between t-J model and \sigma-t-J model.
const double beta = 20.0; // Inverse temperature.

std::vector<int> spinBasis;
std::vector<int> holeBasis;
const int subDim = SpinBasisConstruct(spinBasis);
const int holeDim = HoleBasisConstruct(holeBasis); 
const int dim = subDim*holeDim;

/*
BasisConstruct() function is to construct spin basis with conserved total Sz and return the dimension of the Hilbert space.
*/
int SpinBasisConstruct(std::vector<int> &vec) {
        int max = int(pow(2, numSiteX*numSiteY-1)); // Possible maxmium dimension of the Hilbert space for the sub-space in terms of pure spin configuration with the hole fixed.
        int dim = 0;
        for (int i = 0; i < max; ++i) {
            std::bitset<numSite> b(i);
            if (b.count() == (unsigned int)(zMag+(numSite-numHole)/2.0)) {
                vec.push_back(i);
                dim++;
                }
            }
        return dim;
        }

int HoleBasisConstruct(std::vector<int> &vec) {
        int dim = 0;
        for (int j = 0; j < numSiteY; ++j) {
            for (int i = 0; i < numSiteX; ++i) {
                vec.push_back(j*numSiteX+i);
                dim++;
                }
            }
        return dim;
        }
