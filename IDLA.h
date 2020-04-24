//
// Created by Matthew Nicoletti on 2020-04-18.
//

#ifndef SIMULATIONS_IDLA_H
#define SIMULATIONS_IDLA_H




#include <random>
#include <vector>
#include <complex>


using namespace std;


class IDLA {

public:

    vector<int> graph;
    int X;
    int H;
    // this is for testing purposes
    int k;

    int Y;

    int RED = 1;
    int EMPTY = 0;

    // one species = RED

    int bottomEmptyY;    // The lowest y coordinate of empty square
    int topRedY;    // The highest y coordinate of red


    vector<complex<double>> M;  // the martingale

    vector<complex<double>> Mvals;  // the value of M_{y0.N^2} for each sample run

    vector<int> numRed;   // the number of Red for every y.

    vector<int> fluctuations;  // the y coordinate of x = 0 as a timeseries. The
    // first non Blue coordinate from bottom.

    vector<double> correlations;  // correlations[i] gives E[(h[0] -
    // E[h[0])]*(h[i] - E[h[i]])], where h[i] is the
    // height of the interface at x = i.

    vector<complex<double>> Fourier;

    vector<complex<double>> Dphivals;  // list of samples of Dphi --
    // one at a fixed time, for each separate run of IDLA

    IDLA(int N, int H, int k);

    void initializeGraph();

    int h(int x);  // height of the interface at x.
    int getColor(int x, int y);
    int setColor(int x, int y, int color);

    void setTopRedY();
    void setBottomEmptyY();

    void MarkovChain(int num_samples, int interval);

    void printGraph();

    // This function is for easier visualization. It computes a new cylinder of dimension X/scale x Y/scale
    // where each matrix element corresponds to a scale*scale sized block of the original cylinder, and is
    // equal to either "1", "." or "0" depending on the proportion of 1 and 0 in the block.
//    void printDensityGraph(int scale);

    vector<int> randomWalk(int startX, int startY, int startColor,
                           std::uniform_real_distribution<double>& dist,
                           std::mt19937& mt);
};


#endif //SIMULATIONS_IDLA_H
