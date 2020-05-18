//
// Created by Matthew Nicoletti on 2020-04-17.
//

#ifndef SIMULATIONS_GRAPH_H
#define SIMULATIONS_GRAPH_H

#include <random>
#include <vector>
#include <complex>




using namespace std;

class CylinderGraph {
public:

    vector<int> graph;
    int X;
    int Y;

    // this is for testing purposes
    int k;

    // The time at which we want approximately measure (D, phi) at
    int T0;

    int RED = 1;
    int BLUE = 0;

    // top = RED
    // bottom = BLUE

    int topBlueY;    // The highest y coordinate of blue
    int bottomRedY;  // the lowest y coordinate of red


    vector<complex<double>> M;  // the martingale
    vector<double> MQ; // values of |M(t+1) - Mt|^2
    vector<complex<double>> Ylist;  // // the value of e^(2 pi k/N (T - T0)/N)M(T)
    vector<double> YQ; // values of |Y(t+1) - Yt|^2
    vector<complex<double>> Dphilist; // trajectory of Dphi over time
    vector<double> DphiQ; // |Dphi(t+1) - Dphi(t)|^2


    vector<int> numRed;   // the number of Red for every y.
    vector<int> numBlue;  // the number of Blue for every y.

    vector<int> fluctuations;  // the y coordinate of x = 0 as a timeseries. The
    // first non Blue coordinate from bottom.

    vector<double> correlations;  // correlations[i] gives E[(h[0] -
    // E[h[0])]*(h[i] - E[h[i]])], where h[i] is the
    // height of the interface at x = i.

    vector<vector<double>> dists;


    CylinderGraph(int N, int M, int k, int T0);

    void initializeGraph();

    int h(int x);  // height of the interface at x.
    int getColor(int x, int y);
    int setColor(int x, int y, int color);

    void setTopBlueY();
    void setBottomRedY();

    void MarkovChain(int num_samples, int interval);

    double computeReturnProb(int M, int x);

    void computeReturnDists();

    void printGraph();

    // This function is for easier visualization. It computes a new cylinder of dimension X/scale x Y/scale
    // where each matrix element corresponds to a scale*scale sized block of the original cylinder, and is
    // equal to either "1", "." or "0" depending on the proportion of 1 and 0 in the block.
    void printDensityGraph(int scale);

    vector<int> randomWalk(int startX, int startY, int startColor,
                           std::uniform_real_distribution<double>& dist,
                           std::mt19937& mt);
};


#endif //SIMULATIONS_GRAPH_H
