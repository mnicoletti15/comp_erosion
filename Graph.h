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

    int RED = 1;
    int BLUE = 0;

    // top = RED
    // bottom = BLUE

    int topBlueY;    // The highest y coordinate of blue
    int bottomRedY;  // the lowest y coordinate of red

    vector<complex<double>> M;  // the martingale

    vector<int> numRed;   // the number of Red for every y.
    vector<int> numBlue;  // the number of Blue for every y.

    vector<int> fluctuations;  // the y coordinate of x = 0 as a timeseries. The
    // first non Blue coordinate from bottom.

    vector<double> correlations;  // correlations[i] gives E[(h[0] -
    // E[h[0])]*(h[i] - E[h[i]])], where h[i] is the
    // height of the interface at x = i.

    CylinderGraph(int N, int M);

    void initializeGraph();

    int h(int x);  // height of the interface at x.
    int getColor(int x, int y);
    int setColor(int x, int y, int color);

    void setTopBlueY();
    void setBottomRedY();

    void MarkovChain(int num_samples, int interval);

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
