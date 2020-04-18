//
// Created by Matthew Nicoletti on 2020-04-17.
//


#include "Graph.h"

#include <chrono>
#include <iostream>
#include <fstream>




using namespace std;

// x is expected to be within 0 and X-1
// y is expect to be within 0 and Y-1

// assumes X is power of 2
// assumes Y is power of 2

CylinderGraph::CylinderGraph(int X, int Y) {
    CylinderGraph::X = X;
    CylinderGraph::Y = Y;
    CylinderGraph::topBlueY = Y - 1;
    CylinderGraph::bottomRedY = 0;
}

void CylinderGraph::initializeGraph() {



    graph = vector<int>(X * Y, BLUE);
    numRed = vector<int>(Y, 0);
    numBlue = vector<int>(Y, X);
    correlations = vector<double>(X,0);

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(0, 1);



    // start with nonsense
    /*
    for (int y = 0; y < Y; y++) {
        for (int x = 0; x < X; x++) {
            setColor(x, y, (x + y) % 2);

        }
    }

    */

  // start with limit shape

    CylinderGraph::topBlueY = Y/2 - 1;
    CylinderGraph::bottomRedY = Y/2;

    for(int y = 0; y < Y; y++) {
        for (int x = 0; x < X; x++) {
            if(y >= Y/2){ setColor(x, y, RED); }
        }
    }

    setTopBlueY();
    setBottomRedY();


}

void CylinderGraph::setTopBlueY() {
    if ((topBlueY + 1 <= Y - 1) and (numBlue[topBlueY + 1] != 0)) {
        topBlueY++;
    }  // assumes that topBlueY can only increase by one.
    while ((topBlueY > 0) and (numBlue[topBlueY] == 0)) {
        topBlueY--;
    }
}

void CylinderGraph::setBottomRedY() {
    if ((bottomRedY - 1 >= 0) and (numRed[bottomRedY - 1] != 0)) {
        bottomRedY--;
    }  // assumes that bottomRedY can only decrease by one.
    while ((bottomRedY < Y - 1) and (numRed[bottomRedY] == 0)) {
        bottomRedY++;
    }
}

int CylinderGraph::getColor(int x, int y) {
    if ((y < 0) or (y >= Y)) {
        return -1;
    }
    return graph[y * X + ((x % X) + X) % X];
}

int CylinderGraph::setColor(int x, int y, int color) {
    if ((y < 0) or (y >= Y)) {
        return -1;
    }
    if (color == BLUE) {
        numBlue[y]++;
        numRed[y]--;
    } else {
        numBlue[y]--;
        numRed[y]++;
    }

    graph[y * X + x] = color;
    return 1;
}

std::vector<int> CylinderGraph::randomWalk(
        int startX, int startY, int startColor,
        std::uniform_real_distribution<double>& dist, std::mt19937& mt) {
    int x = startX;
    int y = startY;

    char next_neighbor;
    double rand_neighbor;
    while (graph[y * X + x] == startColor) {
        rand_neighbor = dist(mt);
        if (y == Y - 1) {
            if (rand_neighbor <= 0.3333) {
                next_neighbor = 'L';
            } else if (rand_neighbor <= 0.6666) {
                next_neighbor = 'D';
            } else {
                next_neighbor = 'R';
            }
        } else if (y == 0) {
            if (rand_neighbor <= 0.3333) {
                next_neighbor = 'L';
            } else if (rand_neighbor <= 0.6666) {
                next_neighbor = 'U';
            } else {
                next_neighbor = 'R';
            }
        } else {
            if (rand_neighbor <= 0.25) {
                next_neighbor = 'L';
            } else if (rand_neighbor <= 0.5) {
                next_neighbor = 'U';
            } else if (rand_neighbor <= 0.75) {
                next_neighbor = 'D';
            } else {
                next_neighbor = 'R';
            }
        }

        if (next_neighbor == 'U') {
            y++;
        }
        if (next_neighbor == 'D') {
            y--;
        }
        if (next_neighbor == 'R') {
            x = (x + 1) % X;
        }
        if (next_neighbor == 'L') {
            x = (x - 1 + X) % X;
        }
    }
    return {x, y};
}

void CylinderGraph::MarkovChain(int num_samples, int interval) {
    const double pi = std::acos(-1);
    const std::complex<double> I(0, 1);
    double N = (double) X;

    int k = 25;
    double q = acosh(2 - cos(2 * pi * k / X));


    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(0, X - 1);

    std::random_device rdRW;
    std::mt19937 mtRW(rd());
    std::uniform_real_distribution<double> distRW(0, 1);

    int x;
    std::vector<int> bdest, rdest;  // TODO make point
    int err;
    complex<double> Mt;


    // initialize martingale
    M.push_back(0);

    for (int sample = 0; sample < num_samples; sample++) {
        std::cout<<"sample "<<sample<<"..."<<std::endl;
        for (int i = 0; i < interval; i++) {

            // bottom
            x = dist(mt);
            bdest = randomWalk(x, max(bottomRedY - 1, 0), BLUE, distRW, mtRW);
            err = setColor(bdest[0], bdest[1], BLUE);
            if (err == -1) {
                break;
            }




            // top
            x = dist(mt);
            rdest = randomWalk(x, min(topBlueY + 1, Y - 1), RED, distRW, mtRW);
            err = setColor(rdest[0], rdest[1], RED);
            if (err == -1) {
                break;
            }

            // update martingale
            Mt = M[sample * interval + i];

            Mt = Mt + exp(I * (2 * pi * bdest[0] * k) / N) * exp(q * (bdest[1] - Y / 2)) / ((double) X)-
                          exp(I * (2 * pi * rdest[0] * k) / N) * exp(q * (rdest[1] - Y / 2)) / ((double) X);

            if (std::abs(real(Mt)) < 0.0001)
                Mt.real(0);
            if (std::abs(imag(Mt)) < 0.0001)
                Mt.imag(0);

            M.emplace_back( Mt );

            setTopBlueY();
            setBottomRedY();
        }
        std::cout<<"done."<<std::endl;

        // printDensityGraph(32);

        fluctuations.push_back(h(0));

        for (int x = 0; x < X; x++) {
            correlations[x] += (h(x) - Y / 2.0) * (h(0) - Y / 2.0) / num_samples;
        }
    }
}

int CylinderGraph::h(int x) {
    int red_interface = bottomRedY;
    while ((red_interface < Y - 1) and (graph[X * red_interface + x] == BLUE)) {
        red_interface++;
    }

    return red_interface;
}

void CylinderGraph::printDensityGraph(int scale) {
    for (int yb = 0; yb < Y / scale; yb++) {
        for (int xb = 0; xb < X / scale; xb++) {
            double s = 0;
            for (int y = 0; y < scale; y++) {
                for (int x = 0; x < scale; x++) {
                    s += getColor(xb * scale + x, yb * scale + y);
                }
            }
            s /= scale * scale;

            if (s < .1)
                std::cout << 1 << " ";
            else if (s > .9)
                std::cout << 0 << " ";
            else
                std::cout << "."
                          << " ";
        }
        std::cout << std::endl;
    }
}

void CylinderGraph::printGraph() {
    for (int y = Y - 1; y >= 0; y--) {
        for (int x = 0; x < X; x++) {
            int c = getColor(x, y);
            std::cout << c;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


