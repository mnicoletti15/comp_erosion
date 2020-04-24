//
// Created by Matthew Nicoletti on 2020-04-18.
//

#include "IDLA.h"

#include <chrono>
#include <iostream>
#include <fstream>




using namespace std;

// x is expected to be within 0 and X-1
// y is expect to be within 0 and Y-1

// assumes X is power of 2
// assumes Y is power of 2

IDLA::IDLA(int X, int H, int k) {
    IDLA::X = X;
    IDLA::H = H;
    IDLA::k = k;
    IDLA::topRedY = 0;
    IDLA::bottomEmptyY = 0;
    IDLA::correlations = vector<double>(X,0);
    Fourier = vector<complex<double>>();
}

void IDLA::initializeGraph() {

    Y = 2 * X;
    graph = vector<int>(X * Y, EMPTY);
    numRed = vector<int>(Y, 0);


    topRedY = 0;
    bottomEmptyY = 0;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(0, 1);




}

void IDLA::setBottomEmptyY() {
    while (numRed[bottomEmptyY] == X) {
        bottomEmptyY++;
    }
}

void IDLA::setTopRedY() {
    //assumes it is already correct in call to this function
    if (numRed[topRedY + 1] != 0) {
        topRedY++;
    }
}

int IDLA::getColor(int x, int y) {
    if ((y < 0) or (y >= Y)) {
        return -1;
    }
    return graph[y * X + ((x % X) + X) % X];
}

int IDLA::setColor(int x, int y, int color) {
    if ((y < 0) or (y >= Y)) {
        return -1;
    }
    if (color != RED) {
        return -1;
    } else {
        numRed[y]++;
    }

    graph[y * X + x] = color;
    return 1;
}

std::vector<int> IDLA::randomWalk(
        int startX, int startY, int startColor,
        std::uniform_real_distribution<double>& dist, std::mt19937& mt) {
    int x = startX;
    int y = startY;

    char next_neighbor;
    double rand_neighbor;
    while (graph[y * X + x] == RED) {
        rand_neighbor = dist(mt);
        if (y == 0) {
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

void IDLA::MarkovChain(int num_samples, int interval) {

    /*
     * This function takes in a number of samples, and an interval size
     * Runs IDLA for time T = 'interval' in each sample
     * at the end, measures values of relevant quantities at time T
     * then resets each quantity to 0, and repeats.
     * does thus num_samples number of times times.
     */



    const double pi = std::acos(-1);
    const std::complex<double> I(0, 1);
    double N = (double) X;

    double q = acosh(2 - cos(2 * pi * k / X));

    double Eh = 0;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(0, X - 1);

    std::random_device rdRW;
    std::mt19937 mtRW(rd());
    std::uniform_real_distribution<double> distRW(0, 1);

    int x;
    std::vector<int> dest;  // TODO make point
    int err;

    // This list is meant to be the actual martingale
    M.push_back(0);

    for (int sample = 0; sample < num_samples; sample++) {
        // re-initialize graph
        initializeGraph();


        // initialize martingale
        complex<double> Mt = 0;
        complex<double> fourier = 0;
        complex<double> Dphi = 0;

        // notification
        std::cout<<"sample "<<sample<<"..."<<std::endl;
        for (int i = 0; i < interval; i++) {

            // bottom
            x = dist(mt);
            dest = randomWalk(x, max(bottomEmptyY - 1, 0), RED, distRW, mtRW);

//            // get rid of current contribution to fourier mode at dest[0]
//            fourier -= (h(dest[0]) - Eh) * exp(I * (2 * pi * dest[0] * k) / N);

//            int newh = max(dest[1], h(dest[0]));

            err = setColor(dest[0], dest[1], RED);
            if (err == -1) {
                break;
            }


            // update expected height
            Eh = (double) i / N;



            // update martingale

            Mt = Mt + exp(I * (2 * pi * dest[0] * k) / N) * exp(q * (dest[1] - H)) / ((double) X);




//            fourier += -1.0 / N * fourier + (newh-Eh)  * exp(I * (2 * pi * dest[0] * k) / N);




//            if (i % X == N/2) {
//
//                fourier = 0;
//                for (int x = 0; x < X; x++) {
//                    fourier += exp(I * (2 * pi * x * k) / N) * (double) h(x);
//                }
//
//                Fourier.emplace_back(fourier);
//
//
//
//                for (int x = 0; x < X; x++) {
//                    correlations[x] += ((double) (h(x) - Eh) * (h(0) - Eh)) / ((double) (2 * N));
//                }
//            }








            setBottomEmptyY();
            setTopRedY();
        }


        fourier = 0;
        for (x = 0; x < X; x++) {
            fourier += exp(I * (2 * pi * x * k) / N) * (double) h(x);
        }
        Fourier.push_back(fourier);


        Dphi = 0;
        for (int y = bottomEmptyY; y <= topRedY; y++) {
            for (int x = 0; x < X; x++) {
                if (getColor(x, y) == RED) {
                    Dphi += exp(I * (2 * pi * x * k) / N)  / (double) X;
                }
            }
        }
        Dphivals.push_back(Dphi);

//                M.emplace_back(Mt);



        if (std::abs(real(Mt)) < 0.0001)
            Mt.real(0);
        if (std::abs(imag(Mt)) < 0.0001)
            Mt.imag(0);
        Mvals.emplace_back( Mt );



        std::cout<<"done."<<std::endl;

        // printDensityGraph(32);

//        fluctuations.push_back(h(0));

        for (int x = 0; x < X; x++) {
            correlations[x] += ((double) (h(x) - H) * (h(0) - H)) / ((double) num_samples);
        }
//        for (int x = 0; x < X; x++) {
//            cout << correlations[x] << ", ";
//        }
//        cout << endl;
    }
}

int IDLA::h(int x) {
    int red_interface = bottomEmptyY;
    while ((red_interface < Y - 1) and (graph[X * red_interface + x] == RED)) {
        red_interface++;
    }

    return red_interface;
}


void IDLA::printGraph() {
    for (int y = Y - 1; y >= 0; y--) {
        for (int x = 0; x < X; x++) {
            int c = getColor(x, y);
            std::cout << c;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}