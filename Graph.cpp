//
// Created by Matthew Nicoletti on 2020-04-17.
//


#include "Graph.h"

#include <chrono>
#include <iostream>
#include <fstream>




using namespace std;

const double pi = std::acos(-1);

// x is expected to be within 0 and X-1
// y is expect to be within 0 and Y-1

// assumes X is power of 2
// assumes Y is power of 2

CylinderGraph::CylinderGraph(int X, int Y, int k, int T0) {
    CylinderGraph::X = X;
    CylinderGraph::Y = Y;
    CylinderGraph::k = k;
    CylinderGraph::topBlueY = Y - 1;
    CylinderGraph::bottomRedY = 0;
    CylinderGraph::T0 = T0;
    CylinderGraph::dists = vector<vector<double>>(X, vector<double>(X, 0));

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

//   start with limit shape

    CylinderGraph::topBlueY = Y/2 - 1;
    CylinderGraph::bottomRedY = Y/2;
    for(int y = 0; y < Y; y++) {
        for (int x = 0; x < X; x++) {
            if (y >= Y/2) {
                setColor(x, y, RED);
            }
        }
    }

//    // start with sine wave
//    CylinderGraph::topBlueY = Y/2 - 1;
//    CylinderGraph::bottomRedY = Y/2;
//    for(int y = 0; y < Y; y++) {
//        for (int x = 0; x < X; x++) {
//            if (y >= Y/2 + 3 * cos(2 * pi * k * x / X)) {
//                setColor(x, y, RED);
//            }
//        }
//    }




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
    int i;
    int M;

    vector<double> d;

    char next_neighbor;
    double rand_neighbor;
    double rand01;
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
            if (startColor == RED and y >= topBlueY + 1 and y - Y/2 < Y / 4) {
                M = 3 * Y / 2 - 2 * y - 1;
                d = dists[M];
                rand01 = dist(mt);
                i = 0;
                double sum = 0;
                while (rand01 > sum) {
                    sum = sum + d[i];
                    i++;
                }
                x = (x + i) % X;
            } else {
                y++;
            }
        }
        if (next_neighbor == 'D') {
            if (startColor == BLUE and y <= bottomRedY-1 and Y/2-y <= Y / 4) {
                M = y;
                d = dists[M];
                rand01 = dist(mt);
                i = 0;
                double sum = 0;
                while (rand01 > sum) {
                    sum = sum + d[i];
                    i++;
                }
                x = (x + i) % X;

            }
            else {
                y--;
            }

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
    const std::complex<double> I(0, 1);
    double N = (double) X;


    double q = N *  acosh(2 - cos(2 * pi * k / X));

    const double Ablue = exp(q/2)/(exp(q/2) + exp(((-1 - N/2) * q)/N));
    const double Bblue = 1 - exp(q/2)/(exp(q/2) + exp(((-1 - N/2) * q)/N));

    const double Ared = 1/(1 + exp(q/2 + ((-1 + N/2) * q)/N));
    const double Bred = 1 - 1/(1 + exp(q/2 + ((-1 + N/2) * q)/N));

    const double lambda = ((-1 + exp(2 * pi * k)) * 2 * pi * k)/(1 + exp(2 * pi * k));

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(0, X - 1);

    std::random_device rdRW;
    std::mt19937 mtRW(rd());
    std::uniform_real_distribution<double> distRW(0, 1);

    int x;
    std::vector<int> bdest, rdest;  // TODO make point
    int err;


    //current value of M
//    complex<double> Mt = 0;
    complex<double> mtemp;

    //The martingale weighted by e^t -- should converge to e^t x_t, x_t o.u. process
//    complex<double> Yt = 0;
    complex<double> ytemp = 0;

    //current value of Dphi
    complex<double> Dphi = 0;
    complex<double> difftemp;

//    // not starting at limit shape
//    for(int y = 0; y < Y; y++) {
//        for (int x = 0; x < X; x++) {
//            if (getColor(x, y) == BLUE) {
////                Mt += cos((2 * pi * x * k) / N) * exp(q * (y - Y / 2))/ N ;
//                Mt += cos((2 * pi * x * k) / N) * cosh(q * (y - Y / 2))/ N ;
//                if (y >= Y/2) {
//                    Dphi += cos((2 * pi * x * k) / N) / N;
//                }
//            } else if (getColor(x, y) == RED and y < Y /2) {
//                Dphi -= cos((2 * pi * x * k) / N) / N;
//            }
//        }
//    }




    //quadratic variation of M, Dphi respectively
//    double Q = 0;
//    double phiQ = 0;
//    double YQt = 0;
    double diff4 = 0;

    //time
    int T = 0;


    for (int sample = 0; sample < num_samples; sample++) {

//      this is just temporary -- only for many samples of same time
        // re-initialize graph
        initializeGraph();
//
//        // initialize martingale
//        complex<double> Mt = 0;
//        // not starting at limit shape
//        for(int y = 0; y < Y; y++) {
//            for (int x = 0; x < X; x++) {
//                if (getColor(x, y) == BLUE) {
//                    Mt += cos((2 * pi * x * k) / N) * exp(q * (y - Y / 2))/ N ;
//                }
//            }
//        }
//        cout << Mt << endl;

//
        complex<double> Mt = 0;
        double Q = 0;
        complex<double> Dphi = 0;
        double phiQ = 0;
        complex<double> Yt = 0;
        double YQt = 0;

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




//          Y(t+1) - Y(t) = \frac{-1}{1 + e^{q_k}} e^{2 \pi i k X / N}(-e^{q_k} \exp(q_k Y / N) - \exp(-q_k Y / N))
//                 - \frac{-1}{1 + e^{-q_k}} e^{2 \pi i k X / N}(-e^{-q_k} \exp(q_k Y / N) - \exp(-q_k Y / N))
//          Fix boundary condition so it isn't for shifted half integer lattice



            ytemp = Yt;

            T = i;



            Yt = Yt +
                 (    cos((2 * pi * bdest[0] * k) / N) * (Ablue * exp(q * (bdest[1] - Y / 2) / N)
                                                          + Bblue * exp(-q * (bdest[1] - Y / 2) / N) )
                      - cos(2 * pi * rdest[0] * k / N) *
                        (Ared * exp(q * (rdest[1] - Y / 2) / N) + Bred * exp(-q * (rdest[1] - Y / 2) / N ))  )
                 * exp(2 * lambda * (T - T0) / (N*N)) / N;

            //quad var of Y
            YQt += abs(Yt - ytemp) * abs(Yt - ytemp);


            mtemp = Mt * exp(2 * lambda * (T-1 - T0) / (N*N));
//
//
            Mt = Mt + ( cos((2 * pi * bdest[0] * k) / N) -
                       cos((2 * pi * rdest[0] * k) / N) ) / N;



//
            if (T == 4 * 256) {
                M.emplace_back( Mt );
            }
//
            Q += abs(Mt * exp(2 * lambda * (T - T0) / (N*N)) - mtemp) * abs(Mt * exp(2 * lambda * (T - T0) / (N*N)) - mtemp);

//            // temporary just for testing
//            diff4 += abs(Mt - mtemp) * abs(Mt - mtemp) * abs(Mt - mtemp) * abs(Mt - mtemp);
//
//
            difftemp = Dphi;

            // Update the value of Dphi, add new value to list

            Dphi = Dphi + cos((2 * pi * bdest[0] * k) / N) / N -
                   cos((2 * pi * rdest[0] * k) / N)  / N;


//            // update quad variation for Dphi
            phiQ += abs(Dphi - difftemp) * abs(Dphi - difftemp);






            setTopBlueY();
            setBottomRedY();
        }


        Dphilist.push_back(Dphi);
        DphiQ.emplace_back( phiQ );

        Ylist.emplace_back( Yt );
        YQ.emplace_back(YQt);

//        M.emplace_back( Mt );
        MQ.emplace_back( Q );
//

//
//
//

//
//         temporarily using the following two lists
//          not for trajectory, but for
//         independent samples of value at fixed time




        std::cout<<"done."<<std::endl;

        // printDensityGraph(32);

//        fluctuations.push_back(h(0));
//
//        for (int x = 0; x < X; x++) {
//            correlations[x] += (h(x) - Y / 2.0) * (h(0) - Y / 2.0) / num_samples;
//        }
    }
    cout << diff4 << endl;
}

int CylinderGraph::h(int x) {
    int red_interface = bottomRedY;
    while ((red_interface < Y - 1) and (graph[X * red_interface + x] == BLUE)) {
        red_interface++;
    }

    return red_interface;
}

double CylinderGraph::computeReturnProb(int M, int x) {
    const double pi = std::acos(-1);
    double N = (double) X;
    double Z = 1.0/4 * (2 * M + 1) * N;
    double S = 0;
    for (int k = 0; k < N; k++) {
        for (int m = 0; m < M; m++) {
            S += 1/(4 - 2 * cos(2 * pi * k / N) - 2 * cos(pi * (m + 1.0/2)/(M + 1.0/2))) *
                    cos(2 * pi * k * x / N) *
                            (- cos(2 * pi * (m + 1.0/2)/(M + 1.0/2)) + 1)/2;
        }
    }
    S = S / Z;
    return S;
}

void CylinderGraph::computeReturnDists() {
    for (int M = X / 4; M < Y / 2; M++) {
        for (int x = 0; x < X; x++) {
            dists[M][x] = computeReturnProb(M, x);
        }
    }
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


