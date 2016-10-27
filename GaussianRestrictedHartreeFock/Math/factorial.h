#pragma once

inline double factorial(int n) {
    double result = 1.0;
    for (int m = 1; m < n+1; m++) {
        result *= m;
    }
    return result;
}
