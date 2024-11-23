#ifndef UTILS_H

#include <list>
#include <utility>

const double eps = 1e-6;
const double M = 1e6;
const double inf = 1e10;

inline bool isl(double a, double b) {
	return a < b - eps;
}

inline bool iseq(double a, double b) {
	return std::abs(a - b) < eps;
}

struct Node {
	std::list<std::pair<int, int>> sep, tog;
};

#endif UTILS_H