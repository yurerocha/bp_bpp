#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include <list>
#include <utility>
#include "combo.h"

const double eps = 1e-6;
const double M = 1e6;

inline bool isl(double a, double b) {
	return a < b - eps;
}

inline bool isg(double a, double b) {
	return a > b + eps;
}

inline bool iseq(double a, double b) {
	return std::abs(a - b) < eps;
}

inline bool isgeq(double a, double b) {
	return isg(a, b) || iseq(a, b);
}

struct Node {
	std::list<std::pair<int, int>> sep, tog;
	bool isRoot;
};

class Timer {
public:
	/**
	 * @brief Starts the timer.
	 */
	void start() { begin = std::chrono::steady_clock::now(); }

	/**
	 * @brief Returns the elapsed time since the start of the timer.
	 */
	double count() const {
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> et{end - begin};
		return et.count();
	}

private:
	std::chrono::steady_clock::time_point begin;
};

struct Combo {
	item *pItems;
	double maxProfit;

	// ~Combo() {
	// 	if (pItems) {
	// 		delete pItems;
	// 	}
	// }
};

#endif // Utils