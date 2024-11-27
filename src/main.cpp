#include <ilcplex/ilocplex.h>
#include <cstdlib>
#include <memory>
#include <utility>
#include <vector>
#include <set>
#include <sstream>
#include "BPP.h"
#include "Data.h"
#include "Utils.h"

using namespace std;

int main(int argc, char **argv) {
	if (argc != 2) {
		std::cout << "Usage:\n./bnp input/instance" << std::endl;
		return -1;
	}

	IloEnv env;

	auto pData = std::make_shared<Data>();
    pData->readData(argv[1]);

	// std::unique_ptr<Model> pBPP(pData, env);
	auto pBPP = std::make_unique<BPP>(pData, env);

	list<Node> tree;
	tree.push_back(Node());
	double bestLowerBound = inf;

	// Begin branch and price
	int it = 0;
	while (!tree.empty()) {
		cout << "It:" << it++ << endl;
		// DFS:
		auto node = tree.back();
		auto itNode = tree.end();
		--itNode;
		cout << "New branch" << endl;

		pBPP->updateBounds(node, 0.0);

		pBPP->solve();

		cout << "Initial lower bound: " << pBPP->getRmpObjValue() << endl;

		cout << "Initial solution: ";
		pBPP->printSol();
		pBPP->printBins();

		if (isl(pBPP->getRmpObjValue(), bestLowerBound)) {
			bestLowerBound = pBPP->getRmpObjValue();
		} else {
			pBPP->updateBounds(node, 1.0);
			// Update tree
			tree.erase(itNode);
			cout << "tree size:" << tree.size() << endl;
			continue;
		}

		cout << "Start column generation" << endl;
		// Begin column generation 
		bool isColAdded = false;
		while (true) {
			if (pBPP->solvePricingProblem(node, env)) {
				isColAdded = true;
			} else {
				break;
			}
		} // End column generation

		if (isColAdded) {
			// Check if the master problem has been changed since the last it
			auto [i, j] = pBPP->computeBranchingItems();
			if (i >= 0) {
				Node n1(node), n2(node);
				n1.sep.push_back(std::make_pair(i, j));
				n2.tog.push_back(std::make_pair(i, j));
				tree.push_back(n1);
				tree.push_back(n2);
			}
		}
		pBPP->updateBounds(node, 1.0);
		// Update tree

		cout << "tree size:" << tree.size() << endl;

		tree.erase(itNode);
	} // End branch and price

	std::cout << "Bins used: "  << bestLowerBound << std::endl;

	env.end();

	return 0;
}