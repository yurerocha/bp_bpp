#include <cstdlib>
#include <ilcplex/ilocplex.h>
#include <iomanip>
#include <memory>
#include <utility>
#include <vector>
#include <set>
#include <sstream>
#include "../include/BPP.h"
#include "../include/Data.h"
#include "../include/Utils.h"

using namespace std;

int main(int argc, char **argv) {
	if (argc != 2) {
		std::cout << "Usage:\n./bnp input/instance" << std::endl;
		return -1;
	}

	IloEnv env;

	auto pData = std::make_shared<Data>();
    pData->readData(argv[1]);

	auto pBPP = std::make_unique<BPP>(pData, env);

	Timer timer;

	list<Node> tree;
	Node n;
	n.sep = {};
	n.tog = {};
	n.isRoot = true;
	tree.push_back(n);

	int it = 0;
	timer.start();
	while (!tree.empty()) { // Begin branch and price
		// cout << "It:" << it++ << endl;
		// DFS:
		auto& rNode = tree.back();
		auto itNode = tree.end();
		--itNode;

		auto b = pBPP->solve(rNode);
		if (b.first >= 0) {
			Node ns, nj;
			ns = rNode; 
			nj = rNode;
			ns.isRoot = false;
			nj.isRoot = false;
			ns.sep.push_back(b);
			nj.tog.push_back(b);
			tree.push_back(ns);
			tree.push_back(nj);
		}
		// Update tree
		tree.erase(itNode);
		// cout << "tree size:" << tree.size() << endl;
	} // End branch and price

	env.end();

	// std::cout << argv[1] << std::endl;
	// std::cout << "Bins used: "  << pBPP->getBestObjValue() << std::endl;
	// std::cout << "Elapsed time (s):" << timer.count() << std::endl;

	std::cout << argv[1] << " " 
			  << pBPP->getBestIntObjValue() << " "
			  << std::fixed << std::setprecision(3) 
			  << timer.count() << std::endl;

	return 0;
}