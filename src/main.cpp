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

	Timer timer;

	list<Node> tree;
	Node n;
	n.sep = {};
	n.tog = {};
	tree.push_back(n);

	// Begin branch and price
	int it = 0;
	timer.start();
	while (!tree.empty()) {
		// cout << "It:" << it++ << endl;
		// DFS:
		auto& node = tree.back();
		auto itNode = tree.end();
		--itNode;
		// cout << "New branch" << endl;

		auto b = pBPP->solve(node);
		if (b.first >= 0) {
			Node ns, nj;
			ns = node; 
			nj = node;
			ns.sep.push_back(b);
			nj.tog.push_back(b);
			tree.push_back(nj);
			tree.push_back(ns);
		}
		// Update tree
		tree.erase(itNode);
		// cout << "tree size:" << tree.size() << endl;
	} // End branch and price

	env.end();

	std::cout << argv[1] << std::endl;
	std::cout << "Bins used: "  << pBPP->getBestIntObjValue() << std::endl;
	std::cout << "Elapsed time (s):" << timer.count() << std::endl;

	return 0;
}