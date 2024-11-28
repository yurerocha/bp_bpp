#include "BPP.h"

BPP::BPP(const std::shared_ptr<Data> pData, IloEnv& rEnv)
    : mpData(pData), mMasterModel(rEnv), 
      mLambdas(rEnv, mpData->getNbItems(), 0, IloInfinity),
      mObj(rEnv), mPartitionConstr(rEnv), 
	  mBestIntObj(std::numeric_limits<int>::max()) {
	for (int i = 0; i < mpData->getNbItems(); i++) {
		std::ostringstream oss;
        oss << "y" << i;

		mLambdas[i].setName(oss.str().c_str());
		mObj += M * mLambdas[i];

		mPartitionConstr.add(mLambdas[i] == 1);
		
		mItems.push_back({i});
	}

	mMasterModel.add(mPartitionConstr);

	mMasterObj = IloMinimize(rEnv, mObj);
	mMasterModel.add(mMasterObj);

	mMaster = IloCplex(mMasterModel);

    // Disables CPLEX log
	mMaster.setOut(rEnv.getNullStream());
}

std::pair<int, int> BPP::solve(const Node& rNode) {
	updateBounds(rNode);

	// Solve restricted master problem
	mMaster.solve();

	IloEnv env;
	// cout << "Initial lower bound: " << pBPP->getRmpObjValue() << endl;

	// cout << "Initial solution: ";
	// pBPP->printSol();
	// pBPP->printBins();

	// std::cout << "Bins used: "  << pBPP->getBestIntObjValue() << std::endl;
	// double bestLowerBound = inf;

	// if (isl(pBPP->getRmpObjValue(), bestLowerBound)) {
	// 	bestLowerBound = pBPP->getRmpObjValue();
	// } else {
	// 	// Update tree
	// 	tree.erase(itNode);
	// 	// cout << "tree size:" << tree.size() << endl;
	// 	continue;
	// }

	// cout << "Start column generation" << endl;
	// Begin column generation 
	bool isColAdded = false;
	while (true) {
		// Get the dual variables
		auto pi = getDuals(env);

		// for (size_t i = 0; i < mpData->getNbItems(); i++) {
		// 	std::cout << "Dual variable of constraint " 
		// 			  << i << " = " << pi[i] << std::endl;
		// }

		// Build and solve the pricing problem
		IloModel pricingModel(env);

		IloNumVarArray x(env, mpData->getNbItems(), 0, 1, IloNumVar::Bool);
		IloExpr sumPacked(env);
		IloExpr pricingObj(env, 1);

		for (int i = 0; i < mpData->getNbItems(); i++) {
			std::ostringstream oss;
			oss << "x" << i;

			x[i].setName(oss.str().c_str());

			pricingObj -= pi[i] * x[i];

			sumPacked += mpData->getItemWeight(i) * x[i];
		}
		addPricingConstrs(rNode, env, pricingModel, x);
		pricingModel.add(IloMinimize(env, pricingObj));
		pricingModel.add(sumPacked <= mpData->getBinCapacity());

		IloCplex pricingProblem(pricingModel);

		// Disables CPLEX log
		pricingProblem.setOut(env.getNullStream());

		pricingProblem.solve();

		// std::cout << "Pricing obj:" << pricingProblem.getObjValue() << std::endl;
		if (isl(pricingProblem.getObjValue(), 0.0)) {
			isColAdded = true;
			// std::cout << "Reduced cost is equal to " 
			// 		  << pricingProblem.getObjValue() 
			// 		  << ", which is less than 0..." << std::endl;

			IloNumArray enteringCol(env, mpData->getNbItems());

			pricingProblem.getValues(enteringCol, x);

			insertColumn(enteringCol);

			// std::cout << "Solving the RMP again..." << std::endl;

			mMaster.solve();

			// enteringCol.clear();
            // enteringCol.end();
            // pricingProblem.clear();
            // pricingProblem.end();
		} else {
			// std::cout << "No column with negative reduced costs found. "
			// 		  << "The current basis is optimal" << std::endl;
			// std::cout << "Final master problem: " << std::endl;
			// system("cat model.lp");
			// pricingProblem.clear();
            // pricingProblem.end();
			break;
		}
	} // End column generation

	std::pair<int, int> b = {-1, -1};
	if (isColAdded) {
		// Check if the master problem has been changed since the last it
		b = computeBranchingItems();
	}

	bool isSolInt = true;
	for (int j = 0; j < mLambdas.getSize(); ++j) {
		double l = mMaster.getValue(mLambdas[j]);
        if (!(iseq(l, 0.0) || iseq(l, 1.0))) {
			isSolInt = false;
			break;
        }
    }

	if (isSolInt && isl(mMaster.getObjValue(), mBestIntObj)) {
		mBestIntObj = mMaster.getObjValue();
	}

	return b;
}

void BPP::updateBounds(const Node& rNode) {
	double ub = 0.0;
	// Set the bounds for the master problem
	for (int i = 0; i < mLambdas.getSize(); ++i) {
		// Restore the bounds for the following iterations
		mLambdas[i].setUB(1.0);
		// Force s.first and s.second to be in separate bins
		for (auto s : rNode.sep) {
			if (mItems[i].contains(s.first) && mItems[i].contains(s.second)) {
				mLambdas[i].setUB(ub);
			}
		}
		// Force t.first and t.second to be together in the same bin
		for (auto t : rNode.tog) {
			if (mItems[i].contains(t.first) != mItems[i].contains(t.second)) {
				mLambdas[i].setUB(ub);
			}
		}
	}
}

void BPP::addPricingConstrs(const Node& rNode, 
                                IloEnv& rEnv,
							    IloModel& rPricingModel,
                                IloNumVarArray& x) {
	// Set the bounds for the master problem
	for (int i = 0; i < mpData->getNbItems(); ++i) {
		// Force s.first and s.second to be in separate bins
		for (auto s : rNode.sep) {
			rPricingModel.add(x[s.first] + x[s.second] <= 1);
		}
		// Force t.first and t.second to be together in the same bin
		for (auto t : rNode.tog) {
			rPricingModel.add(x[t.first] == x[t.second]);
		}
	}
}

std::pair<int, int> BPP::computeBranchingItems() const {
	std::vector<std::vector<double>> z(mpData->getNbItems(), 
                                std::vector<double>(mpData->getNbItems(), 0.0));
	double m = std::numeric_limits<double>::max();
	std::pair<int, int> b(-1, 0);
	for (int i = 0; i < mpData->getNbItems(); ++i) {
		for (int j = i + 1; j < mpData->getNbItems(); ++j) {
			for (int k = mpData->getNbItems(); k < mItems.size(); ++k) {
				if (mItems[k].contains(i) && mItems[k].contains(j)) {
					z[i][j] += mMaster.getValue(mLambdas[k]);
				}
			}

			if (isl(std::abs(z[i][j] - 0.5), m)) {
				m = abs(z[i][j] - 0.5);
				// std::cout << "m:" << m << std::endl;
				b = {i, j};
			}
		}
	}

	return b;
}

IloNumArray BPP::getDuals(IloEnv& rEnv) const {
    IloNumArray pi(rEnv, mpData->getNbItems());
    mMaster.getDuals(pi, mPartitionConstr);

    return pi;
}

void BPP::insertColumn(IloNumArray& rCol) {
    std::set<int> newItemsCol;
    // std::cout << std::endl << "Entering column:" << std::endl;
    for (int i = 0; i < mpData->getNbItems(); ++i) {
        if (isg(rCol[i], 0.5)) {
            // std::cout << 1 << std::endl;
			// rCol[i] = 1.0;
            newItemsCol.insert(i);
        } 
		// else {
        //     std::cout << 0 << std::endl;
		// 	rCol[i] = 0.0;
        // }
    }
    mItems.push_back(newItemsCol);

    // Add the column to the master problem
    // (the cost of the new variable is always 1)
    std::ostringstream oss;
    oss << "y" << getNbLambda() + 1;

    // Essa linha é bruxaria pura; adiciona a variável na fo com custo 1 e já 
    // adiciona a restrição de partição
    IloNumVar newLambda(mMasterObj(1) + mPartitionConstr(rCol), 0, IloInfinity);
    newLambda.setName(oss.str().c_str());

    mLambdas.add(newLambda);
}

void BPP::printSol() const {
    for (int j = 0; j < mLambdas.getSize(); ++j) {
        std::cout << mMaster.getValue(mLambdas[j]) << " ";
    }
    std::cout << std::endl;
}

void BPP::printBins() const {
    for (int j = 0; j < mLambdas.getSize(); ++j) {
        if (!iseq(mMaster.getValue(mLambdas[j]), 0.0)) {
            std::cout << "Bin " << j << ":" << std::endl;
            for (auto v : mItems[j]) {
                std::cout << v << " ";
            }
            std::cout << std::endl;
        }
    }
}