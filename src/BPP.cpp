#include "BPP.h"

BPP::BPP(const std::shared_ptr<Data> pData, IloEnv& rEnv)
    : mpData(pData), mMasterModel(rEnv), 
      mLambdas(rEnv, mpData->getNbItems(), 0, IloInfinity),
      mObj(rEnv), mPartitionConstr(rEnv) {
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

	mRmp = IloCplex(mMasterModel);

    // Disables CPLEX log
	mRmp.setOut(rEnv.getNullStream());
}

bool BPP::solvePricingProblem(const Node& rNode, IloEnv& rEnv) {
	// Get the dual variables
	auto pi = getDuals(rEnv);

	for (size_t i = 0; i < mpData->getNbItems(); i++) {
		std::cout << "Dual variable of constraint " 
				  << i << " = " << pi[i] << std::endl;
	}

	// Build and solve the pricing problem
	IloModel pricingModel(rEnv);

	IloNumVarArray x(rEnv, mpData->getNbItems(), 0, 1, IloNumVar::Bool);
	IloExpr sum_w_packed(rEnv);
	IloExpr pricingObj(rEnv, 1);

	for (int i = 0; i < mpData->getNbItems(); i++) {
		std::ostringstream oss;
		oss << "x" << i;

		x[i].setName(oss.str().c_str());

		pricingObj -= pi[i] * x[i];

		sum_w_packed += mpData->getItemWeight(i) * x[i];
	}
	auto constrs = computeConstraints(rNode, rEnv, x);
	pricingModel.add(constrs);
	pricingModel.add(IloMinimize(rEnv, pricingObj));
	pricingModel.add(sum_w_packed <= mpData->getBinCapacity());

	IloCplex pricingProblem(pricingModel);

	// Disables CPLEX log
	pricingProblem.setOut(rEnv.getNullStream());

	pricingProblem.solve();

	std::cout << "Pricing obj:" << pricingProblem.getObjValue() << std::endl;
	if (isl(pricingProblem.getObjValue(), 0.0)) {
		std::cout << "Reduced cost is equal to " 
				  << pricingProblem.getObjValue() 
				  << ", which is less than 0..." << std::endl;

		IloNumArray enteringCol(rEnv, mpData->getNbItems());

		pricingProblem.getValues(x, enteringCol);

		insertColumn(enteringCol);

		std::cout << "Solving the RMP again..." << std::endl;

		solve();
	} else {
		std::cout << "No column with negative reduced costs found. "
				  << "The current basis is optimal" << std::endl;
		std::cout << "Final master problem: " << std::endl;
		// system("cat model.lp");
		return false;
	}

	return true;
}

void BPP::updateBounds(const Node& rNode, bool newUB) {
	// Set the bounds for the master problem
	for (int i = 0; i < mLambdas.getSize(); ++i) {
		// Force s.first and s.second to be in separate bins
		for (auto s : rNode.sep) {
			if (mItems[i].contains(s.first) && mItems[i].contains(s.second)) {
				mLambdas[i].setUB(newUB);
			}
		}
		// Force t.first and t.second to be together in the same bin
		for (auto t : rNode.tog) {
			if (mItems[i].contains(t.first) != mItems[i].contains(t.second)) {
				mLambdas[i].setUB(newUB);
			}
		}
	}
}

IloRangeArray BPP::computeConstraints(const Node& rNode, 
                                        IloEnv& rEnv, 
                                        IloNumVarArray& x) const {
	IloRangeArray constrs(rEnv);
	// Set the bounds for the master problem
	for (int i = 0; i < mpData->getNbItems(); ++i) {
		// Force s.first and s.second to be in separate bins
		for (auto s : rNode.sep) {
			if (mItems[i].contains(s.first) && mItems[i].contains(s.second)) {
				constrs.add(x[i] <= 0);
			}
		}
		// Force t.first and t.second to be together in the same bin
		for (auto t : rNode.tog) {
			if (mItems[i].contains(t.first) != mItems[i].contains(t.second)) {
				constrs.add(x[i] <= 0);
			}
		}
	}

	return constrs;
}

std::pair<int, int> BPP::computeBranchingItems() const {
	std::vector<std::vector<double>> z(mpData->getNbItems(), 
                                std::vector<double>(mpData->getNbItems(), 0.0));
	for (int i = 0; i < mpData->getNbItems(); ++i) {
		for (int j = i + 1; j < mpData->getNbItems(); ++j) {
			for (int k = 0; k < mItems.size(); ++k) {
				if (mItems[k].contains(i) && mItems[k].contains(j)) {
					z[i][j] += mRmp.getValue(mLambdas[k]);
				}
			}
		}
	}

	double m = 1;
	int k = -1, l;
	for (int i = 0; i < mpData->getNbItems(); i++) {
		for (int j = i + 1; j < mpData->getNbItems(); j++) {
			if (abs(z[i][j] - 0.5) < m) {
				m = abs(z[i][j] - 0.5);
				k = i;
				l = j;
			}
		}
	}

	return std::make_pair(k, l);
}

IloNumArray BPP::getDuals(IloEnv& rEnv) const {
    IloNumArray pi(rEnv, mpData->getNbItems());
    mRmp.getDuals(pi, mPartitionConstr);

    return pi;
}

void BPP::insertColumn(const IloNumArray& rCol) {
    std::set<int> newItemsCol;
    std::cout << std::endl << "Entering column:" << std::endl;
    for (int i = 0; i < mpData->getNbItems(); ++i) {
        if (isl(rCol[i], 0.5)) {
            std::cout << 0 << std::endl;
        } else {
            std::cout << 1 << std::endl;
            newItemsCol.insert(i);
        }
    }
    std::cout << std::endl;
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
        std::cout << mRmp.getValue(mLambdas[j]) << " ";
    }
    std::cout << std::endl;
}

void BPP::printBins() const {
    for (int j = 0; j < mLambdas.getSize(); ++j) {
        if (iseq(mRmp.getValue(mLambdas[j]), 1.0)) {
            std::cout << "Bin " << j << ":" << std::endl;
            for (auto v : mItems[j]) {
                std::cout << v << " ";
            }
            std::cout << std::endl;
        }
    }
}