#ifndef BIN_PACKING_PROBLEM_H
#define BIN_PACKING_PROBLEM_H

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <memory>
#include <limits>
#include <sstream>
#include <unordered_set>
#include <utility>
#include <vector>
#include "combo.h"
#include "Data.h"
#include "Utils.h"

class BPP {
public:
    BPP(const std::shared_ptr<Data> pData, IloEnv& rEnv);

    /**
     * @brief Solve the restricted master problem.
     */
    std::pair<int, int> solve(const Node& rNode);

    /**
     * @brief Update the bounds of the lambda variables in the master problem.
     */
    void updateBounds(const Node& rNode);
    /**
     * @brief Insert constraints to force the branchin in the pricing problem.
     */
    void addPricingConstrs(const Node& rNode,
                           IloModel& rPricingModel,
						   IloBoolVarArray& x);

    /**
     * @brief Compute the items over which the branching is performed.
     * 
     * The first element of the returned pair is equal to -1 if branching is not
     * required.
     */
    std::pair<int, int> computeBranchingItems();

    IloNumArray getDuals(IloEnv& rEnv) const;

    double getBestIntObjValue() const { return mBestIntObj; }
    double getRmpObjValue() const { return mMaster.getObjValue(); }

    int getNbLambda() const { return mLambdas.getSize(); }

    void insertColumn(IloNumArray& rCol);

    Combo runCombo(double lowerBound, 
				   double upperBound, 
				   const IloNumArray& pi);

    void printSol() const;
    void printBins() const;

private:
    std::shared_ptr<Data> mpData;
    std::vector<std::unordered_set<int>> mItems;

	IloModel mMasterModel;
	IloNumVarArray mLambdas;
	IloExpr mObj;
	IloRangeArray mPartitionConstr;
	IloObjective mMasterObj;
	IloCplex mMaster;

    double mBestIntObj;
};

#endif