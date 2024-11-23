#include <ilcplex/ilocplex.h>
#include <cstdlib>
#include <utility>
#include <vector>
#include <set>
#include "Utils.h"

using namespace std;

void updateBounds(const Node& node,
				  int lambda_counter,
				  const vector<set<int>>& items, 
				  IloNumVarArray& lambda, 
				  bool newUB);
IloRangeArray computeConstraints(const Node& node,
								 int n,
								 const vector<set<int>>& items,
								 IloEnv& env, 
							     IloNumVarArray& x);
std::pair<int, int> computeBranchingItems(int n,
										  const vector<set<int>>& items,
										  const IloCplex& rmp,
										  const IloNumVarArray& lambda);

int main() {
	vector<int> weight = {2, 1, 3, 3, 5};
	vector<set<int>> items;
	int capacity = 7;
	int n = weight.size();

	IloEnv env;
	IloModel master_model(env);

	IloNumVarArray lambda(env, n, 0, IloInfinity);

	IloExpr sum_obj(env);
	IloRangeArray partition_constraint(env);

	for (int i = 0; i < n; i++) {
		char var_name[50];
		sprintf(var_name, "y%d", i);

		lambda[i].setName(var_name);
		sum_obj += M * lambda[i];

		partition_constraint.add(lambda[i] == 1);
		
		items.push_back({i});
	}

	master_model.add(partition_constraint);

	IloObjective master_objective = IloMinimize(env, sum_obj);
	master_model.add(master_objective);

	IloCplex rmp(master_model);

	rmp.setOut(env.getNullStream()); // disables CPLEX log

	// rmp.solve();

	// cout << "Initial lower bound: " << rmp.getObjValue() << endl;

	// cout << "Initial solution: ";
	// for (size_t j = 0; j < lambda.getSize(); j++) {
	// 	cout << rmp.getValue(lambda[j]) << " ";
	// }
	// cout << endl;

	list<Node> tree;
	tree.push_back(Node());
	double best_lower_bound = inf;

	int lambda_counter = n;
	// Begin branch and price
	int it = 0;
	while (!tree.empty()) {
		cout << "It:" << it++ << endl;
		// DFS:
		auto node = tree.back();
		auto itNode = tree.end();
		--itNode;
		cout << "New branch" << endl;

		// TODO: merge do seguinte algoritmo com o laço da linha 69 
		// Setar os limites de acordo com sep e tog
		// Os limites devem ser ajustados em ambos os problemas
		// -Na relaxação: seta-se os bounds
		// -No pricing: adiciona-se uma restrição
		// Resolve a GC novamente, com essa nova configuração
		// Armazena novamente as colunas adicionadas
		// Calcula o z
		// Pega a mais fracionária e adiciona novas regras de sep e tog
		// Adiciona o nó na árvore
		// Setar os limites de volta aos originais

		// z é recalculado em cada iteração
		// o vector com itens vinculados a cada lambda é o mesmo para todos os nós
		// o problema mestre é o mesmo, mas os pricings são calculados em cada 
		// iteração

		updateBounds(node, lambda_counter, items, lambda, 0.0);

		rmp.solve();

		cout << "Initial lower bound: " << rmp.getObjValue() << endl;

		cout << "Initial solution: ";
		for (size_t j = 0; j < lambda.getSize(); j++) {
			cout << rmp.getValue(lambda[j]) << " ";
		}
		cout << endl;
		for (size_t j = 0; j < lambda.getSize(); j++) {
			if (iseq(rmp.getValue(lambda[j]), 1.0)) {
				cout << "Bin " << j << ":" << endl;
				for (auto i : items[j]) {
					cout << i << " ";
				}
				cout << endl;
			}
		}
		getchar();
		if (isl(rmp.getObjValue(), best_lower_bound)) {
			best_lower_bound = rmp.getObjValue();
		} else {
			updateBounds(node, lambda_counter, items, lambda, 1.0);
			// Update tree

			tree.erase(itNode);
			cout << "tree size:" << tree.size() << endl;
			continue;
		}


		cout << "Start column generation" << endl;
		// Begin column generation 
		bool isColAdded = false;
		while (true) {
			// Get the dual variables
			IloNumArray pi(env, n);

			rmp.getDuals(pi, partition_constraint);

			for (size_t i = 0; i < n; i++) {
				cout << "Dual variable of constraint " << i << " = " << pi[i] << endl;
			}

			// Build and solve the pricing problem
			IloModel pricing_model(env);

			IloNumVarArray x(env, n, 0, 1, IloNumVar::Bool);
			IloExpr sum_w_packed(env);
			IloExpr sum_obj_pricing(env, 1);

			for (int i = 0; i < n; i++) {
				char var_name[50];
				sprintf(var_name, "x%d", i);

				lambda[i].setName(var_name);
				sum_obj_pricing -= pi[i] * x[i];

				sum_w_packed += weight[i] * x[i];
			}
			auto constrs = computeConstraints(node, n, items, env, x);
			pricing_model.add(constrs);
			pricing_model.add(IloMinimize(env, sum_obj_pricing));
			pricing_model.add(sum_w_packed <= capacity);

			IloCplex pricing_problem(pricing_model);

			pricing_problem.setOut(env.getNullStream()); // disables CPLEX log

			pricing_problem.solve();

			cout << "Pricing obj:" << pricing_problem.getObjValue() << endl;
			if (isl(pricing_problem.getObjValue(), 0.0)) {
				isColAdded = true;
				
				cout << "Reduced cost is equal to " << pricing_problem.getObjValue() << ", which is less than 0..." << endl;

				IloNumArray entering_col(env, n);

				pricing_problem.getValues(x, entering_col);

				set<int> col;
				cout << endl << "Entering column:" << endl;
				for (size_t i = 0; i < n; i++) {
					if (isl(entering_col[i], 0.5)) {
						cout << 0 << endl;
					} else {
						cout << 1 << endl;
						col.insert(i);
					}
				}
				cout << endl;
				items.push_back(col);

				// Add the column to the master problem
				// (the cost of the new variable is always 1)
				char var_name[50];
				sprintf(var_name, "y%d", lambda_counter++);
				// Essa linha é bruxaria pura; adiciona a variável na fo com custo 1 e já adiciona a restrição de partição
				IloNumVar new_lambda(master_objective(1) + partition_constraint(entering_col), 0, IloInfinity);
				new_lambda.setName(var_name);

				lambda.add(new_lambda);

				cout << "Solving the RMP again..." << endl;

				rmp.solve();
			} else {
				cout << "No column with negative reduced costs found. The current basis is optimal" << endl;
				cout << "Final master problem: " << endl;
				system("cat model.lp");
				break;
			}
		} // End column generation
		if (isColAdded) {
			// Check if the master problem has been changed since the last it
			auto [i, j] = computeBranchingItems(n, items, rmp, lambda);
			if (i >= 0) {
				Node n1(node), n2(node);
				n1.sep.push_back(std::make_pair(i, j));
				n2.tog.push_back(std::make_pair(i, j));
				tree.push_back(n1);
				tree.push_back(n2);
			}
		}
		updateBounds(node, lambda_counter, items, lambda, 1.0);
		// Update tree

		cout << "tree size:" << tree.size() << endl;

		tree.erase(itNode);
	} // End branch and price

	for (int i = 0; i < items.size(); ++i) {
		cout << "lambda:" << i << " items:";
		for (auto v : items[i]) {
			cout << v << " ";
		}
		cout << endl;
	}

	cout << endl;
	cout << "Forcing items 1 and 2 to be separated in the master (for branch-and-price only): " << endl;
	// 0 1 2 3 4 5 6 7 8 9 10 11
	//                           
	// 1 0 0 0 0 1 1 1 0 1  0  0
	// 0 1 0 0 0 1 1 0 0 0  1  1
	// 0 0 1 0 0 1 0 1 1 0  0  1
	// 0 0 0 1 0 0 1 0 1 0  0  1
	// 0 0 0 0 1 0 0 0 0 1  1  0

	// itens 1 and 2 are together only on columns 5 and 11
	lambda[11].setUB(0.0);
	lambda[5].setUB(0.0);

	// to allow them again:
	// lambda[5].setUB(IloInfinity);
	// lambda[11].setUB(IloInfinity);

	rmp.solve();

	for (size_t j = 0; j < lambda.getSize(); j++)
	{
		cout << rmp.getValue(lambda[j]) << " ";
	}
	cout << endl;

	env.end();

	return 0;
}

/**
 * @brief Update the bounds of the lambda variables in the master problem.
 */
void updateBounds(const Node& node,
				  int lambda_counter,
				  const vector<set<int>>& items, 
				  IloNumVarArray& lambda, 
				  bool newUB) {
	// Set the bounds for the master problem
	for (int i = 0; i < lambda_counter; ++i) {
		// Force s.first and s.second to be in separate bins
		for (auto s : node.sep) {
			if (items[i].contains(s.first) && items[i].contains(s.second)) {
				lambda[i].setUB(newUB);
			}
		}
		// Force t.first and t.second to be together in the same bin
		for (auto t : node.tog) {
			if (items[i].contains(t.first) != items[i].contains(t.second)) {
				lambda[i].setUB(newUB);
			}
		}
	}
}

/**
 * @brief Insert constraints to force the branchin in the pricing problem.
 */
IloRangeArray computeConstraints(const Node& node,
								 int n,
								 const vector<set<int>>& items,
								 IloEnv& env, 
							     IloNumVarArray& x) {
	IloRangeArray constrs(env);
	// Set the bounds for the master problem
	for (int i = 0; i < n; ++i) {
		// Force s.first and s.second to be in separate bins
		for (auto s : node.sep) {
			if (items[i].contains(s.first) && items[i].contains(s.second)) {
				constrs.add(x[i] <= 0);
			}
		}
		// Force t.first and t.second to be together in the same bin
		for (auto t : node.tog) {
			if (items[i].contains(t.first) != items[i].contains(t.second)) {
				constrs.add(x[i] <= 0);
			}
		}
	}

	return constrs;
}

/**
 * @brief Compute the items over which the branching is performed.
 * 
 * The first element of the returned pair is equal to -1 if branching is not
 * required.
 */
std::pair<int, int> computeBranchingItems(int n,
										  const vector<set<int>>& items,
										  const IloCplex& rmp,
										  const IloNumVarArray& lambda) {
	vector<vector<double>> z(n, vector<double>(n, 0.0));
	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			for (int k = 0; k < items.size(); ++k) {
				if (items[k].contains(i) && items[k].contains(j)) {
					z[i][j] += rmp.getValue(lambda[k]);
				}
			}
		}
	}

	double m = 1;
	int k = -1, l;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (abs(z[i][j] - 0.5) < m) {
				m = abs(z[i][j] - 0.5);
				k = i;
				l = j;
			}
		}
	}

	return std::make_pair(k, l);
}