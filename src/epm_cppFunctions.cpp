// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

// return index of x in integer vec
// [[Rcpp::export(name = c_which_int, rng = false)]]
int c_which_int(std::vector<int> vec, int x) {
	int nx = vec.size();
	for (int i = 0; i < nx; i++) {
		if (vec[i] == x) {
			return i;
		}
	}
	return -1;
}

// return index of x in character vec
// [[Rcpp::export(name = c_which_char, rng = false)]]
int c_which_char(std::vector<std::string> vec, std::string x) {
	int nx = vec.size();
	for (int i = 0; i < nx; i++) {
		if (vec[i] == x) {
			return i;
		}
	}
	return -1;
}

// return intersect of two integer vectors
// [[Rcpp::export(name = intersect_int, rng = false)]]
std::vector<int> intersect_int(std::vector<int> vec1, std::vector<int> vec2) {

	// intersect(vec1, vec2)
	std::vector<int> out;
	for (int k = 0; k < vec1.size(); k++) { 
		if (std::find(vec2.begin(), vec2.end(), vec1[k]) != vec2.end()) {
			out.push_back(vec1[k]);
		}
	}

	return out;
}

// setdiff for integers
// [[Rcpp::export(name = setdiff_int, rng = false)]]
std::vector<int> setdiff_int(std::vector<int> vec1, std::vector<int> vec2) {

	// setdiff(vec1, vec2)
	std::vector<int> out;

	sort(vec1.begin(), vec1.end());
	sort(vec2.begin(), vec2.end());
	set_difference(
		vec1.begin(),
	    vec1.end(),
	    vec2.begin(),
	    vec2.end(),
	    back_inserter(out));
	
	return out;
}



// return intersect of two species vectors
// [[Rcpp::export(name = getComponentA, rng = false)]]
std::vector<std::string> getComponentA(std::vector<std::string> commI, std::vector<std::string> commJ) {

	// intersect(cellI, cellJ)
	std::vector<std::string> a;
	// for each species in cellI, is it present in cellJ?
	for (int k = 0; k < commI.size(); k++) { 
		if (std::find(commJ.begin(), commJ.end(), commI[k]) != commJ.end()) {
			a.push_back(commI[k]);
		}
	}

	return a;
}


// return species in commJ but not commI
// [[Rcpp::export(name = getComponentB, rng = false)]]
std::vector<std::string> getComponentB(std::vector<std::string> commI, std::vector<std::string> commJ) {

	// setdiff(cellI, cellJ)
	std::vector<std::string> diffIJ;

	sort(commI.begin(), commI.end());
	sort(commJ.begin(), commJ.end());
	set_difference(
		commI.begin(),
	    commI.end(),
	    commJ.begin(),
	    commJ.end(),
	    back_inserter(diffIJ));
	
	return diffIJ;
}


// return species in commI but not commJ
// [[Rcpp::export(name = getComponentC, rng = false)]]
std::vector<std::string> getComponentC(std::vector<std::string> commI, std::vector<std::string> commJ) {

	// setdiff(cellJ, cellI)
	std::vector<std::string> diffJI;

	sort(commI.begin(), commI.end());
	sort(commJ.begin(), commJ.end());
	set_difference(
		commJ.begin(),
	    commJ.end(),
	    commI.begin(),
	    commI.end(),
	    back_inserter(diffJI));
	
	return diffJI;
}


// return edge indices
// [[Rcpp::export(name = getRootToTipEdges, rng = false)]]
List getRootToTipEdges(List phylo) {

	// extract components from phylo list
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
	NumericVector edge1a = edge(_, 0);
	NumericVector edge2a = edge(_, 1);
  
	std::vector<int> edge1 = as< std::vector<int> >(edge1a);
	std::vector<int> edge2 = as< std::vector<int> >(edge2a);
	int rootnode = tipLabels.size() + 1;

	List out = tipLabels.size();
	
	for (int i = 0; i < tipLabels.size(); i++) {
		std::vector<int> nodes;
		int childnode = i + 1;
		while (childnode != rootnode) {
			int parentnode = edge1[c_which_int(edge2, childnode)];
			nodes.push_back(childnode);
			childnode = parentnode;
		}

		std::vector<int> edgesInd(nodes.size());
		for (int j = 0; j < nodes.size(); j++) {
			edgesInd[j] = c_which_int(edge2, nodes[j]);
		}

		out[i] = edgesInd;
	}


	return out;
}



// from vector of values, create pairwise matrix and take the mean of the by-row minimums.
// [[Rcpp::export(name = meanNNdist, rng = false)]]
double meanNNdist(NumericVector input) {

	int n = input.size();

	if (n == 1) {
		return 0;
	} else {
	
		// get all pairwise distances
		NumericVector minVals(n);
		for (int i = 0; i < n; i++) {
			NumericVector vec(n, NumericVector::get_na());
			for (int j = 0; j < n; j++) {
				if (i != j) {
					vec[j] = std::abs(input[i] - input[j]);
				}
			}
			minVals[i] = min(na_omit(vec));
		}

		return mean(na_omit(minVals));
	}
}

// from vector of values, create pairwise matrix and take the minimum of the by-row minimums.
// [[Rcpp::export(name = minNNdist, rng = false)]]
double minNNdist(NumericVector input) {

	int n = input.size();

	if (n == 1) {
		return 0;
	} else {
		
		// get all pairwise distances
		NumericVector minVals(n);
		for (int i = 0; i < n; i++) {
			NumericVector vec(n, NumericVector::get_na());
			for (int j = 0; j < n; j++) {
				if (i != j) {
					vec[j] = std::abs(input[i] - input[j]);
				}
			}
			minVals[i] = min(na_omit(vec));
		}

		return min(na_omit(minVals));
	}
}

// from vector of values, create pairwise matrix and take the variance of the by-row minimums.
// [[Rcpp::export(name = minNNdist, rng = false)]]
double varNNdist(NumericVector input) {

	int n = input.size();

	if (n == 1) {
		return 0;
	} else {
		
		// get all pairwise distances
		NumericVector minVals(n);
		for (int i = 0; i < n; i++) {
			NumericVector vec(n, NumericVector::get_na());
			for (int j = 0; j < n; j++) {
				if (i != j) {
					vec[j] = std::abs(input[i] - input[j]);
				}
			}
			minVals[i] = min(na_omit(vec));
		}

		return var(na_omit(minVals));
	}
}


// return either mean or median value of trait per cell
// [[Rcpp::export(name = cellAvg, rng = false)]]
NumericVector cellAvg(List input, NumericVector trait, String stat) {
		
	int n = input.size();
	NumericVector out(n);
	
	for (int i = 0; i < n; i++) {
		
		CharacterVector ind = as< CharacterVector >(input[i]);
		if (ind[0] != "NA") {
			NumericVector vals = trait[ind];
		
			if (stat == "mean") {
				out[i] = double(mean(vals));
			} else if (stat == "median") {
				out[i] = double(median(vals));				
			} else if (stat == "variance") {
				out[i] = double(var(vals));
			} else if (stat == "mean_NN_dist") {
				out[i] = double(meanNNdist(vals));
			} else if (stat == "min_NN_dist") {
				out[i] = double(minNNdist(vals));
			} else if (stat == "range") {
				out[i] = double(max(vals) - min(vals));
			} else if (stat == "evenness") {
				out[i] = double(varNNdist(vals));
			}
		} else {
			out[i] = NA_REAL;
		}
	}
	
	return out;	
}


// given a list of character vectors, and a character vector, 
// return the list of vectors, intersected with the second vector
// [[Rcpp::export(name = intersectList, rng = false)]]
List intersectList(List input, StringVector vec) {

	int n = input.size();
	List out(n);

	for (int i = 0; i < n; i++) {

		Rcpp::checkUserInterrupt();		
		StringVector sp = as< StringVector > (input[i]);

		if (all(!is_na(sp))) {

			StringVector res = intersect(sp, vec);
	
			if (res.size() > 0) {
				out[i] = res;
			} else {
				out[i] = NA_REAL;
			}
		} else {
			out[i] = NA_REAL;
		}
	}

	return out;

}


// for each species in vec, count how many cells it is found in
// [[Rcpp::export(name = countCells, rng = false)]]
NumericVector countCells(List cellList, StringVector vec) {
	
	std::vector<std::string> uniqueSp = as< std::vector<std::string> >(vec);
	std::vector<int> out(uniqueSp.size());
	
	for (int i = 0; i < cellList.size(); i++) {
		
		std::vector<std::string> cell = as< std::vector<std::string> >(cellList[i]);
		
		for (int j = 0; j < uniqueSp.size(); j++) {
		
			if (std::find(cell.begin(), cell.end(), uniqueSp[j]) != cell.end()) {
				out[j] = out[j] + 1;
			}
		}
	}
	
	return wrap(out);
}


// Calculate taxonomic turnover for all pairwise communities
// This will be 1 - sorensen
// [[Rcpp::export(name = calcPairwiseTaxonomicSorensen, rng = false)]]
NumericMatrix calcPairwiseTaxonomicSorensen(List allComm, String component) {

	int nComm = allComm.size();
	NumericMatrix out(nComm, nComm);

	for (int i = 0; i < nComm; i++) {
		for (int j = 0; j < nComm; j++) {
			if (i <= j) {

				Rcpp::checkUserInterrupt();
				std::vector<std::string> commI = as< std::vector<std::string> >(allComm[i]);
				std::vector<std::string> commJ = as< std::vector<std::string> >(allComm[j]);

				// if a community is empty, there is no distance
				if ((commI[0] == "empty") | (commJ[0] == "empty")) {
					out(i,j) = -1.0;
					out(j,i) = -1.0;
				} else {

					std::vector<std::string> a = getComponentA(commI, commJ);

					// if the intersect is the same as the sizes of both I and J, then
					// the communities are identical, so 0 turnover
					// if (commI.size() == commJ.size() && commI.size() == a.size()) {
					// 	out(i,j) = 0.0;
					// 	out(j,i) = 0.0;

					// // if intersect is of length 0, then there is complete turnover, so 1.0
					// } else if (a.size() == 0) {
					// 	out(i,j) = 1.0;
					// 	out(j,i) = 1.0;
					// } else {

						std::vector<std::string> b = getComponentB(commI, commJ);
						std::vector<std::string> c = getComponentC(commI, commJ);

						double nA = double(a.size());
						double nB = double(b.size());
						double nC = double(c.size());

						// Simpson's beta diversity index
						// cellVec[j] = 1.0 - double(a.size()) / (double(a.size()) + std::min(double(b.size()), double(c.size())));
						
						// Sorenson metric
						if (component == "turnover") {
							// calculate beta SIM
							out(i,j) = std::min(nB, nC) / (nA + std::min(nB, nC));
							
						} else if (component == "nestedness") {
							// calculate beta SNE
							out(i,j) = ((std::max(nB, nC) - std::min(nB, nC)) / (2 * nA + nB + nC)) * (nA / (nA + std::min(nB, nC)));


						} else if (component == "full") {
							// calculate beta SOR
							out(i,j) = (nB + nC) / (2 * nA + nB + nC);
						}

						out(j,i) = out(i,j);
					
						// Jaccard metric
						//cellVec[j] = 1.0 - double(a.size()) / (double(a.size()) + double(b.size() + double(c.size())));
					//}
				}
			}
		}
	}
	return out;
}


// given:
// a list of species names: a,
// the full vector of species names: tipLabels,
// a list of edge indices that trace the path from the root to each tip: spEdges,
// find the union of branch indices for all species in a.
// [[Rcpp::export(name = uniqueBranchesForSet, rng = false)]]
std::vector<int> uniqueBranchesForSet(std::vector<std::string> a, std::vector<std::string> tipLabels, List spEdges) {

	//std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);

	std::vector<int> out;

	for (int i = 0; i < a.size(); i++) {
		
		int spInd = c_which_char(tipLabels, a[i]);
		std::vector<int> branchInd = as< std::vector<int> >(spEdges[spInd]);

		for (int j = 0; j < branchInd.size(); j++) {
			if (std::find(out.begin(), out.end(), branchInd[j]) == out.end()) {
				out.push_back(branchInd[j]);
			}
		}
	}

	return out;

}
	


// Calculate phylogenetic beta diversity for each pair of communities in list
// [[Rcpp::export(name = calcPairwisePhylosor2, rng = false)]]
NumericMatrix calcPairwisePhylosor2(List allComm, List phylo, String component) {

	// extract relevant info from input tree for function FaithPD_branchIndices
	// 		tip labels, branch lengths, and tipward nodes vector
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	std::vector<double> edgeLengths = as< std::vector<double> >(phylo["edge.length"]);

	// get leaves for all nodes as well as tip edges
//	List nodeLeaves = getLeavesForNodes(phylo);
	List spEdges = getRootToTipEdges(phylo);


	int nComm = allComm.size();
	NumericMatrix out(nComm, nComm);

	// Rcout << "Starting pairwise matrix..." << std::endl;
	Progress p(nComm, true);

	for (int i = 0; i < nComm; i++) {
		// Rcout << "row " << i << std::endl;
		p.increment(); 

		for (int j = 0; j < nComm; j++) {
			if (i <= j) {

				// Rcout << "\tcol " << j << std::endl;

				Rcpp::checkUserInterrupt();
				std::vector<std::string> commI = as< std::vector<std::string> >(allComm[i]);
				std::vector<std::string> commJ = as< std::vector<std::string> >(allComm[j]);

				if ((commI[0] == "empty") | (commJ[0] == "empty")) {
					out(i,j) = -1.0;
					out(j,i) = -1.0;
				} else {

					// get intersection of commI and J
					// std::vector<std::string> a = getComponentA(commI, commJ);

					double a_pd = 0;
					double b_pd = 0;
					double c_pd = 0;

					// get branches that contribute to commI PD
					std::vector<int> edgesCommI = uniqueBranchesForSet(commI, tipLabels, spEdges);

					// get branches that contribute to commJ PD
					std::vector<int> edgesCommJ = uniqueBranchesForSet(commJ, tipLabels, spEdges);

					// identify branches that are part of both commI and commJ
					std::vector<int> edgesCommIJ = intersect_int(edgesCommI, edgesCommJ);

					// identify branches that are part of commI PD, but not commJ PD
					std::vector<int> edgesInotJ = setdiff_int(edgesCommI, edgesCommJ);

					// identify branches that are part of commJ PD, but not commI PD
					std::vector<int> edgesJnotI = setdiff_int(edgesCommJ, edgesCommI);

					// sum the branch lengths for a
					for (int k = 0; k < edgesCommIJ.size(); k++) {
						a_pd = a_pd + edgeLengths[edgesCommIJ[k]];
					}

					// sum the branch lengths for b
					for (int k = 0; k < edgesInotJ.size(); k++) {
						b_pd = b_pd + edgeLengths[edgesInotJ[k]];
					}

					// sum the branch lengths for c
					for (int k = 0; k < edgesJnotI.size(); k++) {
						c_pd = c_pd + edgeLengths[edgesJnotI[k]];
					}

					// Rcout << "a_pd " << a_pd << std::endl;
					// Rcout << "b_pd " << b_pd << std::endl;
					// Rcout << "c_pd " << c_pd << std::endl;

					if (component == "turnover") {
						// calculate beta SIM
						out(i,j) = std::min(b_pd, c_pd) / (a_pd + std::min(b_pd, c_pd));

					} else if (component == "nestedness") {
						// calculate beta SNE
						out(i,j) = ((std::max(b_pd, c_pd) - std::min(b_pd, c_pd)) / (2 * a_pd + b_pd + c_pd)) * (a_pd / (a_pd + std::min(b_pd, c_pd)));

					} else if (component == "full") {
						// calculate beta SOR
						out(i,j) = (b_pd + c_pd) / (2 * a_pd + b_pd + c_pd);
					}
					
					out(j,i) = out(i,j);
				}
			}
		}
	}
	return out;
}




