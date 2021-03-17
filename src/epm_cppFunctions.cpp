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


// return a list of all nodes, containing tiplabels that are descendant from each
// [[Rcpp::export(name = getLeavesForNodes, rng = false)]]
List getLeavesForNodes(List phylo) {

	// extract components from phylo list
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
	NumericVector edge1a = edge(_, 0);
	NumericVector edge2a = edge(_, 1);
  
	std::vector<int> edge1 = as< std::vector<int> >(edge1a);
	std::vector<int> edge2 = as< std::vector<int> >(edge2a);
	int rootnode = tipLabels.size() + 1;
	int nodemax = *std::max_element(std::begin(edge2), std::end(edge2));
	std::vector<int> allnodes(nodemax);
	std::iota(allnodes.begin(), allnodes.end(), 1);

	// create list of tips, containing root-to-tip nodes
	List rootToTipNodes(tipLabels.size());
	for (int i = 0; i < tipLabels.size(); i++) {
		std::vector<int> nodes;
		int childnode = i + 1;
		while (childnode != rootnode) {
			int parentnode = edge1[c_which_int(edge2, childnode)];
			nodes.push_back(childnode);
			childnode = parentnode;
		}
		nodes.push_back(rootnode);
		rootToTipNodes[i] = nodes;
	}

	// create list of all nodes, and for each node, list the descendant taxa
	List nodeLeaves(allnodes.size());
	for (int i = 0; i < allnodes.size(); i++) {
		std::vector<std::string> foundLeaves;
		for (int j = 0; j < rootToTipNodes.size(); j++) {
			std::vector<int> tmp = as< std::vector<int> >(rootToTipNodes[j]);
			if (std::find(tmp.begin(), tmp.end(), allnodes[i]) != tmp.end()) {
				foundLeaves.push_back(tipLabels[j]);
			}
		}
		nodeLeaves[i] = foundLeaves;
	}


	return nodeLeaves;
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


// function to receive the nodeLeaves list and a set of tiplabels, to return the MRCA
// [[Rcpp::export(name = getMRCA_from_nodeLeaves, rng = false)]]
int getMRCA_from_nodeLeaves(List nodeLeaves, std::vector<std::string> taxa) {

//	std::vector<std::string> taxa = as< std::vector<std::string> >(tip);

	// for each node, does it contain all target taxa downstream?
	std::vector<int> commonNodes;
	for (int i = 0; i < nodeLeaves.size(); i++) {
		std::vector<std::string> node = as< std::vector<std::string> >(nodeLeaves[i]);
		int counter = 0;

		for (int j = 0; j < taxa.size(); j++) {
			if (std::find(node.begin(), node.end(), taxa[j]) != node.end()) {
				counter = counter + 1;
			}
		}

		// if node contained all taxa downstream, then counter should be equal to number of taxa
		if (counter == taxa.size()) {
			commonNodes.push_back(i + 1);
		}
	}

	// most recent common ancestor will be max node number
	// if just one tip, then it will be the minimum

	int mrca;

	if (taxa.size() > 1) {
	 	mrca = *std::max_element(std::begin(commonNodes), std::end(commonNodes));
	} else {
		mrca = *std::min_element(std::begin(commonNodes), std::end(commonNodes));
	}

	return mrca;
}



// [[Rcpp::export(name = FaithPD_branchIndices, rng = false)]]
std::vector<int> FaithPD_branchIndices(std::vector<std::string> a, List phylo, List nodeLeaves, List spEdges, bool includeRoot = false) {
	
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	std::vector<double> edgeLengths = as< std::vector<double> >(phylo["edge.length"]);

	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
	NumericVector edge2a = edge(_, 1);
  
	std::vector<int> edge2 = as< std::vector<int> >(edge2a);

	std::vector<int> branchIndices;

	// get mrca node of taxon set a
	if (a.size() > 1) {
		int mrca = getMRCA_from_nodeLeaves(nodeLeaves, a);

		// get mrca index in edge2 (= edge index)
		int mrcaInd = c_which_int(edge2, mrca);
		
		// get index of species labels and get branch indices from spEdges
		// only keep branch if branch index is greater than MRCA index
		for (int i = 0; i < a.size(); i++) {
			
			int tmpInd = c_which_char(tipLabels, a[i]);
			std::vector<int> branchInd = as< std::vector<int> >(spEdges[tmpInd]);
			
			for (int j = 0; j < branchInd.size(); j++) {
				if (std::find(branchIndices.begin(), branchIndices.end(), branchInd[j]) == branchIndices.end()) {
					if (branchInd[j] > mrcaInd) {
						branchIndices.push_back(branchInd[j]);
					}
				}
			}
		}

	} else {
		// if 1 species, return terminal branch length
		int spInd = c_which_char(tipLabels, a[0]);
		branchIndices.push_back(c_which_int(edge2, spInd + 1));
	}

	if (includeRoot) {
		// get the index of one of the species
		int spInd = c_which_char(tipLabels, a[0]);

		// get root-to-tip edge indices for this species
		std::vector<int> branchInd = as< std::vector<int> >(spEdges[spInd]);

		// add these indices to edge indices is not already present
		for (int i = 0; i < branchInd.size(); i++) {
			if (std::find(branchIndices.begin(), branchIndices.end(), branchInd[i]) == branchIndices.end()) {
				branchIndices.push_back(branchInd[i]);
			}
		}
	}

	return branchIndices;
}

// from vector of values, create pairwise matrix and take minimum
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

// from vector of values, create pairwise matrix and take minimum
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




// function to determine the number of geographic cells for every branch in phylogeny
// [[Rcpp::export(name = phyloBranchRanges, rng = false)]]
List phyloBranchRanges(List phylo, List speciesList, List tipEdges) {

	// extract components from phylo list
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	std::vector<double> edgeLengths = as< std::vector<double> >(phylo["edge.length"]);
	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
	NumericVector edge1a = edge(_, 0);
	NumericVector edge2a = edge(_, 1);
  
	std::vector<int> edge1 = as< std::vector<int> >(edge1a);
	std::vector<int> edge2 = as< std::vector<int> >(edge2a);

	std::vector<int> allnodes = edge2;

	// tipEdges are edge indices. Convert to child node number
	// tipNodes and tipEdges are in order of tiplabels
	List tipNodes(tipEdges.size());
	for (int i = 0; i < tipEdges.size(); i++) {
		std::vector<int> spTipEdges = as< std::vector<int> >(tipEdges[i]);
		for (int j = 0; j < spTipEdges.size(); j++) {
			spTipEdges[j] = edge2[spTipEdges[j]];
		}
		tipNodes[i] = spTipEdges;
	}

	Rcpp::checkUserInterrupt();

	// For each node number, find which species have it listed
	// speciesFromNode is in order of edge2 nodes
	List speciesFromNode(allnodes.size());
	for (int i = 0; i < allnodes.size(); i++) {
		Rcpp::checkUserInterrupt();
		std::vector<int> isPresent;
		if (allnodes[i] <= tipLabels.size()) {
			isPresent.push_back(allnodes[i]);
		} else {
			for (int j = 0; j < tipNodes.size(); j++) { 
				std::vector<int> nodesForSp = as< std::vector<int> >(tipNodes[j]);
				for (int k = 0; k < nodesForSp.size(); k++) { 
					if (nodesForSp[k] == allnodes[i]) {
						isPresent.push_back(j+1);
						break;
					}
				}
			}
		}
		speciesFromNode[i] = isPresent;
	}

//	convert tip indices to names
//	allLeaves is in order of edge2 nodes
	List allLeaves(speciesFromNode.size());
	for (int i = 0; i < speciesFromNode.size(); i++) {
		Rcpp::checkUserInterrupt();
		std::vector<int> tipIndices = as< std::vector<int> >(speciesFromNode[i]);
		std::vector<std::string> tipNames(tipIndices.size());
		for (int j = 0; j < tipIndices.size(); j++) { 
			tipNames[j] = tipLabels[tipIndices[j] - 1];
		}
		allLeaves[i] = tipNames;
	}

	// for each branch (= child node), count the number of cells that contain any of the
	// tip taxa
	std::vector<int> branchCellCount(allLeaves.size());
	for (int i = 0; i < allLeaves.size(); i++) {
		Rcpp::checkUserInterrupt();
		int counter = 0;
		std::vector<std::string> nodeSp = as< std::vector<std::string> >(allLeaves[i]);
		for (int j = 0; j < speciesList.size(); j++) {
			int counter2 = 0;
			std::vector<std::string> cellSp = as< std::vector<std::string> >(speciesList[j]);
			for (int k = 0; k < nodeSp.size(); k++) {
				if (std::find(cellSp.begin(), cellSp.end(), nodeSp[k]) != cellSp.end()) {
					counter2 = counter2 + 1;
				}
			}
			if (counter2 > 0) {
				counter = counter + 1;
			}
		}
		branchCellCount[i] = counter;
	}

	List out(2);
	out[0] = edgeLengths;
	out[1] = branchCellCount;

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
				if (commI[0] == "empty" | commJ[0] == "empty") {
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
// [[Rcpp::export(name = calcPairwisePhylosor, rng = false)]]
NumericMatrix calcPairwisePhylosor(List allComm, List phylo, String component) {

	// extract relevant info from input tree for function FaithPD_branchIndices
	// 		tip labels, branch lengths, and tipward nodes vector
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	std::vector<double> edgeLengths = as< std::vector<double> >(phylo["edge.length"]);
	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
	NumericVector edge2a = edge(_, 1);
	std::vector<int> edge2 = as< std::vector<int> >(edge2a);


	// get leaves for all nodes as well as tip edges
	List nodeLeaves = getLeavesForNodes(phylo);
	List spEdges = getRootToTipEdges(phylo);


	int nComm = allComm.size();
	NumericMatrix out(nComm, nComm);

	// Rcout << "Starting pairwise matrix..." << std::endl;
//	Progress p(nComm, true);

	for (int i = 0; i < nComm; i++) {
		// Rcout << "row " << i << std::endl;
//		p.increment(); 

		for (int j = 0; j < nComm; j++) {
			if (i <= j) {

				// Rcout << "\tcol " << j << std::endl;

				Rcpp::checkUserInterrupt();
				std::vector<std::string> commI = as< std::vector<std::string> >(allComm[i]);
				std::vector<std::string> commJ = as< std::vector<std::string> >(allComm[j]);

				if (commI[0] == "empty" | commJ[0] == "empty") {
					out(i,j) = -1.0;
					out(j,i) = -1.0;
				} else {

					// get intersection of commI and J
					// std::vector<std::string> a = getComponentA(commI, commJ);

					double a_pd = 0;
					double b_pd = 0;
					double c_pd = 0;

					// get branches that contribute to commI PD
					std::vector<int> edgesCommI = FaithPD_branchIndices(commI, phylo, nodeLeaves, spEdges, true);

					// get branches that contribute to commJ PD
					std::vector<int> edgesCommJ = FaithPD_branchIndices(commJ, phylo, nodeLeaves, spEdges, true);

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

				if (commI[0] == "empty" | commJ[0] == "empty") {
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




