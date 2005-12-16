/*
  The Broad Institute
  SOFTWARE COPYRIGHT NOTICE AGREEMENT
  This software and its documentation are copyright (2003-2006) by the
  Broad Institute/Massachusetts Institute of Technology. All rights are
  reserved.

  This software is supplied without any warranty or guaranteed support
  whatsoever. Neither the Broad Institute nor MIT can be responsible for its
  use, misuse, or functionality.
*/


package edu.mit.broad.marker.permutation;

/**
 *@author    Joshua Gould
 */
public class UnbalancedCompletePermuter implements  Permuter {
	CombinationGenerator combinationGenerator;
	int size;
	int numClassZero;
	
	public UnbalancedCompletePermuter(int size, int numClassZero) {
		this.size = size;
		this.numClassZero = numClassZero;
		combinationGenerator = new CombinationGenerator(size, size - numClassZero);
	}

	public java.math.BigInteger getTotal() {
		return combinationGenerator.getTotal();
	}


	public int[] next() {
		int[] classOneIndices = combinationGenerator.getNext(); // choose class one indices
		int[] values = new int[size];
		for(int i = 0; i < classOneIndices.length; i++) {
			values[classOneIndices[i]] = 1;
		}
		return values;
	}
	
}
