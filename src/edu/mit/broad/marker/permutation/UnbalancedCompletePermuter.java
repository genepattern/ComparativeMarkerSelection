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
		int[] classOneIndices = combinationGenerator.getNext();
		int[] values = new int[size];
		for(int i = 0; i < numClassZero; i++) {
			values[classOneIndices[i]] = 1;
		}
		return values;
	}
	
}
