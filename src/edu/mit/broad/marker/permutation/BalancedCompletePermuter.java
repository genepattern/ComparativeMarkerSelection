package edu.mit.broad.marker.permutation;

import java.util.*;

public class BalancedCompletePermuter implements Permuter {
	Random rand;
	int[] assignments;
	int permutations;
	
	public BalancedCompletePermuter(int i, int j, int[] assignments) {
		this.assignments = assignments;
		rand = new Random(3);
		permutations = (int) getNumBalancedPermutations();
	}

	public int[] next() {
		throw new RuntimeException();
		
	}
	
	public int getTotal() {
		return permutations;	
	}
	
	private long getNumBalancedPermutations() {
		return 0;
		/*int numClassZero = classVector.getIndices(0).length;

		return cern.jet.math.Arithmetic.longFactorial(classVector.size()) / cern.jet.math.Arithmetic.longFactorial(
				numClassZero / 2);
				*/
	}
	
  
}
