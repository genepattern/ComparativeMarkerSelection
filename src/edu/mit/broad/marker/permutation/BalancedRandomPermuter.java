package edu.mit.broad.marker.permutation;

import java.util.*;

public class BalancedRandomPermuter implements Permuter {

	long[] values;
	int[] classZeroIndices, classOneIndices;
	cern.jet.random.engine.MersenneTwister random;


	public BalancedRandomPermuter(int[] classZeroIndices, int[] classOneIndices) {

		this.classZeroIndices = classZeroIndices;
		this.classOneIndices = classOneIndices;
		values = new long[classZeroIndices.length / 2];
		random = new cern.jet.random.engine.MersenneTwister();
	}


	public int[] next() {
		int[] newAssignments = new int[classZeroIndices.length * 2];

		cern.jet.random.sampling.RandomSampler.sample(values.length, classZeroIndices.length, values.length, 0, values, 0, random);// choose half the indices from class 0

		// gives indices of classZeroIndices that will be class 1
		for(int i = 0; i < values.length; i++) {
			int classOneIndex = (int) values[i];
			newAssignments[classZeroIndices[classOneIndex]] = 1;
		}

		cern.jet.random.sampling.RandomSampler.sample(values.length, classZeroIndices.length, values.length, 0, values, 0, random);// choose half the indices from class 1

		for(int i = 0; i < values.length; i++) {
			int classOneIndex = (int) values[i];
			newAssignments[classOneIndices[classOneIndex]] = 1;
		}
		return newAssignments;
	}


	public static void main(String[] args) {
		int[] classZeroIndices = {0, 1, 2, 3};
		int[] classOneIndices = {4, 5, 6, 7};
		BalancedRandomPermuter p = new BalancedRandomPermuter(classZeroIndices, classOneIndices);
		for(int i = 0; i < 10; i++) {
			edu.mit.broad.marker.Util.print(p.next());
		}

	}

}

