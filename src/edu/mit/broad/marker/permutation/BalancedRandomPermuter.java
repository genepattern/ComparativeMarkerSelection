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

import java.util.*;

public class BalancedRandomPermuter implements Permuter {

	long[] values;
	int[] classZeroIndices, classOneIndices;
	cern.jet.random.engine.MersenneTwister random;


	public BalancedRandomPermuter(int[] classZeroIndices, int[] classOneIndices, int seed) {

		this.classZeroIndices = classZeroIndices;
		this.classOneIndices = classOneIndices;
		values = new long[classZeroIndices.length / 2];
		random = new cern.jet.random.engine.MersenneTwister(seed);
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

}

