package edu.mit.broad.marker.permutation;

import java.util.*;

public class UnbalancedRandomPermuter implements Permuter {
	
	cern.jet.random.engine.MersenneTwister random;
	long[] values;
	int size;
	int numClassOne;
   
	public UnbalancedRandomPermuter(int size, int numClassOne, int seed) {
		random = new cern.jet.random.engine.MersenneTwister(seed);
		this.size = size;
		this.numClassOne = numClassOne;
		values = new long[numClassOne];
	}

	public int[] next() {
		int[] assignments = new int[size];
		
		cern.jet.random.sampling.RandomSampler.sample(numClassOne, size, numClassOne, 0, values, 0, random); 
		
		for(int i = 0; i < values.length; i++) {
			int classOneIndex = (int) values[i];
			assignments[classOneIndex] = 1;
		}
		return assignments;
	}
	
	static String toString(int[] a) {
		String s = "";
		for(int i = 0; i < a.length; i++) {
			s+= a[i];	
		}
		return s;
	}
}
