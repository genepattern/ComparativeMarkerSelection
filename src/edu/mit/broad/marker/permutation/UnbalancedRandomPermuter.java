package edu.mit.broad.marker.permutation;

import java.util.*;

public class UnbalancedRandomPermuter implements Permuter {
	Random rand;
	int[] assignments;
	
	public UnbalancedRandomPermuter(int[] assignments) {
		this.assignments = assignments;
		rand = new Random(3);
	}

	public int[] next() {
		int N = assignments.length;
		return gc_permute(assignments);
		
	}
	
    int[] gc_permute(int[] _aArr) {
		int[] aArr = (int[]) _aArr.clone();
      int len=aArr.length;
      int r;
      int tmpInd;
      for(int i=len; i>1; i--) {
         r=(int) (rand.nextFloat() * i);
         tmpInd=aArr[r];
         aArr[r]=aArr[i - 1];
         aArr[i - 1]=tmpInd;
      }
      return aArr;
   }
}
