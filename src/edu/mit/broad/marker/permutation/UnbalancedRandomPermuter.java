package edu.mit.broad.marker.permutation;

import java.util.*;

public class UnbalancedRandomPermuter implements Permuter {
	
	cern.jet.random.engine.MersenneTwister random;
	long[] values;
	int size;
	int numClassOne;
	public UnbalancedRandomPermuter(int size, int numClassOne) {
		random = new cern.jet.random.engine.MersenneTwister();
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
	
	public static void main(String[] args) {
		UnbalancedRandomPermuter perm = new UnbalancedRandomPermuter(4, 1);
		HashMap map = new HashMap();
		for(int i = 0; i < 100; i++) {
			String p = toString(perm.next());
			Integer count = (Integer) map.get(p);
			if(count == null) {
				count = new Integer(1);	
			} else {
				count = new Integer(count.intValue()+1);
			}
			map.put(p, count);
		}
		System.out.println(map);
	}
	
  /*  int[] gc_permute(int[] _aArr) {
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
   }*/
}
