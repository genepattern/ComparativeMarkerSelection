package edu.mit.broad.marker.permutation;



public class BalancedCompletePermuter implements Permuter {
	CombinationGenerator comb1;
	CombinationGenerator comb2;
	int[] comb1Indices;
	int[] classZeroIndices;
	int[] classOneIndices;
	
	public BalancedCompletePermuter(int[] classZeroIndices, int[] classOneIndices, int seed) {
		this.classZeroIndices = classZeroIndices;
		this.classOneIndices = classOneIndices;
		comb1 = new CombinationGenerator(classOneIndices.length, classOneIndices.length/2);
		comb2 = new CombinationGenerator(classOneIndices.length, classOneIndices.length/2);
		comb1Indices = comb1.getNext();
	}

	public int[] next() {
		int[] values = new int[classOneIndices.length*2];
		//edu.mit.broad.marker.Util.print(comb1Indices);
		for(int i = 0; i < comb1Indices.length; i++) {
			values[classOneIndices[comb1Indices[i]]] = 1;
		}
		
		int[] indices = comb2.getNext();
		for(int i = 0; i < indices.length; i++) {
			values[classZeroIndices[indices[i]]] = 1;
		}
		if(!comb2.hasMore()) {
			comb2.reset();
			if(comb1.hasMore()) {
				comb1Indices = comb1.getNext();
			}
		}
		return values;
	}
	
	
	public java.math.BigInteger getTotal() {
		return comb1.getTotal().multiply(comb2.getTotal());
	}
  
}
