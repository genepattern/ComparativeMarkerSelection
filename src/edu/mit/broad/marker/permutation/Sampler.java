package edu.mit.broad.marker.permutation;
import java.util.Random;

import gnu.trove.*;
import java.math.BigInteger;

/**
 *@author    Joshua Gould
 */
public class Sampler {
	static Random rand = new Random(3);
	
	/*
	    public int median(double[] values) {
	    int n = values.length;
	    int half <- (n + 1)/2
	    if (n%2 == 1) {
	    return Arrays.sort(x, partial = half)[half];
	    }
	    else {
	    sum(Arrays.sort(x, partial = c(half, half + 1))[c(half, half +
	    1)])/2
	    }
	    }
	  */
	/**
	 *  Fills the values array with an ordered subset of numbers from the interval
	 *  from - to
	 *
	 *@param  from         Description of the Parameter
	 *@param  to           Description of the Parameter
	 *@param  numToChoose  Description of the Parameter
	 *@param  values       Description of the Parameter
	 */
	public static void choose(int from, int to, int numToChoose, long[] values) {
		cern.jet.random.sampling.RandomSampler.sample(numToChoose, to + 1, numToChoose, from, values, 0, null);
	}

	public double sum(double[] values) {
		double sum = 0;
		for(int i = 0, length = values.length; i < length; i++) {
			sum += values[i];	
		}
		return sum;
	}


	/**
	 *  Returns the p-th permutation of the sequence [0,1,...,N-1]
	 *
	 *@param  p  Description of the Parameter
	 *@param  N  Description of the Parameter
	 *@return    Description of the Return Value
	 */
	public static int[] permutation(long p, int N) {
		return cern.colt.GenericPermuting.permutation(p, N);
	}

	/**
     * Shuffle the elements of the list using the specified random
     * number generator.
     *
     * @param rand a <code>Random</code> value
     */
    public static int[] shuffle(int[] data) {
		 int[] copy = (int[]) data.clone();
        for (int i = copy.length; i-- > 1;) {
            swap(copy, i, rand.nextInt(i));
        }
		  return copy;
    }
	 
    /**
     * Swap the values at offsets <tt>i</tt> and <tt>j</tt>.
     *
     * @param i an offset into the data array
     * @param j an offset into the data array
     */
    private static final void swap(int[] _data, int i, int j) {
        int tmp = _data[i];
        _data[i] = _data[j];
        _data[j] = tmp;
    }

	/**
	 *  Returns a random permutation of the sequence [0,1,...,N-1]
	 *
	 *@param  N  Description of the Parameter
	 *@return    Description of the Return Value
	 */
	public static int[] randomPermutation(int N) {
		long p = 10l; //random.nextLong();
		return cern.colt.GenericPermuting.permutation(p, N);
	}


	static void print(int[] a) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i] + " ");
		}
		System.out.println();
	}


	

}
