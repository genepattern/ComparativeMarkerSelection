package edu.mit.broad.marker.permutation;
import java.util.*;

/**
 *@author    Joshua Gould
 */
public class PermutationTest {

	

	static String toString(int[] a) {
		StringBuffer buf = new StringBuffer(a.length); 
		for(int i = 0; i < a.length; i++) {
			buf.append(a[i]);
		}
		return buf.toString();
	}


	static Map doTest(Permuter permuter, int runs) {
		HashMap map = new HashMap();
		for(int i = 0; i < runs; i++) {
			String p = toString(permuter.next());
			Integer count = (Integer) map.get(p);
			if(count == null) {
				count = new Integer(1);
			} else {
				count = new Integer(count.intValue() + 1);
			}
			map.put(p, count);
		}
		return map;
	}


	static void test1() {
		UnbalancedRandomPermuter permuter = new UnbalancedRandomPermuter(6, 3);
		System.out.println(doTest(permuter, 1000));
		//6!/3!3!=20
		// {010101=45, 010011=55, 101010=41, 000111=57, 101100=46, 100101=43, 110010=55, 011010=56, 111000=44, 001011=45, 011100=50, 010110=65, 100011=57, 100110=46, 110001=39, 001101=53, 110100=47, 001110=47, 011001=62, 101001=47}
	}
	
	static void test2() {
		int[] classZeroIndices = {0,1,2};
		int[] classOneIndices = {3,4,5};
		BalancedRandomPermuter permuter = new BalancedRandomPermuter(classZeroIndices, classOneIndices);
		System.out.println(doTest(permuter, 1000));
		// 3!/1!2!
		// {001010=121, 001001=109, 001100=97, 100100=108, 010010=125, 100010=111, 010001=86, 010100=112, 100001=131}
	}
	
	
	static void test3() {
		int[] classZeroIndices = {0,1,2};
		int[] classOneIndices = {3,4,5};
		BalancedCompletePermuter permuter = new BalancedCompletePermuter(classZeroIndices, classOneIndices);
		int total = permuter.getTotal().intValue();
		System.out.println(doTest(permuter, total));
		// 3!/1!2!
		//{001001=1, 001010=1, 001100=1, 100100=1, 010010=1, 010001=1, 100010=1, 010100=1, 100001=1}
	}
	
	static void test4() {
		UnbalancedCompletePermuter permuter = new UnbalancedCompletePermuter(6, 3);
		int total = permuter.getTotal().intValue();
		System.out.println(doTest(permuter, total));
		// {010101=1, 000111=1, 010011=1, 101010=1, 101100=1, 100101=1, 110010=1, 011010=1, 001011=1, 111000=1, 011100=1, 010110=1, 100110=1, 100011=1, 110001=1, 001101=1, 110100=1, 001110=1, 011001=1, 101001=1}
	}


	


	public static void main(String[] args) {
		test4();
	}
}
