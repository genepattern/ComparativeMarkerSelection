package edu.mit.broad.marker.permutation;
import java.util.*;

/**
 * @author     Joshua Gould
 * @created    October 4, 2004
 */
public class PermutationTest {



	public static void main(String[] args) {
		test2();
	}



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
		int[] classZeroIndices = {0, 1, 2, 3};
		int[] classOneIndices = {4, 5, 6, 7};
		BalancedRandomPermuter permuter = new BalancedRandomPermuter(classZeroIndices, classOneIndices);
		System.out.println(doTest(permuter, 1000));
		// {00111010=24, 01101001=35, 11001010=34, 10101001=21, 10100101=20, 01010110=30, 01011001=30, 11001100=25, 00110101=34, 01011100=18, 10010011=29, 01010101=35, 10011001=22, 11000011=24, 10010101=19, 10011100=21, 10011010=30, 00110110=34, 10100011=28, 01011010=28, 11000110=30, 01101100=32, 01100011=26, 00111100=24, 11000101=30, 11001001=30, 00111001=21, 01010011=26, 01100101=38, 10101100=31, 00110011=33, 10010110=25, 01100110=25, 10100110=25, 10101010=31, 01101010=32}

	}


	static void test3() {
		int[] classZeroIndices = {0, 1, 2, 3};
		int[] classOneIndices = {4, 5, 6, 7};
		BalancedCompletePermuter permuter = new BalancedCompletePermuter(classZeroIndices, classOneIndices);
		int total = permuter.getTotal().intValue();
		System.out.println(doTest(permuter, total));
		// 4!/2!2!*4!/2!2!
		//	00111010=1, 01101001=1, 11001010=1, 10100101=1, 10101001=1, 01010110=1, 01011001=1, 11001100=1, 00110101=1, 01011100=1, 10010011=1, 01010101=1, 11000011=1, 10011001=1, 10010101=1, 10011010=1, 10011100=1, 00110110=1, 10100011=1, 01011010=1, 11000110=1, 01101100=1, 01100011=1, 00111100=1, 11000101=1, 11001001=1, 01010011=1, 00111001=1, 01100101=1, 10101100=1, 00110011=1, 10010110=1, 01100110=1, 10100110=1, 01101010=1, 10101010=1}

	}


	static void test4() {
		UnbalancedCompletePermuter permuter = new UnbalancedCompletePermuter(6, 3);
		int total = permuter.getTotal().intValue();
		System.out.println(doTest(permuter, total));
		// {010101=1, 000111=1, 010011=1, 101010=1, 101100=1, 100101=1, 110010=1, 011010=1, 001011=1, 111000=1, 011100=1, 010110=1, 100110=1, 100011=1, 110001=1, 001101=1, 110100=1, 001110=1, 011001=1, 101001=1}
	}
}
