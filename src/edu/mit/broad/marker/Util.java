package edu.mit.broad.marker;

/**
 *@author    Joshua Gould
 */
public class Util {
	public static void print(double[] a) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i] + " ");
		}
		System.out.println();
	}
	
	static public void print(int[] a) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i] + "\t");
		}
		System.out.println();
	}
	
	static public void print(long[] a) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i] + "\t");
		}
		System.out.println();
	}


	public static void print(double[][] a) {
		for(int i = 0; i < a.length; i++) {
			for(int j = 0; j < a[0].length; j++) {
				System.out.print(a[i][j] + "\t");
			}
			System.out.println();
		}
	}


	public static void print(double[][] a, int col) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i][col] + "\t");
		}
		System.out.println();
		System.out.println();
	}
}
