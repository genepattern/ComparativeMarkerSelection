package edu.mit.broad.marker;
import java.io.*;
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


	public static void print(int[] a) {
		for(int i = 0; i < a.length; i++) {
			System.out.print(a[i] + "\t");
		}
		System.out.println();
	}


	public static void print(long[] a) {
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


	public static void print(double[][] a, String file) {
		PrintWriter temp = null;
		try {
			temp = new PrintWriter(new FileWriter(file));
			for(int i = 0; i < a.length; i++) {
				for(int j = 0; j < a[0].length; j++) {
					if(j > 0) {
						temp.print("\t");
					}
					temp.print(a[i][j]);
				}
				temp.println();
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		} finally {
			if(temp != null) {
				temp.close();
			}
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
