package edu.mit.broad.marker;
import edu.mit.broad.dataobj.Dataset;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
/**
 *@author    Joshua Gould
 */
public class Util {
	public static int ASCENDING = 0;
	public static int DESCENDING = 1;
	public static int ABSOLUTE = 2;

	private static final double kMinVariancePercent = .20;

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


	public static double[] sort(final double[] values, int order) {
		if(order == 0) {
			cern.colt.GenericSorting.quickSort(0, values.length, new AscendingComparator(values), new DataSwapper(values));
		} else if(order == 1) {
			cern.colt.GenericSorting.quickSort(0, values.length, new DescendingComparator(values), new DataSwapper(values));
		} else if(order == 2) {
			cern.colt.GenericSorting.quickSort(0, values.length, new AbsoluteDescendingComparator(values), new DataSwapper(values));
		} else {
			throw new IllegalArgumentException("Invalid order");
		}

		return values;
	}


	/**
	 *@param  order   one of ASCENDING, DESCENDING, ABSOLUTE
	 *@param  values  Description of the Parameter
	 *@return         The indices
	 */
	public static int[] getIndices(final double[] values, int order) {
		final int[] indices = new int[values.length];
		for(int i = 0; i < indices.length; i++) {
			indices[i] = i;
		}
		if(order == 0) {
			cern.colt.GenericSorting.quickSort(0, values.length, new AscendingIndexComparator(values, indices), new IndexSwapper(indices));
		} else if(order == 1) {
			cern.colt.GenericSorting.quickSort(0, values.length, new DescendingIndexComparator(values, indices), new IndexSwapper(indices));
		} else if(order == 2) {
			cern.colt.GenericSorting.quickSort(0, values.length, new AbsoluteDescendingIndexComparator(values, indices), new IndexSwapper(indices));
		} else {
			throw new IllegalArgumentException("Invalid order");
		}

		return indices;
	}


	public static void main(String[] args) {
		double[] d = {7, 5, 10, 8};
		int[] indices = getIndices(d, 1);
		print(indices);
		print(rank(indices));
	}


	/**
	 *  Computes the rank using the given index array
	 *
	 *@param  indices  Description of the Parameter
	 *@return          Description of the Return Value
	 */
	public static int[] rank(int[] indices) {
		int[] rank = new int[indices.length];
		for(int j = 0; j < indices.length; j++) {
			rank[indices[j]] = j + 1;
		}
		return rank;
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


	public static double mean(Dataset dataset, int[] indices, int row) {

		double sum = 0;

		for(int j = 0; j < indices.length; j++) {
			sum += dataset.get(row, indices[j]);
		}

		return sum / indices.length;
	}


	public static double median(Dataset dataset, int[] indices, int row) {
		double[] data = new double[indices.length];

		for(int j = 0; j < indices.length; j++) {
			data[j] = dataset.get(row, indices[j]);
		}

		java.util.Arrays.sort(data);
		int half = indices.length / 2;
		return data[half];
	}


	public static void ttest(Dataset dataset,
			int[] classOneIndices,
			int[] class2Indices,
			double[] scores, boolean fixStdDev) {

		int rows = dataset.getRowDimension();

		for(int i = 0; i < rows; i++) {

			double class1Mean = Util.mean(dataset, classOneIndices, i);
			double class2Mean = Util.mean(dataset, class2Indices, i);
			double class1Std = Util.standardDeviation(dataset, classOneIndices,
					i, class1Mean);
			double class2Std = Util.standardDeviation(dataset, class2Indices, i,
					class2Mean);
					if(fixStdDev) {
						class1Std = fixStdDev(class1Std, class1Mean);
						class2Std = fixStdDev(class2Std, class2Mean);
					}
			double denom = Math.sqrt((class1Std * class1Std) + (class2Std * class2Std));
			double Sxi = (class1Mean - class2Mean) / denom;
			scores[i] = Sxi;
		}
	}


	public static double fixStdDev(double s, double m) {
		double absMean = Math.abs(m);
		double minS = kMinVariancePercent * absMean;
		if(minS < s) {
			minS = s;
		}

		if(minS == 0) {
			minS = .1;
		}
		return minS;
	}


	public static void snr(Dataset dataset,
			int[] classOneIndices,
			int[] class2Indices,
			double[] scores, boolean fixStdDev) {

		int rows = dataset.getRowDimension();

		for(int i = 0; i < rows; i++) {

			double class1Mean = Util.mean(dataset, classOneIndices, i);
			double class2Mean = Util.mean(dataset, class2Indices, i);
			double class1Std = Util.standardDeviation(dataset, classOneIndices,
					i, class1Mean);
			double class2Std = Util.standardDeviation(dataset, class2Indices, i,
					class2Mean);
			if(fixStdDev) {
						class1Std = fixStdDev(class1Std, class1Mean);
						class2Std = fixStdDev(class2Std, class2Mean);
					}
					
			double Sxi = (class1Mean - class2Mean) / (class1Std +
					class2Std);
			scores[i] = Sxi;
		}
	}


	public static double signal2Noise(Dataset dataset, int[] classOneIndices,
			int[] class2Indices, int row) {

		double mean1 = mean(dataset, classOneIndices, row);
		double mean2 = mean(dataset, class2Indices, row);
		double std1 = standardDeviation(dataset, classOneIndices, row, mean1);
		double std2 = standardDeviation(dataset, class2Indices, row, mean2);

		return (mean1 - mean2) / (std1 + std2);
	}


	public static double standardDeviation(Dataset dataset, int[] indices,
			int row, double mean) {

		double sum = 0;

		for(int j = 0; j < indices.length; j++) {

			double x = dataset.get(row, indices[j]);
			double diff = x - mean;
			diff = diff * diff;
			sum += diff;
		}

		double variance = sum / (indices.length - 1);

		return Math.sqrt(variance);
	}


	/**
	 *@author    Joshua Gould
	 */
	static class AscendingIndexComparator implements cern.colt.function.IntComparator {
		double[] values;
		int[] indices;


		public AscendingIndexComparator(double[] values, int[] indices) {
			this.values = values;
			this.indices = indices;
		}


		public int compare(int o1, int o2) {
			if(values[indices[o1]] < values[indices[o2]]) {
				return -1;
			} else if(values[indices[o1]] > values[indices[o2]]) {
				return 1;
			}
			return 0;
		}
	}


	/**
	 *@author    Joshua Gould
	 */
	static class DescendingIndexComparator implements cern.colt.function.IntComparator {
		double[] values;
		int[] indices;


		public DescendingIndexComparator(double[] values, int[] indices) {
			this.values = values;
			this.indices = indices;
		}


		public int compare(int o1, int o2) {
			if(values[indices[o1]] < values[indices[o2]]) {
				return 1;
			} else if(values[indices[o1]] > values[indices[o2]]) {
				return -1;
			}
			return 0;
		}
	}


	/**
	 *@author    Joshua Gould
	 */
	static class AbsoluteDescendingIndexComparator implements cern.colt.function.IntComparator {
		double[] values;
		int[] indices;


		public AbsoluteDescendingIndexComparator(double[] values, int[] indices) {
			this.values = values;
			this.indices = indices;
		}


		public int compare(int o1, int o2) {
			if(Math.abs(values[indices[o1]]) < Math.abs(values[indices[o2]])) {
				return 1;
			} else if(Math.abs(values[indices[o1]]) > Math.abs(values[indices[o2]])) {
				return -1;
			}
			return 0;
		}
	}


	/**
	 *@author    Joshua Gould
	 */
	static class IndexSwapper implements cern.colt.Swapper {
		int[] indices;


		public IndexSwapper(int[] indices) {
			this.indices = indices;
		}


		public void swap(int o1, int o2) {
			int ind1 = indices[o1];
			indices[o1] = indices[o2];
			indices[o2] = ind1;
		}
	}


	/**
	 *@author    Joshua Gould
	 */
	static class DataSwapper implements cern.colt.Swapper {
		double[] values;


		public DataSwapper(double[] values) {
			this.values = values;
		}


		public void swap(int o1, int o2) {
			double tmp = values[o1];
			values[o1] = values[o2];
			values[o2] = tmp;
		}
	}


	/**
	 *@author    Joshua Gould
	 */
	static class AscendingComparator implements cern.colt.function.IntComparator {
		double[] values;



		public AscendingComparator(double[] values) {
			this.values = values;

		}


		public int compare(int o1, int o2) {
			if(values[o1] < values[o2]) {
				return -1;
			} else if(values[o1] > values[o2]) {
				return 1;
			}
			return 0;
		}
	}


	/**
	 *@author    Joshua Gould
	 */
	static class DescendingComparator implements cern.colt.function.IntComparator {
		double[] values;


		public DescendingComparator(double[] values) {
			this.values = values;

		}


		public int compare(int o1, int o2) {
			if(values[o1] < values[o2]) {
				return 1;
			} else if(values[o1] > values[o2]) {
				return -1;
			}
			return 0;
		}
	}


	/**
	 *@author    Joshua Gould
	 */
	static class AbsoluteDescendingComparator implements cern.colt.function.IntComparator {
		double[] values;


		public AbsoluteDescendingComparator(double[] values) {
			this.values = values;

		}


		public int compare(int o1, int o2) {
			if(Math.abs(values[o1]) < Math.abs(values[o2])) {
				return 1;
			} else if(Math.abs(values[o1]) > Math.abs(values[o2])) {
				return -1;
			}
			return 0;
		}
	}

}
