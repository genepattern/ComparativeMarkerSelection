package edu.mit.broad.marker;
import java.io.*;
import java.io.FileReader;

import java.util.Arrays;
import java.util.Random;
import java.util.StringTokenizer;

import edu.mit.broad.dataobj.*;
import edu.mit.broad.dataobj.microarray.*;
import edu.mit.broad.gp.*;
import edu.mit.broad.io.microarray.*;

import edu.mit.broad.marker.permutation.*;

/**
 * @author     Joshua Gould
 * @created    September 17, 2004
 */
public class MarkerSelection {

	/**  input expression values */
	Dataset dataset;

	/**  input classes */
	ClassVector classVector;
	/**  number of features to select */
	int numFeatures;
	/**  number of permutations to perform */
	int numPermutations;
	final static int CLASS_ZERO_GREATER_THAN_CLASS_ONE = 0;
	final static int CLASS_ZERO_LESS_THAN_CLASS_ONE = 1;
	final static int TWO_SIDED = 2;

	/**  whether to do one-sided >, one-sided <, or two-sided test */
	int side;

	int metric;
	final static int T_TEST = 0;
	final static int SNR = 1;

	/**  whether permutations are balanced */
	boolean balanced;

	/**  whether to compute all possible permutations */
	boolean complete;
	/**  whether to read permutations from a file */
	boolean readPermutations = false;
	BufferedReader testClassPermutationsReader;
	int[][] testAssignments;
	String outputFileName;


	public MarkerSelection(Dataset _dataset, ClassVector _classVector,
			int _numFeatures, int _numPermutations, int _side,
			String _outputFileName, boolean _balanced,
			boolean complete) {
		this.dataset = _dataset;
		this.classVector = _classVector;
		this.numFeatures = _numFeatures;
		this.numPermutations = _numPermutations;
		this.side = _side;
		this.balanced = _balanced;
		this.outputFileName = _outputFileName;
		this.complete = complete;

		String permutationsFile = System.getProperty("permutations");

		if(permutationsFile != null) {
			readPermutations = true;
			readPermutations(permutationsFile);
		}

		run();

		if(readPermutations) {

			if(testClassPermutationsReader != null) {

				try {
					testClassPermutationsReader.close();
				} catch(Exception e) {
				}
			}
		}
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

		Arrays.sort(data);
		int half = indices.length / 2;
		return data[half];
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


	public static void main(String[] args)
			 throws Exception {

		String datasetFile = args[0];
		String clsFile = args[1];
		int _numFeatures = Integer.parseInt(args[2]);
		int _numPermutations = Integer.parseInt(args[3]);
		int side = Integer.parseInt(args[4]);
		String outputFileName = args[5];
		boolean balanced = Boolean.valueOf(args[6]).booleanValue();
		boolean complete = Boolean.valueOf(args[7]).booleanValue();
		ExpressionDataReader reader = GPUtil.getExpressionReader(datasetFile);
		ExpressionDataImpl e = (ExpressionDataImpl) reader.read(datasetFile);
		new MarkerSelection(e.getExpressionMatrix(),
				GPUtil.getClassVector(clsFile), _numFeatures,
				_numPermutations, side, outputFileName, balanced,
				complete);
	}


	int[] invert(ClassVector cv) {
		int[] assignments = new int[cv.size()];
		int[] classZeroIndices = cv.getIndices(0);

		for(int i = 0, length = classZeroIndices.length; i < length; i++) {
			assignments[classZeroIndices[i]] = 1;
		}
		return assignments;
	}


	void run() {
		int[] classZeroIndices = null;
		int[] classOneIndices = null;

		if(side == CLASS_ZERO_LESS_THAN_CLASS_ONE) {
			classZeroIndices = classVector.getIndices(1);
			classOneIndices = classVector.getIndices(0);
		} else {
			classZeroIndices = classVector.getIndices(0);
			classOneIndices = classVector.getIndices(1);
		}

		if(balanced) {

			if(classZeroIndices.length != classOneIndices.length) {
				GPUtil.exit(
						"The number of items in each class must be equal for balanced permutations.");
			}

			if((classZeroIndices.length % 2) != 0) {
				GPUtil.exit(
						"The number of items in each class 0 must be an even number.");
			}

			if((classOneIndices.length % 2) != 0) {
				GPUtil.exit(
						"The number of items in each class 1 must be an even number.");
			}
		}

		RankFeatureSelector datasetFeatureSelector = new RankFeatureSelector(numFeatures,
				dataset,
				classZeroIndices,
				classOneIndices, side == TWO_SIDED);
		Dataset topFeaturesDataset = dataset.slice(
				datasetFeatureSelector.sortedIndices,
				null);

		Permuter permuter = null;

		if(!complete && balanced) {
			permuter = new BalancedRandomPermuter(classZeroIndices, classOneIndices);
		} else if(!complete && !balanced) {
			int[] classAssignments = null;
			if(side == CLASS_ZERO_LESS_THAN_CLASS_ONE) {
				classAssignments = invert(classVector);
			} else {
				classAssignments = classVector.getAssignments();
			}
			permuter = new UnbalancedRandomPermuter(classAssignments);
		} else if(complete && !balanced) {
			permuter = new UnbalancedCompletePermuter(classVector.size(),
					classZeroIndices.length);
			java.math.BigInteger totalPermutations = ((UnbalancedCompletePermuter) (permuter)).getTotal();
			if((totalPermutations.compareTo(new java.math.BigInteger("" + Integer.MAX_VALUE))) == 1) {
				GPUtil.exit("Number of permutations exceeds maximum of " + Integer.MAX_VALUE);
			}
			numPermutations  = totalPermutations.intValue();
		} else if(complete && balanced) {
			GPUtil.exit("Algorithm not yet implemented for complete and balanced permutations.");
		}
		double[][] rankBasedScores = new double[numFeatures][numPermutations];
		double[][] geneSpecificScores = new double[numFeatures][numPermutations];

		for(int perm = 0; perm < numPermutations; perm++) {

			int[] permutedAssignments = null;
			int[] permutedClassZeroIndices = null;
			int[] permutedClassOneIndices = null;

			if(readPermutations) {
				permutedAssignments = testAssignments[perm];

				ClassVector temp = new ClassVector(permutedAssignments,
						new int[]{0, 1});

				if(side == CLASS_ZERO_LESS_THAN_CLASS_ONE) {
					permutedClassZeroIndices = temp.getIndices(1);
					permutedClassOneIndices = temp.getIndices(0);
				} else {
					permutedClassZeroIndices = temp.getIndices(0);
					permutedClassOneIndices = temp.getIndices(1);
				}
			} else {
				permutedAssignments = permuter.next();
				ClassVector temp = new ClassVector(permutedAssignments,
						new int[]{0, 1});

				permutedClassZeroIndices = temp.getIndices(0);
				permutedClassOneIndices = temp.getIndices(1);

			}

			RankFeatureSelector rankFeatureSelector = new RankFeatureSelector(
					numFeatures,
					dataset,
					permutedClassZeroIndices,
					permutedClassOneIndices,
					side == TWO_SIDED);

			for(int feature = 0; feature < numFeatures; feature++) {
				rankBasedScores[feature][perm] = rankFeatureSelector.Sx[feature];
			}

			new GeneSpecificFeatureSelector(topFeaturesDataset,
					permutedClassZeroIndices,
					permutedClassOneIndices,
					geneSpecificScores, perm);
		}

		double[] rankBasedPValues = null;
		double[] geneSpecificPValues = null;

		if(side != TWO_SIDED) {
			rankBasedPValues = computeTwoSidedPValues(
					datasetFeatureSelector.Sx,
					rankBasedScores);
			geneSpecificPValues = computeOneSidedPValues(
					datasetFeatureSelector.Sx,
					geneSpecificScores);
		} else {
			rankBasedPValues = computeTwoSidedPValues(
					datasetFeatureSelector.Sx,
					rankBasedScores);
			geneSpecificPValues = computeTwoSidedPValues(
					datasetFeatureSelector.Sx,
					geneSpecificScores);
		}

		double[] fwer = new double[numFeatures];

		for(int feature = 0; feature < numFeatures; feature++) {

			int fwerCount = 0;
			double actualScore = Math.abs(datasetFeatureSelector.Sx[feature]);

			for(int i = 0; i < numFeatures; i++) {

				boolean iterationGreaterThanS = false;

				for(int j = numPermutations - 1; j >= 0; j--) {

					if(Math.abs(geneSpecificScores[i][j]) > actualScore) {
						iterationGreaterThanS = true;

						break;
					}
				}

				if(iterationGreaterThanS) {
					fwerCount++;
				}
			}

			fwer[feature] = (double) (fwerCount) / numPermutations;
		}

		PrintWriter pw = null;
		java.util.Vector features = new java.util.Vector(numFeatures);

		try {
			pw = new PrintWriter(new FileWriter(outputFileName));

			for(int feature = 0; feature < numFeatures; feature++) {
				features.add(topFeaturesDataset.getRowName(feature));
				pw.println(topFeaturesDataset.getRowName(feature) + "\t" +
						rankBasedPValues[feature] + "\t" +
						geneSpecificPValues[feature] + "\t " +
						fwer[feature]);
			}
		} catch(Exception e) {
			GPUtil.exit("An error occurred while saving the output file.", e);
		}

		if(pw != null) {

			try {
				pw.close();
			} catch(Exception x) {
			}
		}

		new MarkerSelectionFrame(features, geneSpecificPValues,
				rankBasedPValues, fwer);
	}


	void readPermutations(String fileName) {

		try {
			testClassPermutationsReader = new BufferedReader(new FileReader(
					fileName));
			testAssignments = new int[numPermutations][];

			for(int i = 0; i < numPermutations; i++) {

				String s = testClassPermutationsReader.readLine();
				StringTokenizer st = new StringTokenizer(s, "\t");
				int cols = dataset.getColumnDimension();
				int[] a = new int[cols];

				for(int j = 0; j < cols; j++) {
					a[j] = Integer.parseInt(st.nextToken());
				}

				testAssignments[i] = a;
			}
		} catch(Exception e) {
			GPUtil.exit(
					"An error occurred while reading the permutations file.",
					e);
		}
	}


	double[] computeOneSidedPValues(double[] actualScores,
			double[][] permutedScores) {

		double[] oneSidedPValues = new double[numFeatures];

		for(int i = 0; i < numFeatures; i++) {

			for(int j = 0; j < numPermutations; j++) {

				if(permutedScores[i][j] >= actualScores[i]) {
					oneSidedPValues[i] = oneSidedPValues[i] + 1;
				}
			}

			oneSidedPValues[i] = oneSidedPValues[i] / (double) numPermutations;
		}

		return oneSidedPValues;
	}


	double[] computeTwoSidedPValues(double[] actualScores,
			double[][] permutedScores) {

		double[] twoSidedPValues = new double[numFeatures];

		for(int i = 0; i < numFeatures; i++) {

			for(int j = 0; j < numPermutations; j++) {

				if(Math.abs(permutedScores[i][j]) >= Math.abs(actualScores[i])) {
					twoSidedPValues[i] = twoSidedPValues[i] + 1;
				}
			}

			twoSidedPValues[i] = twoSidedPValues[i] / (double) numPermutations;
		}

		return twoSidedPValues;
	}


	public long getNumUnBalancedPermutations() {

		int numClassZero = classVector.getIndices(0).length;
		int numClassOne = classVector.getIndices(1).length;

		return cern.jet.math.Arithmetic.longFactorial(classVector.size()) / (cern.jet.math.Arithmetic.longFactorial(
				numClassZero) * cern.jet.math.Arithmetic.longFactorial(
				numClassOne));
	}


	/**
	 * @author     Joshua Gould
	 * @created    September 17, 2004
	 */
	static class RankFeatureSelector {

		int[] sortedIndices;
		double[] Sx;


		public RankFeatureSelector(int numFeatures, Dataset dataset,
				int[] classOneIndices, int[] class2Indices) {
			this(numFeatures, dataset, classOneIndices, class2Indices, false);
		}


		public RankFeatureSelector(int numFeatures, Dataset dataset,
				int[] classOneIndices, int[] class2Indices,
				boolean absolute) {
			sortedIndices = new int[numFeatures];
			Sx = new double[numFeatures];
			Arrays.fill(Sx, -Double.MAX_VALUE);

			for(int i = 0, rows = dataset.getRowDimension(); i < rows; i++) {

				double class1Mean = mean(dataset, classOneIndices, i);
				double class2Mean = mean(dataset, class2Indices, i);
				double class1Std = standardDeviation(dataset, classOneIndices,
						i, class1Mean);
				double class2Std = standardDeviation(dataset, class2Indices, i,
						class2Mean);
				double Sxi = (class1Mean - class2Mean) / (class1Std +
						class2Std);

				if(absolute) {
					Sxi = Math.abs(Sxi);
				}

				int insertionIndex = -1;

				for(int j = 0; j < numFeatures; j++) {

					if(Sxi > Sx[j]) {
						insertionIndex = j;

						break;
					}
				}

				if(insertionIndex != -1) {

					// shift everything to the right
					System.arraycopy(Sx, insertionIndex, Sx,
							insertionIndex + 1,
							numFeatures - insertionIndex - 1);
					Sx[insertionIndex] = Sxi;
					System.arraycopy(sortedIndices, insertionIndex,
							sortedIndices, insertionIndex + 1,
							numFeatures - insertionIndex - 1);
					sortedIndices[insertionIndex] = i;
				}
			}
		}
	}


	/**
	 * @author     Joshua Gould
	 * @created    September 17, 2004
	 */
	static class GeneSpecificFeatureSelector {

		double[] Sx;


		public GeneSpecificFeatureSelector(Dataset dataset,
				int[] classOneIndices,
				int[] class2Indices,
				double[][] scores,
				int permutationIndex) {

			int rows = dataset.getRowDimension();
			Sx = new double[rows];

			for(int i = 0; i < rows; i++) {

				double class1Mean = mean(dataset, classOneIndices, i);
				double class2Mean = mean(dataset, class2Indices, i);
				double class1Std = standardDeviation(dataset, classOneIndices,
						i, class1Mean);
				double class2Std = standardDeviation(dataset, class2Indices, i,
						class2Mean);
				double Sxi = (class1Mean - class2Mean) / (class1Std +
						class2Std);
				scores[i][permutationIndex] = Sxi;
			}
		}
	}
}

