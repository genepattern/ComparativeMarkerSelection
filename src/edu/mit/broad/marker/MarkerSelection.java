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
 *@author     Joshua Gould
 *@created    September 17, 2004
 */
public class MarkerSelection {

	/**  input expression values */
	Dataset dataset;

	/**  input classes */
	ClassVector classVector;

	/**  number of permutations to perform */
	int numPermutations;
	final static int CLASS_ZERO_GREATER_THAN_CLASS_ONE = 0;
	final static int CLASS_ZERO_LESS_THAN_CLASS_ONE = 1;
	final static int TWO_SIDED = 2;

	/**  whether to do one-sided >, one-sided <, or two-sided test */
	int testDirection;

	int metric;
	final static int T_TEST = 0;
	final static int SNR = 1;

	/**  whether permutations are balanced */
	boolean balanced;

	/**  whether to compute all possible permutations */
	boolean complete;
	/**  whether to read permutations from a file */
	boolean readPermutations = false;
	/**  reader when reading permutations from file */
	BufferedReader testClassPermutationsReader;

	/**  name of output file */
	String outputFileName;

	/**  unpermuted scores */
	double[] scores;

	/**  number of rows in input data */
	int N;

	/**  Whether to fix the standard deviation, as is done in GeneCluster */
	boolean fixStdev;
	StatisticalMeasure statisticalMeasure;


	public MarkerSelection(Dataset _dataset, ClassVector _classVector,
			int _numPermutations, int _side,
			String _outputFileName, boolean _balanced,
			boolean complete, boolean fixStdev, int metric, String permutationsFile) {
		this.dataset = _dataset;
		this.classVector = _classVector;

		GPUtil.checkDimensions(dataset, classVector);

		this.numPermutations = _numPermutations;
		this.testDirection = _side;
		this.balanced = _balanced;
		this.outputFileName = _outputFileName;
		this.complete = complete;
		this.N = dataset.getRowDimension();
		this.fixStdev = fixStdev;
		this.metric = metric;
		if(permutationsFile != null) {
			readPermutations = true;
			initFile(permutationsFile);
		}
		if(metric == T_TEST) {
			statisticalMeasure = new TTest();
		} else if(metric == SNR) {
			statisticalMeasure = new SNR();
		} else {
			GPUtil.exit("Unknown test statistic");
		}
		computePValues();

		if(readPermutations) {

			if(testClassPermutationsReader != null) {

				try {
					testClassPermutationsReader.close();
				} catch(Exception e) {
				}
			}
		}
	}


	public static void main(String[] args)
			 throws Exception {

		String datasetFile = args[0];
		String clsFile = args[1];
		int _numPermutations = Integer.parseInt(args[2]);
		int testDirection = Integer.parseInt(args[3]);
		String outputFileName = args[4];
		boolean balanced = Boolean.valueOf(args[5]).booleanValue();
		boolean complete = Boolean.valueOf(args[6]).booleanValue();
		boolean fixStdev = Boolean.valueOf(args[7]).booleanValue();
		int metric = Integer.parseInt(args[8]);
		String permutationsFile = null;
		if(args.length == 10) {
			permutationsFile = args[9];
		}
		ExpressionDataReader reader = GPUtil.getExpressionReader(datasetFile);
		ExpressionDataImpl e = (ExpressionDataImpl) reader.read(datasetFile);
		new MarkerSelection(e.getExpressionMatrix(),
				GPUtil.getClassVector(clsFile),
				_numPermutations, testDirection, outputFileName, balanced,
				complete, fixStdev, metric, permutationsFile);
	}


	void computePValues() {
		int[] classZeroIndices = classVector.getIndices(0);
		int[] classOneIndices = classVector.getIndices(1);

		if(balanced) {
			if(classZeroIndices.length != classOneIndices.length) {
				GPUtil.exit(
						"The number of items in each class must be equal for balanced permutations.");
			}
			if((classZeroIndices.length % 2) != 0) {
				GPUtil.exit(
						"The number of items in class 0 must be an even number for balanced permutations.");
			}
			if((classOneIndices.length % 2) != 0) {
				GPUtil.exit(
						"The number of items in class 1 must be an even number for balanced permutations.");
			}
		}

		scores = new double[N];

		statisticalMeasure.compute(dataset,
				classZeroIndices,
				classOneIndices, scores, fixStdev);

		Permuter permuter = null;

		if(!complete && balanced) {
			permuter = new BalancedRandomPermuter(classZeroIndices, classOneIndices);
		} else if(!complete && !balanced) {
			permuter = new UnbalancedRandomPermuter(classVector.size(), classVector.getIndices(1).length);
		} else if(complete && !balanced) {
			permuter = new UnbalancedCompletePermuter(classVector.size(),
					classZeroIndices.length);
			java.math.BigInteger totalPermutations = ((UnbalancedCompletePermuter) (permuter)).getTotal();
			if((totalPermutations.compareTo(new java.math.BigInteger("" + Integer.MAX_VALUE))) == 1) {
				GPUtil.exit("Number of permutations exceeds maximum of " + Integer.MAX_VALUE);
			}
			numPermutations = totalPermutations.intValue();
		} else if(complete && balanced) {
			permuter = new BalancedCompletePermuter(classZeroIndices, classOneIndices);
			java.math.BigInteger totalPermutations = ((BalancedCompletePermuter) (permuter)).getTotal();
			if((totalPermutations.compareTo(new java.math.BigInteger("" + Integer.MAX_VALUE))) == 1) {
				GPUtil.exit("Number of permutations exceeds maximum of " + Integer.MAX_VALUE);
			}
		}

		int[] descendingIndices = Util.getIndices(scores, Util.DESCENDING);
		double[] absoluteScores = (double[]) scores.clone();
		for(int i = 0; i < N; i++) {
			absoluteScores[i] = Math.abs(absoluteScores[i]);
		}
		int[] descendingAbsIndices = Util.getIndices(absoluteScores, Util.DESCENDING);

		double[] rankBasedPValues = new double[N];
		double[] fwer = new double[N];
		double[] fpr = new double[N];
		double[] geneSpecificPValues = new double[N];
		double[] permutedScores = new double[N];
		int[] levels = {0, 1};

		for(int perm = 0; perm < numPermutations; perm++) {
			int[] permutedAssignments = null;
			int[] permutedClassZeroIndices = null;
			int[] permutedClassOneIndices = null;

			if(readPermutations) {
				permutedAssignments = nextPermutation();
				ClassVector temp = new ClassVector(permutedAssignments,
						levels);
				permutedClassZeroIndices = temp.getIndices(0);
				permutedClassOneIndices = temp.getIndices(1);

			} else {
				permutedAssignments = permuter.next();
				ClassVector temp = new ClassVector(permutedAssignments,
						levels);

				permutedClassZeroIndices = temp.getIndices(0);
				permutedClassOneIndices = temp.getIndices(1);

			}

			statisticalMeasure.compute(dataset,
					permutedClassZeroIndices,
					permutedClassOneIndices, permutedScores, fixStdev); // compute scores using permuted class labels

			for(int i = 0; i < N; i++) {
				double score = scores[i];
				if(testDirection == TWO_SIDED || testDirection == CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
					if(permutedScores[i] >= score) {
						geneSpecificPValues[i] += 1.0;
					}
				} else {
					if(permutedScores[i] <= score) {
						geneSpecificPValues[i] += 1.0;
					}
				}

			}

			Util.sort(permutedScores, Util.DESCENDING);

			for(int i = 0; i < N; i++) {
				double score = scores[descendingIndices[i]];

				if(testDirection == TWO_SIDED) {
					if(score >= 0) {
						if(permutedScores[i] >= score) {
							rankBasedPValues[i] += 1.0;
						}
					} else {
						if(permutedScores[i] <= score) {
							rankBasedPValues[i] += 1.0;
						}
					}
				} else if(testDirection == CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
					if(permutedScores[i] >= score) {
						rankBasedPValues[i] += 1.0;
					}
				} else {
					if(permutedScores[i] <= score) {
						rankBasedPValues[i] += 1.0;
					}
				}
			}

			for(int i = 0; i < N; i++) {
				permutedScores[i] = Math.abs(permutedScores[i]);
			}
			Util.sort(permutedScores, Util.DESCENDING);

			int j = 0;
			int count = 0;
			for(int i = 0; i < N; i++) {
				double score = absoluteScores[descendingAbsIndices[i]];
				while(j < N && score < permutedScores[j]) {
					count++;
					j++;
				}
				fpr[descendingAbsIndices[i]] += count;

				if(count > 0) {
					fwer[descendingAbsIndices[i]]++;
				}
			}
		}

		for(int i = 0; i < N; i++) {
			fpr[i] /= N;
			fpr[i] /= numPermutations;
			geneSpecificPValues[i] /= numPermutations;
			if(testDirection == TWO_SIDED) {
				geneSpecificPValues[i] = 2.0 * Math.min(geneSpecificPValues[i], 1.0 - geneSpecificPValues[i]);
			}
			fwer[i] /= numPermutations;
			rankBasedPValues[i] /= numPermutations;
		}

		double[] fdr = new double[N];
		int[] pValueIndices = Util.getIndices(fpr, Util.ASCENDING);
		int[] ranks = Util.rank(pValueIndices);

		// check for ties
		for(int i = pValueIndices.length - 1; i > 0; i--) {
			double bigPValue = fpr[pValueIndices[i]];
			double smallPValue = fpr[pValueIndices[i - 1]];
			if(bigPValue == smallPValue) {
				ranks[pValueIndices[i - 1]] = ranks[pValueIndices[i]];
			}
		}

		for(int i = 0; i < N; i++) {
			int index = i;
			int rank = ranks[index];
			double p = fpr[index];
			fdr[i] = (p * N) / rank;
		}

		int[] _ranks = null;

		if(testDirection == CLASS_ZERO_GREATER_THAN_CLASS_ONE || testDirection == CLASS_ZERO_LESS_THAN_CLASS_ONE) {
			_ranks = Util.rank(descendingIndices);
		} else if(testDirection == TWO_SIDED) {
			_ranks = Util.rank(descendingAbsIndices);
		}

		PrintWriter pw = null;

		try {
			pw = new PrintWriter(new FileWriter(outputFileName));
			pw.println("ODF 1.0");
			pw.println("HeaderLines=10");
			pw.println("COLUMN_NAMES:Rank\tFeature\tScore\tGene Specific P Value\tFPR\tFWER\tRank Based P Value\tFDR\tBonferroni");
			pw.println("COLUMN_TYPES:int\tString\tdouble\tdouble\tdouble\tdouble\tdouble\tdouble\tdouble");
			pw.println("Model=Comparative Marker Selection");
			pw.println("Permutations=" + numPermutations);
			pw.println("Balanced=" + balanced);
			pw.println("Complete=" + complete);
			if(testDirection == CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
				pw.println("Test=Class 0");
			} else if(testDirection == CLASS_ZERO_LESS_THAN_CLASS_ONE) {
				pw.println("Test=Class 1");
			} else {
				pw.println("Test=2 Sided");
			}
			if(metric == SNR) {
				pw.println("Statistical Measure=SNR");
			} else {
				pw.println("Statistical Measure=T-Test");
			}

			pw.println("Fix Standard Deviation=" + fixStdev);
			pw.println("DataLines=" + N);


			for(int i = 0; i < N; i++) {
				int index = descendingIndices[i];

				int rank = _ranks[index];
				if(testDirection == CLASS_ZERO_LESS_THAN_CLASS_ONE) {
					rank = N - rank + 1;
				}
				double bonferroni = Math.min(geneSpecificPValues[index] * N, 1.0);
				pw.println(rank + "\t" + dataset.getRowName(index) + "\t" +
						scores[index] + "\t" + geneSpecificPValues[index] + "\t " + fpr[index] + "\t" + fwer[index] + "\t" + rankBasedPValues[i] + "\t" + fdr[index] + "\t" + bonferroni);
			}

		} catch(Exception e) {
			GPUtil.exit("An error occurred while saving the output file.", e);
		} finally {

			if(pw != null) {

				try {
					pw.close();
				} catch(Exception x) {
				}
			}
		}

	}



	/**
	 *  Reads the next permutation from a file
	 *
	 *@return    the next permutation
	 */
	int[] nextPermutation() {
		try {
			String s = testClassPermutationsReader.readLine();
			StringTokenizer st = new StringTokenizer(s, " \t");
			int cols = dataset.getColumnDimension();
			int[] a = new int[cols];

			for(int j = 0; j < cols; j++) {
				a[j] = Integer.parseInt(st.nextToken());
			}

			return a;
		} catch(Exception e) {
			GPUtil.exit("An error occurred while reading the permutations file.", e);
			return null;
		}

	}


	/**
	 *  Opens the permutations file for reading
	 *
	 *@param  fileName  the file name
	 */
	void initFile(String fileName) {

		try {
			testClassPermutationsReader = new BufferedReader(new FileReader(
					fileName));
		} catch(Exception e) {
			GPUtil.exit(
					"An error occurred while reading the permutations file.",
					e);
		}
	}


	/**
	 *@author     Joshua Gould
	 *@created    September 29, 2004
	 */
	static class TTest implements StatisticalMeasure {
		public void compute(Dataset dataset,
				int[] classZeroIndices,
				int[] classOneIndices, double[] scores, boolean fixStdev) {
			Util.ttest(dataset,
					classZeroIndices,
					classOneIndices, scores, fixStdev);

		}
	}


	/**
	 *@author     Joshua Gould
	 *@created    September 29, 2004
	 */
	static class SNR implements StatisticalMeasure {
		public void compute(Dataset dataset,
				int[] classZeroIndices,
				int[] classOneIndices, double[] scores, boolean fixStdev) {
			Util.snr(dataset,
					classZeroIndices,
					classOneIndices, scores, fixStdev);

		}
	}


	/**
	 *@author     Joshua Gould
	 *@created    September 29, 2004
	 */
	static interface StatisticalMeasure {
		public void compute(Dataset dataset,
				int[] classZeroIndices,
				int[] classOneIndices, double[] scores, boolean fixStdev);
	}

}

