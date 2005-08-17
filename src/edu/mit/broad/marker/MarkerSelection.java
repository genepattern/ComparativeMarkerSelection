package edu.mit.broad.marker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;

import org.genepattern.data.expr.ExpressionData;
import org.genepattern.data.matrix.ClassVector;
import org.genepattern.data.matrix.DoubleMatrix2D;
import org.genepattern.io.expr.ExpressionDataCreator;
import org.genepattern.io.expr.IExpressionDataReader;
import org.genepattern.module.AnalysisUtil;
import org.genepattern.stats.ITestStatistic;
import org.genepattern.stats.Sorting;
import org.genepattern.stats.TestStatistics;
import org.genepattern.ioutil.Util;

import edu.mit.broad.marker.permutation.BalancedCompletePermuter;
import edu.mit.broad.marker.permutation.BalancedRandomPermuter;
import edu.mit.broad.marker.permutation.Permuter;
import edu.mit.broad.marker.permutation.UnbalancedCompletePermuter;
import edu.mit.broad.marker.permutation.UnbalancedRandomCovariatePermuter;
import edu.mit.broad.marker.permutation.UnbalancedRandomPermuter;

/**
 * @author     Joshua Gould
 */
public class MarkerSelection {

	/** input expression values */
	private DoubleMatrix2D dataset;

	/** data array stored in <tt>dataset</tt> */
	private double[][] dataArray;
	
	/** input classes */
	private ClassVector classVector;

	/** total number of permutations to perform */
	private int numPermutations;

	/** whether to do one-sided >, one-sided <, or two-sided test */
	private int testDirection;

	private int metric;

	/**
	 * minimum standard deviation when metric==T_TEST_MIN_STD or SNR_MIN_STD or
	 * SNR_MEDIAN_MIN_STD
	 */
	private double minStd;

	/** whether permutations are balanced */
	private boolean balanced;

	/** whether to compute all possible permutations */
	private boolean complete;

	/** whether to read permutations from a file */
	private boolean readPermutations = false;

	/** reader when reading permutations from file */
	private BufferedReader testClassPermutationsReader;

	/** name of output file */
	private String outputFileName;

	/** unpermuted scores */
	private double[] scores;

	/** number of rows in input data */
	private int numFeatures;

	private ITestStatistic statisticalMeasure;

	/** Input dataset file name */
	private String datasetFile;

	/** Input cls file name */
	private String clsFile;
	
	/** Input confounding variable cls file name */
	private String confoundingClsFile;

	/** Covariate assignments */
	private ClassVector covariate;

	/** Seed for generating permutations */
	private int seed;

	private static Debugger debugger;

	/** Whether a seed for used for generating permutations */
	private boolean seedUsed = false;

	/** Family-wise error rate, size=numFeatures */
	private double[] fwer;

	/** size=numFeatures */
	private double[] featureSpecificPValues;

	/** Temporary array used in calculation of permuted scores */
	private double[] permutedScores;

	/** Temporary array used in calculation of permuted score */
	private double[] monotonicPermutedScores;

	/** size=numFeatures */
	private double[] maxT;

	/** Absolute value of the test statistic, size=numFeatures */
	private double[] absoluteScores;

	/** descending sorted indices into absoluteScores */
	private int[] descendingAbsIndices;

	private Permuter permuter;

	/** Whether to try to trim the number of features to permute */
	private boolean speedUp;

	/** Whether to smooth p values */
	private boolean smoothPValues;

	public MarkerSelection(DoubleMatrix2D dataset, String datasetFile, 
			ClassVector classVector, String clsFile,
			int _numPermutations, int _side, String _outputFileName,
			boolean _balanced, boolean complete, int metric,
			double minStd, int seed, ClassVector confoundingClassVector,
			String confoundingClsFile, boolean speedup, boolean smoothPValues) {
		this.dataset = dataset;
		this.dataArray = dataset.getArray();
		this.datasetFile = datasetFile;
		this.classVector = AnalysisUtil.readClassVector(clsFile);
		AnalysisUtil.checkDimensions(dataset, classVector);
		if (classVector.getClassCount() != 2) {
			AnalysisUtil.exit("Class file must contain 2 classes.");
		}
		this.clsFile = clsFile;
		this.numPermutations = _numPermutations;
		this.testDirection = _side;
		this.outputFileName = _outputFileName;
		this.balanced = _balanced;
		this.complete = complete;
		this.metric = metric;
		this.minStd = minStd;
		this.seed = seed;
		this.covariate = confoundingClassVector;
		this.confoundingClsFile = confoundingClsFile;
		this.speedUp = speedup;
		this.smoothPValues = smoothPValues;
		if(complete && speedup) {
			AnalysisUtil.exit("Speedup option can only be used when performing random permutations.");
		}
		
		this.numFeatures = dataset.getRowCount();
		
		String permutationsFile = System
				.getProperty("edu.mit.broad.marker.perm");
		if (permutationsFile != null) {
			System.out.println("reading permutations from file "
					+ permutationsFile);
			readPermutations = true;
			initFile(permutationsFile);
		}
		if (metric == Constants.T_TEST) {
			statisticalMeasure = new TestStatistics.TTest();
		} else if (metric == Constants.T_TEST_MEDIAN) {
			statisticalMeasure = new TestStatistics.TTestMedian();
		} else if (metric == Constants.T_TEST_MIN_STD) {
			if (minStd <= 0) {
				AnalysisUtil
						.exit("Minimum standard deviation must be greater than zero.");
			}
			statisticalMeasure = new TestStatistics.TTestMinStd(minStd);
		}  else if (metric == Constants.T_TEST_MEDIAN_MIN_STD) {
			if (minStd <= 0) {
				AnalysisUtil
						.exit("Minimum standard deviation must be greater than zero.");
			}
			statisticalMeasure = new TestStatistics.TTestMedianMinStd(minStd);
					
		} else if (metric == Constants.SNR) {
			statisticalMeasure = new TestStatistics.SNR();
		} else if (metric == Constants.SNR_MEDIAN) {
			statisticalMeasure = new TestStatistics.SNRMedian();
		
		} else if (metric == Constants.SNR_MIN_STD) {
			if (minStd <= 0) {
				AnalysisUtil
						.exit("Minimum standard deviation must be greater than zero.");
			}
			statisticalMeasure = new TestStatistics.SNRMinStd(minStd);
		} else if (metric == Constants.SNR_MEDIAN_MIN_STD) {
			if (minStd <= 0) {
				AnalysisUtil
						.exit("Minimum standard deviation must be greater than zero.");
			}
			statisticalMeasure = new TestStatistics.SNRMedianMinStd(minStd);
			
		} else {
			AnalysisUtil.exit("Unknown test statistic");
		}

		if (testDirection != Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE
				&& testDirection != Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE
				&& testDirection != Constants.TWO_SIDED) {
			AnalysisUtil.exit("Unknown test direction.");
		}

		computePValues();

		if (readPermutations) {

			if (testClassPermutationsReader != null) {

				try {
					testClassPermutationsReader.close();
				} catch (Exception e) {
				}
			}
		}
	}

	public static void main(String[] args) {
		try {
			run(args);
		} catch (Throwable t) {
			AnalysisUtil.exit("An error occurred while running the algorithm.");
		}
	}

	static void run(String[] args) {
		String datasetFile = args[0];
		String clsFile = args[1];
		int _numPermutations = 0;
		try {
			_numPermutations = Integer.parseInt(args[2]);
		} catch (NumberFormatException nfe) {
			AnalysisUtil.exit("number of permutations is not an integer");
		}
		int testDirection = Integer.parseInt(args[3]);
		String outputFileName = args[4];
		boolean balanced = Boolean.valueOf(args[5]).booleanValue();
		boolean complete = Boolean.valueOf(args[6]).booleanValue();
		boolean fixStdev = Boolean.valueOf(args[7]).booleanValue();
		int metric = Integer.parseInt(args[8]);
		int seed = 0;
		try {
			seed = Integer.parseInt(args[9]);
		} catch (NumberFormatException nfe) {
			AnalysisUtil.exit("random seed is not an integer");
		}
		boolean speedup = Boolean.valueOf(args[10]).booleanValue();
		boolean smoothPValues = Boolean.valueOf(args[11]).booleanValue();
		
		double minStd = -1;
		String confoundingClsFile = null;
		for (int i = 12; i < args.length; i++) {
			String arg = args[i].substring(0, 2);
			String value = args[i].substring(2, args[i].length());
			if (value.equals("")) {
				continue;
			}
			if (arg.equals("-m")) {
				try {
					minStd = Double.parseDouble(value);
				} catch (NumberFormatException nfe) {
					AnalysisUtil
							.exit("minimum standard deviation is not a number");
				}
			} else if (arg.equals("-c")) {
				confoundingClsFile = value;
			}
		}

		ClassVector classVector = AnalysisUtil.readClassVector(clsFile);
		IExpressionDataReader reader = AnalysisUtil
				.getExpressionReader(datasetFile);
		ExpressionData expressionData = AnalysisUtil
				.readExpressionData(reader, datasetFile);

		DoubleMatrix2D dataset = expressionData.getExpressionMatrix();
		
		ClassVector confoundingClassVector = null;
		if (confoundingClsFile != null) {
			confoundingClassVector = AnalysisUtil.readClassVector(confoundingClsFile);
			/*  ClassVector[] cv = new edu.mit.broad.internal.MultiClassReader().read(clsFile, expressionData);
			 this.classVector = cv[0];
			 if(cv.length > 1) {
			 covariate = cv[1];
			 for(int i = 2; i < cv.length; i++) {
			 covariate = covariate.union(cv[i]);
			 }
			 }
			 */
		}


		if(classVector.getClassCount() > 2) {
			ClassVector[] oneVersusAll = classVector.getOneVersusAll();
			String baseOutputFileName = Util.getBaseFileName(outputFileName);
			for(int i = 0; i < classVector.getClassCount(); i++) {
				String tempOutputFileName = baseOutputFileName + "." + 
					classVector.getClassName(i) + ".vs.All.odf";
				new MarkerSelection(dataset, datasetFile, oneVersusAll[i], 
				clsFile, _numPermutations, testDirection, tempOutputFileName, 
				balanced, complete, metric, minStd, seed, 
				confoundingClassVector, confoundingClsFile, speedup, smoothPValues);
			}
		} else {
			new MarkerSelection(dataset, datasetFile, classVector, clsFile,
				_numPermutations, testDirection, outputFileName, balanced, 
				complete, metric, minStd, seed, 
				confoundingClassVector, confoundingClsFile, speedup, smoothPValues);
		}
	}

	final void computePValues() {
		int[] classZeroIndices = classVector.getIndices(0);
		int[] classOneIndices = classVector.getIndices(1);
		if (balanced) {
			if (classZeroIndices.length != classOneIndices.length) {
				AnalysisUtil
						.exit("The number of items in each class must be equal for balanced permutations.");
			}
			if ((classZeroIndices.length % 2) != 0) {
				AnalysisUtil
						.exit("The number of items in class 0 must be an even number for balanced permutations.");
			}
			if ((classOneIndices.length % 2) != 0) {
				AnalysisUtil
						.exit("The number of items in class 1 must be an even number for balanced permutations.");
			}
		}
		scores = new double[numFeatures];
		
		statisticalMeasure.compute(dataArray, classZeroIndices, classOneIndices,
				scores, null, -1); // index of observed scores, greatest to smallest
		int[] descendingIndices = Sorting.index(scores, Sorting.DESCENDING);
		
		if (covariate != null) {
			seedUsed = true;
			permuter = new UnbalancedRandomCovariatePermuter(classVector,
					covariate, seed);
			if (complete || balanced) {
				AnalysisUtil
						.exit("Covariate permutations not yet implemented for complete or balanced permutations.");
			}
		} else if (!complete && balanced) {
			seedUsed = true;
			permuter = new BalancedRandomPermuter(classZeroIndices,
					classOneIndices, seed);

		} else if (!complete && !balanced) {
			seedUsed = true;
			permuter = new UnbalancedRandomPermuter(classVector.size(),
					classVector.getIndices(1).length, seed);
		} else if (complete && !balanced) {
			seedUsed = false;
			permuter = new UnbalancedCompletePermuter(classVector.size(),
					classZeroIndices.length);
			java.math.BigInteger totalPermutations = ((UnbalancedCompletePermuter) (permuter))
					.getTotal();
			if ((totalPermutations.compareTo(new java.math.BigInteger(""
					+ Integer.MAX_VALUE))) == 1) {
				AnalysisUtil.exit("Number of permutations exceeds maximum of "
						+ Integer.MAX_VALUE);
			}
			numPermutations = totalPermutations.intValue();
		} else if (complete && balanced) {
			seedUsed = false;
			permuter = new BalancedCompletePermuter(classZeroIndices,
					classOneIndices);

			java.math.BigInteger totalPermutations = ((BalancedCompletePermuter) (permuter))
					.getTotal();
			if ((totalPermutations.compareTo(new java.math.BigInteger(""
					+ Integer.MAX_VALUE))) == 1) {
				AnalysisUtil.exit("Number of permutations exceeds maximum of "
						+ Integer.MAX_VALUE);
			}
			numPermutations = totalPermutations.intValue();
		}

		featureSpecificPValues = new double[numFeatures];
		permutedScores = new double[numFeatures];
		absoluteScores = (double[]) scores.clone();
		for (int i = 0; i < numFeatures; i++) {
				absoluteScores[i] = Math.abs(absoluteScores[i]);
		}
		descendingAbsIndices = Sorting
					.index(absoluteScores, Sorting.DESCENDING);
					
		if (speedUp) {
			long allowableCost = numFeatures*numPermutations;
			int[] permutationsPerFeature = new int[numFeatures];
			
			int numFeaturesToPermute = numFeatures;
			int permutationsToPerform = 100;
			double removalThreshold = 0.1; // FIXME
			int[] indicesToPermute = new int[numFeatures];
			for(int i = 0; i < numFeatures; i++) {
				indicesToPermute[i] = i;
			}
			int lengthOfIndicesToPermute = numFeatures;
			
			while (allowableCost > 0 && numFeaturesToPermute > 0) {
				System.out.println("performing " + permutationsToPerform);
				permute(permutationsToPerform, indicesToPermute, lengthOfIndicesToPermute);
				lengthOfIndicesToPermute = 0;
				
				for (int i = 0; i < numFeaturesToPermute; i++) {
					int index = indicesToPermute[i]; // index into feature p and score array
					permutationsPerFeature[index] += permutationsToPerform;
					double k = featureSpecificPValues[index];
					int N = permutationsPerFeature[index];
					double p;
					if (smoothPValues) {
						p = (k + 1) / (N + 2);
					} else {
						p = k / N;
					}
					if (testDirection == Constants.TWO_SIDED) {
						p = 2.0 * Math.min(p, 1.0 - p);
					}

					double shape1 = k + 1;
					double shape2 = N - k + 1;
					JSci.maths.statistics.BetaDistribution beta = new JSci.maths.statistics.BetaDistribution(
							shape1, shape2);
					double lowerBound = beta.inverse(0.025);
					double upperBound = beta.inverse(0.975);
				//	System.out.println("lowerBound " + lowerBound + " upperBound " + upperBound + " p " + p);
					
					double obsThreshold = 2 * (upperBound - lowerBound)
							/ (upperBound + lowerBound);
					if (obsThreshold < removalThreshold) { // exclude marker
						featureSpecificPValues[index] = p;
					} else {
						indicesToPermute[lengthOfIndicesToPermute] = index;
						lengthOfIndicesToPermute++;
					}
				}
				allowableCost -= numFeaturesToPermute*permutationsToPerform;
				numFeaturesToPermute = lengthOfIndicesToPermute;
				permutationsToPerform = 100; // FIXME
				
				long newCost = allowableCost - numFeaturesToPermute*permutationsToPerform;
				if(newCost < 0) {
					permutationsToPerform = (int)(allowableCost/numFeaturesToPermute);
				}
			}
			
			System.out.println("Features left " + numFeaturesToPermute);
			for (int i = 0; i < numFeaturesToPermute; i++) { // calculate p value for features that are left
				int index = indicesToPermute[i]; 
				System.out.println(dataset.getRowName(index));
				double k = featureSpecificPValues[index];
				int N = permutationsPerFeature[index];
				double p;
				if (smoothPValues) {
					p = (k + 1) / (N + 2);
				} else {
					p = k / N;
				}
				if (testDirection == Constants.TWO_SIDED) {
					p = 2.0 * Math.min(p, 1.0 - p);
				}
				featureSpecificPValues[index] = p;
			}
			
		} else {
			monotonicPermutedScores = new double[numFeatures];
			maxT = new double[numFeatures];
			fwer = new double[numFeatures];
			permute(numPermutations, null, -1);
			
			// free temporary storage
			
			monotonicPermutedScores = null;
			
		}
		absoluteScores = null;
		permutedScores = null; 
		
		if (!speedUp) {
			for (int i = 0; i < numFeatures; i++) {
				// fpr[i] /= numFeatures;
				// fpr[i] /= numPermutations;
				if (smoothPValues) {
					featureSpecificPValues[i] = (featureSpecificPValues[i] + 1)
							/ (numPermutations + 2);
				} else {
					featureSpecificPValues[i] /= numPermutations;
				}
				if (testDirection == Constants.TWO_SIDED) {
					featureSpecificPValues[i] = 2.0 * Math.min(
							featureSpecificPValues[i],
							1.0 - featureSpecificPValues[i]);
				}

				fwer[i] /= numPermutations;
				// rankBasedPValues[i] /= numPermutations;
				maxT[i] /= numPermutations;

			}
		}

		double[] fdr = new double[numFeatures];
		int[] pValueIndices = Sorting.index(featureSpecificPValues,
				Sorting.ASCENDING);
		int[] ranks = Sorting.rank(pValueIndices);

		// check for ties
		for (int i = pValueIndices.length - 1; i > 0; i--) {
			double bigPValue = featureSpecificPValues[pValueIndices[i]];
			double smallPValue = featureSpecificPValues[pValueIndices[i - 1]];
			if (bigPValue == smallPValue) {
				ranks[pValueIndices[i - 1]] = ranks[pValueIndices[i]];
			}
		}

		for (int i = 0; i < numFeatures; i++) {
			int rank = ranks[i];
			double p = featureSpecificPValues[i];
			fdr[i] = (p * numFeatures) / rank;
		}

		// ensure fdr is monotonically decreasing
		int[] pIndices = Sorting.index(featureSpecificPValues,
				Sorting.DESCENDING);
		for (int i = 0; i < pIndices.length - 1; i++) {
			int highIndex = pIndices[i];
			int lowIndex = pIndices[i + 1];
			fdr[lowIndex] = Math.min(fdr[lowIndex], fdr[highIndex]);
		}

		for (int i = 0; i < numFeatures; i++) {
			fdr[i] = Math.min(fdr[i], 1);
		}

		if (!speedUp) {
			for (int i = 1; i < numFeatures; i++) { // ensure maxT is monotonic
				maxT[descendingAbsIndices[i]] = Math.max(
						maxT[descendingAbsIndices[i]],
						maxT[descendingAbsIndices[i - 1]]);
			}
		}
		int[] _ranks = null;

		if (testDirection == Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE
				|| testDirection == Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE) {
			_ranks = Sorting.rank(descendingIndices);
		} else if (testDirection == Constants.TWO_SIDED) {
			_ranks = Sorting.rank(descendingAbsIndices);
		}

		// write p values to temporary file
		PrintWriter tempFileWriter = null;
		String tempFileName = "pvalues.txt";
		new File(tempFileName).deleteOnExit();
		try {
			tempFileWriter = new PrintWriter(new FileWriter(tempFileName));
			for (int i = 0; i < numFeatures; i++) {
				tempFileWriter
						.println(featureSpecificPValues[descendingIndices[i]]);
			}
		} catch (IOException ioe) {
			AnalysisUtil.exit(
					"An error occurred while saving the output file.", ioe);
		} finally {
			if (tempFileWriter != null) {
				tempFileWriter.close();
			}
		}

		edu.mit.broad.marker.qvalue.QValue.qvalue(tempFileName);
		String[] qvalues = new String[numFeatures];
		int qvalueHeaderLines = 4;
		String[] qvalueHeaders = new String[qvalueHeaderLines];
		BufferedReader br = null;
		String qvalueOutputFileName = "qvalues.txt"; // produced by R code above
		new File(qvalueOutputFileName).deleteOnExit();
		
		try {// read in results from qvalue
			br = new BufferedReader(new FileReader(qvalueOutputFileName));
			for (int i = 0; i < qvalueHeaderLines; i++) {
				qvalueHeaders[i] = br.readLine();
			}
			for (int i = 0; i < numFeatures; i++) {
				qvalues[i] = br.readLine();
			}
		} catch (IOException ioe) {
			AnalysisUtil.exit(
					"An error occurred while saving the output file.");
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException x) {
				}
			}
		}

		PrintWriter pw = null;

		try {
			if (!outputFileName.toLowerCase().endsWith(".odf")) {
				outputFileName += ".odf";
			}
			pw = new PrintWriter(new FileWriter(outputFileName));
			pw.println("ODF 1.0");
			int numHeaderLines = 18;
			if(seedUsed) {
				numHeaderLines++;
			}
			if(covariate!=null) {
				numHeaderLines++;
			}
			
			pw.println("HeaderLines=" + numHeaderLines);

			pw
					.println("COLUMN_NAMES:Rank\tFeature\tScore\tFeature P\tFWER\tFDR(BH)\tBonferroni\tQ Value\tmaxT");
			pw.println("Model=Comparative Marker Selection");
			pw.println("Dataset File=" + AnalysisUtil.getFileName(datasetFile));
			pw.println("Class File=" + AnalysisUtil.getFileName(clsFile));
			if(covariate!=null) {
				pw.println("Confounding Class File=" + 
					AnalysisUtil.getFileName(confoundingClsFile));
			}
			pw.println("Permutations=" + numPermutations);
			pw.println("Balanced=" + balanced);
			pw.println("Complete=" + complete);
			if (testDirection == Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
				pw.println("Test Direction=Class 0");
			} else if (testDirection == Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE) {
				pw.println("Test Direction=Class 1");
			} else {
				pw.println("Test Direction=2 Sided");
			}
			pw.println("Class 0=" + classVector.getClassName(0));
			pw.println("Class 1=" + classVector.getClassName(1));
			pw.println("Test Statistic=" + statisticalMeasure.toString());
			
			for (int i = 0; i < qvalueHeaderLines; i++) {
				pw.println(qvalueHeaders[i]);
			}
			if (seedUsed) {
				pw.println("Random Seed=" + seed);
			}
			
			pw.println("Speedup=" + speedUp);
			pw.println("Smooth p-values=" + smoothPValues);
			pw.println("DataLines=" + numFeatures);
	
			for (int i = 0; i < numFeatures; i++) {
				int index = descendingIndices[i];

				int rank = _ranks[index];
				if (testDirection == Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE) {
					rank = numFeatures - rank + 1;
				}
				double bonferroni = Math.min(featureSpecificPValues[index] * numFeatures,
						1.0);
				pw.print(rank);
				pw.print("\t");
				pw.print(dataset.getRowName(index));
				pw.print("\t");
				pw.print(scores[index]);
				pw.print("\t");
				pw.print(featureSpecificPValues[index]);
				pw.print("\t");
				if(!speedUp) {
					pw.print(fwer[index]);
				} else {
					pw.print("NaN");
				}
				pw.print("\t");
				pw.print(fdr[index]);
				pw.print("\t");
				pw.print(bonferroni);
				pw.print("\t");
				pw.print(qvalues[i]);
				pw.print("\t");
				if(!speedUp) {
					pw.println(maxT[index]);
				} else {
					pw.println("NaN");
				}
			}

		} catch (Exception e) {
			AnalysisUtil.exit(
					"An error occurred while saving the output file.", e);
		} finally {

			if (pw != null) {

				try {
					pw.close();
				} catch (Exception x) {
				}
			}
		}

		if (debugger != null) {
			debugger.print();
		}
	}

	private final void permute(int numPermutations, int[] rowIndices, int lengthOfRowIndices) {
		for (int perm = 0; perm < numPermutations; perm++) {
			int[] permutedClassZeroIndices = null;
			int[] permutedClassOneIndices = null;
			int[] permutedAssignments = null;
			if (readPermutations) {
				permutedAssignments = nextPermutation();
			} else {
				permutedAssignments = permuter.next();
			}

			if (debugger != null) {
				debugger.addAssignment(permutedAssignments);
			}
			java.util.List zeroIndices = new java.util.ArrayList();
			java.util.List oneIndices = new java.util.ArrayList();
			for (int i = 0, length = permutedAssignments.length; i < length; i++) {
				if (permutedAssignments[i] == 0) {
					zeroIndices.add(new Integer(i));
				} else {
					oneIndices.add(new Integer(i));
				}
			}
			permutedClassZeroIndices = new int[zeroIndices.size()];
			for (int i = 0, length = permutedClassZeroIndices.length; i < length; i++) {
				permutedClassZeroIndices[i] = ((Integer) zeroIndices.get(i))
						.intValue();
			}

			permutedClassOneIndices = new int[oneIndices.size()];
			for (int i = 0, length = permutedClassOneIndices.length; i < length; i++) {
				permutedClassOneIndices[i] = ((Integer) oneIndices.get(i))
						.intValue();
			}

			statisticalMeasure.compute(dataArray, permutedClassZeroIndices,
					permutedClassOneIndices, permutedScores, rowIndices, 
					lengthOfRowIndices);// compute scores using permuted class labels

			int length = rowIndices!=null?rowIndices.length:numFeatures;
			for (int i = 0; i < length; i++) {
				int index = rowIndices!=null?rowIndices[i]:i;
				double score = scores[index];
				if (testDirection == Constants.TWO_SIDED
						|| testDirection == Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
					if (permutedScores[index] >= score) {
						featureSpecificPValues[index] += 1.0;
					}
				} else {
					if (permutedScores[index] <= score) {
						featureSpecificPValues[index] += 1.0;
					}
				}

			}

			/*Sorting.sort(permutedScores, Sorting.DESCENDING);

			 for(int i = 0; i < numFeatures; i++) {
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
			 }*/

			if (!speedUp) {
				for (int i = 0; i < numFeatures; i++) {
					permutedScores[i] = Math.abs(permutedScores[i]);
					monotonicPermutedScores[i] = permutedScores[i];
				}

				for (int i = numFeatures - 2; i >= 0; i--) { // make permuted
																// scores
																// monotonicically
																// non-decreasing
																// for maxT
					int index = descendingAbsIndices[i];
					int indexPlusOne = descendingAbsIndices[i + 1];

					if (monotonicPermutedScores[indexPlusOne] > monotonicPermutedScores[index]) {
						monotonicPermutedScores[index] = monotonicPermutedScores[indexPlusOne];
					}

				}

				for (int i = 0; i < numFeatures; i++) {
					if (monotonicPermutedScores[i] >= absoluteScores[i]) {
						maxT[i] += 1.0;
					}
				}

				Sorting.sort(permutedScores, Sorting.DESCENDING);
				int j = 0;
				int count = 0;
				for (int i = 0; i < numFeatures; i++) {
					double score = absoluteScores[descendingAbsIndices[i]];

					while (j < numFeatures && score < permutedScores[j]) {
						count++;
						j++;
					}
					// fpr[descendingAbsIndices[i]] += count;

					if (count > 0) {
						fwer[descendingAbsIndices[i]]++;
					}
				}
			}

		}
	}

	/**
	 *  Reads the next permutation from a file
	 *
	 * @return    the next permutation
	 */
	final int[] nextPermutation() {
		try {
			String s = testClassPermutationsReader.readLine();
			StringTokenizer st = new StringTokenizer(s, " \t");
			int cols = dataset.getColumnCount();
			int[] a = new int[cols];

			for (int j = 0; j < cols; j++) {
				a[j] = Integer.parseInt(st.nextToken());
			}

			return a;
		} catch (Exception e) {
			AnalysisUtil
					.exit(
							"An error occurred while reading the permutations file.",
							e);
			return null;
		}

	}

	/**
	 *  Opens the permutations file for reading
	 *
	 * @param  fileName  the file name
	 */
	void initFile(String fileName) {

		try {
			testClassPermutationsReader = new BufferedReader(new FileReader(
					fileName));
		} catch (Exception e) {
			AnalysisUtil
					.exit(
							"An error occurred while reading the permutations file.",
							e);
		}
	}

	static class Debugger {
		java.util.Map assignment2Occurences = new java.util.HashMap();

		public static String toString(int[] a) {
			StringBuffer buf = new StringBuffer();
			for (int i = 0; i < a.length; i++) {
				buf.append(a[i]);
			}
			return buf.toString();
		}

		public void addAssignment(int[] a) {
			String s = toString(a);
			Integer occurences = (Integer) assignment2Occurences.get(s);
			if (occurences == null) {
				occurences = new Integer(0);
			}
			assignment2Occurences
					.put(s, new Integer(occurences.intValue() + 1));
		}

		public void print() {
			System.out.println("Unique permutations "
					+ assignment2Occurences.size());
			System.out.println("permutations");
			System.out.println(assignment2Occurences);
		}
	}

	static {
		if ("true".equalsIgnoreCase(System
				.getProperty("edu.mit.broad.marker.debug"))) {
			debugger = new Debugger();
		}
	}

}
