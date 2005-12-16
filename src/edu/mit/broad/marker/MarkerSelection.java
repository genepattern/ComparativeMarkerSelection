/*
  The Broad Institute
  SOFTWARE COPYRIGHT NOTICE AGREEMENT
  This software and its documentation are copyright (2003-2006) by the
  Broad Institute/Massachusetts Institute of Technology. All rights are
  reserved.

  This software is supplied without any warranty or guaranteed support
  whatsoever. Neither the Broad Institute nor MIT can be responsible for its
  use, misuse, or functionality.
*/


package edu.mit.broad.marker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.StringTokenizer;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.TTestImpl;
import org.apache.commons.math.stat.regression.SimpleRegression;
import org.genepattern.data.expr.ExpressionData;
import org.genepattern.data.matrix.ClassVector;
import org.genepattern.io.OdfWriter;
import org.genepattern.io.Util;
import org.genepattern.io.expr.IExpressionDataReader;
import org.genepattern.module.AnalysisUtil;
import org.genepattern.stats.Function;
import org.genepattern.stats.ITestStatistic;
import org.genepattern.stats.Sorting;
import org.genepattern.stats.TestStatistics;
import org.genepattern.stats.ZeroFinder;

import edu.mit.broad.marker.permutation.BalancedCompletePermuter;
import edu.mit.broad.marker.permutation.BalancedRandomPermuter;
import edu.mit.broad.marker.permutation.Permuter;
import edu.mit.broad.marker.permutation.UnbalancedCompletePermuter;
import edu.mit.broad.marker.permutation.UnbalancedRandomCovariatePermuter;
import edu.mit.broad.marker.permutation.UnbalancedRandomPermuter;
import edu.mit.broad.marker.qvalue.QValue;

/**
 * @author Joshua Gould
 */
public class MarkerSelection {
	private static boolean printStackTraces;

	private static boolean printPermutations;

	private static Debugger debugger;

	/** input expression values */
	private ExpressionData expressionData;

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

	/** Whether a seed was used for generating permutations */
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

	/** Lower confidence interval for feature p value */
	private double[] lowerBound;

	/** Upper confidence interval for feature p value */
	private double[] upperBound;

	/** Whether to try to trim the number of features to permute */
	private boolean significanceBooster;

	/** Keeps track of number of permutations performed per feature */
	private int[] permutationsPerFeature;

	/**
	 * Threshold for removing features when <tt>removeFeatures</tt> is
	 * <tt>true</tt>
	 */
	private double theta = 0.25;

	/** Whether to smooth p values */
	private boolean smoothPValues;

	private double gamma = 0;

	private boolean estimateGamma = false;

	/**
	 * if <tt>true</tt> don't calculate fwer, maxT even when significance
	 * booster is off
	 */
	private static boolean testGamma = false;

	private int[] featureIndicesToPermute;

	private int lengthOfIndicesToPermute;

	/** whether to calculate asymptotic p-values */
	private boolean asymptotic = false;

	/**
	 * Creates a new instance
	 * 
	 * @param _dataset
	 * @param _datasetFile
	 * @param _classVector
	 * @param _clsFile
	 * @param _numPermutations
	 * @param _side
	 * @param _outputFileName
	 * @param _balanced
	 * @param _complete
	 * @param _metric
	 * @param _minStd
	 * @param _seed
	 * @param _confoundingClassVector
	 * @param _confoundingClsFile
	 * @param _removeFeatures
	 * @param _theta
	 * @param _smoothPValues
	 * @param _gamma
	 *            if gamma >= 0 set gamma manually, otherwise estimate it
	 */
	public MarkerSelection(ExpressionData _dataset, String _datasetFile,
			ClassVector _classVector, String _clsFile, int _numPermutations,
			int _side, String _outputFileName, boolean _balanced,
			boolean _complete, int _metric, double _minStd, int _seed,
			ClassVector _confoundingClassVector, String _confoundingClsFile,
			boolean _removeFeatures, double _theta, boolean _smoothPValues,
			double _gamma) {
		this.expressionData = _dataset;
		this.dataArray = _dataset.getArray();
		this.datasetFile = _datasetFile;
		this.classVector = _classVector;
		AnalysisUtil.checkDimensions(_dataset, _classVector);
		if (_classVector.getClassCount() != 2) {
			AnalysisUtil.exit("Class file must contain 2 classes.");
		}
		this.clsFile = _clsFile;
		this.numPermutations = _numPermutations;
		this.testDirection = _side;
		this.outputFileName = _outputFileName;
		this.balanced = _balanced;
		this.complete = _complete;
		this.metric = _metric;
		this.minStd = _minStd;
		this.seed = _seed;
		this.covariate = _confoundingClassVector;
		this.confoundingClsFile = _confoundingClsFile;
		this.significanceBooster = _removeFeatures;
		this.theta = _theta;
		this.smoothPValues = _smoothPValues;
		if (significanceBooster) {
			this.gamma = _gamma;
			this.estimateGamma = gamma < 0;
		}

		if (_complete && _removeFeatures) {
			AnalysisUtil
					.exit("Speedup option can only be used when performing random permutations.");
		}
		if (_complete && _smoothPValues) {
			System.out
					.println("Smooth p-values set to false. Smoothing p-values disabled when performing all possible permutations.");
			this.smoothPValues = false;
		}
		if (significanceBooster && !smoothPValues) {
			System.out
					.println("Smooth p-values set to true. Smooth p-values must be true when using significance booster.");
			this.smoothPValues = true;
		}

		this.numFeatures = expressionData.getRowCount();

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
		} else if (metric == Constants.T_TEST_MEDIAN_MIN_STD) {
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
		// debug = true;
		// debugger = new Debugger();

		try {
			run(args);
		} catch (Throwable t) {
			if (printStackTraces) {
				t.printStackTrace();
			}
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
			AnalysisUtil.exit("Number of permutations is not an integer.");
		}
		int testDirection = Integer.parseInt(args[3]);
		String outputFileName = args[4];
		boolean balanced = Boolean.valueOf(args[5]).booleanValue();
		boolean complete = Boolean.valueOf(args[6]).booleanValue();
		boolean fixStdev = Boolean.valueOf(args[7]).booleanValue(); // ignored
		int metric = Integer.parseInt(args[8]);
		int seed = 0;
		try {
			seed = Integer.parseInt(args[9]);
		} catch (NumberFormatException nfe) {
			AnalysisUtil.exit("Random seed is not an integer.");
		}
		boolean removeFeatures = Boolean.valueOf(args[10]).booleanValue();
		boolean smoothPValues = Boolean.valueOf(args[11]).booleanValue();

		double minStd = -1;
		String confoundingClsFile = null;
		double theta = -1;
		boolean performAllPairs = true;
		double gamma = -1;

		for (int i = 12; i < args.length; i++) {
			String arg = args[i].substring(0, 2);
			String value = args[i].substring(2, args[i].length());
			if (value.equals("")) {
				continue;
			}
			if (arg.equals("-m")) {
				try {
					minStd = Double.parseDouble(value);
					if (minStd <= 0) {
						AnalysisUtil
								.exit("Minimum standard deviation must be > 0.");
					}
				} catch (NumberFormatException nfe) {
					AnalysisUtil
							.exit("Minimum standard deviation is not a number.");
				}
			} else if (arg.equals("-c")) {
				confoundingClsFile = value;
			} else if (arg.equals("-t")) {
				try {
					theta = Double.parseDouble(value);
					if (theta < 0) {
						AnalysisUtil.exit("Theta must be >= 0.");
					}
				} catch (NumberFormatException nfe) {
					AnalysisUtil.exit("Theta is not a number.");
				}
			} else if (arg.equals("-p")) {
				performAllPairs = value.equalsIgnoreCase("all pairs");
			} else if (arg.equals("-z")) {
				MarkerSelection.testGamma = true; // don't compute fwer, etc
			} else if (arg.equals("-g")) {
				try {
					gamma = Double.parseDouble(value);
					if (gamma < 0) {
						AnalysisUtil.exit("Gamma must be >= 0.");
					}
				} catch (NumberFormatException nfe) {
					AnalysisUtil.exit("Gamma is not a number.");
				}
			}
		}
		ClassVector classVector = AnalysisUtil.readClassVector(clsFile);
		IExpressionDataReader reader = AnalysisUtil
				.getExpressionReader(datasetFile);
		ExpressionData expressionData = AnalysisUtil.readExpressionData(reader,
				datasetFile);

		ClassVector confoundingClassVector = null;
		if (confoundingClsFile != null) {
			confoundingClassVector = AnalysisUtil
					.readClassVector(confoundingClsFile);
			/*
			 * ClassVector[] cv = new
			 * edu.mit.broad.internal.MultiClassReader().read(clsFile,
			 * expressionData); this.classVector = cv[0]; if(cv.length > 1) {
			 * covariate = cv[1]; for(int i = 2; i < cv.length; i++) { covariate =
			 * covariate.union(cv[i]); } }
			 */
		}

		if (classVector.getClassCount() > 2) {
			if (performAllPairs) {
				ClassVector[] allPairs = classVector.getAllPairs();

				String baseOutputFileName = Util
						.getBaseFileName(outputFileName);
				int n = classVector.getClassCount();
				int allPairsIndex = 0;
				for (int i = 0; i < n; i++) {
					for (int j = i + 1; j < n; j++) {
						int[] iIndices = classVector.getIndices(i);
						int[] jIndices = classVector.getIndices(j);
						int[] indices = new int[iIndices.length
								+ jIndices.length];
						for (int k = 0; k < iIndices.length; k++) {
							indices[k] = iIndices[k];
						}
						for (int k = 0; k < jIndices.length; k++) {
							indices[k + iIndices.length] = jIndices[k];
						}
						ExpressionData pairwiseDataset = expressionData.slice(
								null, indices);
						ClassVector pairwiseConfounder = null;
						if (confoundingClassVector != null) {
							pairwiseConfounder = confoundingClassVector
									.slice(indices);
						}
						String tempOutputFileName = baseOutputFileName + "."
								+ classVector.getClassName(i) + ".vs."
								+ classVector.getClassName(j) + ".odf";
						new MarkerSelection(pairwiseDataset, datasetFile,
								allPairs[allPairsIndex++], clsFile,
								_numPermutations, testDirection,
								tempOutputFileName, balanced, complete, metric,
								minStd, seed, pairwiseConfounder,
								confoundingClsFile, removeFeatures, theta,
								smoothPValues, gamma);

					}
				}

			} else {
				ClassVector[] oneVersusAll = classVector.getOneVersusAll();
				String baseOutputFileName = Util
						.getBaseFileName(outputFileName);
				for (int i = 0; i < classVector.getClassCount(); i++) {
					String tempOutputFileName = baseOutputFileName + "."
							+ classVector.getClassName(i) + ".vs.Rest.odf";
					new MarkerSelection(expressionData, datasetFile,
							oneVersusAll[i], clsFile, _numPermutations,
							testDirection, tempOutputFileName, balanced,
							complete, metric, minStd, seed,
							confoundingClassVector, confoundingClsFile,
							removeFeatures, theta, smoothPValues, gamma);
				}
			}

		} else {
			new MarkerSelection(expressionData, datasetFile, classVector,
					clsFile, _numPermutations, testDirection, outputFileName,
					balanced, complete, metric, minStd, seed,
					confoundingClassVector, confoundingClsFile, removeFeatures,
					theta, smoothPValues, gamma);
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
		lowerBound = new double[numFeatures];
		upperBound = new double[numFeatures];

		statisticalMeasure.compute(dataArray, classZeroIndices,
				classOneIndices, scores, null, -1);

		// index of observed scores, greatest to smallest
		int[] descendingIndices = Sorting.index(scores, Sorting.DESCENDING);

		if (numPermutations <= 0) {
			// calculate asymptotic p-values
			asymptotic = true;
			statisticalMeasure = new TestStatistics.TTest();
			// if (!(statisticalMeasure instanceof TestStatistics.TTest)
			// && !(statisticalMeasure instanceof TestStatistics.TTestMedian)
			// && !(statisticalMeasure instanceof TestStatistics.TTestMinStd)
			// && !(statisticalMeasure instanceof
			// TestStatistics.TTestMedianMinStd)) {
			// AnalysisUtil
			// .exit("Selected test statistic not supported when calculating
			// asymptotic p-values.");

			// }
		} else if (covariate != null) {
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
						+ Integer.MAX_VALUE + ".");
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
						+ Integer.MAX_VALUE + ".");
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

		if (!significanceBooster) {
			monotonicPermutedScores = new double[numFeatures];
			maxT = new double[numFeatures];
			fwer = new double[numFeatures];
		} else {
			if (estimateGamma) {

				long start = System.currentTimeMillis();
				for (int i = 0; i < 5; i++) {
					gamma = estimateGamma2(); // throw away
				}

				int tries = 50;
				double[] gammas = new double[tries];
				for (int i = 0; i < tries; i++) {
					gammas[i] = estimateGamma2();
				}
				// System.out.println(Util.toString(gammas));
				Arrays.sort(gammas);
				gamma = (gammas[gammas.length / 2] + gammas[gammas.length / 2 - 1]) / 2.0;
				long end = System.currentTimeMillis();
				double t = (end - start) / 1000.0;
				// System.out.println("Time to compute gamma=" + t);
			}
		}
		if (asymptotic) {
			TTestImpl test = new TTestImpl();
			double[] classZeroValues = new double[classZeroIndices.length];
			double[] classOneValues = new double[classOneIndices.length];

			for (int i = 0; i < numFeatures; i++) {
				double[] row = this.dataArray[i];
				for (int j = 0, length = classZeroIndices.length; j < length; j++) {
					classZeroValues[j] = row[classZeroIndices[j]];
				}
				for (int j = 0, length = classOneIndices.length; j < length; j++) {
					classOneValues[j] = row[classOneIndices[j]];
				}
				try {
					featureSpecificPValues[i] = test.tTest(classZeroValues,
							classOneValues);
					if (testDirection != Constants.TWO_SIDED) {
						featureSpecificPValues[i] /= 2;
					}
				} catch (IllegalArgumentException e) {
					AnalysisUtil
							.exit("An error occurred while computing the p-value for "
									+ expressionData.getRowName(i));
				} catch (MathException e) {
					AnalysisUtil
							.exit("An error occurred while computing the p-value for "
									+ expressionData.getRowName(i));
				}
			}
		} else {
			permute();
		}

		// free temporary storage
		monotonicPermutedScores = null;
		absoluteScores = null;
		permutedScores = null;

		double[] foldChange = new double[numFeatures];
		for (int i = 0; i < numFeatures; i++) {
			double meanClassZero = 0;
			double[] row = dataArray[i];
			for (int j = 0, length = classZeroIndices.length; j < length; j++) {
				meanClassZero += row[classZeroIndices[j]];
			}
			meanClassZero = meanClassZero / classZeroIndices.length;

			double meanClassOne = 0;
			for (int j = 0, length = classOneIndices.length; j < length; j++) {
				meanClassOne += row[classOneIndices[j]];
			}
			meanClassOne = meanClassOne / classOneIndices.length;
			if (meanClassOne == 0) {
				meanClassOne = 0.01;
			}
			foldChange[i] = meanClassZero / meanClassOne;
		}
		if (!asymptotic) {
			for (int i = 0; i < numFeatures; i++) {
				// fpr[i] /= numFeatures;
				// fpr[i] /= numPermutations;
				int N = permutationsPerFeature[i];
				double k = featureSpecificPValues[i];
				double p;
				if (smoothPValues) {
					p = (k + 1) / (N + 2);
				} else {
					p = k / N;
				}
				if (testDirection == Constants.TWO_SIDED) {
					double oneMinusP = 1.0 - p;
					if (oneMinusP < p) {
						p = oneMinusP;
					}
					p *= 2;

					if (p == 0) {
						// ensure not degenerate case where profile is
						// completely
						// flat
						// TODO handle cases where profile is flat (but not
						// completely)
						double[] row = dataArray[i];
						double val = row[0];
						boolean flat = true;
						for (int j = 1, cols = row.length; j < cols && flat; j++) {
							if (row[j] != val) {
								flat = false;
							}
						}
						if (flat) {
							p = 1;
						}
					}
				}
				featureSpecificPValues[i] = p;

				if (testDirection == Constants.TWO_SIDED) {
					k = 2 * Math.min(k, N - k);
				}

				if (complete) {
					lowerBound[i] = p;
					upperBound[i] = p;
				} else {
					double shape1 = k + 1;
					double shape2 = N - k + 1;
					final JSci.maths.statistics.BetaDistribution betaDist = new JSci.maths.statistics.BetaDistribution(
							shape1, shape2);
					double plow = betaDist.inverse(0.025);
					double phigh = betaDist.inverse(0.975);

					Function function = new Function() {
						public double evaluate(double x) {
							double d = Math.min(1, 0.95 + betaDist
									.cumulative(x));
							return betaDist.probability(x)
									- betaDist.probability(betaDist.inverse(d));
						}
					};

					if (k == 0) {
						plow = 0;
						phigh = betaDist.inverse(0.05);
					} else if (k == N) {
						plow = betaDist.inverse(0.95);
						phigh = 1;
					} else if (k < (N / 2)) {
						double accuracy = Math.min(p / 1000, 0.000001);
						// search between 0 and plow
						plow = ZeroFinder
								.bisection(0, plow, accuracy, function);
						phigh = betaDist.inverse(0.95 + betaDist
								.cumulative(plow));
					} else if (k > (N / 2)) {
						double accuracy = Math.min(p / 1000, 0.000001);
						// search between plow and beta.inverse(0.05)
						plow = ZeroFinder.bisection(plow, betaDist
								.inverse(0.05), accuracy, function);
						phigh = betaDist.inverse(0.95 + betaDist
								.cumulative(plow));
					} else {
						plow = betaDist.inverse(0.025);
						phigh = betaDist.inverse(0.975);
					}
					lowerBound[i] = plow;
					upperBound[i] = phigh;
				}

				if (!significanceBooster) {
					fwer[i] /= N;
					// rankBasedPValues[i] /= N;
					maxT[i] /= N;
				}
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

		if (!significanceBooster) {
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
		if (!printStackTraces) {
			new File(tempFileName).deleteOnExit();
		}
		try {
			tempFileWriter = new PrintWriter(new FileWriter(tempFileName));
			for (int i = 0; i < numFeatures; i++) {
				tempFileWriter
						.println(featureSpecificPValues[descendingIndices[i]]);
			}
		} catch (IOException ioe) {
			if (printStackTraces) {
				ioe.printStackTrace();
			}
			AnalysisUtil
					.exit("An error occurred while saving the output file.");
		} finally {
			if (tempFileWriter != null) {
				tempFileWriter.close();
			}
		}

		boolean qvalueSuccess = QValue.qvalue(tempFileName,
				printStackTraces);
		String[] qvalues = new String[numFeatures];
		int qvalueHeaderLines = 4;
		String[] qvalueHeaders = new String[qvalueHeaderLines];
		BufferedReader br = null;
		File qvalueOutputFile = new File("qvalues.txt"); // produced by R code
		// above
		
		
		if (qvalueOutputFile.exists()) {
			qvalueOutputFile.deleteOnExit();
			try {// read in results from qvalue
				br = new BufferedReader(new FileReader(qvalueOutputFile));
				for (int i = 0; i < qvalueHeaderLines; i++) {
					qvalueHeaders[i] = br.readLine();
				}
				for (int i = 0; i < numFeatures; i++) {
					qvalues[i] = br.readLine();
				}
			} catch (IOException ioe) {
				if (printStackTraces) {
					ioe.printStackTrace();
				}
				qvalueHeaders = new String[] {};
				Arrays.fill(qvalues, "NaN");
				System.err
						.println("An error occurred while computing the q-value-continuing anyway");
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (IOException x) {
					}
				}
			}
		} else {
			qvalueHeaders = new String[] {};
			Arrays.fill(qvalues, "NaN");
		}

		OdfWriter pw = null;
		try {
			String[] columnNames = null;
			String[] columnTypes = null;
			String[] rowDescriptions = expressionData.getRowDescriptions();

			if (!significanceBooster) {
				if (asymptotic) {
					if (rowDescriptions != null) {
						columnNames = new String[] { "Rank", "Feature",
								"Description", "Score", "Feature P", "FDR(BH)",
								"Q Value", "Bonferroni", "Fold Change" };
						columnTypes = new String[] { "int", "String", "String",
								"double", "double", "double", "double",
								"double", "double" };
					} else {
						columnNames = new String[] { "Rank", "Feature",
								"Score", "Feature P", "FDR(BH)", "Q Value",
								"Bonferroni", "Fold Change" };
						columnTypes = new String[] { "int", "String", "double",
								"double", "double", "double", "double",
								"double" };
					}
				} else if (complete) {
					if (rowDescriptions != null) {
						columnNames = new String[] { "Rank", "Feature",
								"Description", "Score", "Feature P", "FDR(BH)",
								"Q Value", "Bonferroni", "maxT", "FWER",
								"Fold Change" };
						columnTypes = new String[] { "int", "String", "String",
								"double", "double", "double", "double",
								"double", "double", "double", "double" };
					} else {
						columnNames = new String[] { "Rank", "Feature",
								"Score", "Feature P", "FDR(BH)", "Q Value",
								"Bonferroni", "maxT", "FWER", "Fold Change" };
						columnTypes = new String[] { "int", "String", "double",
								"double", "double", "double", "double",
								"double", "double", "double" };
					}
				} else {
					if (rowDescriptions != null) {
						columnNames = new String[] { "Rank", "Feature",
								"Description", "Score", "Feature P",
								"Feature P Low", "Feature P High", "FDR(BH)",
								"Q Value", "Bonferroni", "maxT", "FWER",
								"Fold Change" };
						columnTypes = new String[] { "int", "String", "String",
								"double", "double", "double", "double",
								"double", "double", "double", "double",
								"double", "double" };
					} else {
						columnNames = new String[] { "Rank", "Feature",
								"Score", "Feature P", "Feature P Low",
								"Feature P High", "FDR(BH)", "Q Value",
								"Bonferroni", "maxT", "FWER", "Fold Change" };
						columnTypes = new String[] { "int", "String", "double",
								"double", "double", "double", "double",
								"double", "double", "double", "double",
								"double" };
					}
				}

			} else {
				if (rowDescriptions != null) {
					columnNames = new String[] { "Rank", "Feature",
							"Description", "Score", "Feature P",
							"Feature P Low", "Feature P High", "FDR(BH)",
							"Q Value", "Bonferroni", "Fold Change",
							"Permutations", "Active" };
					columnTypes = new String[] { "int", "String", "String",
							"double", "double", "double", "double", "double",
							"double", "double", "double", "int", "boolean" };
				} else {
					columnNames = new String[] { "Rank", "Feature", "Score",
							"Feature P", "Feature P Low", "Feature P High",
							"FDR(BH)", "Q Value", "Bonferroni", "Fold Change",
							"Permutations", "Active" };
					columnTypes = new String[] { "int", "String", "double",
							"double", "double", "double", "double", "double",
							"double", "double", "int", "boolean" };
				}
			}
			pw = new OdfWriter(outputFileName, columnNames,
					"Comparative Marker Selection", numFeatures, true);
			pw.setColumnTypes(columnTypes);
			pw.addHeader("Dataset File", AnalysisUtil.getFileName(datasetFile));
			pw.addHeader("Class File", AnalysisUtil.getFileName(clsFile));
			if (covariate != null) {
				pw.addHeader("Confounding Class File", AnalysisUtil
						.getFileName(confoundingClsFile));
			}

			if (!asymptotic) {
				pw.addHeader("Permutations", numPermutations);
			}

			if (!asymptotic) {
				pw.addHeader("Balanced", String.valueOf(balanced));
				pw.addHeader("Complete", String.valueOf(complete));
			}
			if (testDirection == Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
				pw.addHeader("Test Direction", "Class 0");
			} else if (testDirection == Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE) {
				pw.addHeader("Test Direction", "Class 1");
			} else {
				pw.addHeader("Test Direction", "2 Sided");
			}
			pw.addHeader("Class 0", classVector.getClassName(0));
			pw.addHeader("Class 1", classVector.getClassName(1));
			pw.addHeader("Test Statistic", statisticalMeasure.toString());

			for (int i = 0; i < qvalueHeaders.length; i++) {
				String[] tokens = qvalueHeaders[i].split("=");
				pw.addHeader(tokens[0], tokens[1]);
			}
			if (seedUsed) {
				pw.addHeader("Random Seed", seed);
			}

			if (!asymptotic) {
				pw.addHeader("Significance Booster", String
						.valueOf(significanceBooster));
			}
			boolean[] active = null;
			if (significanceBooster) {
				pw.addHeader("Theta", theta);
				pw.addHeader("Gamma", gamma);
				active = new boolean[numFeatures];
				for (int i = 0; i < lengthOfIndicesToPermute; i++) {
					active[featureIndicesToPermute[i]] = true;
				}
			}
			if (!asymptotic) {
				pw.addHeader("Smooth p-values", String.valueOf(smoothPValues));
			}
			pw.printHeader();

			for (int i = 0; i < numFeatures; i++) {
				int index = descendingIndices[i];

				int rank = _ranks[index];
				if (testDirection == Constants.CLASS_ZERO_LESS_THAN_CLASS_ONE) {
					rank = numFeatures - rank + 1;
				}
				double bonferroni = Math.min(featureSpecificPValues[index]
						* numFeatures, 1.0);
				pw.print(rank);
				pw.print("\t");
				pw.print(expressionData.getRowName(index));
				pw.print("\t");
				if (rowDescriptions != null) {
					String s = rowDescriptions[index];
					if (s == null) {
						s = "";
					}
					pw.print(s);
					pw.print("\t");
				}
				pw.print(scores[index]);
				pw.print("\t");
				pw.print(featureSpecificPValues[index]);
				if (!complete && !asymptotic) {
					pw.print("\t");
					pw.print(lowerBound[index]);
					pw.print("\t");
					pw.print(upperBound[index]);
				}
				pw.print("\t");
				pw.print(fdr[index]);
				pw.print("\t");
				pw.print(qvalues[i]);
				pw.print("\t");
				pw.print(bonferroni);
				if (!significanceBooster && !asymptotic) {
					pw.print("\t");
					pw.print(maxT[index]);
					pw.print("\t");
					pw.print(fwer[index]);
				}

				pw.print("\t");
				pw.print(foldChange[index]);

				if (significanceBooster) {
					pw.print("\t");
					pw.print(permutationsPerFeature[index]);
					pw.print("\t");
					pw.print(active[index] ? "true" : "false");
				}
				pw.println();
			}

		} catch (Exception e) {
			if (printStackTraces) {
				e.printStackTrace();
			}
			AnalysisUtil
					.exit("An error occurred while saving the output file.");
		} finally {
			if (pw != null) {
				try {
					pw.close();
				} catch (Exception x) {
				}
			}
		}

		if (printPermutations) {
			debugger.print();
		}
	}

	private final double estimateGamma2() {
		Permuter permuter = null;
		int[] classZeroIndices = classVector.getIndices(0);
		int[] classOneIndices = classVector.getIndices(1);

		if (permuter instanceof UnbalancedRandomCovariatePermuter) {
			new UnbalancedRandomCovariatePermuter(classVector, covariate, seed);
		} else if (permuter instanceof BalancedRandomPermuter) {
			permuter = new BalancedRandomPermuter(classZeroIndices,
					classOneIndices, seed);

		} else {
			permuter = new UnbalancedRandomPermuter(classVector.size(),
					classVector.getIndices(1).length, seed);
		}
		long gammaStart = System.currentTimeMillis();

		int[] permutedClassZeroIndices = null;
		int[] permutedClassOneIndices = null;
		int tries = 1000;
		for (int test = 0; test < tries; test++) {
			int[] permutedAssignments = permuter.next();

			if (printPermutations) {
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
		}
		long gammaEnd = System.currentTimeMillis();
		double gammaElapsed = gammaEnd - gammaStart;

		double[] dummyArray = new double[1];
		int dummyInt = 0;
		double[] dummyRow = (double[]) dataArray[0].clone();
		double[][] dummyData = new double[1][];
		dummyData[0] = dummyRow;

		long permStart = System.currentTimeMillis();
		int[] featureIndicesToPermute = { 0 };

		for (int test = 0; test < tries; test++) {
			// calculate time for 1 feature
			statisticalMeasure.compute(dummyData, permutedClassZeroIndices,
					permutedClassOneIndices, permutedScores,
					featureIndicesToPermute, 1);
			int index = featureIndicesToPermute != null ? featureIndicesToPermute[0]
					: 0;
			double score = scores[index];
			if (testDirection == Constants.TWO_SIDED
					|| testDirection == Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
				if (permutedScores[index] >= score) {
					dummyArray[index] += 1.0;
				}
			} else {
				if (permutedScores[index] <= score) {
					dummyArray[index] += 1.0;
				}
			}

			if (significanceBooster) {
				dummyArray[index]++;
				double k = dummyArray[index];
				int N = (int) dummyArray[index];
				if (testDirection == Constants.TWO_SIDED) {
					k = Math.min(k, N - k);
				}
				double d = 2 * 1.96 * Math.sqrt((N + 1 - k)
						/ ((N + 3) * (k + 1)));
				if (d >= 0) { // theta include marker
					dummyArray[0] = index;
					dummyInt++;
				}
			}
		}
		long permEnd = System.currentTimeMillis();
		double permElapsed = permEnd - permStart;
		double permTimeForOne = permElapsed / tries;
		double gammaTimeForOne = gammaElapsed / tries;
		double gamma = gammaTimeForOne / permTimeForOne;
		return gamma;
	}

	private final double estimateGamma() {
		Permuter permuter = null;
		int[] classZeroIndices = classVector.getIndices(0);
		int[] classOneIndices = classVector.getIndices(1);

		if (permuter instanceof UnbalancedRandomCovariatePermuter) {
			new UnbalancedRandomCovariatePermuter(classVector, covariate, seed);
		} else if (permuter instanceof BalancedRandomPermuter) {
			permuter = new BalancedRandomPermuter(classZeroIndices,
					classOneIndices, seed);

		} else {
			permuter = new UnbalancedRandomPermuter(classVector.size(),
					classVector.getIndices(1).length, seed);
		}

		int[] numGenesArray = { 1, 2, 3, 4, 10, 20 };
		int tries = 1000; // num permuations to perform
		double[] times = new double[numGenesArray.length];
		for (int n = 0; n < numGenesArray.length; n++) {
			int numGenes = numGenesArray[n];

			double bank = Double.MAX_VALUE;
			double[][] tempDataArray = new double[numGenes][];
			for (int i = 0; i < numGenes; i++) {
				tempDataArray[i] = dataArray[i];
			}

			double[] permutedScores = new double[numGenes];
			double[] tempArray = new double[numGenes];
			int[] tempArray2 = new int[numGenes];
			int[] featureIndicesToPermute = new int[numGenes];
			for (int i = 0; i < numGenes; i++) {
				featureIndicesToPermute[i] = i;
			}
			int lengthOfIndicesToPermute = numGenes;

			long start = System.currentTimeMillis();

			for (int test = 0; test < tries; test++) {
				int[] permutedClassZeroIndices = null;
				int[] permutedClassOneIndices = null;
				int[] permutedAssignments = null;
				if (readPermutations) {
					permutedAssignments = nextPermutation();
				} else {
					permutedAssignments = permuter.next();
				}

				if (printPermutations) {
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

				statisticalMeasure.compute(tempDataArray,
						permutedClassZeroIndices, permutedClassOneIndices,
						permutedScores, featureIndicesToPermute,
						lengthOfIndicesToPermute);// compute
				// scores
				// using
				// permuted
				// class
				// labels

				int length = featureIndicesToPermute != null ? lengthOfIndicesToPermute
						: numFeatures;
				lengthOfIndicesToPermute = 0;
				if (length == 0) {
					bank = 0;
				}

				for (int i = 0; i < length && bank > 0; i++) {
					bank--;
					int index = featureIndicesToPermute != null ? featureIndicesToPermute[i]
							: i;
					double score = scores[index];
					if (testDirection == Constants.TWO_SIDED
							|| testDirection == Constants.CLASS_ZERO_GREATER_THAN_CLASS_ONE) {
						if (permutedScores[index] >= score) {
							tempArray[index] += 1.0;
						}
					} else {
						if (permutedScores[index] <= score) {
							tempArray[index] += 1.0;
						}
					}

					if (significanceBooster) {
						tempArray2[index]++;
						double k = tempArray[index];
						int N = tempArray2[index];
						if (testDirection == Constants.TWO_SIDED) {
							k = Math.min(k, N - k);
						}
						double d = 2 * 1.96 * Math.sqrt((N + 1 - k)
								/ ((N + 3) * (k + 1)));
						if (d >= 0) { // theta include marker
							featureIndicesToPermute[lengthOfIndicesToPermute] = index;
							lengthOfIndicesToPermute++;
						}
					}
				}

			}
			long end = System.currentTimeMillis();
			double elapsed = (end - start);
			elapsed /= tries;
			times[n] = elapsed;
		}

		SimpleRegression regression = new SimpleRegression();
		for (int i = 0; i < times.length; i++) {
			System.out.println(numGenesArray[i] + " " + times[i]);
			regression.addData(numGenesArray[i], times[i]);
		}
		System.out.println();
		double gamma = regression.getIntercept() / regression.getSlope();
		if (gamma < 0) {
			gamma = 0;
		}
		return gamma;
	}

	private final void permute() {
		BigDecimal maxLong = BigDecimal.valueOf(Long.MAX_VALUE);
		BigDecimal bigIntBank = new BigDecimal(numFeatures + gamma)
				.multiply(BigDecimal.valueOf(numPermutations));
		if (bigIntBank.compareTo(maxLong) > 0) {
			AnalysisUtil
					.exit("Arithmetic overflow. Try decreasing the number of permutations.");
		}
		long bank = bigIntBank.longValue();

		permutationsPerFeature = new int[numFeatures];
		if (significanceBooster) {
			featureIndicesToPermute = new int[numFeatures];
			for (int i = 0; i < numFeatures; i++) {
				featureIndicesToPermute[i] = i;
			}
		} else {
			Arrays.fill(permutationsPerFeature, numPermutations);
		}

		lengthOfIndicesToPermute = numFeatures;

		while (bank > 0) {
			int[] permutedClassZeroIndices = null;
			int[] permutedClassOneIndices = null;
			int[] permutedAssignments = null;
			if (readPermutations) {
				permutedAssignments = nextPermutation();
			} else {
				permutedAssignments = permuter.next();
			}

			if (printPermutations) {
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
					permutedClassOneIndices, permutedScores,
					featureIndicesToPermute, lengthOfIndicesToPermute);// compute
			// scores
			// using
			// permuted
			// class
			// labels

			int length = featureIndicesToPermute != null ? lengthOfIndicesToPermute
					: numFeatures;
			lengthOfIndicesToPermute = 0;
			if (length == 0) {
				bank = 0;
			}
			for (int i = 0; i < length && bank > 0; i++) {
				bank--;
				int index = featureIndicesToPermute != null ? featureIndicesToPermute[i]
						: i;
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

				if (significanceBooster) {
					permutationsPerFeature[index]++;
					double k = featureSpecificPValues[index];
					int N = permutationsPerFeature[index];
					if (testDirection == Constants.TWO_SIDED) {
						k = Math.min(k, N - k);
					}
					double d = 2 * 1.96 * Math.sqrt((N + 1 - k)
							/ ((N + 3) * (k + 1)));
					if (d >= theta) { // include marker
						featureIndicesToPermute[lengthOfIndicesToPermute] = index;
						lengthOfIndicesToPermute++;
					}
				}
			}
			bank -= gamma;

			/*
			 * Sorting.sort(permutedScores, Sorting.DESCENDING);
			 * 
			 * for(int i = 0; i < numFeatures; i++) { double score =
			 * scores[descendingIndices[i]];
			 * 
			 * if(testDirection == TWO_SIDED) { if(score >= 0) {
			 * if(permutedScores[i] >= score) { rankBasedPValues[i] += 1.0; } }
			 * else { if(permutedScores[i] <= score) { rankBasedPValues[i] +=
			 * 1.0; } } } else if(testDirection ==
			 * CLASS_ZERO_GREATER_THAN_CLASS_ONE) { if(permutedScores[i] >=
			 * score) { rankBasedPValues[i] += 1.0; } } else {
			 * if(permutedScores[i] <= score) { rankBasedPValues[i] += 1.0; } } }
			 */

			if (!testGamma && !significanceBooster) {
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
	 * Reads the next permutation from a file
	 * 
	 * @return the next permutation
	 */
	final int[] nextPermutation() {
		try {
			String s = testClassPermutationsReader.readLine();
			StringTokenizer st = new StringTokenizer(s, " \t");
			int cols = expressionData.getColumnCount();
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
	 * Opens the permutations file for reading
	 * 
	 * @param fileName
	 *            the file name
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
}
