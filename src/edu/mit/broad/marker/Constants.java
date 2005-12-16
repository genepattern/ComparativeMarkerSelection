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

public class Constants {
	
	private Constants(){}

	/** Test direction */	
	final static int CLASS_ZERO_GREATER_THAN_CLASS_ONE = 0;

	final static int CLASS_ZERO_LESS_THAN_CLASS_ONE = 1;

	final static int TWO_SIDED = 2;
	
	/** Test statistic */

	final static int T_TEST = 0;
	
	final static int T_TEST_MEDIAN = 2;
	
	final static int T_TEST_MIN_STD = 4;

	final static int T_TEST_MEDIAN_MIN_STD = 7;
	

	final static int SNR = 1;

	final static int SNR_MEDIAN = 3;
	
	final static int SNR_MIN_STD = 5;

	final static int SNR_MEDIAN_MIN_STD = 6;
	
}
