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


package edu.mit.broad.marker.qvalue;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * @author Joshua Gould
 */
public class QValue {
	private final static String R_FILE = "qvalue.R";

	public static boolean qvalue(String markerInputFile, boolean printRErrors) {
		try {
			String[] args = { "'" + markerInputFile + "'" };
			edu.mit.broad.marker.r.R.run(new File(System.getProperty("libdir"),
					R_FILE).getCanonicalPath(), "compute.qvalue", args, System
					.getProperty("R"), printRErrors, printRErrors);
			return true;
		} catch (Exception ioe) {
			System.err
					.println("An error occurred while computing q value-continuing anyway");
			return false;
		}
	}
}
