package edu.mit.broad.marker.qvalue;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
/**
 *@author    Joshua Gould
 */
public class QValue {
	private final static String R_FILE = "qvalue.R";


	public static void qvalue(String markerInputFile) {
		try {
			String[] args = {"'" + markerInputFile + "'"};
			edu.mit.broad.marker.r.R.run(new File(System.getProperty("libdir"), R_FILE).getCanonicalPath(), "compute.qvalue", args, System.getProperty("R"), false, false);
		} catch(Exception ioe) {
			System.err.println("Error computing q value");
		} 
	}
}
