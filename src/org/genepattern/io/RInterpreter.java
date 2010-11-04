package org.genepattern.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * RInterpreter is a simple interface that translates a command line into something that the R interpreter can handle.
 * The goal is to start R running, have it read in a script, and then begin execution of a particular method in the
 * script with specified arguments. R is not natively capable of doing this. So the RunR class invokes R, then feeds it
 * a series of commands via the standard input stream and copies the results to its own stdout and stderr output
 * streams.
 *
 * @author Jim Lerner
 * @author Joshua Gould
 */

public class RInterpreter {
    private Process process;
    private OutputStream stdin;
    private Thread outputReader;
    private Thread errorReader;

    /**
     * Invoke the R interpreter, create a few lines of input to feed to the stdin input stream of R, and spawn two
     * threads that copy stdout and stderr from R to this process' version of the same.
     *
     * @param args
     *            The command line arguments
     * @param pathToRHome
     *            Description of the Parameter
     * @author Jim Lerner
     * @throws Exception
     */
    public RInterpreter(String pathToRHome) throws Exception {
	this(createCommandLine(pathToRHome));
    }

    RInterpreter(String[] commandLine) throws Exception {

	String[] rFlags = new String[] { "--no-save", "--quiet", "--slave", "--no-restore" };
	String[] fullCommandLine = new String[commandLine.length + rFlags.length];
	System.arraycopy(commandLine, 0, fullCommandLine, 0, commandLine.length);

	System.arraycopy(rFlags, 0, fullCommandLine, commandLine.length, rFlags.length);
	process = Runtime.getRuntime().exec(fullCommandLine, null, null);
	// create threads to read from the command's stdout and stderr streams
	outputReader = streamCopier(process.getInputStream(), System.out);
	errorReader = streamCopier(process.getErrorStream(), System.err);

	// drain the output and error streams
	outputReader.start();
	errorReader.start();

	stdin = process.getOutputStream();

    }

    public void quit() throws Exception {
	sendCmd(stdin, "q(save=\"no\")\n");
	stdin.close();
	// wait for all output before attempting to send it back to the client
	outputReader.join();
	errorReader.join();

	// the process will be dead by now
	process.waitFor();
    }

    public void read(InputStream is) throws Exception {
	BufferedReader br = new BufferedReader(new InputStreamReader(is));
	String s = null;
	while ((s = br.readLine()) != null) {
	    sendCmd(stdin, s + "\n");
	}

    }

    public void runMethod(String methodName, String[] args) throws Exception {
	sendCmd(stdin, "output <- " + methodName);
	sendCmd(stdin, "(");

	for (int i = 0; i < args.length; i++) {
	    if (i > 0) {
		sendCmd(stdin, ", ");
	    }

	    sendCmd(stdin, "'" + fixPath(args[i]) + "'");
	}
	sendCmd(stdin, ")\n");
    }

    public void source(File file) throws Exception {
	BufferedReader br = new BufferedReader(new FileReader(file));
	String s = null;
	while ((s = br.readLine()) != null) {
	    sendCmd(stdin, s + "\n");
	}
	br.close();
    }

    public void source(String text) throws Exception {
	sendCmd(stdin, text);
    }

    /**
     * launch a new thread (an instance of this) that will copy stdout or stderr to a PrintStream
     *
     */
    static class StreamCopier extends Thread {
	InputStream is = null;
	PrintStream os = null;

	/**
	 * Constructor for the StreamCopier object
	 *
	 * @param is
	 *            Description of the Parameter
	 * @param os
	 *            Description of the Parameter
	 */
	public StreamCopier(InputStream is, PrintStream os) {
	    this.is = is;
	    this.os = os;
	    this.setDaemon(true);
	}

	/**
	 * Main processing method for the StreamCopier object
	 */

	@Override
	public void run() {
	    BufferedReader in = new BufferedReader(new InputStreamReader(is));
	    String line;

	    try {
		while ((line = in.readLine()) != null) {
		    os.print(line);
		    os.flush();
		    // show it to the user ASAP
		}

	    } catch (IOException ioe) {
		System.err.println(ioe + " while reading from process stream");
	    }
	}
    }

    public static void main(String[] args) throws Exception {
	String rHome = args[0];
	RInterpreter r = new RInterpreter(rHome);
	String rSourceFile = args[1];
	String method = args[2];
	r.source(new File(rSourceFile));
	String[] methodArgs = new String[args.length - 3];
	System.arraycopy(args, 3, methodArgs, 0, methodArgs.length);
	r.runMethod(method, methodArgs);
	r.quit();
    }

    /**
     * convert Windows path separators to Unix, which R prefers!
     *
     * @param path
     *            path to convert to Unix format
     * @return String path with delimiters replaced
     * @author Jim Lerner
     */
    protected static String fixPath(String path) {
	return path.replace('\\', '/');
    }

    /**
     * write a string to stdin of the R process
     *
     * @param command
     *            string to send to R
     * @param stdin
     *            Description of the Parameter
     * @exception IOException
     *                Description of the Exception
     * @author Jim Lerner
     */
    protected static void sendCmd(OutputStream stdin, String command) throws IOException {
	stdin.write(command.getBytes());
	stdin.flush();
    }

    /**
     * copies one of the output streams from R to this process' output stream
     *
     * @param is
     *            InputStream to read from (from R)
     * @param os
     *            PrintStream to write to (stdout of this process)
     * @return Description of the Return Value
     * @exception IOException
     *                Description of the Exception
     * @author Jim Lerner
     */
    protected static Thread streamCopier(InputStream is, PrintStream os) {
	return new StreamCopier(is, os);
    }

    private static String[] createCommandLine(String pathToRHome) {
	String[] commandLine = null;
	boolean runningOnWindows = System.getProperty("os.name").startsWith("Windows");
	if (runningOnWindows) {
	    if (pathToRHome == null) {
		// assume Rterm is in path
		commandLine = new String[] { "cmd", "/c", "Rterm" };
	    } else {
		commandLine = new String[] { "cmd", "/c", pathToRHome + "\\bin\\Rterm" };
	    }
	} else {
	    if (pathToRHome == null) {
		// assume R is in path
		commandLine = new String[] { "R" };
	    } else {
		commandLine = new String[] { pathToRHome + "/bin/R" };
	    }
	}
	return commandLine;
    }
}
