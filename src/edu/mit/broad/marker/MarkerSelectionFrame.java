package edu.mit.broad.marker;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Toolkit;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableModel;
import java.util.StringTokenizer;
import java.util.Vector;
import java.io.BufferedReader;
import java.io.FileReader;
import gnu.trove.*;
import edu.mit.broad.modules.genelist.plot.*;
import edu.mit.broad.modules.genelist.gc.*;
import edu.mit.broad.modules.genelist.table.*;

/**
 *@author     Joshua Gould
 *@created    September 17, 2004
 */
public class MarkerSelectionFrame extends JFrame {
	GPTable table;
	GPPlot plot;
	PlotMenu plotMenu;
	GeneCruiserMenu geneCruiserMenu;
	GeneCruiserModel geneCruiserModel;
	Vector features;
	double[] geneSpecificPValues;
	double[] rankBasedPValues;
	double[] fwer;
	double[] fdr;
	double[] fpr;
	int N;
	final static Class[] COLUMN_CLASSES = {String.class, Double.class, Double.class, Double.class, Double.class, Double.class};
	final static String[] COLUMN_NAMES = {"Feature", "Feature Specific P Value", "Rank Based P Value", "FWER", "FPR", "FDR (BH)"};


	public MarkerSelectionFrame(Vector features, double[] rankBasedPValues, double[] geneSpecificPValues, double[] fwer, double[] fpr, double[] fdr) {
		this.features = features;
		this.geneSpecificPValues = geneSpecificPValues;
		this.rankBasedPValues = rankBasedPValues;
		this.fwer = fwer;
		this.fdr = fdr;
		this.fpr = fpr;
		plot = new MarkerSelectionPlot();
		plot.setMarksStyle("points", 0);
		plot.setTitle("Markers");
		plot.setXLabel("Feature Specific P Value");
		plot.setYLabel("Rank Based P Value");
		plot.addLegend(0, "P Values");
		N = features.size();
		double max = 0;
		for(int i = 0; i < N; i++) {
			plot.addPoint(0, geneSpecificPValues[i], rankBasedPValues[i], false);
			max = Math.max(max, geneSpecificPValues[i]);
			max = Math.max(max, rankBasedPValues[i]);
		}
		plot.setXRange(0, max);
		plot.setYRange(0, max);

		geneCruiserModel = new GeneCruiserModel(new MyTableModel(), features);
		geneCruiserMenu = new GeneCruiserMenu(geneCruiserModel, this);
		table = new GPTable(geneCruiserModel);

		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setTitle("Marker Selection");
		Container c = getContentPane();
		c.setLayout(new BorderLayout());
		JPanel tablePanel = new JPanel(new BorderLayout());
		/*
		    String[] summaryKeys={
		    "Select Method", "Distance Function", "Neighbors", "Class Estimate",
		    "User Significance Level", "Permutations", "Distinct Permutations"
		    };
		    String[] summaryValues={
		    parser.getMarkerSelectionMethod(), parser.getDistanceFunction(),
		    parser.getNumNeighbors(), parser.getClassEstimate(),
		    parser.getUserSigLevel(), parser.getNumPermutations(),
		    parser.getNumDistinctPermutations()
		    };
		    tablePanel.add(new SummaryPanel(summaryKeys, summaryValues),
		    BorderLayout.NORTH);
		  */
		JScrollPane sp = new JScrollPane(table);
		tablePanel.add(sp, BorderLayout.CENTER);

		final JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
				plot, tablePanel);
		splitPane.setDividerLocation(plot.getPreferredSize().height);
		c.add(splitPane, BorderLayout.CENTER);
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		plotMenu = new PlotMenu(plot);
		menuBar.add(plotMenu);
		JMenu viewMenu = new JMenu("View");
		menuBar.add(viewMenu);
		menuBar.add(geneCruiserMenu);
		Dimension size = Toolkit.getDefaultToolkit().getScreenSize();
		setSize((int) (size.width * .9), (int) (size.height * .9));
		show();
	}


	public static void main(String[] args) {
		String inputFile = args[0];
		Vector features = new Vector();
		TDoubleArrayList rankBasedPValues = new TDoubleArrayList();
		TDoubleArrayList geneSpecificPValues = new TDoubleArrayList();
		TDoubleArrayList fwer = new TDoubleArrayList();
		TDoubleArrayList fpr = new TDoubleArrayList();
		TDoubleArrayList fdr = new TDoubleArrayList();
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(inputFile));
			String s = null;
			while((s = br.readLine()) != null) {
				StringTokenizer st = new StringTokenizer(s, "\t");
				features.add(st.nextToken());
				rankBasedPValues.add(Double.parseDouble(st.nextToken()));
				geneSpecificPValues.add(Double.parseDouble(st.nextToken()));
				fwer.add(Double.parseDouble(st.nextToken()));
				fpr.add(Double.parseDouble(st.nextToken()));
				fdr.add(Double.parseDouble(st.nextToken()));
			}
		} catch(Exception e) {
			edu.mit.broad.gp.GPUtil.exit("An error occurred while reading the input file.", e);
		}
		if(br != null) {
			try {
				br.close();
			} catch(Exception x){}
		}
		new MarkerSelectionFrame(features, rankBasedPValues.toNativeArray(), geneSpecificPValues.toNativeArray(), fwer.toNativeArray(), fpr.toNativeArray(), fdr.toNativeArray());
	}


	/**
	 *@author    Joshua Gould
	 */
	static class MarkerSelectionPlot extends GPPlot {
		public MarkerSelectionPlot() {
			_padding = 0;
		}


		protected void _drawPlot(java.awt.Graphics graphics, boolean clearfirst,
				java.awt.Rectangle drawRectangle) {
			super._drawPlot(graphics, clearfirst, drawRectangle);
			graphics.setColor(Color.black);

			int ybottom = (int) yToPix(_yMin);
			int ytop = (int) yToPix(_yMax);
			int xleft = (int) xToPix(_xMin);
			int xright = (int) xToPix(_xMax);
			graphics.drawLine(xleft, ybottom, xright, ytop);
		}
	}


	/**
	 *@author     Joshua Gould
	 *@created    September 17, 2004
	 */
	class MyTableModel implements TableModel {

		public void addTableModelListener(TableModelListener l) {

		}


		public void removeTableModelListener(TableModelListener l) {

		}


		public void setValueAt(Object aValue, int rowIndex, int columnIndex) {

		}


		public int getRowCount() {
			return N;
		}


		public int getColumnCount() {
			return COLUMN_NAMES.length;
		}


		public String getColumnName(int columnIndex) {
			return COLUMN_NAMES[columnIndex];
		}


		public Class getColumnClass(int columnIndex) {
			return COLUMN_CLASSES[columnIndex];
		}


		public boolean isCellEditable(int rowIndex, int columnIndex) {
			return false;
		}


		public Object getValueAt(int rowIndex, int columnIndex) {
			if(columnIndex == 0) {
				return features.get(rowIndex);
			} else if(columnIndex == 1) {
				return new Double(geneSpecificPValues[rowIndex]);
			} else if(columnIndex == 2) {
				return new Double(rankBasedPValues[rowIndex]);
			} else if(columnIndex == 3) {
				return new Double(fwer[rowIndex]);
			} else if(columnIndex == 4) {
				return new Double(fpr[rowIndex]);
			} else {
				return new Double(fdr[rowIndex]);
			}
		}

	}

}

