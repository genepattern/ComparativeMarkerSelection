package edu.mit.broad.marker;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.util.Vector;
import edu.mit.broad.modules.genelist.plot.*;
import edu.mit.broad.modules.genelist.gc.*;
import edu.mit.broad.modules.genelist.table.*;

/**
 * @author     Joshua Gould
 * @created    September 17, 2004
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
	int N;
	final static Class[] COLUMN_CLASSES = {String.class, Double.class, Double.class, Double.class, Double.class};
	final static String[] COLUMN_NAMES = {"Feature", "Feature Specific P Value", "Rank Based P Value", "FWER", "FDR (BH)"};


	public MarkerSelectionFrame(Vector features, double[] geneSpecificPValues, double[] rankBasedPValues, double[] fwer, double[] fdr) {
		this.features = features;
		this.geneSpecificPValues = geneSpecificPValues;
		this.rankBasedPValues = rankBasedPValues;
		this.fwer = fwer;
		this.fdr = fdr;

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
	 * @author     Joshua Gould
	 * @created    September 17, 2004
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
			} else {
				return new Double(fdr[rowIndex]);	
			}
		}

	}

}

