package edu.mit.broad.marker;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.*;
import javax.swing.JFrame;
import javax.swing.*;
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
import ptolemy.plot.*;

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
	double[] scores;
	/**  number of features to view */
	int numFeatures;
	/**  total number of data points */
	int N;
	final static Class[] COLUMN_CLASSES = {String.class, Double.class, Double.class, Double.class, Double.class, Double.class, Double.class};
	final static String[] COLUMN_NAMES = {"Feature", "Score", "Feature Specific P Value", "FPR", "FWER", "Rank Based P Value", "FDR (BH)"};

	HistPanel fprHistPanel, featureSpecificHistPanel;
	JSplitPane splitPane;


	public MarkerSelectionFrame(Vector features, double[] scores, double[] geneSpecificPValues, double[] fpr, double[] fwer, double[] rankBasedPValues, double[] fdr) {
		this.features = features;
		this.scores = scores;
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
		numFeatures = N;
	//	numFeatures = Math.min(N, 100);
		double max = 0;

		fprHistPanel = new HistPanel("FPR Histogram", "P Value", "Occurences", 0.05, fpr);

		featureSpecificHistPanel = new HistPanel("Feature Specific Histogram", "P Value", "Occurences", 0.05, geneSpecificPValues);

		geneCruiserModel = new GeneCruiserModel(new MyTableModel(), features);
		table = new GPTable(geneCruiserModel);
		final de.qfs.lib.gui.SortedTableHelper helper = new de.qfs.lib.gui.SortedTableHelper(table);
		
		helper.getTableModelSorter().setSortColumn(1);
		helper.getTableModelSorter().setSortAscending(false);
		helper.prepareTable();
		
		for(int i = 0; i < numFeatures; i++) {
			int sortedRow = helper.getSortedTableModel().getMappedRow(i);
			plot.addPoint(0, geneSpecificPValues[sortedRow], rankBasedPValues[sortedRow], false);
			max = Math.max(max, geneSpecificPValues[sortedRow]);
			max = Math.max(max, rankBasedPValues[sortedRow]);
		}

		plot.setXRange(0, max);
		plot.setYRange(0, max);

		
		geneCruiserMenu = new GeneCruiserMenu(geneCruiserModel, this);

		

		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setTitle("Marker Selection");
		Container c = getContentPane();
		c.setLayout(new BorderLayout());
		JPanel tablePanel = new JPanel(new BorderLayout());
		JScrollPane sp = new JScrollPane(table);
		tablePanel.add(sp, BorderLayout.CENTER);

		splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
				plot, tablePanel);
		splitPane.setDividerLocation(plot.getPreferredSize().height);
		c.add(splitPane, BorderLayout.CENTER);
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		plotMenu = new PlotMenu(plot);
		menuBar.add(plotMenu);
		JMenu viewMenu = new JMenu("View");
		final JMenuItem featureSpecific_vs_rankBasedMenuItem = new JMenuItem("Feature Specific vs Rank Based");
		final JMenuItem fprHistMenuItem = new JMenuItem("FPR Histogram");
		final JMenuItem fspHistMenuItem = new JMenuItem("Feature Specific Histogram");
		
		featureSpecific_vs_rankBasedMenuItem.addActionListener(
			new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					featureSpecific_vs_rankBasedMenuItem.setSelected(true);
					fspHistMenuItem.setSelected(false);
					fprHistMenuItem.setSelected(false);
					setPlot(plot);
				}
			});
		viewMenu.add(featureSpecific_vs_rankBasedMenuItem);

		
		fspHistMenuItem.addActionListener(
			new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					featureSpecific_vs_rankBasedMenuItem.setSelected(false);
					fspHistMenuItem.setSelected(true);
					fprHistMenuItem.setSelected(false);
					setPlot(featureSpecificHistPanel);
				}
			});
		viewMenu.add(fspHistMenuItem);

		
		fprHistMenuItem.addActionListener(
			new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					featureSpecific_vs_rankBasedMenuItem.setSelected(false);
					fspHistMenuItem.setSelected(false);
					fprHistMenuItem.setSelected(true);
					setPlot(fprHistPanel);

				}
			});
		viewMenu.add(fprHistMenuItem);
		menuBar.add(viewMenu);

		menuBar.add(geneCruiserMenu);
		Dimension size = Toolkit.getDefaultToolkit().getScreenSize();
		setSize((int) (size.width * .7), (int) (size.height * .7));
		show();
	}


	public static void main(String[] args) {
		String inputFile = args[0];
		Vector features = new Vector();
		TDoubleArrayList scores = new TDoubleArrayList();
		TDoubleArrayList geneSpecificPValues = new TDoubleArrayList();
		TDoubleArrayList fpr = new TDoubleArrayList();
		TDoubleArrayList fwer = new TDoubleArrayList();
		TDoubleArrayList rankBasedPValues = new TDoubleArrayList();
		TDoubleArrayList fdr = new TDoubleArrayList();

		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(inputFile));
			String s = null;
			while((s = br.readLine()) != null) {
				StringTokenizer st = new StringTokenizer(s, "\t");
				features.add(st.nextToken());
				scores.add(Double.parseDouble(st.nextToken()));
				geneSpecificPValues.add(Double.parseDouble(st.nextToken()));
				fpr.add(Double.parseDouble(st.nextToken()));
				fwer.add(Double.parseDouble(st.nextToken()));
				rankBasedPValues.add(Double.parseDouble(st.nextToken()));
				fdr.add(Double.parseDouble(st.nextToken()));
			}
		} catch(Exception e) {
			edu.mit.broad.gp.GPUtil.exit("An error occurred while reading the input file.", e);
		}
		if(br != null) {
			try {
				br.close();
			} catch(Exception x) {}
		}

		new MarkerSelectionFrame(features, scores.toNativeArray(), geneSpecificPValues.toNativeArray(), fpr.toNativeArray(), fwer.toNativeArray(), rankBasedPValues.toNativeArray(), fdr.toNativeArray());
	}


	void setPlot(JPanel p) {
		if(p instanceof HistPanel) {
			HistPanel panel = (HistPanel) p;
			plotMenu.setPlot(panel.hist);
			splitPane.setLeftComponent(panel);
			splitPane.setDividerLocation(panel.getPreferredSize().height);
		} else {
			PlotBox plot = (PlotBox) p;
			plotMenu.setPlot(plot);
			splitPane.setLeftComponent(plot);
			splitPane.setDividerLocation(plot.getPreferredSize().height);
		}
	}


	/**
	 *@author    Joshua Gould
	 */
	class HistPanel extends JPanel {
		GPHistogram hist = new GPHistogram();


		public HistPanel(String title, String xlabel, String ylabel, double width, final double[] data) {
			hist.setTitle(title);
			hist.setYLabel(ylabel);
			hist.setXLabel(xlabel);
			hist.setDiscrete(false);
			hist.setBinWidth(width);
			for(int i = 0; i < data.length; i++) {
				hist.append(data[i]);
			}
			setLayout(new BorderLayout());
			add(hist, BorderLayout.CENTER);
			JPanel temp2 = new JPanel();
			JLabel l = new JLabel("Bin Width:");
			temp2.add(l);
			final JTextField tf = new JTextField(20);
			tf.setText("" + width);
			tf.addActionListener(
				new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						try {
							double width = Double.parseDouble(tf.getText().trim());
							if(width <= 0 || width > 1) {
								JOptionPane.showMessageDialog(MarkerSelectionFrame.this, "Bin width must be between 0 and 1.");
								return;
							}
							HistPanel panel = new HistPanel(hist.getTitle(), hist.getXLabel(), hist.getYLabel(), width, data);
							setPlot(panel);
						} catch(NumberFormatException nfe) {
							JOptionPane.showMessageDialog(MarkerSelectionFrame.this, "Bin width is not a number");
						}
					}
				});

			temp2.add(tf);
			add(temp2, BorderLayout.SOUTH);
		}

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
			return numFeatures;
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
				return new Double(scores[rowIndex]);
			} else if(columnIndex == 2) {
				return new Double(geneSpecificPValues[rowIndex]);
			} else if(columnIndex == 3) {
				return new Double(rankBasedPValues[rowIndex]);
			} else if(columnIndex == 4) {
				return new Double(fwer[rowIndex]);
			} else if(columnIndex == 5) {
				return new Double(fpr[rowIndex]);
			} else {
				return new Double(fdr[rowIndex]);
			}
		}

	}

}

