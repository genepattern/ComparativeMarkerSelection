/*
 *    The Broad Institute
 *    SOFTWARE COPYRIGHT NOTICE AGREEMENT
 *    This software and its documentation are copyright (2003-2006) by the
 *    Broad Institute/Massachusetts Institute of Technology. All rights are
 *    reserved.
 *
 *    This software is supplied without any warranty or guaranteed support
 *    whatsoever. Neither the Broad Institute nor MIT can be responsible for its
 *    use, misuse, or functionality.
 */

    package org.genepattern.matrix;

    import gnu.trove.TObjectIntHashMap;

    import java.io.PrintStream;
    import java.util.ArrayList;
    import java.util.Arrays;
    import java.util.HashMap;
    import java.util.Iterator;
    import java.util.List;
    import java.util.Map;

    import Jama.Matrix;

    /**
     * Default implementation of Dataset interface.
     *
     * @author Joshua Gould
     */
    public class DefaultDataset implements Dataset {
        private Matrix matrix;

        private Map<String, double[][]> matrices = new HashMap<String, double[][]>();

        private MetaData rowMetaData;

        private MetaData columnMetaData;

        private TObjectIntHashMap rowNameToRowIndexMap;

        private TObjectIntHashMap columnNameToColumnIndexMap;

        private String[] rowNames;

        private String[] columnNames;

        private boolean allowDuplicateNames;

        private List<String> seriesNames = new ArrayList<String>();

        /**
         * Creates a new dataset
         *
         * @param data
         *            The data.
         */
        public DefaultDataset(double[][] data) {
        int rows = data.length;
        int columns = data[0].length;
        matrix = new Matrix(data);
        rowNameToRowIndexMap = new TObjectIntHashMap(rows);
        columnNameToColumnIndexMap = new TObjectIntHashMap(columns);
        this.rowNames = new String[rows];
        this.columnNames = new String[columns];
        fillInRows(0);
        fillInColumns(0);
        rowMetaData = new MetaData(rows);
        columnMetaData = new MetaData(columns);
        }

        /**
         * Creates a new matrix
         *
         * @param data
         *            The data.
         * @param rowNames
         *            the row names
         * @param columnNames
         *            the column names
         */
        public DefaultDataset(double[][] data, String[] rowNames, String[] columnNames) {
        this(data, rowNames, columnNames, (String[]) null, (String[]) null);
        }

        public DefaultDataset(double[][] data, String[] rowNames, String[] columnNames, MetaData _rowMetaData,
            MetaData _columnMetaData) {
        this(data, rowNames, columnNames, _rowMetaData, _columnMetaData, false);
        }

        public DefaultDataset(double[][] data, String[] rowNames, String[] columnNames, MetaData _rowMetaData,
            MetaData _columnMetaData, boolean allowDuplicateNames) {
        this.allowDuplicateNames = allowDuplicateNames;
        int rows = data.length;
        if (rows == 0) {
            throw new IllegalArgumentException("Number of rows must be greater than 0");
        }
        int columns = data[0].length;
        if (columns == 0) {
            throw new IllegalArgumentException("Number of columns must be greater than 0");
        }
        this.matrix = new Matrix(data);
        this.rowNameToRowIndexMap = new TObjectIntHashMap(rows);
        this.columnNameToColumnIndexMap = new TObjectIntHashMap(columns);
        if (columnNames.length != columns) {
            throw new IllegalArgumentException("Length of column names must be equal to number of columns in data.");
        }
        if (rowNames.length != rows) {
            throw new IllegalArgumentException("Length of row names must be equal to number of rows in data.");
        }
        this.rowNames = new String[rows];
        this.columnNames = new String[columns];
        setRowNames(Arrays.asList(rowNames));
        setColumnNames(Arrays.asList(columnNames));
        this.rowMetaData = _rowMetaData == null ? new MetaData(rows) : _rowMetaData;
        this.columnMetaData = _columnMetaData == null ? new MetaData(columns) : _columnMetaData;
        }

        /**
         * Creates a new matrix
         *
         * @param data
         *            The data.
         * @param rowNames
         *            the row names
         * @param columnNames
         *            the column names
         */
        public DefaultDataset(double[][] data, String[] rowNames, String[] columnNames, String[] rowDescriptions,
            String[] columnDescriptions) {
        int rows = data.length;
        if (rows == 0) {
            throw new IllegalArgumentException("Number of rows must be greater than 0");
        }
        int columns = data[0].length;
        if (columns == 0) {
            throw new IllegalArgumentException("Number of columns must be greater than 0");
        }
        this.matrix = new Matrix(data);
        this.rowNameToRowIndexMap = new TObjectIntHashMap(rows);
        this.columnNameToColumnIndexMap = new TObjectIntHashMap(columns);
        if (columnNames.length != columns) {
            throw new IllegalArgumentException("Length of column names (" + columnNames.length
                + ") must be equal to number of columns (" + columns + ") in data.");
        }
        if (rowNames.length != rows) {
            throw new IllegalArgumentException("Length of row names (" + rowNames.length + ") must be equal to number of rows ("
                + rows + ") in data.");
        }
        this.rowNames = new String[rows];
        this.columnNames = new String[columns];
        setRowNames(Arrays.asList(rowNames));
        setColumnNames(Arrays.asList(columnNames));
        rowMetaData = new MetaData(rows);
        columnMetaData = new MetaData(columns);
        if (rowDescriptions != null) {
            rowMetaData.setMetaData(DatasetConstants.DESCRIPTION_KEY, rowDescriptions);
        }
        if (columnDescriptions != null) {
            columnMetaData.setMetaData(DatasetConstants.DESCRIPTION_KEY, columnDescriptions);
        }
        }

        /**
         * Creates a new matrix
         *
         * @param rows
         *            The number of rows.
         * @param columns
         *            The number of columns.
         */
        public DefaultDataset(int rows, int columns) {
        if (rows <= 0) {
            throw new IllegalArgumentException("Number of rows must be greater than 0");
        }
        if (columns <= 0) {
            throw new IllegalArgumentException("Number of columns must be greater than 0");
        }
        matrix = new Matrix(rows, columns);
        rowNameToRowIndexMap = new TObjectIntHashMap(rows);
        columnNameToColumnIndexMap = new TObjectIntHashMap(columns);
        rowNames = new String[rows];
        columnNames = new String[columns];
        fillInRows(0);
        fillInColumns(0);
        rowMetaData = new MetaData(rows);
        columnMetaData = new MetaData(columns);
        }

        /** Creates a new uninitalized matrix with 0 rows and 0 columns */

        private DefaultDataset() {
        }

        public void addMatrix(String name, double[][] m) {
        matrices.put(name, m);
        seriesNames.add(name);
        }

        /**
         * Make a deep copy of this matrix
         *
         * @return the copy.
         */
        public DefaultDataset copy() {
        Matrix matrixCopy = matrix.copy();
        DefaultDataset doubleMatrixCopy = new DefaultDataset();
        doubleMatrixCopy.matrix = matrixCopy;
        doubleMatrixCopy.rowNames = this.rowNames.clone();
        doubleMatrixCopy.columnNames = this.columnNames.clone();
        doubleMatrixCopy.columnNameToColumnIndexMap = (TObjectIntHashMap) this.columnNameToColumnIndexMap.clone();
        doubleMatrixCopy.rowNameToRowIndexMap = (TObjectIntHashMap) this.rowNameToRowIndexMap.clone();
        return doubleMatrixCopy;
        }

        /**
         * DoubleMatrix2D determinant
         *
         * @return The determinant.
         */
        public double det() {
        return this.matrix.det();
        }

        /**
         * Gets the underlying double[][] array
         *
         * @return the array
         */
        public double[][] getArray() {
        return matrix.getArray();
        }

        public int getColumnCount() {
        return matrix.getColumnDimension();
        }

        public String getColumnDescription(int columnIndex) {
        return columnMetaData.getMetaData(DatasetConstants.DESCRIPTION_KEY)[columnIndex];
        }

        public String[] getColumnDescriptions() {
        return columnMetaData.getArray(DatasetConstants.DESCRIPTION_KEY);
        }

        /**
         * Gets the column index for the column name .
         *
         * @param columnName
         *            the column name.
         * @return the column index, or -1 if the column name is not contained in this matrix
         */
        public int getColumnIndex(String columnName) {
        if (columnNameToColumnIndexMap.containsKey(columnName)) {
            return columnNameToColumnIndexMap.get(columnName);
        }
        return -1;
        }

        public String getColumnMetadata(int columnIndex, String name) {
        return columnMetaData.contains(name) ? columnMetaData.getArray(name)[columnIndex] : null;
        }

        public int getColumnMetadataCount() {
        return columnMetaData.size();
        }

        public String getColumnMetadataName(int index) {
        return columnMetaData.getNames()[index];
        }

        public String getColumnName(int columnIndex) {
        return columnNames[columnIndex];
        }

        /**
         * Allocates a new array containing the column names.
         *
         * @return the column names.
         */
        public String[] getColumnNames() {
        return columnNames.clone();
        }

        /**
         * Gets the underlying array at the given row
         *
         * @param row
         *            Row index
         * @return the row array
         */
        public double[] getRow(int row) {
        return matrix.getArray()[row];
        }

        public int getRowCount() {
        return matrix.getRowDimension();
        }

        public String getRowDescription(int rowIndex) {
        return rowMetaData.getMetaData(DatasetConstants.DESCRIPTION_KEY)[rowIndex];
        }

        public String[] getRowDescriptions() {
        return rowMetaData.getArray(DatasetConstants.DESCRIPTION_KEY);
        }

        /**
         * Gets the row index for the row name .
         *
         * @param rowName
         *            the row name.
         * @return the row index, or -1 if the row name is not contained in this matrix
         */
        public int getRowIndex(String rowName) {
        if (rowNameToRowIndexMap.containsKey(rowName)) {
            return rowNameToRowIndexMap.get(rowName);
        }
        return -1;
        }

        public String getRowMetadata(int rowIndex, String name) {
        return rowMetaData.contains(name) ? rowMetaData.getArray(name)[rowIndex] : null;
        }

        public int getRowMetadataCount() {
        return rowMetaData.size();
        }

        public String getRowMetadataName(int index) {
        return rowMetaData.getNames()[index];
        }

        public String getRowName(int rowIndex) {
        return rowNames[rowIndex];
        }

        /**
         * Allocates a new array contains the row names
         *
         * @return The row names.
         */
        public String[] getRowNames() {
        return rowNames.clone();
        }

        public int getSeriesCount() {
        return matrices.size();
        }

        public String getSeriesName(int series) {
        return seriesNames.get(series);
        }

        /**
         * Gets a single element
         *
         * @param row
         *            Row index.
         * @param column
         *            Column index.
         * @return The value at A[row,column]
         */
        public double getValue(int row, int column) {
        return matrix.get(row, column);
        }

        public double getValue(int rowIndex, int columnIndex, String name) {
        double[][] m = matrices.get(name);
        return m != null ? m[rowIndex][columnIndex] : null;
        }

        public boolean isAllowDuplicateNames() {
        return allowDuplicateNames;
        }

        /**
         * C = A - B
         *
         * @param B
         *            another matrix
         * @return A - B
         */
        public DefaultDataset minus(DefaultDataset B) {
        DefaultDataset C = new DefaultDataset();
        C.matrix = matrix.minus(B.matrix);
        C.rowNames = this.rowNames.clone();
        C.columnNames = this.columnNames.clone();
        C.columnNameToColumnIndexMap = (TObjectIntHashMap) this.columnNameToColumnIndexMap.clone();
        C.rowNameToRowIndexMap = (TObjectIntHashMap) this.rowNameToRowIndexMap.clone();
        return C;
        }

        /**
         * C = A + B
         *
         * @param B
         *            another matrix
         * @return A + B
         */
        public DefaultDataset plus(DefaultDataset B) {
        DefaultDataset C = new DefaultDataset();
        C.matrix = matrix.plus(B.matrix);
        C.rowNames = this.rowNames.clone();
        C.columnNames = this.columnNames.clone();
        C.columnNameToColumnIndexMap = (TObjectIntHashMap) this.columnNameToColumnIndexMap.clone();
        C.rowNameToRowIndexMap = (TObjectIntHashMap) this.rowNameToRowIndexMap.clone();
        return C;
        }

        /**
         * Prints this matrix in delimitted format using the default number format.
         *
         * @param ps
         *            the print stream
         * @param delimiter
         *            the delimiter between columns
         * @see #print(java.io.PrintStream,String,java.text.NumberFormat)
         */
        public void print(java.io.PrintStream ps, String delimiter) {
        print(ps, delimiter, java.text.NumberFormat.getInstance());
        }

        /**
         * Prints this matrix in delimitted format.
         *
         * @param ps
         *            the print stream
         * @param delimiter
         *            the delimiter between columns
         * @param nf
         *            the formatter
         */
        public void print(java.io.PrintStream ps, String delimiter, java.text.NumberFormat nf) {
        for (int j = 0, columns = getColumnCount(); j < columns; j++) {
            ps.print(delimiter);
            ps.print(columnNames[j]);
        }
        int columns = getColumnCount();
        for (int i = 0, rows = getRowCount(); i < rows; i++) {
            ps.println();
            ps.print(rowNames[i]);
            ps.print(delimiter);
            for (int j = 0; j < columns - 1; j++) {
            ps.print(nf.format(matrix.get(i, j)));
            ps.print(delimiter);
            }
            ps.print(nf.format(matrix.get(i, columns - 1)));// don't print the
            // delimmiter after
            // the last column
        }
        }

        /**
         * Prints the data in this dataset to the given stream.
         *
         * @param ps
         *            The stream.
         */
        public void print(PrintStream ps) {
        for (int i = 0, rows = getRowCount(); i < rows; i++) {
            ps.print(getRowName(i));
            ps.print(", ");
            for (int j = 0, columns = getColumnCount(); j < columns; j++) {
            if (j > 0) {
                ps.print(", ");
            }
            ps.print(getValue(i, j));
            }
            ps.print("\n");
        }
        }

        /**
         * Returns the effective numerical rank, obtained from Singular Value Decomposition.
         *
         * @return The rank.
         */
        public int rank() {
        return this.matrix.rank();
        }

        public void setAllowDuplicateNames(boolean allowDuplicateNames) {
        this.allowDuplicateNames = allowDuplicateNames;
        }

        /**
         * Sets the column descriptions
         *
         * @param descs
         * @throws IllegalArgumentException
         *             if desc.length != getColumnDimension()
         */
        public void setColumnDescriptions(String[] descs) {
        if (descs != null && descs.length != getColumnCount()) {
            throw new IllegalArgumentException("Length of descriptions must be equal to the number of columns.");
        }
        columnMetaData.setMetaData(DatasetConstants.DESCRIPTION_KEY, descs);
        }

        /**
         * Sets the column name at the specified index
         *
         * @param columnIndex
         *            The column index.
         * @param name
         *            The new column name value
         */
        public void setColumnName(int columnIndex, String name) {
        if (name == null) {
            throw new IllegalArgumentException("Null column names are not allowed.");
        }
        if (!allowDuplicateNames && columnNameToColumnIndexMap.containsKey(name)
            && columnNameToColumnIndexMap.get(name) != columnIndex) {
            throw new IllegalArgumentException("Duplicate column names are not allowed:" + name);
        }
        columnNameToColumnIndexMap.put(name, columnIndex);
        columnNames[columnIndex] = name;
        }

        public void setRowDescription(int rowIndex, String value) {
        rowMetaData.setMetaData(rowIndex, DatasetConstants.DESCRIPTION_KEY, value);
        }

        /**
         * Sets the row name at the specified index
         *
         * @param rowIndex
         *            The row index
         * @param name
         *            The new row name value
         */
        public void setRowName(int rowIndex, String name) {
        if (name == null) {
            throw new IllegalArgumentException("Null row names are not allowed.");
        }
        if (!allowDuplicateNames && rowNameToRowIndexMap.containsKey(name) && rowNameToRowIndexMap.get(name) != rowIndex) {
            throw new IllegalArgumentException("Duplicate row names are not allowed:" + name);
        }
        rowNameToRowIndexMap.put(name, rowIndex);
        rowNames[rowIndex] = name;
        }

        /**
         * Sets a single element.
         *
         * @param value
         *            A(row,column).
         * @param rowIndex
         *            The row index.
         * @param columnIndex
         *            The column index.
         */

        public void setValue(int rowIndex, int columnIndex, double value) {
        matrix.set(rowIndex, columnIndex, value);
        }

        /**
         * Sets the value at the given row and column.
         *
         * @param rowIndex
         *            the row index.
         * @param columnIndex
         *            the column index.
         * @param seriesname
         *            the series name.
         * @throws IndexOutOfBoundsException
         *             if <tt>rowIndex</tt> or <tt>columnIndex</tt> are out of range.
         */
        public void setValue(int rowIndex, int columnIndex, String seriesname, double value) {
        double[][] m = matrices.get(seriesname);
        if (m != null) {
            m[rowIndex][columnIndex] = value;
        }

        }

        /**
         * Sets a single element.
         *
         * @param value
         *            A(row,column).
         * @param rowName
         *            The row name.
         * @param columnName
         *            The column name.
         */

        public void setValue(String rowName, String columnName, double value) {
        if (!rowNameToRowIndexMap.containsKey(rowName)) {
            throw new IllegalArgumentException("row name " + rowName + " not found.");
        }
        if (!columnNameToColumnIndexMap.containsKey(columnName)) {
            throw new IllegalArgumentException("column name " + columnName + " not found.");
        }
        matrix.set(rowNameToRowIndexMap.get(rowName), columnNameToColumnIndexMap.get(columnName), value);
        }

        public DefaultDataset slice(int[] rowIndices, int[] columnIndices) {
        if (rowIndices == null) {
            rowIndices = new int[this.getRowCount()];
            for (int i = 0; i < this.getRowCount(); i++) {
            rowIndices[i] = i;
            }
        }
        if (columnIndices == null) {
            columnIndices = new int[this.getColumnCount()];
            for (int j = 0; j < this.getColumnCount(); j++) {
            columnIndices[j] = j;
            }
        }
        double[][] sData = matrix.getMatrix(rowIndices, columnIndices).getArray();
        String[] sRowKeys = new String[rowIndices.length];
        for (int i = 0; i < rowIndices.length; i++) {
            sRowKeys[i] = rowNames[rowIndices[i]];
        }

        String[] sColumnKeys = new String[columnIndices.length];
        for (int j = 0; j < columnIndices.length; j++) {
            sColumnKeys[j] = columnNames[columnIndices[j]];
        }
        MetaData sRowMetaData = rowMetaData.slice(rowIndices);
        MetaData sColumnMetaData = columnMetaData.slice(columnIndices);
        DefaultDataset d = new DefaultDataset(sData, sRowKeys, sColumnKeys, sRowMetaData, sColumnMetaData);
        for (Iterator<String> it = matrices.keySet().iterator(); it.hasNext();) {
            String name = it.next();
            double[][] m = matrices.get(name);
            double[][] slicedData = new double[rowIndices.length][columnIndices.length];
            for (int i = 0, rows = rowIndices.length; i < rows; i++) {
            for (int j = 0, cols = columnIndices.length; j < cols; j++) {
                slicedData[i][j] = m[rowIndices[i]][columnIndices[j]];
            }
            }
            d.addMatrix(name, slicedData);
        }
        return d;
        }

        /**
         * Linear algebraic matrix multiplication, A * B
         *
         * @param B
         *            another matrix
         * @return DoubleMatrix2D product, A * B
         */
        public DefaultDataset times(DefaultDataset B) {
        Matrix product = matrix.times(B.matrix);
        DefaultDataset C = new DefaultDataset();
        C.matrix = product;
        C.rowNames = this.rowNames.clone();
        C.columnNames = B.columnNames.clone();
        C.columnNameToColumnIndexMap = (TObjectIntHashMap) B.columnNameToColumnIndexMap.clone();
        C.rowNameToRowIndexMap = (TObjectIntHashMap) this.rowNameToRowIndexMap.clone();
        return C;
        }

        /**
         * Multiply a matrix by a scalar, C = s*A
         *
         * @param d
         *            a scalar
         * @return s*A
         */
        public DefaultDataset times(double d) {
        Matrix scaledMatrix = matrix.times(d);
        DefaultDataset C = new DefaultDataset();
        C.matrix = scaledMatrix;
        C.rowNames = this.rowNames.clone();
        C.columnNames = this.columnNames.clone();
        C.columnNameToColumnIndexMap = (TObjectIntHashMap) this.columnNameToColumnIndexMap.clone();
        C.rowNameToRowIndexMap = (TObjectIntHashMap) this.rowNameToRowIndexMap.clone();
        return C;
        }

        /**
         * Computes the sum of the diagonal elements of matrix A; Sum(A[i,i]).
         *
         * @return The trace.
         */
        public double trace() {
        return this.matrix.trace();
        }

        /**
         * matrix transpose.
         *
         * @return A'
         */
        public DefaultDataset transpose() {
        DefaultDataset transpose = new DefaultDataset();
        transpose.rowNames = this.columnNames.clone();
        transpose.columnNames = this.rowNames.clone();
        transpose.columnNameToColumnIndexMap = (TObjectIntHashMap) rowNameToRowIndexMap.clone();
        transpose.rowNameToRowIndexMap = (TObjectIntHashMap) columnNameToColumnIndexMap.clone();
        transpose.matrix = this.matrix.transpose();
        return transpose;
        }

        /**
         * Constructs and returns a new slice view representing the rows of the given column. The returned view is backed by
         * this dataset, so changes in the returned view are reflected in this dataset, and vice-versa.
         *
         * @return The list.
         */
        public DoubleList viewColumn(int columnIndex) {
        final int _columnIndex = columnIndex;
        return new DoubleList() {

            public double get(int index) {
            return getValue(index, _columnIndex);
            }

            public void set(int index, double value) {
            setValue(index, _columnIndex, value);
            }

            public int size() {
            return getRowCount();
            }

        };
        }

        /**
         * Constructs and returns a new slice view representing the columns of the given row. The returned view is backed by
         * this dataset, so changes in the returned view are reflected in this dataset, and vice-versa.
         *
         * @return The list.
         */
        public DoubleList viewRow(int rowIndex) {
        final int _rowIndex = rowIndex;
        return new DoubleList() {

            public double get(int index) {
            return getValue(_rowIndex, index);
            }

            public void set(int index, double value) {
            setValue(_rowIndex, index, value);
            }

            public int size() {
            return getColumnCount();
            }

        };
        }

        /**
         * Fills in names that the user did not specify.
         *
         * @param columnIndex
         *            The column index to start filling in columns.
         */
        private void fillInColumns(int columnIndex) {
        for (int i = columnIndex, columns = getColumnCount(); i < columns; i++) {
            String name = String.valueOf(i + 1);
            if (!allowDuplicateNames && columnNameToColumnIndexMap.containsKey(name)) {
            throw new IllegalArgumentException("Duplicate column names are not allowed:" + name);
            }
            columnNameToColumnIndexMap.put(name, i);
            columnNames[i] = name;
        }
        }

        /**
         * Fills in row names that the user did not specify.
         *
         * @param rowIndex
         *            The row index to start filling in rows.
         */
        private void fillInRows(int rowIndex) {
        for (int i = rowIndex, rows = getRowCount(); i < rows; i++) {
            String name = String.valueOf(i + 1);
            if (!allowDuplicateNames && rowNameToRowIndexMap.containsKey(name)) {
            throw new IllegalArgumentException("Duplicate row names are not allowed:" + name);
            }
            rowNameToRowIndexMap.put(name, i);
            rowNames[i] = name;
        }
        }

        /**
         * Sets the column names to the specified value. Duplicate column names are not allowed. If the length of s is less
         * than getColumnDimension(), the remaining names will be filled in automatically.
         *
         * @param s
         *            The list containing the names
         */
        private void setColumnNames(List<String> s) {
        if (s.size() > getColumnCount()) {
            throw new IllegalArgumentException("Invalid column names length. getColumnDimension():" + getColumnCount()
                + " column names length:" + s.size());
        }
        for (int i = 0; i < s.size(); i++) {
            String name = s.get(i);
            if (!allowDuplicateNames && columnNameToColumnIndexMap.containsKey(name) && columnNameToColumnIndexMap.get(name) != i) {
            throw new IllegalArgumentException("Duplicate column names are not allowed:" + name);
            }
            if (name == null) {
            throw new IllegalArgumentException("Null column names are not allowed.");
            }
            columnNameToColumnIndexMap.put(name, i);
            columnNames[i] = name;
        }
        }

        /**
         * Sets the row names to the specified value. Duplicate row names are not allowed. If the length of s is less than
         * getRowDimension(), the remaining row names will be filled in automatically
         *
         * @param s
         *            The list containing the row names
         */
        private void setRowNames(List<String> s) {
        if (s.size() > getRowCount()) {
            throw new IllegalArgumentException("Invalid row names length. getRowDimension():" + getRowCount()
                + " row names length:" + s.size());
        }
        for (int i = 0, size = s.size(); i < size; i++) {
            String rowName = s.get(i);
            if (!allowDuplicateNames && rowNameToRowIndexMap.containsKey(rowName) && rowNameToRowIndexMap.get(rowName) != i) {
            throw new IllegalArgumentException("Duplicate row names are not allowed:" + rowName);
            }
            if (rowName == null) {
            throw new IllegalArgumentException("Null row names are not allowed.");
            }
            rowNameToRowIndexMap.put(rowName, i);
            rowNames[i] = rowName;
        }
        }
    }

