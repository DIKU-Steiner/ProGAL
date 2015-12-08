package ProGAL.math;

import java.util.Arrays;
import java.util.Random;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;

public class Matrix {
	protected double[][] coords;
	protected int M,N;

	/** Construct an M by N matrix with zeros in all entries. */
	public Matrix(int M, int N){
		coords = new double[M][N];
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) coords[i][j] = 0;
		this.M = M;
		this.N = N;
	}

	/** 
	 * Construct a matrix from an array-of-array-of-double. The array should
	 * index rows first and columns next. The results are undefined if the 
	 * number of elements in each row-array are different.   
	 */
	public Matrix(double[][] coords){
		this.M = coords.length;
		this.N = M>0?coords[0].length:0;
		this.coords = coords;
	}
	
	public Matrix (int M, int N, int seed) {
		coords = new double[M][N];
		Random r = new Random(seed);
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) coords[i][j] = r.nextDouble();
		this.M = M;
		this.N = N;		
	}

	public void set(int r, int c, double v){
		coords[r][c] = v;
	}

	public double get(int i, int j) {
		return coords[i][j];
	}
	
	public double[][] getCoords(){
		double[][] ret = new double[M][N];
		for(int i=0;i<M;i++){
			for(int j=0;j<M;j++){
				ret[i][j] = coords[i][j];
			}
		}
		return ret;
	}

	public int getM(){ return M; }
	public int getN(){ return N; }



	public ProGAL.geomNd.Vector getRow(int r) {
		ProGAL.geomNd.Vector ret = new ProGAL.geomNd.Vector(N);
		for(int c=0;c<N;c++) ret.set(c, coords[r][c]);
		return ret;
	}

	public ProGAL.geomNd.Vector getColumn(int c) {
		ProGAL.geomNd.Vector ret = new ProGAL.geomNd.Vector(M);
		for(int r=0;r<M;r++) ret.set(r, coords[r][c]);
		return ret;
	}

	public void setRow(int r, ProGAL.geomNd.Vector v){
		for(int c=0;c<Math.min(v.getDimensions(), N);c++){
			coords[r][c] = v.get(c);
		}
	}

	public void setColumn(int c, ProGAL.geomNd.Vector v){
		for(int r=0;r<Math.min(v.getDimensions(), M);r++){
			coords[r][c] = v.get(r);
		}
	}


	public Matrix getSubmatrix(int i1, int i2, int j1, int j2) {
		Matrix X = new Matrix(i2-i1, j2-j1);
		for (int i = 0; i < i2-i1; i++)
				for (int j = 0; j < j2-j1; j++) X.coords[i][j] = coords[i1+i][j1+j];
		return X;
	}

	public void setSubmatrix(int i1, int j1, Matrix M) {
		for (int i = 0; i < M.M; i++) 
			for (int j = 0; j < M.N; j++ ) coords[i1+i][j1+j] = M.coords[i][j];
	}


	/** Get the transpose of this matrix */
	public Matrix getTranspose() {
		double[][] newCoords = new double[N][M];
		for(int i=0;i<M;i++) 
			for(int j=0;j<N;j++) newCoords[j][i] = coords[i][j];
		return new Matrix(newCoords);
	}

	/** Get the transpose of this matrix */
	public Matrix transpose(){	return getTranspose();	}


	/** Apply this matrix to the vector v and return the result (this will change v).  */
	public ProGAL.geomNd.Vector multiplyIn(ProGAL.geomNd.Vector v){
		if(N!=v.getDimensions()) throw new Error("Dimensions dont match");
		double[] newCoords = new double[N];
		for(int r=0;r<M;r++){
			for(int c=0;c<N;c++){
				newCoords[r]+=coords[r][c]*v.get(c);
			}
		}
		v.setCoords(newCoords);
		return v;
		//		if(M==3 && N==3){
		//			double[] ret = new double[3];
		//			for(int i=0;i<3;i++) {
		//				for(int j=0;j<3;j++)
		//					ret[i] += v.get(j)*coords[i][j];
		//			}
		//			for(int i=0;i<3;i++) v.set(i, ret[i]);
		//			return v;
		//		}
		//		if( (M==4 || M==3) && N==4){
		//			double[] ret = new double[4];
		//			for(int i=0;i<4;i++) {
		//				for(int j=0;j<4;j++){
		//					ret[i] += coords[i][j]*(j==3?1:v.get(j));
		//				}
		////				ret[i] += coords[i][3];
		//			}
		//			if(Math.abs(ret[3]-1)>Constants.EPSILON)
		//				throw new RuntimeException("Multiplication with non-homogeneous coordinates failed");
		//			for(int i=0;i<3;i++) v.set(i,ret[i]);
		//			
		//			return v;
		//		}
		//		throw new Error("Can only apply 3x3, 3x4 or 4x4 matrices to vectors");
	}

	/** Apply this matrix to the point p and return the result (this will NOT change p). 
	 * This method requires the matrix to be a 3x3, 3x4 or 4x4 matrix. If it is a 
	 * 4x4 matrix the bottom row is assumed to be (0,0,0,1).*/
	public ProGAL.geomNd.Point multiply(ProGAL.geomNd.Point p){
		return multiplyIn(p.clone());
	}

	public Point multiply(Point p){
		return (Point)multiplyIn(p.clone());
	}

	/** Apply this matrix to the point p and return the result (this will change p). 
	 * This method requires the matrix to be a 3x3, 3x4 or 4x4 matrix. If it is a 
	 * 4x4 matrix the bottom row is assumed to be (0,0,0,1).*/
	public ProGAL.geomNd.Point multiplyIn(ProGAL.geomNd.Point p){
		double[] ret = new double[M];
		for(int r=0;r<M;r++) {
			for(int c=0;c<N;c++)
				ret[r] += p.get(c)*coords[r][c];
		}
		for(int r=0;r<M;r++){
			p.set(r, ret[r]);
		}
		return p;
	}

	/** Apply this matrix to the vector v and return the result (without changing v). 
	 * This method requires the matrix to be a 3x3, 3x4 or 4x4 matrix. If it is a 
	 * 4x4 matrix the bottom row is assumed to be (0,0,0,1).*/
	public Vector multiply(Vector v){
		if(M==3 && N==3){
			double[] ret = new double[3];
			for(int i=0;i<3;i++) {
				for(int j=0;j<3;j++)
					ret[i] += v.get(j)*coords[i][j];
			}
			return new Vector(ret[0],ret[1],ret[2]);
		}
		if( (M==4 || M==3) && N==4){//Assumes that the bottom row is (0,0,0,1)
			double[] ret = new double[3];
			for(int i=0;i<3;i++) {
				for(int j=0;j<3;j++)
					ret[i] += v.get(j)*coords[i][j];
				ret[i] += coords[i][3];
			}
			return new Vector(ret[0],ret[1],ret[2]);
		}
		throw new Error("Can only apply 3x3, 3x4 or 4x4 matrices to vectors");
	}

	/** Multiply this matrix to another matrix. This matrix is changed and then returned */
	public Matrix multiplyThis(Matrix m){
		double[][] newCoords = new double[M][N];
		for(int r=0;r<M;r++){
			for(int c=0;c<N;c++){
				newCoords[r][c] = 0; 
				for(int i=0;i<N;i++) newCoords[r][c]+=coords[r][i]*m.coords[i][c];
			}
		}
		this.coords = newCoords;
		return this;
	}

	/** Multiply this matrix to another matrix and return the result. */
	public Matrix multiply(Matrix m) {
		if(N!=m.M)
			throw new Error("Incompatible matrix sizes");
		Matrix ret = new Matrix(M, m.N);
		double[][] newCoords = ret.coords;
		for(int r=0;r<M;r++){
			for(int c=0;c<m.N;c++){
				for(int i=0;i<N;i++) 
					newCoords[r][c]+=coords[r][i]*m.coords[i][c];
			}
		}
		return ret;
	}

	public Matrix multiply_Strassen(Matrix B) {
		if (N != B.M) throw new Error("Incompatible matrix sizes");
		// add appropriate zero rows or columns so both matrices are 2^n x 2^n matrices
		Matrix extA = this.expand();
		Matrix extB = B.expand();
		return extA.multiply_Strassen_regular(extB);
	}
	private Matrix multiply_Strassen_regular(Matrix B)	{
		Matrix C = new Matrix(M, M);
		if (M == 1) {
			C.coords[0][0] = coords[0][0] * B.coords[0][0]; 
			return C;
		}
		// get submatrices A11, A12, A21, A22, B11, B12, B21, B22
		Matrix A11 = getSubmatrix(0, M/2, 0, M/2);
		Matrix A12 = getSubmatrix(0, M/2, M/2, M);
		Matrix A21 = getSubmatrix(M/2, M, 0, M/2);
		Matrix A22 = getSubmatrix(M/2, M, M/2, M);
		Matrix B11 = B.getSubmatrix(0, M/2, 0, M/2);
		Matrix B12 = B.getSubmatrix(0, M/2, M/2, M);
		Matrix B21 = B.getSubmatrix(M/2, M, 0, M/2);
		Matrix B22 = B.getSubmatrix(M/2, M, M/2, M);
				
		// compute M1, M2, M3, M4, M5, M6, M7
		Matrix M1 = A11.add(A22).multiply_Strassen_regular(B11.add(B22));
		Matrix M2 = A21.add(A22).multiply_Strassen_regular(B11);
		Matrix M3 = A11.multiply_Strassen_regular(B12.subtract(B22));
		Matrix M4 = A22.multiply_Strassen_regular(B21.subtract(B11));
		Matrix M5 = A11.add(A12).multiply_Strassen_regular(B22);
		Matrix M6 = A21.subtract(A11).multiply_Strassen_regular(B11.add(B12));
		Matrix M7 = A12.subtract(A22).multiply_Strassen_regular(B21.add(B22));
		
		// compute C
		C.setSubmatrix(0, 0, M1.add(M4).subtract(M5).add(M7));
		C.setSubmatrix(0, M/2, M3.add(M5));
		C.setSubmatrix(M/2, 0, M2.add(M4));
		C.setSubmatrix(M/2, M/2, M1.subtract(M2).add(M3).add(M6));
		
		return C;
	}
	
	/** Add the components of two matrices. The result is stored in a new matrix and then returned. */
	public Matrix add(Matrix m){
		Matrix ret = new Matrix(M, N);
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) ret.set(i, j, coords[i][j]+m.coords[i][j]);
		return ret;
	}

	/** Add the components of two matrices. The result is stored in this matrix and then returned. */
	public Matrix addThis(Matrix m){
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) coords[i][j]+=m.coords[i][j];
		return this;
	}

	/** Subtract the components of two matrices. The result is stored in a new matrix and then returned. */
	public Matrix subtract(Matrix m){
		Matrix ret = new Matrix(M, N);
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) ret.set(i, j, coords[i][j] - m.coords[i][j]);
		return ret;
	}

	/** Subtract the components of two matrices. The result is stored in this matrix and then returned. */
	public Matrix subtractThis(Matrix m){
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) coords[i][j] -= m.coords[i][j];
		return this;
	}

	/** Multiply the components of this matrix by a scalar. The result is stored in a new matrix which is returned. */
	public Matrix multiply(double scalar){
		Matrix ret = new Matrix(M, N);
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) ret.set(i,j,coords[i][j]*scalar);
		return ret;
	}

	/** Multiply the components of this matrix by a scalar. Changes and then returns this */
	public Matrix multiplyThis(double scalar){
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) coords[i][j]*=scalar;
		return this;
	}

//	/** Get the determinant of this matrix. Throws an error if the matrix is not square*/
//	public double determinant(){
//		if(M!=N)
//			throw new Error("Determinant undefined for non-square matrix");
//		if(M==1)
//			return coords[0][0];
//
//		double ret = 0;
//		for(int c=0;c<N;c++){
//			double minorDet = minor(0,c).determinant();
//			if(c%2==0) 	ret+=coords[0][c]*minorDet;
//			else		ret-=coords[0][c]*minorDet;
//		}
//
//		return ret;
//	}

	/** Return the minor, i.e. the matrix that results from removing row r and column c from this matrix. */
	public Matrix minor(int r, int c){
		if (M<2 || N<2) 
			throw new Error("The minor matrix is undefined for "+M+"x"+N+" matrices");

		Matrix ret = new Matrix(M-1, N-1);
		for(int i=0;i<M-1;i++){
			for(int j=0; j<N-1; j++) ret.set(i, j, coords[i >=r? i+1 : i][j >= c? j+1 : j]);
		}
		return ret;
	}
	
	/** Get the determinant of this matrix. Throws an error if the matrix is not square */
	public double determinant(){
		if(M!=N) throw new Error("Determinant undefined for non-square matrix");
		if(M==1) return coords[0][0];
		return new LUDecomposition(this).det();
	}

	
	public boolean isSquare(){ return M==N; }
	
	
	/** extends to the smallest square matrix by adding appropriate zero rows or columns */
	public Matrix extend() {
		Matrix X;
		if (M < N) {
			X = new Matrix(N, N);
			X.setSubmatrix(0, 0, this);
			for (int i = M; i < N; i++)
				for (int j = 0; j < N; j++) X.coords[i][j] = 0.0;
			
		}
		else {
			X = new Matrix(M, M);
			X.setSubmatrix(0,  0, this);
			for (int j = N; j < M; j++) 
				for (int i = 0; i < M; i++) X.coords[i][j] = 0.0;
		}
		return X;
	}
	
	/** extends to the smallest power of 2 square matrix by adding appropriate zero rows and columns */
	public Matrix expand() {
		Matrix X = extend();
		int ext = Functions.roundUpToPowerOf2(X.M);
		Matrix Y = new Matrix(ext, ext);
		Y.setSubmatrix(0, 0, X);
		for (int i = 0; i < X.M; i++) for (int j = X.M; j < ext; j++) Y.coords[i][j] = 0.0;
		for (int i = X.M; i < ext; i++) for (int j = 0; j < ext; j++) Y.coords[i][j] = 0.0;
		return Y;
	}

	/** Return the inverse of this matrix. */
	public Matrix invert(){
		Matrix ret = clone();
		return ret.invertThis();
	}

	/** Invert this matrix (overwrites this and returns it).  */
	public Matrix invertThis(){
		if(!isSquare()) 
			throw new Error("Cant invert non-square matrix ("+M+"x"+N+")");
		Matrix tmp = new Matrix(M,2*N);
		for(int r=0;r<M;r++){
			for(int c=0;c<N;c++) tmp.set(r, c, get(r,c));
			tmp.set(r,N+r, 1);
		}
		tmp.reduceThis();

		for(int r=0;r<M;r++){
			for(int c=0;c<N;c++) set(r, c, tmp.get(r,c+N));
		}
		return this;
	}

	/** Reduce this matrix to row canonical form (reduced row echelon form). A new 
	 * object is returned. */
	public Matrix reduce(){
		Matrix ret = this.clone();
		ret.reduceThis();
		return ret;
	}

	/** Reduce this matrix to row canonical form (reduced row echelon form). This 
	 * matrix is changed and returned.
	 * @hops For n by m matrix: 0 to m*n+(n-1)^2*m
	 * @todo TODO Optimize loops to remove 0-nominator divisions.  */
	public Matrix reduceThis(){
		int rowCount = M;
		int colCount = N;
		int lead = 0;
		for(int r=0;r<rowCount;r++){
			if(lead>=colCount)	return this;
			int i=r;
			while(Math.abs(coords[i][lead])<=Constants.EPSILON){
				i++;
				if(rowCount==i){
					i=r;
					lead++;
					if(colCount==lead) {
						return this;
					}
				}
			}
			//Swap rows i and r
			if(i!=r){
				double[] tmpArr = coords[i];
				coords[i] = coords[r];
				coords[r] = tmpArr;
			}

			//Divide row r by coords[r][lead]
			double tmp = coords[r][lead];
			for(int c=0;c<coords[r].length;c++) coords[r][c]/=tmp;

			for(i=0;i<rowCount;i++){
				if(i!=r){
					//Subtract coords[i][lead] multiplied by row r from row i
					tmp = coords[i][lead];
					for(int c=0;c<coords[i].length;c++) 
						coords[i][c]-=tmp*coords[r][c]; 
				}
			}
			lead++;
		}
		return this;
	}

	/**
	 * Run Gram-Schmidt orthonormalization procedure. Assumes that the columns of this
	 * are linearly independent.
	 * @return this
	 */
	public Matrix orthonormalizeThis(){
		ProGAL.geomNd.Vector[] bs = new ProGAL.geomNd.Vector[N]; 
		bs[0] = this.getColumn(0);
		double[] bSq = new double[N];
		for(int i=1;i<N;i++){
			ProGAL.geomNd.Vector ai = this.getColumn(i);
			bs[i] = ai.clone();
			bSq[i-1] = bs[i-1].dot(bs[i-1]); 
			for(int j=0;j<i;j++) bs[i].addThis(bs[j].multiply(-ai.dot(bs[j])/bSq[j])); 
		}
		bSq[N-1] = bs[N-1].dot(bs[N-1]);
		for(int i=0;i<N;i++){
			if(bSq[i]!=1){
				bs[i].divideThis(Math.sqrt(bSq[i]));
			}

			for(int r=0;r<M;r++){
				coords[i][r] = bs[i].get(r);
			}
		}
		return this;
	}

	public boolean equals(Matrix m){
		if(M!=m.M) return false;
		for(int r=0;r<M;r++)
			if(!Arrays.equals(coords[r], m.coords[r])) return false;
		return true;
	}
	public boolean equals(Object o){
		if(o instanceof Matrix) return equals((Matrix)o);
		return false;
	}

	public EigenvalueDecomposition getEigenvalueDecomposition(){
		return new EigenvalueDecomposition();
	}

	public String toString(){ return toString(4); }

	public String toString(int dec){
		StringBuilder sb = new StringBuilder();
		for(int i=0;i<M;i++) {
			sb.append('|');
			for(int j=0;j<N;j++) {
				sb.append(String.format("% 9."+dec+"f", coords[i][j]));
				sb.append(' ');
			}
			sb.append('|');
			sb.append('\n');
		}
		return sb.toString();
	}

	public void toConsole(){ toConsole(2); 	}

	public void toConsole(int dec){
		System.out.print(toString(dec));
	}

	/** Returns a clone of this matrix. */ 
	public Matrix clone(){
		Matrix ret = new Matrix(M, N);
		for(int r=0;r<M;r++) for(int c=0;c<N;c++)
			ret.coords[r][c] = coords[r][c];
		return ret;
	}


	public static Matrix createRandomMatrix(int n) {
		Matrix ret = new Matrix(n,n);
		for(int i=0;i<n;i++){ 
			for(int j=0;j<n;j++){ 
					ret.coords[i][j] = Randomization.randBetween(0.0, 1.0);
			}
		}
		return ret;
	}
	public static Matrix createIdentityMatrix(int n) {
		Matrix ret = new Matrix(n,n);
		for(int i=0;i<n;i++) ret.coords[i][i] = 1;
		return ret;
	}

	public static Matrix3x3 createColumnMatrix(Vector v1, Vector v2, Vector v3){
		Matrix3x3 ret = new Matrix3x3();
		for(int i=0;i<3;i++) ret.coords[i][0] = v1.get(i);
		for(int i=0;i<3;i++) ret.coords[i][1] = v2.get(i);
		for(int i=0;i<3;i++) ret.coords[i][2] = v3.get(i);
		return ret;
	}
	public static Matrix create4x4ColumnMatrix(Vector v1, Vector v2, Vector v3, Vector v4){
		Matrix ret = new Matrix(4,4);
		for(int i=0;i<3;i++) ret.coords[i][0] = v1.get(i);
		for(int i=0;i<3;i++) ret.coords[i][1] = v2.get(i);
		for(int i=0;i<3;i++) ret.coords[i][2] = v3.get(i);
		for(int i=0;i<3;i++) ret.coords[i][3] = v4.get(i);
		ret.coords[3][3] = 1;
		return ret;
	}
	public static Matrix3x3 createRowMatrix(Vector v1, Vector v2, Vector v3){
		Matrix3x3 ret = new Matrix3x3();
		for(int i=0;i<3;i++) ret.coords[0][i] = v1.get(i);
		for(int i=0;i<3;i++) ret.coords[1][i] = v2.get(i);
		for(int i=0;i<3;i++) ret.coords[2][i] = v3.get(i);
		return ret;
	}

	//Daisy says : The returned rotation matrix changes the size of the vector it is applied to
	public static Matrix createRotationMatrix(double angle, Vector v){
		double ux = v.x();
		double uy = v.y();
		double uz = v.z();
		double cosA = Math.cos(angle);
		double sinA = Math.sin(angle);

		Matrix ret = new Matrix3x3();
		ret.coords[0][0] = ux*ux + cosA*(1.0f-ux*ux);
		ret.coords[1][0] = ux*uy*(1.0f-cosA) + uz*sinA;
		ret.coords[2][0] = uz*ux*(1.0f-cosA) - uy*sinA;

		ret.coords[0][1] = ux*uy*(1.0f-cosA) - uz*sinA;
		ret.coords[1][1] = uy*uy + cosA*(1.0f-uy*uy);
		ret.coords[2][1] = uy*uz*(1.0f-cosA) + ux*sinA;

		ret.coords[0][2] = uz*ux*(1.0f-cosA) + uy*sinA;
		ret.coords[1][2] = uy*uz*(1.0f-cosA) - ux*sinA;
		ret.coords[2][2] = uz*uz + cosA*(1.0f-uz*uz);
		return ret;
	}



	/** Eigenvalues and eigenvectors of a real matrix. If A is symmetric, 
	 * then A = V*D*V' where the eigenvalue matrix D is diagonal and the 
	 * eigenvector matrix V is orthogonal. I.e. A = V.times(D.times(V.transpose())) 
	 * and	V.times(V.transpose()) equals the identity matrix.
	 * 
	 * If A is not symmetric, then the eigenvalue matrix D is block diagonal 
	 * with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues, 
	 * lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The 
	 * columns of V represent the eigenvectors in the sense that A*V = V*D, 
	 * i.e. A.times(V) equals V.times(D).  The matrix V may be badly 
	 * conditioned, or even singular, so the validity of the equation 
	 * A = V*D*inverse(V) depends upon V.cond().
	 */
	public class EigenvalueDecomposition {

		/** Row and column dimension (square matrix). */
		private int n;
		/** Symmetry flag. */
		private boolean issymmetric;
		/** Arrays for internal storage of eigenvalues. */
		private double[] d, e;
		/** Array for internal storage of eigenvectors. */
		private double[][] V;
		/** Array for internal storage of nonsymmetric Hessenberg form. */
		private double[][] H;
		/** Working storage for nonsymmetric algorithm. */
		private double[] ort;

		/** Symmetric Householder reduction to tridiagonal form. */
		private void tred2 () {

			//  This is derived from the Algol procedures tred2 by
			//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
			//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
			//  Fortran subroutine in EISPACK.
			for (int j = 0; j < n; j++) d[j] = V[n-1][j];

			// Householder reduction to tridiagonal form.
			for (int i = n-1; i > 0; i--) {

				// Scale to avoid under/overflow.
				double scale = 0.0;
				double h = 0.0;
				for (int k = 0; k < i; k++) scale = scale + Math.abs(d[k]);
				if (scale == 0.0) {
					e[i] = d[i-1];
					for (int j = 0; j < i; j++) {
						d[j] = V[i-1][j];
						V[i][j] = 0.0;
						V[j][i] = 0.0;
					}
				} else {

					// Generate Householder vector.
					for (int k = 0; k < i; k++) {
						d[k] /= scale;
						h += d[k] * d[k];
					}
					double f = d[i-1];
					double g = Math.sqrt(h);
					if (f > 0) g = -g;
					e[i] = scale * g;
					h = h - f * g;
					d[i-1] = f - g;
					for (int j = 0; j < i; j++) e[j] = 0.0;

					// Apply similarity transformation to remaining columns.
					for (int j = 0; j < i; j++) {
						f = d[j];
						V[j][i] = f;
						g = e[j] + V[j][j] * f;
						for (int k = j+1; k <= i-1; k++) {
							g += V[k][j] * d[k];
							e[k] += V[k][j] * f;
						}
						e[j] = g;
					}
					f = 0.0;
					for (int j = 0; j < i; j++) {
						e[j] /= h;
						f += e[j] * d[j];
					}
					double hh = f / (h + h);
					for (int j = 0; j < i; j++) {
						e[j] -= hh * d[j];
					}
					for (int j = 0; j < i; j++) {
						f = d[j];
						g = e[j];
						for (int k = j; k <= i-1; k++) 
							V[k][j] -= (f * e[k] + g * d[k]);
						d[j] = V[i-1][j];
						V[i][j] = 0.0;
					}
				}
				d[i] = h;
			}

			// Accumulate transformations.
			for (int i = 0; i < n-1; i++) {
				V[n-1][i] = V[i][i];
				V[i][i] = 1.0;
				double h = d[i+1];
				if (h != 0.0) {
					for (int k = 0; k <= i; k++) 
						d[k] = V[k][i+1] / h;
					for (int j = 0; j <= i; j++) {
						double g = 0.0;
						for (int k = 0; k <= i; k++) 
							g += V[k][i+1] * V[k][j];
						for (int k = 0; k <= i; k++) 
							V[k][j] -= g * d[k];
					}
				}
				for (int k = 0; k <= i; k++) V[k][i+1] = 0.0;
			}
			for (int j = 0; j < n; j++) {
				d[j] = V[n-1][j];
				V[n-1][j] = 0.0;
			}
			V[n-1][n-1] = 1.0;
			e[0] = 0.0;
		} 

		/** Symmetric tridiagonal QL algorithm. */
		private void tql2() {

			//  This is derived from the Algol procedures tql2, by
			//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
			//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
			//  Fortran subroutine in EISPACK.

			for (int i = 1; i < n; i++) e[i-1] = e[i];
			e[n-1] = 0.0;

			double f = 0.0;
			double tst1 = 0.0;
			double eps = Math.pow(2.0,-52.0);
			for (int l = 0; l < n; l++) {

				// Find small subdiagonal element
				tst1 = Math.max(tst1,Math.abs(d[l]) + Math.abs(e[l]));
				int m = l;
				while (m < n) {
					if (Math.abs(e[m]) <= eps*tst1) break;
					m++;
				}

				// If m == l, d[l] is an eigenvalue,
				// otherwise, iterate.
				if (m > l) {
					int iter = 0;
					do {
						iter = iter + 1;  // (Could check iteration count here.)

						// Compute implicit shift
						double g = d[l];
						double p = (d[l+1] - g) / (2.0 * e[l]);
						double r = hypot(p,1.0);
						if (p < 0) r = -r;
						d[l] = e[l] / (p + r);
						d[l+1] = e[l] * (p + r);
						double dl1 = d[l+1];
						double h = g - d[l];
						for (int i = l+2; i < n; i++) d[i] -= h;
						f = f + h;

						// Implicit QL transformation.
						p = d[m];
						double c = 1.0;
						double c2 = c;
						double c3 = c;
						double el1 = e[l+1];
						double s = 0.0;
						double s2 = 0.0;
						for (int i = m-1; i >= l; i--) {
							c3 = c2;
							c2 = c;
							s2 = s;
							g = c * e[i];
							h = c * p;
							r = hypot(p,e[i]);
							e[i+1] = s * r;
							s = e[i] / r;
							c = p / r;
							p = c * d[i] - s * g;
							d[i+1] = h + s * (c * g + s * d[i]);

							// Accumulate transformation.
							for (int k = 0; k < n; k++) {
								h = V[k][i+1];
								V[k][i+1] = s * V[k][i] + c * h;
								V[k][i] = c * V[k][i] - s * h;
							}
						}
						p = -s * s2 * c3 * el1 * e[l] / dl1;
						e[l] = s * p;
						d[l] = c * p;

						// Check for convergence.
					} while (Math.abs(e[l]) > eps*tst1);
				}
				d[l] = d[l] + f;
				e[l] = 0.0;
			}

			// Sort eigenvalues and corresponding vectors.
			for (int i = 0; i < n-1; i++) {
				int k = i;
				double p = d[i];
				for (int j = i+1; j < n; j++) {
					if (d[j] < p) {
						k = j;
						p = d[j];
					}
				}
				if (k != i) {
					d[k] = d[i];
					d[i] = p;
					for (int j = 0; j < n; j++) {
						p = V[j][i];
						V[j][i] = V[j][k];
						V[j][k] = p;
					}
				}
			}
		}

		/** Nonsymmetric reduction to Hessenberg form. */
		private void orthes() {

			//  This is derived from the Algol procedures orthes and ortran,
			//  by Martin and Wilkinson, Handbook for Auto. Comp.,
			//  Vol.ii-Linear Algebra, and the corresponding
			//  Fortran subroutines in EISPACK.
			int low = 0;
			int high = n-1;

			for (int m = low+1; m <= high-1; m++) {

				// Scale column.
				double scale = 0.0;
				for (int i = m; i <= high; i++) 
					scale = scale + Math.abs(H[i][m-1]);

				if (scale != 0.0) {

					// Compute Householder transformation.

					double h = 0.0;
					for (int i = high; i >= m; i--) {
						ort[i] = H[i][m-1]/scale;
						h += ort[i] * ort[i];
					}
					double g = Math.sqrt(h);
					if (ort[m] > 0) g = -g;
					h = h - ort[m] * g;
					ort[m] = ort[m] - g;

					// Apply Householder similarity transformation
					// H = (I-u*u'/h)*H*(I-u*u')/h)
					for (int j = m; j < n; j++) {
						double f = 0.0;
						for (int i = high; i >= m; i--) 
							f += ort[i]*H[i][j];
						f = f/h;
						for (int i = m; i <= high; i++) 
							H[i][j] -= f*ort[i];
					}

					for (int i = 0; i <= high; i++) {
						double f = 0.0;
						for (int j = high; j >= m; j--) 
							f += ort[j]*H[i][j];
						f = f/h;
						for (int j = m; j <= high; j++) 
							H[i][j] -= f*ort[j];
					}
					ort[m] = scale*ort[m];
					H[m][m-1] = scale*g;
				}
			}

			// Accumulate transformations (Algol's ortran).
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					V[i][j] = (i == j ? 1.0 : 0.0);
				}
			}

			for (int m = high-1; m >= low+1; m--) {
				if (H[m][m-1] != 0.0) {
					for (int i = m+1; i <= high; i++) 
						ort[i] = H[i][m-1];
					for (int j = m; j <= high; j++) {
						double g = 0.0;
						for (int i = m; i <= high; i++) 
							g += ort[i] * V[i][j];
						// Double division avoids possible underflow
						g = (g / ort[m]) / H[m][m-1];
						for (int i = m; i <= high; i++) 
							V[i][j] += g * ort[i];
					}
				}
			}
		}


		private transient double cdivr, cdivi;
		/** Complex scalar division. */
		private void cdiv(double xr, double xi, double yr, double yi) {
			double r,d;
			if (Math.abs(yr) > Math.abs(yi)) {
				r = yi/yr;
				d = yr + r*yi;
				cdivr = (xr + r*xi)/d;
				cdivi = (xi - r*xr)/d;
			} else {
				r = yr/yi;
				d = yi + r*yr;
				cdivr = (r*xr + xi)/d;
				cdivi = (r*xi - xr)/d;
			}
		}


		/** Nonsymmetric reduction from Hessenberg to real Schur form. */
		private void hqr2() {

			//  This is derived from the Algol procedure hqr2,
			//  by Martin and Wilkinson, Handbook for Auto. Comp.,
			//  Vol.ii-Linear Algebra, and the corresponding
			//  Fortran subroutine in EISPACK.

			// Initialize
			int nn = this.n;
			int n = nn-1;
			int low = 0;
			int high = nn-1;
			double eps = Math.pow(2.0,-52.0);
			double exshift = 0.0;
			double p=0,q=0,r=0,s=0,z=0,t,w,x,y;

			// Store roots isolated by balanc and compute matrix norm
			double norm = 0.0;
			for (int i = 0; i < nn; i++) {
				if (i < low | i > high) {
					d[i] = H[i][i];
					e[i] = 0.0;
				}
				for (int j = Math.max(i-1,0); j < nn; j++) {
					norm = norm + Math.abs(H[i][j]);
				}
			}

			// Outer loop over eigenvalue index
			int iter = 0;
			while (n >= low) {

				// Look for single small sub-diagonal element
				int l = n;
				while (l > low) {
					s = Math.abs(H[l-1][l-1]) + Math.abs(H[l][l]);
					if (s == 0.0) s = norm;
					if (Math.abs(H[l][l-1]) < eps * s) break;
					l--;
				}

				// Check for convergence
				// One root found
				if (l == n) {
					H[n][n] = H[n][n] + exshift;
					d[n] = H[n][n];
					e[n] = 0.0;
					n--;
					iter = 0;

					// Two roots found
				} else if (l == n-1) {
					w = H[n][n-1] * H[n-1][n];
					p = (H[n-1][n-1] - H[n][n]) / 2.0;
					q = p * p + w;
					z = Math.sqrt(Math.abs(q));
					H[n][n] = H[n][n] + exshift;
					H[n-1][n-1] = H[n-1][n-1] + exshift;
					x = H[n][n];

					// Real pair
					if (q >= 0) {
						if (p >= 0)	z = p + z;
						else 		z = p - z;
						d[n-1] = x + z;
						d[n] = d[n-1];
						if (z != 0.0) d[n] = x - w / z;
						e[n-1] = 0.0;
						e[n] = 0.0;
						x = H[n][n-1];
						s = Math.abs(x) + Math.abs(z);
						p = x / s;
						q = z / s;
						r = Math.sqrt(p * p+q * q);
						p = p / r;
						q = q / r;

						// Row modification
						for (int j = n-1; j < nn; j++) {
							z = H[n-1][j];
							H[n-1][j] = q * z + p * H[n][j];
							H[n][j] = q * H[n][j] - p * z;
						}

						// Column modification
						for (int i = 0; i <= n; i++) {
							z = H[i][n-1];
							H[i][n-1] = q * z + p * H[i][n];
							H[i][n] = q * H[i][n] - p * z;
						}

						// Accumulate transformations
						for (int i = low; i <= high; i++) {
							z = V[i][n-1];
							V[i][n-1] = q * z + p * V[i][n];
							V[i][n] = q * V[i][n] - p * z;
						}

						// Complex pair
					} else {
						d[n-1] = x + p;
						d[n] = x + p;
						e[n-1] = z;
						e[n] = -z;
					}
					n = n - 2;
					iter = 0;

					// No convergence yet
				} else {

					// Form shift
					x = H[n][n];
					y = 0.0;
					w = 0.0;
					if (l < n) {
						y = H[n-1][n-1];
						w = H[n][n-1] * H[n-1][n];
					}

					// Wilkinson's original ad hoc shift
					if (iter == 10) {
						exshift += x;
						for (int i = low; i <= n; i++) H[i][i] -= x;

						s = Math.abs(H[n][n-1]) + Math.abs(H[n-1][n-2]);
						x = y = 0.75 * s;
						w = -0.4375 * s * s;
					}

					// MATLAB's new ad hoc shift
					if (iter == 30) {
						s = (y - x) / 2.0;
						s = s * s + w;
						if (s > 0) {
							s = Math.sqrt(s);
							if (y < x) s = -s;
							s = x - w / ((y - x) / 2.0 + s);
							for (int i = low; i <= n; i++) H[i][i] -= s;
							exshift += s;
							x = y = w = 0.964;
						}
					}

					iter = iter + 1;   // (Could check iteration count here.)

					// Look for two consecutive small sub-diagonal elements
					int m = n-2;
					while (m >= l) {
						z = H[m][m];
						r = x - z;
						s = y - z;
						p = (r * s - w) / H[m+1][m] + H[m][m+1];
						q = H[m+1][m+1] - z - r - s;
						r = H[m+2][m+1];
						s = Math.abs(p) + Math.abs(q) + Math.abs(r);
						p = p / s;
						q = q / s;
						r = r / s;
						if (m == l) break;
						if (Math.abs(H[m][m-1]) * (Math.abs(q) + Math.abs(r)) <
								eps * (Math.abs(p) * (Math.abs(H[m-1][m-1]) + Math.abs(z) +
										Math.abs(H[m+1][m+1])))) {
							break;
						}
						m--;
					}

					for (int i = m+2; i <= n; i++) {
						H[i][i-2] = 0.0;
						if (i > m+2) H[i][i-3] = 0.0;
					}

					// Double QR step involving rows l:n and columns m:n
					for (int k = m; k <= n-1; k++) {
						boolean notlast = (k != n-1);
						if (k != m) {
							p = H[k][k-1];
							q = H[k+1][k-1];
							r = (notlast ? H[k+2][k-1] : 0.0);
							x = Math.abs(p) + Math.abs(q) + Math.abs(r);
							if (x != 0.0) {
								p = p / x;
								q = q / x;
								r = r / x;
							}
						}
						if (x == 0.0) break;
						s = Math.sqrt(p * p + q * q + r * r);
						if (p < 0) s = -s;
						if (s != 0) {
							if (k != m) {
								H[k][k-1] = -s * x;
							} else if (l != m) {
								H[k][k-1] = -H[k][k-1];
							}
							p = p + s;
							x = p / s;
							y = q / s;
							z = r / s;
							q = q / p;
							r = r / p;

							// Row modification
							for (int j = k; j < nn; j++) {
								p = H[k][j] + q * H[k+1][j];
								if (notlast) {
									p = p + r * H[k+2][j];
									H[k+2][j] = H[k+2][j] - p * z;
								}
								H[k][j] = H[k][j] - p * x;
								H[k+1][j] = H[k+1][j] - p * y;
							}

							// Column modification
							for (int i = 0; i <= Math.min(n,k+3); i++) {
								p = x * H[i][k] + y * H[i][k+1];
								if (notlast) {
									p = p + z * H[i][k+2];
									H[i][k+2] = H[i][k+2] - p * r;
								}
								H[i][k] = H[i][k] - p;
								H[i][k+1] = H[i][k+1] - p * q;
							}

							// Accumulate transformations
							for (int i = low; i <= high; i++) {
								p = x * V[i][k] + y * V[i][k+1];
								if (notlast) {
									p = p + z * V[i][k+2];
									V[i][k+2] = V[i][k+2] - p * r;
								}
								V[i][k] = V[i][k] - p;
								V[i][k+1] = V[i][k+1] - p * q;
							}
						}  // (s != 0)
					}  // k loop
				}  // check convergence
			}  // while (n >= low)

			// Backsubstitute to find vectors of upper triangular form

			if (norm == 0.0) return;

			for (n = nn-1; n >= 0; n--) {
				p = d[n];
				q = e[n];

				// Real vector
				if (q == 0) {
					int l = n;
					H[n][n] = 1.0;
					for (int i = n-1; i >= 0; i--) {
						w = H[i][i] - p;
						r = 0.0;
						for (int j = l; j <= n; j++) {
							r = r + H[i][j] * H[j][n];
						}
						if (e[i] < 0.0) {
							z = w;
							s = r;
						} else {
							l = i;
							if (e[i] == 0.0) {
								if (w != 0.0) 	H[i][n] = -r / w;
								else			H[i][n] = -r / (eps * norm);

								// Solve real equations
							} else {
								x = H[i][i+1];
								y = H[i+1][i];
								q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
								t = (x * s - z * r) / q;
								H[i][n] = t;
								if (Math.abs(x) > Math.abs(z)) 	H[i+1][n] = (-r - w * t) / x;
								else 							H[i+1][n] = (-s - y * t) / z;
							}

							// Overflow control
							t = Math.abs(H[i][n]);
							if ((eps * t) * t > 1) {
								for (int j = i; j <= n; j++) {
									H[j][n] = H[j][n] / t;
								}
							}
						}
					}

					// Complex vector
				} else if (q < 0) {
					int l = n-1;

					// Last vector component imaginary so matrix is triangular
					if (Math.abs(H[n][n-1]) > Math.abs(H[n-1][n])) {
						H[n-1][n-1] = q / H[n][n-1];
						H[n-1][n] = -(H[n][n] - p) / H[n][n-1];
					} else {
						cdiv(0.0,-H[n-1][n],H[n-1][n-1]-p,q);
						H[n-1][n-1] = cdivr;
						H[n-1][n] = cdivi;
					}
					H[n][n-1] = 0.0;
					H[n][n] = 1.0;
					for (int i = n-2; i >= 0; i--) {
						double ra,sa,vr,vi;
						ra = 0.0;
						sa = 0.0;
						for (int j = l; j <= n; j++) {
							ra = ra + H[i][j] * H[j][n-1];
							sa = sa + H[i][j] * H[j][n];
						}
						w = H[i][i] - p;

						if (e[i] < 0.0) {
							z = w;
							r = ra;
							s = sa;
						} else {
							l = i;
							if (e[i] == 0) {
								cdiv(-ra,-sa,w,q);
								H[i][n-1] = cdivr;
								H[i][n] = cdivi;
							} else {

								// Solve complex equations
								x = H[i][i+1];
								y = H[i+1][i];
								vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
								vi = (d[i] - p) * 2.0 * q;
								if (vr == 0.0 & vi == 0.0) {
									vr = eps * norm * (Math.abs(w) + Math.abs(q) +
											Math.abs(x) + Math.abs(y) + Math.abs(z));
								}
								cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
								H[i][n-1] = cdivr;
								H[i][n] = cdivi;
								if (Math.abs(x) > (Math.abs(z) + Math.abs(q))) {
									H[i+1][n-1] = (-ra - w * H[i][n-1] + q * H[i][n]) / x;
									H[i+1][n] = (-sa - w * H[i][n] - q * H[i][n-1]) / x;
								} else {
									cdiv(-r-y*H[i][n-1],-s-y*H[i][n],z,q);
									H[i+1][n-1] = cdivr;
									H[i+1][n] = cdivi;
								}
							}

							// Overflow control
							t = Math.max(Math.abs(H[i][n-1]),Math.abs(H[i][n]));
							if ((eps * t) * t > 1) {
								for (int j = i; j <= n; j++) {
									H[j][n-1] = H[j][n-1] / t;
									H[j][n] = H[j][n] / t;
								}
							}
						}
					}
				}
			}

			// Vectors of isolated roots
			for (int i = 0; i < nn; i++) {
				if (i < low | i > high) {
					for (int j = i; j < nn; j++) {
						V[i][j] = H[i][j];
					}
				}
			}

			// Back transformation to get eigenvectors of original matrix
			for (int j = nn-1; j >= low; j--) {
				for (int i = low; i <= high; i++) {
					z = 0.0;
					for (int k = low; k <= Math.min(j,high); k++) {
						z = z + V[i][k] * H[k][j];
					}
					V[i][j] = z;
				}
			}
		}

		/** Calculate <code>sqrt(a^2 + b^2)</code> without under/overflow. */
		private double hypot(double a, double b) {
			if (Math.abs(a) > Math.abs(b)) {
				double r = b/a;
				return Math.abs(a)*Math.sqrt(1+r*r);
			} else if (b != 0) {
				double r = a/b;
				return Math.abs(b)*Math.sqrt(1+r*r);
			} else {
				return 0.0;
			}
		}


		/** 
		 * Check for symmetry, then construct the eigenvalue decomposition.
		 * @param A    Square matrix
		 */
		protected EigenvalueDecomposition() {
			double[][] A = coords;//Arg.getArray();
			n = N;//Arg.getColumnDimension();
			V = new double[n][n];
			d = new double[n];
			e = new double[n];

			issymmetric = true;
			for (int j = 0; (j < n) & issymmetric; j++) {
				for (int i = 0; (i < n) & issymmetric; i++) {
					issymmetric = (A[i][j] == A[j][i]);
				}
			}

			if (issymmetric) {
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n; j++) {
						V[i][j] = A[i][j];
					}
				}

				tred2();// Tridiagonalize.
				tql2();	// Diagonalize.

			} else {
				H = new double[n][n];
				ort = new double[n];

				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						H[i][j] = A[i][j];
					}
				}


				orthes();	// Reduce to Hessenberg form.
				hqr2();		// Reduce Hessenberg to real Schur form.
			}
		}

		/** Return the eigenvector matrix (immutable). */
		public Matrix getV() { 					return new ImmutableMatrix(V); 	}

		/** Return the real parts of the eigenvalues real(diag(D)). */
		public double[] getRealEigenvalues() {		return d;	}

		/** Return the imaginary parts of the eigenvalues imag(diag(D)). */
		public double[] getImagEigenvalues() {		return e;	}

		/** Return the block diagonal eigenvalue matrix. */
		public Matrix getD() {
			Matrix X = new Matrix(n,n);
			double[][] D = X.coords;
			for (int i = 0; i < n; i++) {
				D[i][i] = d[i];
				if (e[i] > 0) {
					D[i][i+1] = e[i];
				} else if (e[i] < 0) {
					D[i][i-1] = e[i];
				}
			}
			return X;
		}
	}

	public static class ImmutableMatrix extends Matrix{
		public ImmutableMatrix(double[][] coords) {
			super(coords);
		}
		public void set(int r, int c, double v){throw new RuntimeException("This matrix is immutable");}
		public Matrix multiplyThis(Matrix m){ return multiply(m); }
		public Matrix addThis(Matrix m){ return add(m); }
		public Matrix multiplyThis(double s){ return multiply(s); }
		public Matrix invertThis(){ return invert(); }
	}

	/** Get a submatrix.
	   @param r    Array of row indices.
	   @param i0   Initial column index
	   @param i1   Final column index
	   @return     A(r(:),j0:j1)
	   @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */

	Matrix getMatrix (int[] r, int j0, int j1) {
		Matrix X = new Matrix(r.length,j1-j0+1);
		double[][] B = X.coords;
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = j0; j <= j1; j++) {
					B[i][j-j0] = coords[r[i]][j];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
		return X;
	}

}
