package ProGAL.math;

public class Matrix {
	private double[][] coords;
	public Matrix(int N, int M){
		coords = new double[N][M];
		for(int i=0;i<N;i++) for(int j=0;j<M;j++) coords[i][j] = 0;
	}
	public void set(int r, int c, double v){
		coords[r][c] = v;
	}

	/*public static Matrix createColumnMatrix(Vector v1, Vector v2, Vector v3){
		Matrix ret = new Matrix(3,3);
		ret.coords[0][0] = v1.x;
		ret.coords[1][0] = v1.y;
		ret.coords[2][0] = v1.z;
		ret.coords[0][1] = v2.x;
		ret.coords[1][1] = v2.y;
		ret.coords[2][1] = v2.z;
		ret.coords[0][2] = v3.x;
		ret.coords[1][2] = v3.y;
		ret.coords[2][2] = v3.z;
		return ret;
	}
	public static Matrix create4x4ColumnMatrix(Vector v1, Vector v2, Vector v3, Vector v4){
		Matrix ret = new Matrix(4,4);
		ret.coords[0][0] = v1.x;
		ret.coords[1][0] = v1.y;
		ret.coords[2][0] = v1.z;
		ret.coords[3][0] = 0;
		ret.coords[0][1] = v2.x;
		ret.coords[1][1] = v2.y;
		ret.coords[2][1] = v2.z;
		ret.coords[3][1] = 0;
		ret.coords[0][2] = v3.x;
		ret.coords[1][2] = v3.y;
		ret.coords[2][2] = v3.z;
		ret.coords[3][2] = 0;
		ret.coords[0][3] = v4.x;
		ret.coords[1][3] = v4.y;
		ret.coords[2][3] = v4.z;
		ret.coords[3][3] = 1;
		return ret;
	}
	public static Matrix createRowMatrix(Vector v1, Vector v2, Vector v3){
		Matrix ret = new Matrix(3,3);
		ret.coords[0][0] = v1.x;
		ret.coords[0][1] = v1.y;
		ret.coords[0][2] = v1.z;
		ret.coords[1][0] = v2.x;
		ret.coords[1][1] = v2.y;
		ret.coords[1][2] = v2.z;
		ret.coords[2][0] = v3.x;
		ret.coords[2][1] = v3.y;
		ret.coords[2][2] = v3.z;
		return ret;
	}

	public static Matrix createRotationMatrix(float angle, Vector v){
		float ux = v.x;
		float uy = v.y;
		float uz = v.z;
		float cosA = (float)Math.cos(angle);
		float sinA = (float)Math.sin(angle);

		Matrix ret = new Matrix(3,3);
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
	
	public static Matrix createCovarianceMatrix(List<Vector> points){
		Matrix ret = new Matrix(3,3);
		double[] means = {0,0,0};
		for(int d=0;d<3;d++){
			double mean = 0;
			for(Vector p: points) mean+=p.get(d);
			means[d] = mean/points.size();
		}
		//System.out.println("Mean "+means[0]+" .. "+means[1]+" .. "+means[2]);
		for(int i=0;i<3;i++){
			for(int j=i;j<3;j++){
				ret.coords[i][j] = (float)covariance(points, means, i,j);
				ret.coords[j][i] = ret.coords[i][j];
			}
		}

		return ret;
	}

	private static double covariance(List<Vector> points, double[] means, int dim1, int dim2){
		double ret = 0;
		for(Vector v: points){
			ret+=(v.get(dim1)-means[dim1])*(v.get(dim2)-means[dim2]);
		}

		return ret/points.size();
	}
*/

	/** Apply this matrix to the vector v. The vector v is changed and then returned. */
	/*
	public Vector applyToIn(Vector v){
		assert(N==M);
		if(N==M && N==3){
			float newX = v.x*coords[0][0] + v.y*coords[0][1] + v.z*coords[0][2];
			float newY = v.x*coords[1][0] + v.y*coords[1][1] + v.z*coords[1][2];
			float newZ = v.x*coords[2][0] + v.y*coords[2][1] + v.z*coords[2][2];
			v.x = newX;
			v.y = newY;
			v.z = newZ;
			return v;
		}else if(N==4){
			float newX = v.x*coords[0][0] + v.y*coords[0][1] + v.z*coords[0][2];
			float newY = v.x*coords[1][0] + v.y*coords[1][1] + v.z*coords[1][2];
			float newZ = v.x*coords[2][0] + v.y*coords[2][1] + v.z*coords[2][2];
			v.x = newX+coords[0][3];
			v.y = newY+coords[1][3];
			v.z = newZ+coords[2][3];
			return v;
		}
		return null;
	}
	*/
	/** Apply this matrix to the vector v. A new vector is returned. */
	/*
	public Vector applyTo(Vector p){
		Vector v = p.clone(); 
		assert(N==M);
		if(N==M && N==3){
			float newX = v.x*coords[0][0] + v.y*coords[0][1] + v.z*coords[0][2];
			float newY = v.x*coords[1][0] + v.y*coords[1][1] + v.z*coords[1][2];
			float newZ = v.x*coords[2][0] + v.y*coords[2][1] + v.z*coords[2][2];
			v.x = newX;
			v.y = newY;
			v.z = newZ;
			return v;
		}else if(N==4){
			float newX = v.x*coords[0][0] + v.y*coords[0][1] + v.z*coords[0][2];
			float newY = v.x*coords[1][0] + v.y*coords[1][1] + v.z*coords[1][2];
			float newZ = v.x*coords[2][0] + v.y*coords[2][1] + v.z*coords[2][2];
			v.x = newX+coords[0][3];
			v.y = newY+coords[1][3];
			v.z = newZ+coords[2][3];
			return v;
		}
		return null;
	}
	*/

	/** Apply this matrix to another matrix. This matrix is changed and then returned */
	public Matrix applyToThis(Matrix m){
		double[][] newCoords = new double[coords.length][coords[0].length];
		for(int r=0;r<coords.length;r++){
			for(int c=0;c<coords[0].length;c++){
				newCoords[r][c] = 0; 
				for(int i=0;i<coords[0].length;i++) newCoords[r][c]+=coords[r][i]*m.coords[i][c];
			}
		}
		this.coords = newCoords;
		return this;
	}

	/** Apply this matrix to another matrix. The result is returned */
	public Matrix applyTo(Matrix m) {
		Matrix ret = new Matrix(coords.length, coords[0].length);
		double[][] newCoords = ret.coords;
		for(int r=0;r<coords.length;r++){
			for(int c=0;c<coords[0].length;c++){
				newCoords[r][c] = 0; 
				for(int i=0;i<coords[0].length;i++) newCoords[r][c]+=coords[r][i]*m.coords[i][c];
			}
		}
		return ret;
	}

	/** Add the components of two matrices */
	public Matrix addThis(Matrix m){
		for(int i=0;i<coords.length;i++) for(int j=0;j<coords[0].length;j++) coords[i][j]+=m.coords[i][j];
		return this;
	}
	
	/** Multiply the components of this matrix by a scalar. Changes and then returns this */
	public Matrix multiplyThis(double scalar){
		for(int i=0;i<coords.length;i++) for(int j=0;j<coords[0].length;j++) coords[i][j]*=scalar;
		return this;
	}

	/** Return the inverse of this matrix. Uses 31HOPs for 3x3 matrices.*/
	public Matrix invertIn(){
		if(coords.length==3 && coords[0].length==3){
			double[][] newCoords = new double[coords.length][coords[0].length];
			newCoords[0][0] = coords[1][1]*coords[2][2] - coords[1][2]*coords[2][1];
			newCoords[0][1] = coords[0][2]*coords[2][1] - coords[0][1]*coords[2][2];
			newCoords[0][2] = coords[0][1]*coords[1][2] - coords[0][2]*coords[1][1];
			
			newCoords[1][0] = coords[1][2]*coords[2][0] - coords[1][0]*coords[2][2];
			newCoords[1][1] = coords[0][0]*coords[2][2] - coords[0][2]*coords[2][0];
			newCoords[1][2] = coords[0][2]*coords[1][0] - coords[0][0]*coords[1][2];
			
			newCoords[2][0] = coords[1][0]*coords[2][1] - coords[1][1]*coords[2][0];
			newCoords[2][1] = coords[0][1]*coords[2][0] - coords[0][0]*coords[2][1];
			newCoords[2][2] = coords[0][0]*coords[1][1] - coords[0][1]*coords[1][0];//18HOPs
			
			double det = coords[0][0]*newCoords[0][0]+coords[0][1]*newCoords[1][0]+coords[0][2]*newCoords[2][0]; //3HOPs
			this.coords = newCoords;
			return multiplyThis(1/det);//10HOps
		}

		throw new Error("Inversion not implemented for matrices of size "+coords.length+"x"+coords[0].length);
	}
	public double determinant(){
		if(coords.length==3 && coords[0].length==3){
			double ret = coords[0][0]*(coords[1][1]*coords[2][2]-coords[1][2]*coords[2][1]);
			ret -= 		 coords[0][1]*(coords[1][0]*coords[2][2]-coords[1][2]*coords[2][0]);
			ret += 		 coords[0][2]*(coords[1][0]*coords[2][1]-coords[1][1]*coords[2][0]);
			return ret;
			
		}
		throw new Error("Determinant not implemented for matrices of size "+coords.length+"x"+coords[0].length);
	}

	/*
	public Matrix getEigenDecomposition(){
		double[][] array = new double[M][N];
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) array[i][j] = coords[i][j];
		Jama.Matrix A = new Jama.Matrix(array);
		EigenvalueDecomposition ed = new EigenvalueDecomposition(A);
		Jama.Matrix eigenM = ed.getV();
		Matrix ret = new Matrix(N,M);
		double[] eigenVals = ed.getRealEigenvalues();
		for(int r=0;r<N;r++){
			for(int c=0;c<M;c++){
				ret.set(r, c, (float)(eigenM.get(r, c)*eigenVals[c]));
			}
		}
		return ret;
		
	}
	
	public Vector[] getEigenvectors(){
		double[][] array = new double[M][N];
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) array[i][j] = coords[i][j];
		Jama.Matrix A = new Jama.Matrix(array);
		EigenvalueDecomposition ed = new EigenvalueDecomposition(A);
		Jama.Matrix eigenM = ed.getV();
		Vector[] ret = new Vector[eigenM.getColumnDimension()];
		double[] eigenVals = ed.getRealEigenvalues();
		for(int i=0;i<ret.length;i++){
			ret[i] = new Vector(eigenM.get(0, i), eigenM.get(1, i), eigenM.get(2, i)).timesIn((float)eigenVals[i]);
			//ret[i] = new Vector(eigenM.get(0, i), eigenM.get(1, i), eigenM.get(2, i));
		}
		return ret;
	}
	public double[] getEigenValues(){
		assert M==3&&N==3;
		double a = -1;
		double b = coords[0][0]+coords[1][1]+coords[2][2];
		double c = -coords[0][0]*coords[2][2]-
		coords[1][1]*coords[2][2]-
		coords[0][0]*coords[1][1]+
		coords[1][2]*coords[2][1]+
		coords[0][1]*coords[1][0]+
		coords[0][2]*coords[2][0];
		double d = coords[0][0]*coords[1][1]*coords[2][2]-
		coords[0][0]*coords[1][2]*coords[2][1]+
		coords[0][1]*coords[1][2]*coords[2][0]-
		coords[0][1]*coords[1][0]*coords[2][2]+
		coords[0][2]*coords[1][0]*coords[2][1]-
		coords[0][2]*coords[1][1]*coords[2][0];
		double[] cubePars = {a,b,c,d};
		double[] roots = Polynomial.solve(cubePars);

		if(Math.abs(roots[3])>0.000001){
			return new double[]{roots[0]};
		}else{
			if(Math.abs(roots[1]-roots[2])<0.000001)
				return new double[]{roots[0], roots[1]};
			else
				return new double[]{roots[0], roots[1],roots[2]};
		}
	}
*/

	public String toString(){
		StringBuilder sb = new StringBuilder();
		for(int i=0;i<coords.length;i++) {
			sb.append('|');
			for(int j=0;j<coords[0].length;j++) {
				sb.append(String.format("% 9.2f", coords[i][j]));
				sb.append(' ');
			}
			sb.append('|');
			sb.append('\n');
		}
		return sb.toString();
	}
	public static Matrix createIdentityMatrix(int n) {
		Matrix ret = new Matrix(n,n);
		for(int i=0;i<n;i++) ret.coords[i][i] = 1;
		return ret;
	}
	public double get(int i, int j) {
		return coords[i][j];
	}
	
	public Matrix clone(){
		Matrix ret = new Matrix(coords.length, coords[0].length);
		for(int r=0;r<coords.length;r++) for(int c=0;c<coords[0].length;c++)
			ret.coords[r][c] = coords[r][c];
		return ret;
	}
}
