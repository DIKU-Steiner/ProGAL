package ProGAL.math;

import java.util.Arrays;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Vector;

public class Matrix3x3 extends Matrix{

	public Matrix3x3(){
		super(3,3);
	}

	public Vector getColumn(int c){
		return new Vector(coords[0][c],coords[1][c],coords[2][c]);
	}
	public Vector getRow(int r){
		return new Vector(coords[r][0],coords[r][1],coords[r][2]);
	}

	/**
	 *  Get the determinant of this matrix.
	 *  @hops 9 
	 */
	public double determinant(){
		double ret = coords[0][0]*(coords[1][1]*coords[2][2]-coords[1][2]*coords[2][1]);
		ret -= 		 coords[0][1]*(coords[1][0]*coords[2][2]-coords[1][2]*coords[2][0]);
		ret += 		 coords[0][2]*(coords[1][0]*coords[2][1]-coords[1][1]*coords[2][0]);
		return ret;
	}

	/** Invert this matrix (overwrites this and returns it). 
	 * @hops 31*/
	public Matrix invertThis(){
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
	/**
	 * Get the eigenvectors of a 3x3 matrix sorted by decreasing eigenvalue. If the 
	 * matrix is singular ... TODO: test this.
	 * @return an array containing eigenvectors.
	 * @hops 84
	 */
	public Vector[] getEigenvectors3d(){
		/*	if(coords.length!=3 || coords[0].length!=3) 
			throw new Error("Matrix is "+coords.length+"x"+coords[0].length);

		EigenvalueDecomposition ed = new EigenvalueDecomposition();
		Vector3d[] ret = new Vector3d[3];
		double[][] V = ed.getV();
		ret[0] = new Vector3d(V[0][0],V[1][0],V[2][0]);
		ret[1] = new Vector3d(V[0][1],V[1][1],V[2][1]);
		ret[2] = new Vector3d(V[0][2],V[1][2],V[2][2]);
		return ret;
		 */

		boolean issymmetric = true;
		for (int j = 0; (j < 3) & issymmetric; j++) {
			for (int i = 0; (i < 3) & issymmetric; i++) {
				issymmetric = (coords[i][j] == coords[j][i]);
			}
		}

		if(issymmetric){//Use Cromwells equations
			double c = coords[0][0]*coords[1][1];						//Eqn 7 (1HOP)
			double d = coords[1][2]*coords[2][1];						//Eqn 8 (1HOP)
			double e = coords[0][1]*coords[1][0];						//Eqn 9 (1HOP)
			double f = coords[0][2]*coords[2][0];						//Eqn 10 (1HOP)
			double p = -coords[0][0]-coords[1][1]-coords[2][2]; 		//Eqn 11
			double q = c+(coords[0][0]+coords[1][1])*coords[2][2]-d-e-f;//Eqn 12 (1HOP)
			double r = (e-c)*coords[2][2] + d*coords[0][0] - 2*(coords[0][1]*coords[1][2]*coords[2][0])+f*coords[1][1];//Eqn 13 (6HOPs)
			double pThirds = p/3; // (1HOP)
			double a = q-p*pThirds; 									//Eqn 16 (1HOP)
			double b = (2*pThirds*pThirds-q)*pThirds+r; 				//Eqn 17 (3HOPs)
			double aThirds = a/3; //(1HOP)
			double m = 2*Math.sqrt(-aThirds);							//Eqn 19 (2HOPs)
			double t = Math.acos(b/(aThirds*m))/3; 						//Eqn ? (4HOPs)
			double cosT = Math.cos(t); //(1HOP)
			double sinT = Math.sin(t); //(1HOP)

			//Eigenvalues
			double l1 =  m*cosT - pThirds;								//Eqn 29 (1HOP)
			double l2 = -m*((cosT+Constants.SQRT3*sinT)/2) - pThirds;	//Eqn 30 (3HOPs)
			double l3 = -m*((cosT-Constants.SQRT3*sinT)/2) - pThirds;	//Eqn 31 (3HOPs)


			Matrix3x3 m1 = clone();
			m1.set(0, 0, -l1+m1.get(0, 0));
			m1.set(1, 1, -l1+m1.get(1, 1));
			m1.set(2, 2, -l1+m1.get(2, 2));
			m1.reduceThis();		//18HOps
			m1.toConsole();
			Vector v1 = new Vector(-m1.coords[0][2], -m1.coords[1][2], 1);

			Vector v2=null, v3=null;
			if(l2>Constants.EPSILON){
				Matrix3x3 m2 = clone();
				m2.set(0, 0, -l2+m2.get(0, 0));
				m2.set(1, 1, -l2+m2.get(1, 1));
				m2.set(2, 2, -l2+m2.get(2, 2));
				m2.reduceThis();		//18HOPs
				m2.toConsole();
				v2 = new Vector(-m2.coords[0][2], -m2.coords[1][2], 1);
			}

			if(l3>Constants.EPSILON){
				Matrix3x3 m3 = clone();
				m3.set(0, 0, -l3+m3.get(0, 0));
				m3.set(1, 1, -l3+m3.get(1, 1));
				m3.set(2, 2, -l3+m3.get(2, 2));
				m3.reduceThis();//18HOPs
				m3.toConsole();
				v3 = new Vector(-m3.coords[0][2], -m3.coords[1][2], 1);
			}
		}

		return null;
	}

	public static void main(String[] args){
		PointList points = PointList.genPointsInCube(9);
		for(Point p: points) { for(int i=1;i<3;i++) p.set(i, p.get(i)*0.3); }
		Matrix3x3 m = points.getCovariance();
		System.out.println("m: ");
		m.toConsole();

		m.getEigenvectors3d();

		m.getEigenvalueDecomposition().getV().toConsole();
	}

	public Matrix3x3 clone(){
		Matrix3x3 ret = new Matrix3x3();
		for(int r=0;r<coords.length;r++) for(int c=0;c<coords[0].length;c++)
			ret.coords[r][c] = coords[r][c];
		return ret;
	}
}
