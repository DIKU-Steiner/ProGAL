package ProGAL.math;

import java.util.Arrays;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;

public class Matrix3x3 extends Matrix{

	public Matrix3x3(){
		super(3,3);
	}

	public Matrix3x3(double[][] coords){
		super(coords);
		if(M!=3 || N!=3) throw new RuntimeException("Dimensions dont fit");
	}

	public Vector getColumn(int c){
		return new Vector(coords[0][c],coords[1][c],coords[2][c]);
	}
	public Vector getRow(int r){
		return new Vector(coords[r][0],coords[r][1],coords[r][2]);
	}

	public Matrix3x3 getTranspose(){
		Matrix3x3 ret = clone();
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				double tmp = coords[i][j];
				coords[i][j] = coords[j][i];
				coords[j][i] = tmp;
			}
		}
		return ret;
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
	public Vector[] getEigenvectors(){
		if(coords.length!=3 || coords[0].length!=3) 
			throw new Error("Matrix is "+coords.length+"x"+coords[0].length);

		//		EigenvalueDecomposition ed = new EigenvalueDecomposition();
		//		Vector[] ret = new Vector[3];
		//		double[][] V = ed.getV().coords;
		//		ret[0] = new Vector(V[0][0],V[1][0],V[2][0]).multiplyThis(ed.getRealEigenvalues()[0]);
		//		ret[1] = new Vector(V[0][1],V[1][1],V[2][1]).multiplyThis(ed.getRealEigenvalues()[1]);
		//		ret[2] = new Vector(V[0][2],V[1][2],V[2][2]).multiplyThis(ed.getRealEigenvalues()[2]);
		//		if(true) return ret;


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
				//				m2.toConsole();
				v2 = new Vector(-m2.coords[0][2], -m2.coords[1][2], 1);
			}

			if(l3>Constants.EPSILON){
				Matrix3x3 m3 = clone();
				m3.set(0, 0, -l3+m3.get(0, 0));
				m3.set(1, 1, -l3+m3.get(1, 1));
				m3.set(2, 2, -l3+m3.get(2, 2));
				m3.reduceThis();//18HOPs
				//				m3.toConsole();
				v3 = new Vector(-m3.coords[0][2], -m3.coords[1][2], 1);
			}
			return new Vector[]{v1,v2,v3};
		}

		EigenvalueDecomposition ed = new EigenvalueDecomposition();
		Vector[] ret = new Vector[3];
		double[][] V = ed.getV().coords;
		ret[0] = new Vector(V[0][0],V[1][0],V[2][0]).multiplyThis(ed.getRealEigenvalues()[0]);
		ret[1] = new Vector(V[0][1],V[1][1],V[2][1]).multiplyThis(ed.getRealEigenvalues()[1]);
		ret[2] = new Vector(V[0][2],V[1][2],V[2][2]).multiplyThis(ed.getRealEigenvalues()[2]);
		return ret;
	}

	/* Multiply this matrix with matrix M */
	public Matrix3x3 multiply(Matrix3x3 M) {
		double[] r0 = {coords[0][0]*M.coords[0][0]+coords[0][1]*M.coords[1][0]+coords[0][2]*M.coords[2][0],
				coords[0][0]*M.coords[0][1]+coords[0][1]*M.coords[1][1]+coords[0][2]*M.coords[2][1],
				coords[0][0]*M.coords[0][2]+coords[0][1]*M.coords[1][2]+coords[0][2]*M.coords[2][2]};
		double[] r1 = {coords[1][0]*M.coords[0][0]+coords[1][1]*M.coords[1][0]+coords[1][2]*M.coords[2][0],
				coords[1][0]*M.coords[0][1]+coords[1][1]*M.coords[1][1]+coords[1][2]*M.coords[2][1],
				coords[1][0]*M.coords[0][2]+coords[1][1]*M.coords[1][2]+coords[1][2]*M.coords[2][2]};
		double[] r2 = {coords[2][0]*M.coords[0][0]+coords[2][1]*M.coords[1][0]+coords[2][2]*M.coords[2][0],
				coords[2][0]*M.coords[0][1]+coords[2][1]*M.coords[1][1]+coords[2][2]*M.coords[2][1],
				coords[2][0]*M.coords[0][2]+coords[2][1]*M.coords[1][2]+coords[2][2]*M.coords[2][2]};
		double[][] mat = {r0, r1, r2};
		return new Matrix3x3(mat);
	}

	/**
	 * This routine maps three values (x[0], x[1], x[2]) in the range [0,1]
	 * into a 3x3 rotation matrix, M. Uniformly distributed random variables
	 * x0, x1, and x2 create uniformly distributed random rotation matrices.
	 * To create small uniformly distributed "perturbations", supply
	 * samples in the following ranges
	 * <pre>
	 *     x[0] in [ 0, d ]
	 *     x[1] in [ 0, 1 ]
	 *     x[2] in [ 0, d ]
	 * </pre>
	 * where 0 < d < 1 controls the size of the perturbation.  Any of the
	 * random variables may be stratified (or "jittered") for a slightly more
	 * even distribution.
	 *                                                       *
	 * @author Jim Arvo, 1991
	 * @url https://github.com/erich666/GraphicsGems/blob/master/gemsiii/rand_rotation.c
	 */
	public static Matrix3x3 randRotation(double[] x){
		double theta = x[0] * Constants.PI * 2.0; 	/* Rotation about the pole (Z).      */
		double phi   = x[1] * Constants.PI * 2.0; 	/* For direction of pole deflection. */
		double z     = x[2] * 2.0;        			/* For magnitude of pole deflection. */

		/* Compute a vector V used for distributing points over the sphere  */
		/* via the reflection I - V Transpose(V).  This formulation of V    */
		/* will guarantee that if x[1] and x[2] are uniformly distributed,  */
		/* the reflected points will be uniform on the sphere.  Note that V */
		/* has length sqrt(2) to eliminate the 2 in the Householder matrix. */

		double r  = Math.sqrt( z );
		double Vx = Math.sin( phi ) * r;
		double Vy = Math.cos( phi ) * r;
		double Vz = Math.sqrt( 2.0 - z );

		/* Compute the row vector S = Transpose(V) * R, where R is a simple */
		/* rotation by theta about the z-axis.  No need to compute Sz since */
		/* it's just Vz.                                                    */

		double st = Math.sin( theta );
		double ct = Math.cos( theta );
		double Sx = Vx * ct - Vy * st;
		double Sy = Vx * st + Vy * ct;

		/* Construct the rotation matrix  ( V Transpose(V) - I ) R, which   */
		/* is equivalent to V S - R.                                        */

		Matrix3x3 ret = new Matrix3x3();
		ret.coords[0][0] = -(Vx * Sx - ct);
		ret.coords[0][1] = -(Vx * Sy - st);
		ret.coords[0][2] = -(Vx * Vz);
		
		ret.coords[1][0] = -(Vy * Sx + st);
		ret.coords[1][1] = -(Vy * Sy - ct);
		ret.coords[1][2] = -(Vy * Vz);
		
		ret.coords[2][0] = Vz * Sx;
		ret.coords[2][1] = Vz * Sy;
		ret.coords[2][2] = 1.0 - z;
		return ret;
	}
	
	public static Matrix3x3 randRotation2(double[] x){
		double ax = x[0] * Math.PI * 2;
		double ay = x[1] * Math.PI * 2;
		double az = x[2] * Math.PI * 2;

		Matrix3x3 mx = (Matrix3x3)Matrix3x3.createRotationMatrix(ax, Vector.X);
		Matrix3x3 my = (Matrix3x3)Matrix3x3.createRotationMatrix(ay, Vector.Y);
		Matrix3x3 mz = (Matrix3x3)Matrix3x3.createRotationMatrix(az, Vector.Z);
		return mz.multiply(my).multiply(mx);
	}
	public static Matrix3x3 randRotation3(double scl){
		double theta =  Math.acos(Randomization.randBetween(-1.0, 1.0));
		double phi = Randomization.randBetween(0, 2*Constants.PI);
		Vector rax= new Vector( 
				Math.sin(theta)*Math.cos(phi),
				Math.sin(theta)*Math.sin(phi),
				Math.cos(theta)
				);
		double angle = Randomization.randBetween(0.0, scl);
		return (Matrix3x3)Matrix3x3.createRotationMatrix(angle, rax); 
	}

	public static void main(String[] args){
		
		Point[] ps = new Point[]{
				new Point(1,0,0), 
				new Point(0,1,0), 
				new Point(0,0,1), 
//				new Point(1/Math.sqrt(3),1/Math.sqrt(3),1/Math.sqrt(3)) 
				};
		Vector x = new Vector(0.04,0,0);
		Vector y = new Vector(0,0.04,0);
		Vector z = new Vector(0,0,0.04);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
//		scene.addShape(new Sphere(new Point(0,0,0), 1), new java.awt.Color(100,100,50,100), 64);
		scene.setAxisEnabled(true);
		
		for(int i=0;i<10000;i++){
//			double d1 = 0.05;
//			double d2 = 0.05;
//			Matrix3x3 m = randRotation(new double[]{
//					Randomization.randBetween(-d1/2,d1/2),
//					Randomization.randBetween(0.0,1.0),
//					Randomization.randBetween(0.0,d2)
//			});
//			double[] rand = new double[]{1,1,1};
//			while( Math.sqrt(rand[0]*rand[0]+ rand[1]*rand[1]+ rand[2]*rand[2])>d1 )
//			{
//				rand[0] = Randomization.randBetween(-d1,d1);
//				rand[1] = Randomization.randBetween(-d1,d1);
//				rand[2] = Randomization.randBetween(-d1,d1);
//			}
//			Matrix3x3 m = randRotation2(rand);
			Matrix3x3 m = randRotation3(0.1*Math.PI);
			//			System.out.println(m);
//			for(Point p: ps){
//				Point p_r = m.multiply(p);
//				Vector x_r = m.multiply(x);
//				Vector y_r = m.multiply(y);
//				Vector z_r = m.multiply(z);
//				scene.addShape(new Sphere(p, 0.05), java.awt.Color.BLACK);
//				scene.addShape(new Sphere(p_r, 0.01), java.awt.Color.RED, 3);
//				scene.addShape(new LSS(p_r,p_r.add(x_r), 0.004), java.awt.Color.RED, 3);
//				scene.addShape(new LSS(p_r,p_r.add(y_r), 0.004), java.awt.Color.RED, 3);
//				scene.addShape(new LSS(p_r,p_r.add(z_r), 0.004), java.awt.Color.RED, 3);
//			}
			
			
			double theta =  Math.acos(Randomization.randBetween(-1.0, 1.0));
			double phi = Randomization.randBetween(0, 2*Constants.PI);
			Point rax= new Point( 
					Math.sin(theta)*Math.cos(phi),
					Math.sin(theta)*Math.sin(phi),
					Math.cos(theta)
					);
			scene.addShape(new Sphere(rax, 0.01), java.awt.Color.BLACK, 4);
		}
		
	}

	public Matrix3x3 clone(){
		Matrix3x3 ret = new Matrix3x3();
		for(int r=0;r<coords.length;r++) for(int c=0;c<coords[0].length;c++)
			ret.coords[r][c] = coords[r][c];
		return ret;
	}
}
