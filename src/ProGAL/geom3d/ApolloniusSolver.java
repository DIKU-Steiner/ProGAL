package ProGAL.geom3d;

import java.awt.Color;

import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Matrix;

public class ApolloniusSolver {
	public static void main(String[] args){
		Sphere p0 = new Sphere(new Point(-5.32,-0.92,10.16),1.6);
		Sphere p1 = new Sphere(new Point(-6.41,-1.18,14.34),1.6);
		Sphere p2 = new Sphere(new Point(-4.29, 3.57,10.95),1.6);
		Sphere p3 = new Sphere(new Point(-8.25, 2.21,11.52),1.6);
		Sphere tangent = solveApollonius(p0, p1, p2, p3, 1, 1, 1, 1);
		System.out.println(tangent);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		scene.addShape(p0, Color.WHITE,50);
		scene.addShape(p1, Color.WHITE,50);
		scene.addShape(p2, Color.WHITE,50);
		scene.addShape(p3, Color.WHITE,50);
		scene.addShape(tangent, Color.GREEN.darker(),50);
	}
	/** Solves the Apollonius problem of finding a circle tangent to three other circles in the plane. 
	 * @param c0 One of the spheres in the problem
	 * @param c1 One of the spheres in the problem 
	 * @param c2 One of the spheres in the problem
	 * @param c3 One of the spheres in the problem
	 * @param s0 An indication if the solution should be externally or internally tangent (+1/-1) to c0
	 * @param s1 An indication if the solution should be externally or internally tangent (+1/-1) to c1
	 * @param s2 An indication if the solution should be externally or internally tangent (+1/-1) to c2
	 * @param s3 An indication if the solution should be externally or internally tangent (+1/-1) to c3
	 * @return The solution to the problem of Apollonius. 
	 * @hops 126
	 */
	public static Sphere solveApollonius(Sphere c0, Sphere c1, Sphere c2, Sphere c3, int s0, int s1, int s2, int s3){
		Point[] centers = {c0.getCenter(), c1.getCenter(), c2.getCenter(), c3.getCenter()};
		double[] radii = {c0.getRadius(), c1.getRadius(), c2.getRadius(), c3.getRadius()};
		return solveApollonius(centers, radii, s0, s1, s2, s3);
	}
	
	public static Sphere solveApollonius(Point[] centers, double[] radii, int s0, int s1, int s2, int s3){
		//The method is described in detail in BioRepo/ProGAL/Doc/ProblemOfApollonius/paper.tex
		int[] s = {s0,s1,s2,s3};
		
		//Step 1. Rewrite to linear system
		Matrix A = new Matrix(3,5);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				A.set(  i, j, 2*(centers[i+1].get(j)-centers[0].get(j))  );//1HOp * 3 * 3
			}
			A.set(  i, 3, 2*(s[0]*radii[0]-s[i+1]*radii[i+1])  );//3HOp * 3
			
			double sum = 0;
			for(int j=0;j<3;j++) sum+=centers[i+1].get(j)*centers[i+1].get(j) - centers[0].get(j)*centers[0].get(j);//2HOp * 3 * 3
			sum+=radii[0]*radii[0]-radii[i+1]*radii[i+1];//2HOp * 3
			A.set(i, 4, sum);
		}

		//Step 2. Simplify linear system
		A.reduceThis();// m*n+(n-1)^2*m = 3*5+4*4*3 = 63HOp
		double M = A.get(0, 4);
		double N = -A.get(0, 3);
		double P = A.get(1, 4);
		double Q = -A.get(1, 3);
		double R = A.get(2, 4);
		double S = -A.get(2, 3);
		
		//Step3. Find tangent sphere
		//First find r_s
		double a = N*N+Q*Q+S*S-1;//3HOp
		double b = 2*(  (M-centers[0].get(0))*N + (P-centers[0].get(1))*Q + (R-centers[0].get(2))*S + s[0]*radii[0]  );//5HOp
		double c = 
			(M-centers[0].get(0))*(M-centers[0].get(0)) + 
			(P-centers[0].get(1))*(P-centers[0].get(1)) + 
			(R-centers[0].get(2))*(R-centers[0].get(2)) -
			radii[0]*radii[0]; //4HOp
		double r_s = -(-b+Math.sqrt(b*b-4*a*c))/(2*a);//6HOp
		double x_s = M+N*r_s;//1HOp
		double y_s = P+Q*r_s;//1HOp
		double z_s = R+S*r_s;//1HOp
		
		return new Sphere(new Point(x_s,y_s,z_s), r_s);
	}
	
}
