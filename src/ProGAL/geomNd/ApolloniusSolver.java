package ProGAL.geomNd;

import java.awt.Color;

import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.math.Matrix;

public class ApolloniusSolver {
	public static void main(String args[]){
		HyperSphere p0 = new HyperSphere(new Point(new double[]{-5, 7.37,  5,0}),10);
		HyperSphere p1 = new HyperSphere(new Point(new double[]{ 5, 7.37, -5,0}), 1.52);
		HyperSphere p2 = new HyperSphere(new Point(new double[]{-5, 2.66, -5,0}), 1.7);
		HyperSphere p3 = new HyperSphere(new Point(new double[]{ 5,-2.66,  5,0}), 1.2);
		HyperSphere p4 = new HyperSphere(new Point(new double[]{ 0,0,0, 20}), 0.4);
//		Sphere p0 = new Sphere(new Point(1, 0, 0), 0.4);
//		Sphere p1 = new Sphere(new Point(0, 1, 0), 0.52);
//		Sphere p2 = new Sphere(new Point(0, 0, 1), 0.7);
//		Sphere p3 = new Sphere(new Point(1, 1, 1), 0.2);
		HyperSphere tangent = solveApollonius(new HyperSphere[]{p0,p1,p2,p3,p4}, new int[]{1, 1, 1, 1,1});
		System.out.println(tangent);
		
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		scene.addShape(new ProGAL.geom3d.volumes.Sphere(new ProGAL.geom3d.Point(p0.getCenter().coords), p0.getRadius()), Color.WHITE.darker(),50);
		scene.addShape(new ProGAL.geom3d.volumes.Sphere(new ProGAL.geom3d.Point(p1.getCenter().coords), p1.getRadius()), Color.WHITE.darker(),50);
		scene.addShape(new ProGAL.geom3d.volumes.Sphere(new ProGAL.geom3d.Point(p2.getCenter().coords), p2.getRadius()), Color.WHITE.darker(),50);
		scene.addShape(new ProGAL.geom3d.volumes.Sphere(new ProGAL.geom3d.Point(p3.getCenter().coords), p3.getRadius()), Color.WHITE.darker(),50);
		scene.addShape(new ProGAL.geom3d.volumes.Sphere(new ProGAL.geom3d.Point(tangent.getCenter().coords), tangent.getRadius()), new Color(0,100,0,100),50);
//		scene.addShape(tangent, new Color(0,100,0,100),50);
	}
	/** Solves the Apollonius problem of finding a circle tangent to three other circles in the plane. 
	 * @param c0 One of the spheres in the problem
	 * @param c1 One of the spheres in the problem 
	 * @param c2 One of the spheres in the problem
	 * @param c3 One of the spheres in the problem
	 * @param s0 An indication if the solution should be externally or internally tangent (-1/+1) to c0
	 * @param s1 An indication if the solution should be externally or internally tangent (-1/+1) to c1
	 * @param s2 An indication if the solution should be externally or internally tangent (-1/+1) to c2
	 * @param s3 An indication if the solution should be externally or internally tangent (-1/+1) to c3
	 * @return The solution to the problem of Apollonius. 
	 * @hops d*d+3*d+2*d+2*d*d+2+d^3-d^2+3d+3d+8+d = d^3+2d^2+12d+10
	 */
	public static HyperSphere solveApollonius(HyperSphere[] spheres, int[] s){
		//The method is described in detail in BioRepo/ProGAL/Doc/ProblemOfApollonius/paper.tex
		if(s.length!=spheres.length) throw new RuntimeException("Dimension problem");
		int d = spheres.length-1;
		Point[] centers = new Point[d+1];
		double[] radii = new double[d+1];
		for(int i=0;i<spheres.length;i++){
			centers[i] = spheres[i].getCenter();
			if(centers[i].getDimensions()!=d) throw new RuntimeException("Dimension problem");
			radii[i] = spheres[i].getRadius();
		}
		
		//Step 1. Rewrite to linear system
		Matrix A = new Matrix(d,d+2);
		for(int i=0;i<d;i++){//i: row
			for(int j=0;j<d;j++){ //j: col
				A.set(  i, j, 2*(centers[i+1].get(j)-centers[0].get(j))  );//1HOp * d * d
			}
			A.set(  i, d, 2*(s[0]*radii[0]-s[i+1]*radii[i+1])  );//3HOp * d
			
			double sum = 0;
			for(int j=0;j<d;j++) 
				sum+=centers[i+1].get(j)*centers[i+1].get(j) - centers[0].get(j)*centers[0].get(j);//2HOp * d * d
			sum+=radii[0]*radii[0]-radii[i+1]*radii[i+1];//2HOp * d
			A.set(i, d+1, sum);
		}

		//Step 2. Simplify linear system
		A.reduceThis();// m*n+(n-1)^2*m = d*(d+2)+(d-1)*(d-1)*d = d^3-d^2+3d HOp
		
		//Step3. Find tangent sphere
		//First find r_s
		double a = -1;
		double b = s[0]*radii[0];			//1HOp
		double c = -radii[0]*radii[0];		//1HOp
		for(int i=0;i<d;i++) {
			a+= A.get(i, d)*A.get(i, d); 	//d HOp
			double delta = (A.get(i, d+1)-centers[0].get(i)); 
			b-= delta*A.get(i,d); 			//d HOp
			c+= delta*delta;				//d HOp
		}
		b*=2;
		
		double D = (b*b-4*a*c);		//3HOp
		if( D<0 ) throw new RuntimeException("No solution");
		double Drt = Math.sqrt(D);	//1HOp
		double r_s = (-b+Drt)/(2*a);//2HOp
		r_s = Math.max((-b-Drt)/(2*a), r_s);//2HOp
		double[] coords = new double[d];
		for(int i=0;i<d;i++) coords[i] = A.get(i, d+1)-A.get(i, d)*r_s;	//dHOp
		
		return new HyperSphere(new Point(coords), r_s);
	}
	
}
