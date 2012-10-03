package ProGAL.geom2d;

import java.awt.Color;

import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.math.Matrix;

public class ApolloniusSolver {

	/** Solves the Apollonius problem of finding a circle tangent to three other circles in the plane. 
	 * The method uses approximately 68 heavy operations (multiplication, division, square-roots).
	 * @param c1 One of the circles in the problem 
	 * @param c2 One of the circles in the problem
	 * @param c3 One of the circles in the problem
	 * @param s1 An indication if the solution should be externally or internally tangent (-1/+1) to c1
	 * @param s2 An indication if the solution should be externally or internally tangent (-1/+1) to c2
	 * @param s3 An indication if the solution should be externally or internally tangent (-1/+1) to c3
	 * @return The solution to the problem of Apollonius. 
	 */
	public static Circle solveApollonius(Circle c1, Circle c2, Circle c3, int s1, int s2, int s3){
		Point[] centers = {c1.center, c2.center, c3.center};
		double[] radii = {c1.radius, c2.radius, c3.radius};
		int[] s = {s1,s2,s3};
		return solveApollonius(centers, radii, s);
	}

	/**
	 * @hops 66
	 */
	private static Circle solveApollonius(Point[] centers, double[] radii, int[] s){
		
		//Step 1. Rewrite to linear system
		Matrix A = new Matrix(2,4);
		for(int i=0;i<2;i++){//i: row
			for(int j=0;j<2;j++){ //j: col
				A.set(  i, j, 2*(centers[i+1].get(j)-centers[0].get(j))  );//1HOp * 2 * 2 = 4
			}
			A.set(  i, 2, 2*(s[0]*radii[0]-s[i+1]*radii[i+1])  );//3HOp * 2 = 6

			double sum = 0;
			for(int j=0;j<2;j++) 
				sum+=centers[i+1].get(j)*centers[i+1].get(j) - centers[0].get(j)*centers[0].get(j);//2HOp * 2 * 2
			sum+=radii[0]*radii[0]-radii[i+1]*radii[i+1];//2HOp * 2
			A.set(i, 3, sum);
		}

		//Step 2. Simplify linear system
		A.reduceThis();// m*n+(n-1)^2*m = 2*4+3*3*2 = 26HOp
		double M = A.get(0, 3);
		double N = -A.get(0, 2);
		double P = A.get(1, 3);
		double Q = -A.get(1, 2);

		//Step3. Find tangent sphere
		//First find r_s
		double a = N*N+Q*Q-1;//2HOp
		double b = 2*(  (M-centers[0].get(0))*N + (P-centers[0].get(1))*Q + s[0]*radii[0]  );//4HOp
		double c = 
				(M-centers[0].get(0))*(M-centers[0].get(0)) + 
				(P-centers[0].get(1))*(P-centers[0].get(1)) - 
				radii[0]*radii[0]; //3HOp
		double r_s = (-b+Math.signum(a)*Math.sqrt(b*b-4*a*c))/(2*a);//7HOp
		double x_s = M+N*r_s;//1HOp
		double y_s = P+Q*r_s;//1HOp
		return new Circle(new Point(x_s, y_s), r_s);
	}

	public static void main(String[] args){
		Circle c1 = new Circle(new Point(1,0), 0.3);
		Circle c3 = new Circle(new Point(0,1), 0.5);
		Circle c2 = new Circle(new Point(1,1), 0.4);
		Circle tangent = solveApollonius(c1, c2, c3, -1, -1, -1);
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		scene.addShape(c1);
		scene.addShape(c2);
		scene.addShape(c3);
		scene.addShape(tangent, Color.green.darker(), 0,true);
	}

}
