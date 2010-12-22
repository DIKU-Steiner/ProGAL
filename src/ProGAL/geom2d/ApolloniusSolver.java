package ProGAL.geom2d;

import ProGAL.math.Polynomial;

public class ApolloniusSolver {

	/** Solves the Apollonius problem of finding a circle tangent to three other circles in the plane. 
	 * The method uses approximately 68 heavy operations (multiplication, division, square-roots).
	 * @param c1 One of the circles in the problem 
	 * @param c2 One of the circles in the problem
	 * @param c3 One of the circles in the problem
	 * @param s1 An indication if the solution should be externally or internally tangent (+1/-1) to c1
	 * @param s2 An indication if the solution should be externally or internally tangent (+1/-1) to c2
	 * @param s3 An indication if the solution should be externally or internally tangent (+1/-1) to c3
	 * @return The solution to the problem of Apollonius. 
	 * @hops 68
	 */
	public static Circle solveApollonius(Circle c1, Circle c2, Circle c3, int s1, int s2, int s3){
		double x1 = c1.center.getX();
		double y1 = c1.center.getY();
		double r1 = c1.radius;
		double x2 = c2.center.getX();
		double y2 = c2.center.getY();
		double r2 = c2.radius;
		double x3 = c3.center.getX();
		double y3 = c3.center.getY();
		double r3 = c3.radius;

		double v11 = 2*x2 - 2*x1;
		double v12 = 2*y2 - 2*y1;
		double v13 = x1*x1 - x2*x2 + y1*y1 - y2*y2 - r1*r1 + r2*r2;
		double v14 = 2*s2*r2 - 2*s1*r1;
		//14HOps

		double v21 = 2*x3 - 2*x2;
		double v22 = 2*y3 - 2*y2;
		double v23 = x2*x2 - x3*x3 + y2*y2 - y3*y3 - r2*r2 + r3*r3;
		double v24 = 2*s3*r3 - 2*s2*r2;
		//28HOps
		
		double w12 = v12/v11;
		double w13 = v13/v11;
		double w14 = v14/v11;
		//31HOps
		
		double w22 = v22/v21-w12;
		double w23 = v23/v21-w13;
		double w24 = v24/v21-w14;
		//34HOps
		
		double P = -w23/w22;
		double Q = w24/w22;
		double M = -w12*P-w13;
		double N = w14 - w12*Q;
		//38HOps
		
		double a = N*N + Q*Q - 1;
		double b = 2*M*N - 2*N*x1 + 2*P*Q - 2*Q*y1 + 2*s1*r1;
		double c = x1*x1 + M*M - 2*M*x1 + P*P + y1*y1 - 2*P*y1 - r1*r1;
		//59HOps
		
		double[] quadSols = Polynomial.calcRoots(a,b,c); //7 Hops
		double rs = quadSols[0];
		double xs = M+N*rs;
		double ys = P+Q*rs;
		//68HOps
		
		return new Circle(new Point(xs,ys), rs);
	}
	
}
