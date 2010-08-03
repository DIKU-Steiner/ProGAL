package ProGAL.math;

import static java.lang.Math.sqrt;
import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.cos;
import static java.lang.Math.acos;
import static java.lang.Math.PI;

/** 
 * A representation of a variable sized polynomial. To find the roots of the 
 * parabola f(x) = 2x^2 + 3x - 4, for instance, write
 * <code>
 * 	Polynomial parabola = new Polynomial(new double[]{2,3,4});
 *  double[] roots = parabola.calcRoots();
 * </code>
 * or simply
 * <code>
 * 	double[] roots = Polynomial.calcRoots(2,3,4);
 * </code>
 * Currently this only works for polynomials of degree 2 and 3.
 * @author rfonseca
 */
public class Polynomial {
	public final double[] parameters;
	
	Polynomial(double[] pars){
		parameters = pars;
	}
	
	/**
	 * Find the roots of the polynomial
	 * @return an array containing the roots of the polynomial
	 */
	public double[] calcRoots(){
		return calcRoots(parameters);
	}
	
	/**
	 * Find the roots of the polynomial specified by the parameters
	 * @param parameters an array representing the coefficents of the polynomial
	 * @return an array containing the roots of the polynomial
	 */
	public static double[] calcRoots(double[] parameters){
		if(parameters.length==3) return solveSecondDegree(parameters);
		if(parameters.length==4) return solveThirdDegree(parameters);
		
		return null;
	}

	/**
	 * Find the roots of a second degree polynomial.
	 * @param a coefficient of the squared term
	 * @param b coefficient of the linear term
	 * @param c coefficient of the constant
	 * @return an array containing the roots of the polynomial
	 */
	public static double[] calcRoots(double a, double b, double c){
		return solveSecondDegree(new double[]{a,b,c});
	}

	/**
	 * Find the roots of a third degree polynomial.
	 * @param a coefficient of the cubed term
	 * @param b coefficient of the squared term
	 * @param c coefficient of the linear term
	 * @param d coefficient of the constant term
	 * @return an array containing the roots of the polynomial
	 */
	public static double[] calcRoots(double a, double b, double c, double d){
		return solveThirdDegree(new double[]{a,b,c,d});
	}
	
	/** 
	 * Solves quadratic equation. Uses 3 to 7 HOps.  
	 */
	private static double[] solveSecondDegree(double[] parameters){
		double a = parameters[0];
		double b = parameters[1];
		double c = parameters[2];
		double D = b*b-4*a*c;
		if(D<0) return new double[]{};
		if(D==0) return new double[]{-b/(2*a)};
		else {
			double aa = 2*a;
			double sqD = Math.sqrt(D);
			return new double[]{(-b-sqD)/aa, (-b+sqD)/aa};
		}
	}
	
	private static double[] solveThirdDegree(double[] parameters){
		double a = parameters[1]/parameters[0];
		double b = parameters[2]/parameters[0];
		double c = parameters[3]/parameters[0];
		
		double Q = (3*b-a*a)/9;
		double R = (9*a*b-27*c-2*a*a*a)/54;
	    //Polynomial discriminant
		double D = Q*Q*Q + R*R;
		
		double root1, root2, root3, im;
		
		if(D>=0){		//Complex or duplicate roots
			double S = sign(R+sqrt(D))*pow(abs(R+sqrt(D)),1/3d);
			double T = sign(R-sqrt(D))*pow(abs(R-sqrt(D)),1/3d);
			
			root1 = -(a/3) + S+T;
			root2 = -(a/3) - (S+T)/2;
			root3 = -(a/3) - (S+T)/2;
			im = abs(sqrt(3)*(S-T)/2);
		}else{
			double th = acos(R/sqrt(-Q*Q*Q));
	        
			root1 = 2*sqrt(-Q)*cos(th/3)-a/3;
			root2 = 2*sqrt(-Q)*cos((th+2*PI)/3)-a/3;
			root3 = 2*sqrt(-Q)*cos((th+4*PI)/3)-a/3;
			im = 0;
		}
		
		return new double[]{root1,root2, root3, im};
	}
	
	private static double sign(double n){
		if(n<0) return -1;
		else return 1;
	}
	
}
