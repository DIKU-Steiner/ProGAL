package ProGAL.math;

import static java.lang.Math.sqrt;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.acos;
import static java.lang.Math.PI;

import java.util.Arrays;

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
	 * Find the roots of the polynomial.  The imaginary part is not returned.
	 * @return an array containing the roots of the polynomial
	 */
	public double[] calcRoots(){
		return calcRoots(parameters);
	}
	
	/**
	 * Find the roots of the polynomial specified by the parameters. The imaginary part is not returned.
	 * @param parameters an array representing the coefficents of the polynomial
	 * @return an array containing the roots of the polynomial
	 */
	public static double[] calcRoots(double[] parameters){
		if(parameters.length==3) return solveSecondDegree(parameters);
		if(parameters.length==4) return solveThirdDegree(parameters);
		
		return null;
	}

	/**
	 * Find the roots of a second degree polynomial. The length of the returned array depends on 
	 * how many roots the quadratic equation has. In other words: The imaginary part is removed completely.
	 * @param a coefficient of the squared term
	 * @param b coefficient of the linear term
	 * @param c coefficient of the constant
	 * @return an array containing the roots of the polynomial
	 */
	public static double[] calcRoots(double a, double b, double c){
		return solveSecondDegree(new double[]{a,b,c});
	}

	/**
	 * Find the roots of a third degree polynomial. The length of the returned array depends on 
	 * how many roots the cubic equation has. In other words: The imaginary part is removed completely. 
	 * @param a coefficient of the cubed term
	 * @param b coefficient of the squared term
	 * @param c coefficient of the linear term
	 * @param d coefficient of the constant term
	 * @return an array containing the roots of the polynomial
	 */
	public static double[] calcRoots(double a, double b, double c, double d){
		double[] roots = solveThirdDegree(new double[]{a,b,c,d});
		if(roots.length<4) return roots;
		if(roots[3]>Constants.EPSILON) return new double[]{roots[0]};
		int l = 1;
		for(int i=1;i<3;i++) if(!Double.isNaN(roots[i]) && roots[i-1]!=roots[i]) l++;
		double[] ret = new double[l];
		l=0;
		for(int i=0;i<3;i++) if(!Double.isNaN(roots[i]) && (i==0||roots[i-1]!=roots[i])) ret[l++] = roots[i];
		Arrays.sort(ret);
		return ret;
	}
	
	/** 
	 * Solves quadratic equation. 
	 * @hops 4 to 5.  
	 */
	private static double[] solveSecondDegree(double[] parameters){
		double a = parameters[0];
		double b = parameters[1];
		double c = parameters[2];
		double p = -0.5*b/a;	//2HOps
		double D = p*p-c/a;		//2HOps
//		double D = b*b-4*a*c;
		if(D<0) return new double[]{};
//		if(D==0) return new double[]{-b/(2*a)};
		if(D==0) return new double[]{p};
		else {
//			double aa = 2*a;
			double sqD = Math.sqrt(D);		//1HOp
//			return new double[]{(-b-sqD)/aa, (-b+sqD)/aa};
			return new double[]{p-sqD, p+sqD};
		}
	}
	
	/**
	 * Solves qubic equation.
	 * @hops 3 to 35
	 */
	private static double[] solveThirdDegree(double[] parameters){
		if(parameters[0]==0) return solveSecondDegree(new double[]{parameters[1],parameters[2],parameters[3]});
		double a = parameters[1]/parameters[0];
		double b = parameters[2]/parameters[0];
		double c = parameters[3]/parameters[0];
		
		double Q = (3*b-a*a)/9;
		double R = (9*a*b-27*c-2*a*a*a)/54;
	    //Polynomial discriminant
		double D = Q*Q*Q + R*R;
		double aThirds = a/3;//17HOps so far
		
		double root1, root2, root3, im;
		
		if(D>=0){		//Complex or duplicate roots
			double sqrtD = sqrt(D);
			double S = sign(R+sqrtD)*Math.cbrt(abs(R+sqrtD));
			double T = sign(R-sqrtD)*Math.cbrt(abs(R-sqrtD));//22HOps so far
			
			root1 = -aThirds + S+T;
			root2 = -aThirds - (S+T)/2;
			root3 = root2;
			im = abs(Constants.SQRT3*(S-T)/2);//24HOps so far
		}else {
			double sqrtMQ = sqrt(-Q);
			double th = acos(R/(-Q*sqrtMQ)); //acos(R/sqrt(-Q*Q*Q));//21HOps so far

			root1 = 2*sqrtMQ*cos(th/3)-aThirds;
			root2 = 2*sqrtMQ*cos((th+2*PI)/3)-aThirds;
			root3 = 2*sqrtMQ*cos((th+4*PI)/3)-aThirds;//35HOps so far
			im = 0;
		}
		return new double[]{root1,root2, root3, im};
	}

	private static double sign(double n){	return n<0?-1:1;	}
	
	public static double[]ÊsolveHighDegree(double[] a) {
		double[] der = new double[a.length-1];
		for (int j = 1; j < a.length; j++) {
			der[j] = (j+1)*a[]
		}
		double tolerance = 10^-8;
		int maxIterations = 10^4;
		
		for (int i = 0; i < maxIterations; i++) {
	}
	
}
