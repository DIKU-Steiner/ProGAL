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
	public final double[] coeff;
	private int deg;
	
	Polynomial(double[] pars){
		coeff = pars;
		deg = pars.length-1;
	}
	
    // returns a + b
    public Polynomial plus(Polynomial b) {
    	double[] coeff = new double[Math.max(deg, b.deg)];
        for (int i = 0; i <= deg; i++) coeff[i] = this.coeff[i];
        for (int i = 0; i <= b.deg; i++) coeff[i] += b.coeff[i];
        while ((deg != 0) && (coeff[deg] == 0.0)) deg--;
        return new Polynomial(coeff);
    }

    // returns a - b
    public Polynomial minus(Polynomial b) {
    	double[] coeff = new double[Math.max(deg, b.deg)];
        for (int i = 0; i <= deg; i++) coeff[i] = this.coeff[i];
        for (int i = 0; i <= b.deg; i++) coeff[i] -= b.coeff[i];
        while ((deg != 0) && (coeff[deg] == 0.0)) deg--;
        return new Polynomial(coeff);
    }

    // returns a * b
    public Polynomial times(Polynomial b) {
    	double[] coeff = new double[deg + b.deg];
        for (int i = 0; i <= deg; i++)
            for (int j = 0; j <= b.deg; j++) coeff[i+j] += this.coeff[i] * b.coeff[j];
        return new Polynomial(coeff);
    }

    // returns a(b(x)) 
    public Polynomial compose(Polynomial b) {
    	double[] coeff = new double[deg + b.deg];
        for (int i = deg; i >= 0; i--) {
        	for (int j = b.deg; j >= 0; j--) coeff[i+j] += this.coeff[i]*b.coeff[j];
        }
        return new Polynomial(coeff);
    }

    
    // do a and b represent the same polynomial?
    public boolean equal(Polynomial b) {
        if (deg != b.deg) return false;
        for (int i = deg; i >= 0; i--)
            if (coeff[i] != b.coeff[i]) return false;
        return true;
    }

    // use Horner's method to compute and return the polynomial evaluated at x
    public double evaluate(int x) {
        double p = 0;
        for (int i = deg; i >= 0; i--)
            p = coeff[i] + (x * p);
        return p;
    }

    // differentiate this polynomial and return it
    public Polynomial differentiate() {
    	if ( deg == 0) return new Polynomial(new double[1]);
    	double[] par = new double[deg];
    	for (int i = 1; i < deg; i++) {
    		par[i-1] = i * coeff[i];
    	}
    	return new Polynomial(par);
    }

    // convert to string representation
    public String toString() {
        if (deg ==  0) return "" + coeff[0];
        if (deg ==  1) return coeff[1] + "x + " + coeff[0];
        String s = coeff[deg] + "x^" + deg;
        for (int i = deg-1; i >= 0; i--) {
            if      (coeff[i] == 0) continue;
            else if (coeff[i]  > 0) s = s + " + " + ( coeff[i]);
            else if (coeff[i]  < 0) s = s + " - " + (-coeff[i]);
            if      (i == 1) s = s + "x";
            else if (i >  1) s = s + "x^" + i;
        }
        return s;
    }

    
	/**
	 * Find the roots of the polynomial.  The imaginary part is not returned.
	 * @return an array containing the roots of the polynomial
	 */
	public double[] calcRoots(){
		return calcRoots(coeff);
	}
	
	/**
	 * Find the roots of the polynomial specified by the parameters. The imaginary part is not returned.
	 * @param parameters an array representing the coefficents of the polynomial
	 * @return an array containing the roots of the polynomial
	 */
	public static double[] calcRoots(double[] parameters){
		if (parameters.length-1 == 2) return solveSecondDegree(parameters);
		if (parameters.length-1 == 3) return solveThirdDegree(parameters);
		
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
	
	/* Descarte's Rule of Signs: The number of positive roots is either equal to the number of changes in 
	 * signs of the cooefficient (zero coefficient is assumed to have the same sign as its predecessor) or
	 * is less by a multiple of 2
	 */
	public int numberPositiveRoots() {
		int count = 0;
		for (int i = 0; i < deg; i++) 
			if (coeff[i]*coeff[i+1] < 0.0) count++;
		return count;
	}
	
	/* number of negative roots in P(x) is the same as the number of positive roots in P(-x) */
	public int numberNegativeRoots() {
		int count = deg + 1 - numberPositiveRoots();
		for (int i = 1; i <= deg; i++)
			if (coeff[i] == 0.0) count--;
		return count;
	}
	
	/* polynomium a[n]*x^n + a[n-1]*x^(n-1) + ... + a[1]*x + a[0] */ 
	public static double solveHighDegree(double[] a) {
		boolean convergence = false;
		boolean divisionByZero = false;
		double x = 1;
		double xNext = 1.0;
		double[] d = new double[a.length-1];
		for (int j = 0; j < d.length; j++) d[j] = (j+1)*a[j+1];
		double tolerance = 0.00000001;
		int maxIterations = 100;
		double num = 0.0;
		double denum = 0.0;

		int i = 0;
		while ((i < maxIterations)  && !convergence && !divisionByZero) {
			i++;
			
			num = a[a.length-2] + x*a[a.length-1];
			denum = d[a.length-2];
			for (int j = a.length-3; j >= 0; j--) {
				num = a[j] + x*num;
				denum = d[j] + x*denum;
			}
			if (Math.abs(denum) < Constants.EPSILON) {
				System.out.println("WARNING: denominator is too small");
				divisionByZero = true;
			}
			else {
				xNext = x - num/denum;
				if (Math.abs(xNext - x) < tolerance) convergence = true;
				x = xNext;
			}
		}
		if (divisionByZero) System.out.println("WARNING: division by zero");
		if (convergence) {
			if (Math.abs(num) > Constants.EPSILON) System.out.println("WARNING: convergence not to the root");
			return xNext;
		}
		else {
			System.out.println("WARNING: solution not within specified tolerance");
			return xNext;
		}
	}

	
	public static void main(String[] args){
		double[] a = new double[6];
		a[0] = 6; a[1] = 5; a[2] = 4; a[3] = 3; a[4] = 2; a[5] = 1;
		double root = Polynomial.solveHighDegree(a);
	}
}
