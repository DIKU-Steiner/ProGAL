package ProGAL.math;

import static java.lang.Math.sqrt;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.acos;
import static java.lang.Math.PI;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
	
//	private class PairOfPolynomials {
//		private Polynomial p;
//		private Polynomial q;
//		
//		private PairOfPolynomials(Polynomial p, Polynomial q) {
//			this.p = p;
//			this.q = q;
//		}
//	}
	
	/** 
	 * creates a zero polynomial
	 */
	Polynomial() {
		coeff = new double[1];
		deg = 0;
	}
	
	/**
	 * creates a polynomal
	 * @param pars
	 */
	// Daisy: Made it public from nothing
	public Polynomial(double[] pars){
		coeff = pars;
		deg = pars.length-1;
	}
	
	
	/** returns a copy of the polynomial p
	 * @return copy of the polynomial p
	 */
	public Polynomial clone() {
		double[] coefficients = new double[coeff.length];
		for (int i = 0; i < coefficients.length; i++) coefficients[i] = coeff[i];
		return new Polynomial(coefficients);
	}
	
	/** returns TRUE if this polynomial is a zero polynomial 
	 * 
	 * @return
	 */
	public boolean isZeroPolynomial() { return ((deg == 0) && (coeff[0] == 0.0)); }
	
	public boolean dominates(Polynomial p) {
		if (this.deg > p.deg) return true;
		if (this.deg < p.deg) return false;
		for (int i = deg; i >= 0; i--) {
			if (this.coeff[i] > p.coeff[i]) return true;
			if (this.coeff[i] < p.coeff[i]) return false;
		}
		return false;
	}
	
	public Polynomial leadingTerm() {
		double[] coefficients = new double[deg+1];
		coefficients[deg] = coeff[deg];
		return new Polynomial(coefficients);
	}
	
	// Daisy
	public int getDeg() { return deg; }
	
    // returns a + b
    public Polynomial plus(Polynomial b) {
    	int degree = Math.max(deg,  b.deg);
    	while ((degree != -1) && (degree <= deg) && (degree <= b.deg) && (coeff[degree] == -b.coeff[degree])) degree--;
    	if (degree == -1) return new Polynomial(new double[1]);
   	
    	
    	double[] coeff = new double[degree + 1];
        for (int i = 0; i <= this.deg; i++) 
        	if (i <= degree) coeff[i] = this.coeff[i];
        for (int i = 0; i <= b.deg; i++) 
        	if (i <= degree) coeff[i] += b.coeff[i];
        
        //Daisys version put in corner for now
//    	double[] coeff = new double[Math.max(deg+1, b.deg+1)];
//        for (int i = 0; i <= deg; i++) coeff[i] = this.coeff[i];
//        for (int i = 0; i <= b.deg; i++) {
//        	coeff[i] += b.coeff[i];
//        	if (Math.abs(coeff[i]) <= Constants.EPSILON ) coeff[i] = 0;
//        }
//        while ((deg != 0) && (coeff[deg] == 0.0)) deg--;
        return new Polynomial(coeff);
    }

    // returns a - b
    public Polynomial minus(Polynomial b) {
    	int degree = Math.max(deg, b.deg);
    	while ((degree != -1) && (degree <= deg) && (degree <= b.deg) && (coeff[degree] == b.coeff[degree])) degree--;
    	if (degree == -1) return new Polynomial(new double[1]);
    	double[] coeff = new double[degree + 1];
        for (int i = 0; i <= this.deg; i++) 
        	if (i <= degree) coeff[i] = this.coeff[i];
        for (int i = 0; i <= b.deg; i++) 
        	if (i <= degree) coeff[i] -= b.coeff[i];
        return new Polynomial(coeff);

//		Daisys version put in corner
//    	int maxDeg = Math.max(deg, b.deg);
//    	double[] coeff = new double[maxDeg+1];
//        for (int i = 0; i <= deg; i++) coeff[i] = this.coeff[i];
//        for (int i = 0; i <= b.deg; i++) {
//        	coeff[i] -= b.coeff[i];
//        	if (Math.abs(coeff[i]) <= Constants.EPSILON ) coeff[i] = 0;
//        }
//        int tmpDeg = maxDeg;
//        while ((tmpDeg != 0) && (coeff[tmpDeg] == 0.0)) { tmpDeg--; }
//        if (tmpDeg != maxDeg) {
//        	double[] coeff2 = new double[tmpDeg+1];
//        	for (int i = 0; i <= tmpDeg ; i++) {
//        		coeff2[i] = coeff[i];
//        	}
//        	return new Polynomial(coeff2);
//        } else {
//        	return new Polynomial(coeff);
//        }
    }

    // returns a * b
    public Polynomial times(Polynomial b) {
// 		Pawel put in corner
//    	double[] coeff = new double[deg + b.deg + 1];
    	double[] coeffZero = {0};
    	Polynomial zero = new Polynomial(coeffZero);
    	if (this.equal(zero) || b.equal(zero)) return zero;
    	double[] coeff = new double[deg + b.deg+1];

    	for (int i = 0; i <= deg; i++)
            for (int j = 0; j <= b.deg; j++) coeff[i+j] += this.coeff[i] * b.coeff[j];
        return new Polynomial(coeff);
    }
    
// Daisy put in corner
//    /**
//     * dividing a polynomial by another polynomial of the same or lower degree
//     * @param d
//     * @return
//     */
//    public PairOfPolynomials longDivision(Polynomial d) {
//    	if (d.isZeroPolynomial())  throw new IllegalArgumentException("Divisor is a zero polynomial");
//    	Polynomial q = new Polynomial();
//    	Polynomial r = clone();
//    	while (!r.isZeroPolynomial() && (r.deg >= d.deg)) {
//    		double[] coefficients = new double[r.deg - d.deg + 1];
//    		coefficients[r.deg - d.deg] = r.coeff[r.deg]/d.coeff[d.deg];
//    		Polynomial t = new Polynomial(coefficients);
//    		q = q.plus(t);
//    		r = r.minus(t.times(d));    		
//    	}
//    	return new PairOfPolynomials(q, r);
//    }
//
//    public Polynomial greatestCommonDivisor(Polynomial b) {
//       	if (b.deg == 0) return this;
//       	if (b.dominates(this)) return this.greatestCommonDivisor(b.longDivision(this).q);
//       	else return b.greatestCommonDivisor(this.longDivision(b).q.toMonic());
//    }
//    
//    public ArrayList<Polynomial> squareFreeFactorizationTH() {
//    	ArrayList<Polynomial> p = new ArrayList<Polynomial>();
//    	Polynomial c, cN, d, dN;
//    	c = this.greatestCommonDivisor(this.differentiate());
//    	d = this.longDivision(c).p;
//    	while (c.deg != 0) {
//    		cN = c.greatestCommonDivisor(c.differentiate());
//    		dN = c.longDivision(cN).p;
//    		p.add(d.longDivision(dN).p);
//    		c = cN;
//    		d = dN;
//    	} 
//    	return p;
//    }
//    
//    public ArrayList<Polynomial> squareFreeFactorization() {
//    	Polynomial deriv = this.differentiate();
//    	deriv.makeMonic();
//    	ArrayList<Polynomial> a = new ArrayList<Polynomial>();
//    	Polynomial b;
//    	Polynomial c;
//    	Polynomial d;
//    	a.add(0, this.greatestCommonDivisor(deriv));
//    	b = this.longDivision(a.get(0)).p;
//    	c = deriv.longDivision(a.get(0)).p;
//    	d = c.minus(b).differentiate();
//    	int i = 1;
//    	do {
//    		a.add(i, b.greatestCommonDivisor(d));
//    		b = b.longDivision(a.get(i)).p;
//    		c = d.longDivision(a.get(i)).p;
//    		d = c.minus(b.differentiate());
//    		i++;
//    	} while (b.deg != 0); 
//    	a.remove(0);
//    	return a;
//    }
    
    
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

//    public void makeMonic() {
//    	for (int i = 0; i <= deg; i++) this.coeff[i] = coeff[i]/coeff[deg];
//    }
    
    public Polynomial toMonic() {
    	Polynomial p = this.clone();
       	for (int i = 0; i <= deg; i++) p.coeff[i] = p.coeff[i]/p.coeff[deg];
       	return p;
    }
    
    // differentiate this polynomial and return it
    public Polynomial differentiate() {
    	if (deg == 0) return new Polynomial(new double[1]);
    	double[] par = new double[deg];
    	for (int i = 0; i < deg; i++) {
    		par[i] = (i+1) * coeff[i+1];
//    	for (int i = 1; i <= deg; i++) {
//    		par[i-1] = i * coeff[i];
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
	

	/* 
	 * TODO: Pawel siger .. m�ske ikke f�rdig 
	 */
	private static double[] solveQuadric(double p2, double p1, double p0) {
		double b = p1/p2;
		double c = p0/p2;
		double D = b*b - 4*c;
		if (D > 0.0) {
			double sqrtD = Math.sqrt(D);
			return new double[]{(sqrtD - b)/2, (-sqrtD - b)/2};
		}
		else {
			if (D < 0.0) return new double[]{};
			else return new double[]{-b/2};
		}
	}
	
	// Daisy
	public Double solveFirstDegree() {
		return ((-this.coeff[0])/this.coeff[1]);
	}
	public Double[] solveSecondDegree() { //Copy of solveSecondDegree(double[] parameters)
		double a = this.coeff[2];
		double b = this.coeff[1];
		double c = this.coeff[0];
		double p = -0.5*b/a;
		double D = p*p-c/a;
		if(D<0) return new Double[]{};
		if(D==0) return new Double[]{p};
		else {
			double sqD = Math.sqrt(D);
			return new Double[]{p-sqD, p+sqD};
		}
	}
	
	/** 
	 * Solves quadratic equation. 
	 * @hops 4 to 5.  
	 */
	private static double[] solveSecondDegree(double[] parameters){
		double a = parameters[2];
		double b = parameters[1];
		double c = parameters[0];
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
	 * Solves cubic equation.
	 * @hops 3 to 35
	 * Daisy made it public
	 */
	public static double[] solveThirdDegree(double[] parameters){
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


	private static double[] depressQuarticEquation(double[] parameters) {
		double a = parameters[0];
		double b = parameters[1];
		double c = parameters[2];
		double d = parameters[3];
		double e = parameters[4];
		return new double[]{a, 
						0, 
						-3*b*b/(8*a) + c, 
						b*b*b/(8*a*a) - b*c/(2*a) + d,
						-3*b*b*b*b/(256*a*a*a) + b*b*c/(16*a*a) - b*d/(4*a) + e};
	}
	
	private static double[] solveDepressedQuartic(double[] parameters) {
		double A = parameters[2];
		double B = parameters[3];
		double C = parameters[4];
		double[] thirdDegree = solveThirdDegree(new double[]{-8, 
									  4*A - 24*Math.sqrt(C),
									  8*A*Math.sqrt(C) - 16*C,
									  B*B});
		return thirdDegree;
		
	}
	
	private static double[] solveFourthDegree(double[] parameters){
		if(parameters[0]==0) return solveThirdDegree(new double[]{parameters[1],parameters[2],parameters[3],parameters[4]});

		double a1 = parameters[1]/parameters[0];
		double a2 = parameters[2]/parameters[0];
		double a3 = parameters[3]/parameters[0];
		double a4 = parameters[4]/parameters[0];
		
		double b1 = -a2;
		double b2 = a1*a3 - 4*a4;
		double b3 = 4*a2*a4 - a3*a3 - a1*a1*a4;
		double[] cubicRoots = solveThirdDegree(new double[]{1, b1, b2, b3});
		
		double y1 = cubicRoots[0];
		double d1 = Math.sqrt(a1*a1 - 4*a2 + 4*y1);
		double d2 = Math.sqrt(y1 -4*a4);
		double c1 = (a1 + d1)/2;
		double c2 = (y1*y1 - d2)/2;
		double[] firstPair = solveSecondDegree(new double[]{1, c1, c2});
		c1 = (a1 - d1)/2;
		c2 = (y1*y1 + d2)/2;
		double[] secondPair = solveSecondDegree(new double[]{1, c1, c2});
		return new double[]{firstPair[0], firstPair[1], secondPair[0], secondPair[1]};
	}
	

	public static Double[] solveQuartic(double[] c) {
		DecimalFormat newFormat = new DecimalFormat("#.#########");
		c[0] =  Double.valueOf(newFormat.format(c[0]));
		c[1] =  Double.valueOf(newFormat.format(c[1]));
		c[2] =  Double.valueOf(newFormat.format(c[2]));
		c[3] =  Double.valueOf(newFormat.format(c[3]));
		c[4] =  Double.valueOf(newFormat.format(c[4]));
		Double[] roots = new Double[4];
		double[]  coeffs= new double[ 4 ];
	    double  z, u, v, sub;
	    double  A, B, C, D;
	    double  sq_A, p, q, r;
	    int     i, num = 0;

	    // normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0

	    A = c[ 3 ] / c[ 4 ];
	    B = c[ 2 ] / c[ 4 ];
	    C = c[ 1 ] / c[ 4 ];
	    D = c[ 0 ] / c[ 4 ];

	    //  substitute x = y - A/4 to eliminate cubic term:	x^4 + px^2 + qx + r = 0

	    sq_A = A * A;
	    p = - 3.0/8.0 * sq_A + B;
	    q = 1.0/8.0 * sq_A * A - 1.0/2.0 * A * B + C;
	    r = - 3.0/256.0 * sq_A * sq_A + 1.0/16.0 * sq_A * B - 1.0/4.0 * A*C + D;

	    if (Math.abs(r)==0.0) {
			// no absolute term: y(y^3 + py + q) = 0
			coeffs[ 0 ] = q;
			coeffs[ 1 ] = p;
			coeffs[ 2 ] = 0;
			coeffs[ 3 ] = 1;
	
			Double[] s = solveCubic(coeffs);
			int sLen = 0;
			for (Double d : s) {
				if (d != null) sLen += 1;
			}
			num = sLen;
			roots[ num++ ] = 0.0;
	    } else {
			// solve the resolvent cubic
			coeffs[ 0 ] = 1.0/2.0 * r * p - 1.0/8.0 * q * q;
			coeffs[ 1 ] = - r;
			coeffs[ 2 ] = - 1.0/2.0 * p;
			coeffs[ 3 ] = 1.0;
	
			Double[] s =  solveCubic(coeffs);
	
			// take the one real solution
	
			z = s[ 0 ];
			
			// to build two quadric equations
	
			u = z * z - r;
			v = 2 * z - p;
			if (Math.abs(u)==0.0) {
			    u = 0.0;
			} else if (u > 0) {
			    u = Math.sqrt(u);
			} else {
			    return null;
			}
	
			if (Math.abs(v)==0.0) {
			    v = 0;
			} else if (v > 0) {
			    v = Math.sqrt(v);
			} else {
			    return null;
			}
	
			coeffs[ 0 ] = z - u;
			if (q < 0) {
				coeffs[ 1 ] = -v;
			} else coeffs[ 1 ] = v;
			coeffs[ 2 ] = 1;
			Double[] s1 = solveQuadric(coeffs);
			
			coeffs[ 0 ]= z + u;
			if (q < 0) {
				coeffs[ 1 ] = v;
			} else coeffs[ 1 ] = -v;
			coeffs[ 2 ] = 1;
			Double[] s2 = solveQuadric(coeffs);
			
			int s1Len = 0;
			int s2Len = 0;
			if (s1 != null) {
				for (Double d : s1) {
					if (d != null) s1Len += 1;
				}
			}
			if (s2 != null) {
				for (Double d : s2) {
					if (d != null) s2Len += 1;
				}
			}
			roots = new Double[s1Len+s2Len];
			int k = 0;
			if (s1Len>0) {
				for (Double d : s1) {
					if (d==null) break;
					roots[k] = s1[k];
					k += 1;
				}
			}
			if (s2Len>0) {
				for (Double d : s2) {
					if (d==null) break;
					roots[k] = s2[k%2];
					k += 1;
				}
			}
			num = s1Len+s2Len;
			if (num == 0) return null;
	    }

	    // resubstitute

	    sub = 1.0/4.0 * A;

	    for (i = 0; i < num; ++i){
	    	roots[ i ] -= sub;
	    }

	    return roots;
	}
	
	public static Double[] solveCubic(double[] c) {
		Double[] roots = new Double[3];
		int i, num;
		double sub;
		double A, B, C;
		double sq_A, p, q;
		double cb_p, D;
		
		// normal form: x^3 + Ax^2 + Bx + C = 0

	    A = c[ 2 ] / c[ 3 ];
	    B = c[ 1 ] / c[ 3 ];
	    C = c[ 0 ] / c[ 3 ];

	    //  substitute x = y - A/3 to eliminate quadric term: x^3 +px + q = 0 

	    sq_A = A * A;
	    p = 1.0/3.0 * (- 1.0/3.0 * sq_A + B);
	    q = 1.0/2.0 * (2.0/27.0 * A * sq_A - 1.0/3.0 * A * B + C);
	    
	    // use Cardano's formula

	    cb_p = p * p * p;
	    D = q * q + cb_p;
	    if (Math.abs(D)==0.0) {
			if (Math.abs(q)<Constants.EPSILON) {// one triple solution
			    roots[ 0 ] = 0.0;
			    num = 1;
			} else {// one single and one double solution
			    double u = Math.cbrt(-q);
			    roots[ 0 ] = 2 * u;
			    roots[ 1 ] = - u;
			    num = 2;
			}
	    } else if (D < 0) {
			double phi = 1.0/3.0 * Math.acos(-q / Math.sqrt(-cb_p));
			double t = 2 * Math.sqrt(-p);
	
			roots[ 0 ] =   t * Math.cos(phi);
			roots[ 1 ] = - t * Math.cos(phi + Math.PI / 3.0);
			roots[ 2 ] = - t * Math.cos(phi - Math.PI / 3.0);
			num = 3;
	    } else {// one real solution
			double sqrt_D = Math.sqrt(D);
			double u = Math.cbrt(sqrt_D - q);
			double v = - Math.cbrt(sqrt_D + q);
			roots[ 0 ] = u + v;
			num = 1;
	    }

	    // resubstitute

	    sub = 1.0/3.0 * A;

	    for (i = 0; i < num; ++i)
		roots[ i ] -= sub;

	    return roots;
	}
	
	// From http://tog.acm.org/resources/GraphicsGems/gems/Roots3And4.c
	public static Double[] solveQuadric(double[] c) {
		Double[] roots = new Double[2]; 
		double p, q, D;

	    /* normal form: x^2 + px + q = 0 */

	    p = c[ 1 ] / (2 * c[ 2 ]);
	    q = c[ 0 ] / c[ 2 ];

	    D = p * p - q;

	    if (Math.abs(D)<Constants.EPSILON) {
	    	roots[ 0 ] = - p;
	    	return roots;
	    } else if (D < 0) {
	    	return null;
	    } else if (D > 0) {
	    	double sqrt_D = Math.sqrt(D);

			roots[ 0 ] =   sqrt_D - p;
			roots[ 1 ] = - sqrt_D - p;
			return roots;
	    }
	    return null;
	}
	
	private double gcd(double a, double b){
		while (b!=0) {
			double tmp = Math.abs(a) % Math.abs(b);
			a = Math.abs(b);
			b = tmp;
			if (Math.abs(b)<=Constants.EPSILON) b = 0.0;
		}
		return a;
	}
	public List<Polynomial> squareFreeFact() {
		double gcd1 = gcd(this.coeff[4], this.coeff[3]);
		double gcd2 = gcd(gcd1, this.coeff[2]);
		double gcd3 = gcd(gcd2, this.coeff[1]);
		double gcd4 = gcd(gcd3, this.coeff[0]);
		System.out.println("GCD = "+gcd4);
		System.out.println("Old f = "+this.toString());
		double[] f_newCoeff = {this.coeff[0]/gcd4, this.coeff[1]/gcd4, this.coeff[2]/gcd4, this.coeff[3]/gcd4, this.coeff[4]/gcd4};
		Polynomial f_new = new Polynomial(f_newCoeff);
		System.out.println("New f = "+f_new.toString());
/*		
		double[] depressedPars = new double[5];
		double a = this.coeff[3]/this.coeff[4];
		double b = this.coeff[2]/this.coeff[4];
		double c = this.coeff[1]/this.coeff[4];
		double tmpd = this.coeff[0]/this.coeff[4];
		double p = (8*b-3*a*a)/8;
		double q = (Math.pow(a, 3)-4*a*b+8*c)/8;
		double r = ((-3)*Math.pow(a, 4)+256*tmpd-64*c*a+16*a*a*b)/256;
		depressedPars[0] = r;
		depressedPars[1] = q;
		depressedPars[2] = p;
		depressedPars[4] = 1;
		Polynomial f = new Polynomial(depressedPars);
		
		double[] cubicPars = new double[4];
		cubicPars[0] = (Math.pow(p, 3)/2)-((p*r)/2)-((q*q)/8);
		cubicPars[1] = 2*p*p-r;
		cubicPars[2] = (5/2)*p;
		cubicPars[3] = 1;
//		Polynomial cubic = new Polynomial(cubicPars);
		double[] cubicRoots = solveThirdDegree(cubicPars);
		double y = 0;
		for (double root : cubicRoots) {
			if (p+2*root != 0.0) {
				y = root;
				break;
			}
		}
		System.out.println("y = "+y);
		System.out.println("p = "+p);
		System.out.println("p = p ?"+(p == (8*this.coeff[2]*this.coeff[4]-3*this.coeff[3]*this.coeff[3])/(8*this.coeff[4]*this.coeff[4])));
		System.out.println("q = "+q);
		System.out.println("r = "+r);
		double constant = -(this.coeff[3]/(4*this.coeff[4]));
		double [] xs = new double[4];
		xs[0] = constant + (Math.sqrt(p+2*y)+Math.sqrt(-(3*p+2*y+((2*q)/Math.sqrt(p+2*y)))))/2;
		xs[1] = constant + (-Math.sqrt(p+2*y)+Math.sqrt(-(3*p+2*y-((2*q)/Math.sqrt(p+2*y)))))/2;
		xs[2] = constant + (Math.sqrt(p+2*y)-Math.sqrt(-(3*p+2*y+((2*q)/Math.sqrt(p+2*y)))))/2;
		xs[3] = constant + (-Math.sqrt(p+2*y)-Math.sqrt(-(3*p+2*y-((2*q)/Math.sqrt(p+2*y)))))/2;
		for (int i = 0; i<4; i++) {
			System.out.println("Root = "+xs[i]);
		}*/
		
		
		System.out.println("SquareFree going!");
		List<Polynomial> ret = new ArrayList<Polynomial>();
		Polynomial fPrime = this.differentiate();
		System.out.println("func' = "+fPrime.toString());
		Polynomial d = gcd(this, fPrime);
		System.out.println("d1 = "+d.toString());
		Polynomial e = (this.longDivision(d))[0];
		System.out.println("e1 = "+e.toString());
		double[] coeffOne = {1};
		Polynomial one = new Polynomial(coeffOne);
		Polynomial dPrime;
		Polynomial d_new;
		Polynomial e_new;
		while (!(d.equal(one))) {//d.deg==0 {
			dPrime = d.differentiate();
			System.out.println("d = "+d.toString());
			System.out.println("dPrime = "+dPrime.toString());
			d_new = gcd(d, dPrime);
			System.out.println("d_new = "+d_new.toString());
			e_new = (d.longDivision(d_new)[0]);
			System.out.println("eNew = d/dNew = "+d.toString()+" / "+d_new.toString());
			System.out.println("e/eNew = "+e.toString()+" / "+e_new.toString());
			System.out.println("         "+(e.longDivision(e_new))[0]);
			ret.add((e.longDivision(e_new))[0].makeMonic());
			if (d_new.equal(one)) {
				ret.add(d);
			}
/*			Polynomial b_new;
			Polynomial b = fPrime.longDivision(d.minus(e.differentiate()))[0];
		while (!(e.equal(one))) {
			d = gcd(e, b);
			ret.add(d);
			e_new = e.longDivision(d)[0];
			b_new = b.longDivision(d.minus(e_new.differentiate()))[0];
			d_new = gcd(e, d);
			e_new = e.longDivision(d_new)[0];
			ret.add(e_new);*/
			d = d_new;
			e = e_new;
			double greatCoeff = d.coeff[d.deg];
			for (int i = 0 ; i<d.deg+1 ; i++) {
				d.coeff[i] = d.coeff[i]/greatCoeff;
			}
		}
		if (ret.size() == 0) {
			ret = (this.makeMonic()).distinctFact();
		}
		System.out.println("SFF result = "+ret.toString());
		return ret;
	}
	public static Polynomial gcd(Polynomial a, Polynomial b) {
		if (b.deg>a.deg) {
//			return null;
			Polynomial tmp = a;
			a = b;
			b = tmp;
		}
		Polynomial r; /*= a.longDivision(b, 1);
		a = b;
		b = r;
		System.out.println("a = "+a.toString()+" b = "+b.toString());
		r = a.longDivision(b, 1);
		System.out.println("r = "+r.toString());
		a = b;
		b = r;
		System.out.println("a = "+a.toString()+" b = "+b.toString());
		r = a.longDivision(b, 1);
		a = b;
		b = r;
		System.out.println("a = "+a.toString()+" b = "+b.toString());*/
		while (true) { //!(b.equal(zero))
			System.out.println("gcd Call");
			r = (a.longDivision(b))[1];
			System.out.println("intermediate remaining = "+r.toString());
			a = b;
			b = r;
			if (b.deg==0) {
				if (Math.abs(b.coeff[0])<Constants.EPSILON) break;
			}
		}
		return a.makeMonic();
	}
	
	/* Returns either the quotient (0) or remainder (1) of this/d */
	private Polynomial[] longDivision(Polynomial d) {
		if (d == null) { System.out.println("d is null!"); return null; }
		Polynomial[] ret = new Polynomial[2];
		double[] coeffZero = {0};
		Polynomial q = new Polynomial(coeffZero);
		Polynomial zero = new Polynomial(coeffZero);
		if (d.equal(zero)) { System.out.println("d is zero!"); return null; }
		Polynomial r = this;
		
		int dDeg = d.deg;
		System.out.println("New longDivision round!");
		while (r.deg>=dDeg) {
			Polynomial t = leadDiv(r, d);
			if (t == null) {
				ret[0] = q;
				ret[1] = zero;
				return ret;
			}
			System.out.println("t = "+t.toString());
			q = q.plus(t);
			r = r.minus(t.times(d));
//			System.out.println("intermediate rest = "+r.toString());
			if (r.deg==0) {
				if (Math.abs(r.coeff[0])<Constants.EPSILON) break;
			}
		}
		ret[0] = q;
		ret[1] = r;
		return ret;
	}
	
	private Polynomial leadDiv(Polynomial r, Polynomial d) {
		int rDeg = r.deg;
		int dDeg = d.deg;
		int newDeg = rDeg-dDeg;
		double coef = r.coeff[rDeg]/d.coeff[dDeg];
		if (Math.abs(coef)<Constants.EPSILON) {
			return null;
		}
		double[] pars = new double[newDeg+1]; 
		pars[newDeg] = coef;
		return new Polynomial(pars);
	}
	public Polynomial makeMonic() {
		double div = coeff[deg];
		double[] monicPars = new double[deg+1];
		monicPars[deg] = 1.0;
		for (int i=0 ; i<deg ; i++) {
			monicPars[i] = coeff[i]/div;
			if (Math.abs(monicPars[i]) <= (Constants.EPSILON+1.0) && Math.abs(monicPars[i]) >= (1.0-Constants.EPSILON)) {
				System.out.println("rounding : "+monicPars[i]);
				System.out.println("because    "+Math.abs(monicPars[i])+" <= "+(Constants.EPSILON+1.0));
				monicPars[i] = 1.0;
			}
		}
		return new Polynomial(monicPars);
	}
	private List<Polynomial> distinctFact() {
		List<Polynomial> ret = new ArrayList<Polynomial>();
		Polynomial f = this;
		double[] coeffOne = new double[]{1};
		Polynomial one = new Polynomial(coeffOne);
		int i = 1;
		while (f.getDeg() >= (2*i)) {
			double[] tmpPars = new double[i+1];
			tmpPars[i] = 1;
			tmpPars[1] = -1;
			Polynomial g = gcd(f, new Polynomial(tmpPars));
			if (!g.equal(one)) {
				ret.add(g);
				f = f.longDivision(g)[0];
			}
			i += 1;
		}
		if (!f.equal(one)) {
			ret.add(f);
		}
//		if (ret == null) {
//			ret.add(this);
//		}
		return ret;
	}
	// Daisy ends

	
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
		double[] solution = Polynomial.solveQuadric(3, -2, -5);
		double[] depressedQuartic = Polynomial.depressQuarticEquation(new double[] {1,8,12,2*Math.sqrt(30)-16, 4*Math.sqrt(30)-28});
		Polynomial.solveDepressedQuartic(depressedQuartic);
		Polynomial.solveFourthDegree(new double[]{3.0, 6.0, -123.0, -126.0, 1080.0});
/*		double[] a = new double[6];
		a[0] = 6; a[1] = 5; a[2] = 4; a[3] = 3; a[4] = 2; a[5] = 1;
		double root = Polynomial.solveHighDegree(a);
		*/

		//TODO: Move to test
//		double[] a = new double[12];
//		double[] b = new double[4];
//		double[] b2 = new double[4];
//		double[] c = new double[5];
//		double[] c2 = new double[4];
//		double[] d = new double[7];
//		/*b[0] = -2.5750749675480734E14;
//		b[1]= 1.2142986273055004E10;
//		b[2] = -190866.11854848347;
//		b[3] = 1;*/
//		b[0] = -75;
//		b[1]= 49;
//		b[2] = -11;
//		b[3] = 1;
//		b2[0] = 1;
//		b2[1]= -11;
//		b2[2] = 49;
//		b2[3] = -75;
//		//5.71868E-4x^4 + 0.574504513x^3 - 1.867761297x^2 - 0.758281783x - 0.059311159
//		c[4] = 5.71868E-4;
//		c[3]= 0.574504513;
//		c[2] = - 1.867761297;
//		c[1] = - 0.758281783;
//		c[0] = - 0.059311159;
//		Double[] quarroots = solveQuartic(c);
//		for (double root : quarroots) {
//			System.out.println("Qaud root : "+root);
//		}
//		
//		//c[0] = -228.14885350604553; c[1]= 284.5029410954643; c[2] = -31.343804547815083; c[3] = 1;
//		System.out.println("poly = "+b[3]+"x^3 + "+b[2]+"x^2 + "+b[1]+"x + "+b[0]);
//		//Complex[] cubroots = solveCubic(b);
//		Double[] cubroots = solveCubic(b);
//		for (Double r : cubroots) {
//			System.out.println("cubroot: "+r.toString());
//		}

		/*		double[] pC = new double[3];
		pC[2] = 1; pC[1] = -9; pC[0] = -10;
		double[] qC = new double[2];
		qC[1] = -9; qC[0] = -10;
		Polynomial p = new Polynomial(pC);
		Polynomial q = new Polynomial(qC);
		PairOfPolynomials pair = p.longDivision(q);
		System.out.println("(" + pair.p.toString() +" , " + pair.q + ")");
*/
/*		double[] pC = new double[3];
		pC[2] = 1; pC[1] = 0; pC[0] = 2;
		double[] qC = new double[3];
		qC[2] = 2; qC[1] = 0; qC[0] = 2;
		Polynomial p = new Polynomial(pC).toMonic();
		Polynomial q = new Polynomial(qC).toMonic();
		Polynomial gcd = p.greatestCommonDivisor(q);
		System.out.println(gcd.toString());
*/		
/*		double[] pC = new double[12];
		pC[11] = 1; pC[10] = 4; pC[9] = -9; pC[8] = -46; pC[7] = 23; pC[6] = 192; pC[5] = -7; pC[4] = -358; pC[3] = -24; pC[2] = 304; pC[1] = 16; pC[0] = 96;
		Polynomial p = new Polynomial(pC);
		ArrayList<Polynomial> a = p.squareFreeFactorization();
*/
/*		double[] pC = new double[9];
		pC[8] = 1; pC[7] = 0; pC[6] = 6; pC[5] = 0; pC[4] = 12; pC[3] = 0; pC[2] = 8; pC[1] = 0; pC[0] = 0;
		Polynomial p = new Polynomial(pC);
		ArrayList<Polynomial> a = p.squareFreeFactorizationTH();
*/	}
}
