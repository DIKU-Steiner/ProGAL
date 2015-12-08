package ProGAL.math;

import ProGAL.Function;
 
public class RootFinding {

	/** Root-finding algorithm that repeatedly bisects an interval and then selects a subinterval in which a 
	 * root must lie for further processing. It is a very simple and robust method, but it is also relatively slow. 
	 * Because of this, it is often used to obtain a rough approximation to a solution which is then used as a 
	 * starting point for more rapidly converging methods The method is also called the interval halving method, 
	 * the binary search method, or the dichotomy method.*/
	public static Double bisection(Function f, double a, double b, double tol, int nMax) {
		double fa = f.getValue(a);
		double fb = f.getValue(b);
		double fc;
		if (fa*fb >= 0.0) { 
			System.out.println("Bisection method failed - root is not bracketed.");
			return null;
		}
		double signA = Math.signum(fa);
		double c = 0;
		int n = 0;
		while (++n <= nMax) {
			c = (a + b)/2;
			fc = f.getValue(c);
			System.out.println(n +". c = " + c + ", fc = " + fc);
			if ((Math.abs(fc) < Constants.EPSILON) || ((b - a)/2 < tol)) return c;
			if (Math.signum(fc) == signA)  a = c; else b = c; 
		}
		System.out.println("Bisection method failed - max number of iterations.");
		return c;
	}
		
	/** Root-finding algorithm that uses a succession of roots of secant lines to better approximate 
	 * a root of a function f. The secant method can be thought of as a finite difference approximation 
	 * of Newton's method. However, the method was developed independently of Newton's method, 
	 * and predated the latter by over 3,000 years.*/
	public static Double secant(Function f, double x, double x1, double tol, int nMax) {
		int n = 0;
		double x2 = 0;
		double fx1 = f.getValue(x1);
		double fx = f.getValue(x);
		while (++n <= nMax) {
			if ((Math.abs(fx1) < Constants.EPSILON) || (Math.abs(x1 - x) < tol)) return x;
			x2 = x1 - fx1*(x1-x)/(fx1-fx);
			x  = x1;
			fx = fx1;
			x1 = x2;
			fx1 = f.getValue(x1);
		}
		System.out.println("Bisection method failed - max number of iterations.");
		return x;
	}

	/** Root-finding algorithm that starts with an initial guess which is reasonably close to the true root, 
	 * then the function is approximated by its tangent line (which can be computed using the tools of calculus), 
	 * and one computes the x-intercept of this tangent line (which is easily done with elementary algebra). 
	 * This x-intercept will typically be a better approximation to the function's root than the original guess, 
	 * and the method can be iterated.*/
	public static Double newton(Function f, double x, double x1, double tol, int nMax) {
		Function d = f.getDerivative();
		int n = 0;
		double fx = f.getValue(x);
		double dx;
		while (++n <= nMax) {
			if ((Math.abs(fx) < Constants.EPSILON) || (Math.abs(x1 - x) < tol)) return x;
			dx = fx/d.getValue(x);
			x1 = x;
			x = x - dx;
			fx = f.getValue(x);
		}
	System.out.println("Newton method failed - max number of iterations.");
	return x;
	}
	
	/** Root-finding algorithm that uses quadratic interpolation to approximate the inverse of f. 
	 * This algorithm is rarely used on its own, but it is important because it forms part of the 
	 * popular Brent's method.*/
	public static Double inverseQuadraticInterpolation(Function f, double x, double x1, double x2, double tol, int nMax) {
		int n = 0;
		double x0;
		double fx = f.getValue(x);
		double fx1 = f.getValue(x1);
		double fx2 = f.getValue(x2);
		while (++n <= nMax) {
			if ((Math.abs(fx) < Constants.EPSILON) || (Math.abs(x1 - x) < tol)) return x;
			x0 = x2*fx1*fx/((fx2-fx1)*(fx2-fx)) + x1*fx2*fx/((fx1-fx2)*(fx1-fx)) +x*fx2*fx1/(((fx-fx2)*(fx-fx1)));
			x2 = x1; fx2 = fx1;
			x1 = x;  fx1 = fx;
			x = x0;  fx = f.getValue(x);
		}
		System.out.println("Inverse Quadratic Interpolation method failed - max number of iterations.");
		return x;
	}
	
	public static Double dekker(Function f, double a, double b, double tol, int nMax) {
		double fa = f.getValue(a);
		double fb = f.getValue(b);
		if (fa*fb >= 0.0) { 
			System.out.println("Dekker method failed - root is not bracketed.");
			return null;
		}
		double bOld, fbOld;
		double temp;
		int n = 0;
		while (++n <= nMax) {
			if (Math.abs(fa) < Math.abs(fb)) { temp = a; a = b; b = temp; temp = fa; fa = fb; fb = temp; }
			bOld = b; fbOld = fb;
			if (Math.abs(fb - fa) > Constants.EPSILON) b -= fb*(b-a)/(fb-fa); else b = (a + b)/2;
			fb = f.getValue(b);
			System.out.println(n +". b = " + b + ", fb = " + fb);
			if ((Math.abs(fb) < Constants.EPSILON) || (Math.abs((b - a)/2) < tol)) return b;
			if (fa*fb >= 0.0) { a = bOld; fa = fbOld; }
		}
		System.out.println("Dekker method failed - max number of iterations.");
		return b;
	}
	
	/** Root-finding algorithm combining the bisection method, the secant method and inverse quadratic interpolation. 
	 * It has the reliability of bisection but it can be as quick as some of the less reliable methods. The algorithm 
	 * tries to use the potentially fast-converging secant method or inverse quadratic interpolation if possible, 
	 * but it falls back to the more robust bisection method if necessary.
	 * @param args
	 */
	public static Double brent(Function f, double a, double b, double tol, double delta, int nMax) {
		double fa = f.getValue(a);
		double fb = f.getValue(b);
		double fc = fa;
		if (fa*fb >= 0.0) {
			System.out.println("Brent method failed - root is not bracketed.");
			return null;
		}
		double c = a;
		double d = c;
		boolean mflag = true;
		double temp, s, fs;
		int n = 0;
		while (++n <= nMax) {
			if (Math.abs(fa) < Math.abs(fb)) { temp = a; a = b; b = temp; temp = fa; fa = fb; fb = temp; }
			if ((Math.abs(fc - fa) > Constants.EPSILON) && (Math.abs(fc -fb) > Constants.EPSILON))   
				s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb));
			else s = b - fb*(b-a)/(fb-fa);
			if (((s <= (3*a + b)/4) || s >= b) ||
				(mflag  && ((Math.abs(s-b) >= Math.abs(b-c)/2) || (Math.abs(b-c) < delta))) ||
				(!mflag && ((Math.abs(s-b) >= Math.abs(c-d)/2) || (Math.abs(c-d) < delta)))) {
				s = (a + b)/2;
				mflag = true;
			} 
			else mflag = false;
			fs = f.getValue(s);
			d = c;
			c = b;
			if (fa*fs < 0) { b = s; fb = fs; } else { a = s; fa = fs; }
			if (Math.abs(fa) < Math.abs(fb)) { temp = a; a = b; b = temp; temp = fa; fa = fb; fb = temp; }
			System.out.println(n +". b = " + b + ", fb = " + fb);
			if ((Math.abs(fb) < Constants.EPSILON) || (Math.abs((b - a)/2) < tol)) return b;
		}
		System.out.println("Brent method failed - max number of iterations.");
		return c;
	}
	
	
	public static void main(String[] args) {
		double[] c = new double[4]; 
		c[0] = -2;
		c[1] = -1;
		c[2] = 0;
		c[3] = 1;
		Function f = new Function(c);
		double root = RootFinding.brent(f, 1, 2, 0.000001, 0.00001, 25);
		System.out.println("Root at " + root + " with value " + f.getValue(root));
	}

}
