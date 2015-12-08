package ProGAL.math;
import java.math.*;
public class Trigonometry {

	/** solves equation AsinX + BCosX + C = 0 */
	public static Double[] solveAsinXPlusBcosXplusC(double A, double B, double C) {
		Double[] roots = new Double[2];
		double a = C-B; 
		if (Math.abs(a) < Constants.EPSILON) { /*System.out.println("No solutions");*/ return null; }
		double b = 2 * A;
		double c = B + C;
		double delta = b*b - 4*a*c;
		if (delta < - Constants.EPSILON) { /*System.out.println("No solution");*/ return null; }
		if (delta < Constants.EPSILON) {
			roots[0] = 2*Math.atan(-0.5*b/a);
			roots[1] = null;
		}
		else {
			delta = Math.sqrt(delta);
			roots[0] = Math.atan2(0.5*(delta - b), a);
			if (roots[0] < 0) roots[0] = Constants.TAU + roots[0];
			roots[0] = 2*roots[0];
			if (roots[0] > Constants.TAU) roots[0] = roots[0] - Constants.TAU;
			roots[1] = Math.atan2(-0.5*(delta + b), a);
			if (roots[1] < 0) roots[1] = Constants.TAU + roots[1];
			roots[1] = 2*roots[1];
			if (roots[1] > Constants.TAU) roots[1] = roots[1] - Constants.TAU;
		}
		return roots;
	}
	
	/** solves equation AsinX + BCosX + C = 0 */
	public static Double[] solveBDAsinXPlusBcosXplusC(double A, double B, double C) {
		Double[] roots = new Double[2];
		BigDecimal bdA = BigDecimal.valueOf(A); 
		BigDecimal bdB = BigDecimal.valueOf(B);
		BigDecimal bdC = BigDecimal.valueOf(C);
		BigDecimal a = bdC.subtract(bdB); 
		if (a.abs().compareTo(Constants.EPSILONBD) == -1) { System.out.println("No solution"); return null; }
		BigDecimal b = bdA.multiply(new BigDecimal("2.0"), MathContext.DECIMAL128);
		BigDecimal c = bdB.add(bdC);
		BigDecimal delta = b.pow(2).subtract(a.multiply(c).multiply(new BigDecimal("4.0")));
		if (delta.compareTo(Constants.EPSILONBD.negate()) == -1) { /*System.out.println("No solution");*/ return null; }
		if (delta.compareTo(Constants.EPSILONBD) == -1) {
			BigDecimal d = b.divide(a).multiply(new BigDecimal("-0.5"));
			roots[0] = 2*Math.atan(d.doubleValue());
			roots[1] = null;
//			System.out.println("t = " + roots[0]);
		}
		else {
			delta = BigDecimal.valueOf(Math.sqrt(delta.doubleValue()));
			BigDecimal e = delta.subtract(b).divide(a, 50, RoundingMode.HALF_UP).multiply(new BigDecimal("0.5"));
			roots[0] = 2.0*Math.atan(e.doubleValue()); if (roots[0] < 0) roots[0] = roots[0] + Constants.TAU;
			BigDecimal f = delta.add(b).divide(a, 50, RoundingMode.HALF_UP).multiply(new BigDecimal("-0.5"));
			roots[1] = 2.0*Math.atan(f.doubleValue()); if (roots[1] < 0) roots[1] = roots[1] + Constants.TAU;
//			for (int i = 0; i < 2; i++) System.out.println("t" + i + " = " + Functions.toDeg(roots[i]));
		}
		return roots;
	}

	
	public static void main(String[] args) {
		Trigonometry.solveAsinXPlusBcosXplusC(12, 5, -4);
		Trigonometry.solveAsinXPlusBcosXplusC(1, 1, 0);
	}
}
