package ProGAL;

public class Function {

	private double[] c;
	private int n;
	public Function(double[] c) {
		n = c.length;
		this.c = c;
	}

	public double getValue(double x) {
		double res = c[0];
		double prod = 1.0;
		for (int i = 1; i < n; i++) {
			prod = prod * x;
			res = res + c[i]*prod;
		}
		return res;
	}

	public Function getDerivative() {
		double[] d = new double[c.length-1];
		for (int i = 1; i < c.length; i++) d[i-1] = i*c[i];
		return new Function(d);
	}
}
