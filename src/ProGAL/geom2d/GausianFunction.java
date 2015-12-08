package ProGAL.geom2d;

public class GausianFunction {
	double a;
	double b;
	double c;
	
	public GausianFunction(double a, double b, double c) {
		this.a = a;
		this.b = b;
		this.c = c;
	}
	
	public double compute(double x) {
		return a*Math.exp(-((x-b)*(x-b))/(2*c*c));
	}
}
