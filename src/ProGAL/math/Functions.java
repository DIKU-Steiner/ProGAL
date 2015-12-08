package ProGAL.math;

public class Functions {
	public static double toDeg(double alpha) { return 180.0*alpha/Math.PI; }
	public static double toRad(double alpha) { return Math.PI*alpha/180.0; }
	
	/** returns the smallest power of 2 above n, see http://en.wikipedia.org/wiki/Power_of_two */
	public static int roundUpToPowerOf2(int n) {
		int m = n-1;
		for (int i = 1; i < 32; i=2*i) m = m | m>>i;
		return m+1;
	}

}
