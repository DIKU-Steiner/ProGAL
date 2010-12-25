package ProGAL.geom3d.complex.alphaComplex;

public class Interval {
    private double left;
    private double right;

    public Interval(double left, double right) {
        assert (left <= right);
        this.left  = left;
        this.right = right;
    }

    public double getLeft()  { return left;  }
    public double getRight() { return right; }
    public void setLeft(double left) { this.left = left; }
	public void setRight(double right) { this.right = right; }

	// does this interval a intersect b?
    public boolean intersects(Interval b) {
        Interval a = this;
        if (b.left <= a.right && b.left >= a.left) { return true; }
        if (a.left <= b.right && a.left >= b.left) { return true; }
        return false;
    }

    public String toString() {
        return "[" + left + ", " + right + "]";
    }



    // test client
    public static void main(String[] args) {
        Interval a = new Interval(15, 20);
        Interval b = new Interval(25, 30);
        Interval c = new Interval(10, 40);
        Interval d = new Interval(40, 50);

        System.out.println("a = " + a);
        System.out.println("b = " + b);
        System.out.println("c = " + c);
        System.out.println("d = " + d);

        System.out.println("b intersects a = " + b.intersects(a));
        System.out.println("a intersects b = " + a.intersects(b));
        System.out.println("a intersects c = " + a.intersects(c));
        System.out.println("a intersects d = " + a.intersects(d));
        System.out.println("b intersects c = " + b.intersects(c));
        System.out.println("b intersects d = " + b.intersects(d));
        System.out.println("c intersects d = " + c.intersects(d));
    }

}