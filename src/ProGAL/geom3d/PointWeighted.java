package ProGAL.geom3d;

public class PointWeighted extends Point {
	private static final long serialVersionUID = -7734855058026398183L;
	
	private double weight;
	
	public PointWeighted(double x, double y, double z, double w) {
		super(x, y, z);
		this.weight = w;
	}
	
	public PointWeighted(Point p){
		super(p);
		if(p instanceof PointWeighted) 
			this.weight = ((PointWeighted)p).weight;
		else
			this.weight = 1;
	}

	public double getWeight(){ return weight;	}

	public void setWeight(double w){ this.weight = w; }
	
}
