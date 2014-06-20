package ProGAL.geom2d.viewer;

import ProGAL.geom2d.Point;
import ProGAL.geom2d.Shape;
import ProGAL.geomNd.Vector;

public class TextShape implements Shape {
	private String text;
	private double height;
	private Point pos;
	private double angle;
	
	public TextShape(String text, Point pos, double height){
		this(text, pos, height, 0.0);
	}
	
	public TextShape(String text, Point pos, double height, double angle){
		this.text = text;
		this.height = height;
		this.pos = pos;
		this.angle = angle;
	}
	
	public double getHeight(){ return height; }
	public String getText(){ return text; }
	public Point getPos(){ return pos; }
	public double getAngle(){ return angle; }
	
	@Override
	public Point getCenter() {
		return pos.clone();
	}

	@Override
	public boolean contains(Point p) {
		Vector v = pos.vectorTo(p);
		double dX = v.get(0);
		double dY = v.get(1);
		return (dX>0 && dX<height*text.length()/2.0 && dY>0 && dY<height);
	}

}
