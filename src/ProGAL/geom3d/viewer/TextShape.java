package ProGAL.geom3d.viewer;

import ProGAL.geom3d.*;


public class TextShape implements Shape {
	public String text;
	public Point pos;
	public double height;
	public TextShape(String t, Point p){
		this(t,p,0.1f);
	}
	public TextShape(String t, Point p, double height){
		text = t;
		pos = p;
		this.height = height;
	}
	public Point getCenter() {
		return pos.clone();
	}
}
