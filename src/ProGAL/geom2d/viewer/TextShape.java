package ProGAL.geom2d.viewer;

import ProGAL.geom2d.Point;
import ProGAL.geom2d.Shape;

public class TextShape implements Shape {
	private String text;
	private float height;
	private Point pos;
	
	public TextShape(String text, Point pos, float height){
		this.text = text;
		this.height = height;
		this.pos = pos;
	}
	
	public float getHeight(){ return height; }
	public String getText(){ return text; }
	public Point getPos(){ return pos; }
	
	@Override
	public Point getCenter() {
		return pos.clone();
	}

}
