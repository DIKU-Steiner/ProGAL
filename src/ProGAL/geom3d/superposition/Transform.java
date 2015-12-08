package ProGAL.geom3d.superposition;

import java.util.ArrayList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.math.Matrix;

public class Transform {
	private Matrix rotation;
	private Vector preTranslation, postTranslation;

	public Transform(Matrix rotation, Vector preTranslation, Vector postTranslation){
		this.rotation = rotation;
		this.preTranslation = preTranslation;
		this.postTranslation = postTranslation;
	}

	public Transform(Matrix rotation, Vector translation){
		this.rotation = rotation;
		this.preTranslation = null;
		this.postTranslation = translation;
	}
	
	public List<Point> transform(List<Point> pl) {
		List<Point> ret = new ArrayList<Point>(pl.size());
		for(Point p: pl)
			ret.add(transform(p));
		return ret;
	}
	
	public List<Point> transformIn(List<Point> pl) {
		for(Point p: pl)
			transformIn(p);
		return pl;
	}
	
	public Point transform(Point p){
		Point newPoint = p.clone();
		if(preTranslation!=null) newPoint.addThis(preTranslation);
		rotation.multiplyIn(newPoint);
		newPoint.addThis(postTranslation);
		return newPoint;
	}
	public Point transformIn(Point p){
		if(preTranslation!=null) p.addThis(preTranslation);
		rotation.multiplyIn(p);
		p.addThis(postTranslation);
		return p;
	}

	public String toString(int dec){
		StringBuilder sb = new StringBuilder();
		sb.append("Pre-translation: ");
		sb.append(preTranslation.toString(dec));
		sb.append('\n');
		sb.append("Rotation: \n");
		sb.append(rotation.toString(dec));
		sb.append('\n');
		sb.append("Post-translation: ");
		sb.append(postTranslation.toString(dec));
		return sb.toString();
	}
	
	public void toConsole(int dec) {
		System.out.println(toString(dec));
	}

	public void toConsole(){ toConsole(2); }
}
