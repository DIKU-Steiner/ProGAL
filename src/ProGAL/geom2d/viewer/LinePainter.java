package ProGAL.geom2d.viewer;

import java.awt.Graphics2D;
import java.awt.Rectangle;

import ProGAL.geom2d.Line;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.viewer.J2DScene.ShapeOptions;

class LinePainter implements ShapePainter {

	public void paintShape(ShapeOptions shape, Graphics2D g2d) {

		Line l = (Line)shape.shape;
		
		Rectangle rect = g2d.getClipBounds();
//		double d1,d2;
		Point p1 = shape.transformPoint(new java.awt.Point(rect.x,rect.y));
		Point p2 = shape.transformPoint(new java.awt.Point(rect.x+rect.width,rect.y+rect.height));
		double d1 = l.getDirection().normalize().dot(l.getPoint(0).vectorTo(p1));
		double d2 = l.getDirection().normalize().dot(l.getPoint(0).vectorTo(p2));
		if(d1>d2) {double tmp = d1; d1 = d2; d2=tmp; }
		double extra = p2.distanceSquared(p1);
		d1-=extra/Math.min(1,l.getDirection().length());
		d2+=extra/Math.min(1,l.getDirection().length());
		
		java.awt.Point gPoint0 = shape.transformPoint(l.getPoint(d1));
		java.awt.Point gPoint1 = shape.transformPoint(l.getPoint(d2));
		g2d.drawLine(gPoint0.x, gPoint0.y, gPoint1.x, gPoint1.y);
	}


}
