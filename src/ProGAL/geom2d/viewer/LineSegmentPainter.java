package ProGAL.geom2d.viewer;

import java.awt.Graphics2D;

import ProGAL.geom2d.LineSegment;
import ProGAL.geom2d.viewer.J2DScene.ShapeOptions;

public class LineSegmentPainter implements ShapePainter {

	@Override
	public void paintShape(ShapeOptions shape, Graphics2D g2d) {

		LineSegment ls = (LineSegment)shape.shape;
		java.awt.Point gPoint0 = shape.transformPoint(ls.getA());
		java.awt.Point gPoint1 = shape.transformPoint(ls.getB());
		g2d.drawLine(gPoint0.x, gPoint0.y, gPoint1.x, gPoint1.y);
	}


}
