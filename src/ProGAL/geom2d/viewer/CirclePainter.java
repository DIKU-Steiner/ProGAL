package ProGAL.geom2d.viewer;

import java.awt.BasicStroke;
import java.awt.Graphics2D;
import java.awt.Stroke;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.viewer.J2DScene.ShapeOptions;

class CirclePainter implements ShapePainter {

	public void paintShape(ShapeOptions shape, Graphics2D g2d) {

		Circle circle = (Circle)shape.shape;
		java.awt.Point gPoint = shape.transformPoint(circle.center());
		int diam = (int)(circle.getRadius()*2*shape.getScale());
		g2d.setColor(shape.color);

		Stroke oldStroke = g2d.getStroke();
		g2d.setStroke(new BasicStroke((float)(shape.getScale()*shape.borderWidth)));
		if(shape.fill)
			g2d.fillArc(gPoint.x-diam/2, gPoint.y-diam/2, diam,diam, 0,360);
		else
			g2d.drawArc(gPoint.x-diam/2, gPoint.y-diam/2, diam,diam, 0,360);
		g2d.setStroke(oldStroke);
	}


}
