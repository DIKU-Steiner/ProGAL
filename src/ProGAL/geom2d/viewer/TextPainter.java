package ProGAL.geom2d.viewer;

import java.awt.BasicStroke;
import java.awt.Graphics2D;
import java.awt.Stroke;

import ProGAL.geom2d.viewer.J2DScene.ShapeOptions;

class TextPainter implements ShapePainter {

	public void paintShape(ShapeOptions shape, Graphics2D g2d) {

		TextShape text = (TextShape)shape.shape;
		java.awt.Point gPoint = shape.transformPoint(text.getPos());
		g2d.setFont(g2d.getFont().deriveFont((float)(text.getHeight()*shape.getScale())));
		g2d.setColor(shape.color);

		Stroke oldStroke = g2d.getStroke();
		g2d.setStroke(new BasicStroke((float)(shape.getScale()*shape.borderWidth)));
		g2d.drawString(text.getText(), gPoint.x, gPoint.y);
		g2d.setStroke(oldStroke);
	}


}
