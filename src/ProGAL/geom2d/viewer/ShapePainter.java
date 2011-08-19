package ProGAL.geom2d.viewer;

import java.awt.Graphics2D;

import ProGAL.geom2d.viewer.J2DScene.ShapeOptions;

public interface ShapePainter {
	void paintShape(ShapeOptions shape, Graphics2D g2d);
}
