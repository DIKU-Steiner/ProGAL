package ProGAL.geom2d.viewer;

import java.awt.Graphics2D;

import ProGAL.geom2d.viewer.J2DScene.ShapeOptions;

interface ShapePainter {
	void paintShape(ShapeOptions shape, Graphics2D g2d);
}
