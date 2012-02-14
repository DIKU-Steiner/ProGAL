package ProGAL.geom2d.viewer;

import java.awt.Graphics2D;

import ProGAL.geom2d.Triangle;
import ProGAL.geom2d.viewer.J2DScene.ShapeOptions;

class TrianglePainter implements ShapePainter {

	public void paintShape(ShapeOptions shape, Graphics2D g2d) {

		Triangle tri = (Triangle)shape.shape;
		java.awt.Point p0 = shape.transformPoint(tri.getCorner(0));
		java.awt.Point p1 = shape.transformPoint(tri.getCorner(1));
		java.awt.Point p2 = shape.transformPoint(tri.getCorner(2));
		
		if(shape.fill){
			g2d.fillPolygon(
					new int[]{p0.x, p1.x, p2.x}, 
					new int[]{p0.y, p1.y, p2.y}, 
					3);
		}else{
			g2d.drawPolygon(
					new int[]{p0.x, p1.x, p2.x}, 
					new int[]{p0.y, p1.y, p2.y},
					3);
		}
	}


}
