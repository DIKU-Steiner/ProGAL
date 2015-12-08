package ProGAL.geom2d.viewer;

import java.awt.Graphics2D;

import ProGAL.geom2d.Polygon;
import ProGAL.geom2d.viewer.J2DScene.ShapeOptions;

class PolygonPainter implements ShapePainter {

	public void paintShape(ShapeOptions shape, Graphics2D g2d) {

		Polygon pol = (Polygon)shape.shape;
		int[] xs = new int[pol.size()];
		int[] ys = new int[pol.size()];
		for(int i=0;i<xs.length;i++){
			java.awt.Point p = shape.transformPoint(pol.get(i));
			xs[i] = p.x;
			ys[i] = p.y;
		}
		
		if(shape.fill)	g2d.fillPolygon(xs,ys,xs.length);
		else			g2d.drawPolygon(xs,ys,xs.length);
	}

}
