package ProGAL.geom2d.viewer;

import java.awt.BasicStroke;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.Stroke;

import ProGAL.geom2d.LSC;
import ProGAL.geom2d.viewer.J2DScene.ShapeOptions;

class LSCPainter implements ShapePainter {

	public void paintShape(ShapeOptions shape, Graphics2D g2d) {

		LSC lsc = (LSC)shape.shape;
		java.awt.Point gPoint0 = shape.transformPoint(lsc.getSegment().getA());
		java.awt.Point gPoint1 = shape.transformPoint(lsc.getSegment().getB());
		if(shape.fill){
			g2d.setStroke(new BasicStroke((float)(shape.getScale()*(lsc.getRadius()*2+shape.borderWidth)), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
		}else{
			Stroke str = new CompositeStroke( 
					new BasicStroke((float)(shape.getScale()*lsc.getRadius()*2), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND),
					new BasicStroke((float)(shape.getScale()*shape.borderWidth))  
					);
			g2d.setStroke(str);
		}
		g2d.drawLine(gPoint0.x, gPoint0.y, gPoint1.x, gPoint1.y);
	}


	private static class CompositeStroke implements Stroke {
		private Stroke stroke1, stroke2;

		CompositeStroke( Stroke stroke1, Stroke stroke2 ) {
			this.stroke1 = stroke1;
			this.stroke2 = stroke2;
		}

		public Shape createStrokedShape( Shape shape ) {
			return stroke2.createStrokedShape( stroke1.createStrokedShape( shape ) );
		}
	}
}
