package ProGAL.geom2d.viewer;

import java.awt.event.MouseEvent;

import ProGAL.geom2d.Shape;

public interface ClickListener {
	void shapeClicked(Shape shape, MouseEvent event);
}
