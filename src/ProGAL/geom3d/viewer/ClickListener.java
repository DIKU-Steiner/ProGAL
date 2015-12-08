package ProGAL.geom3d.viewer;

import java.awt.event.MouseEvent;

import ProGAL.geom3d.Shape;

public interface ClickListener {
	void shapeClicked(Shape shape, MouseEvent event);
}
