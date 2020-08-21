package emd_package;
import java.awt.Component;
import java.awt.geom.AffineTransform;
import java.awt.geom.*;
import java.awt.image.AffineTransformOp;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;



 public class picture extends Component {

	 int alpha ;
	 int red ;
	 int green ;
	 int blue;
	 int[][][] temp;
 public void printPixelARGB(int pixel) {
 

 //System.out.println("argb: " + alpha + ", " + red + ", " + green + ", " + blue);
 //System.out.println(pixel);
}

private void marchThroughImage(BufferedImage image) {
int w = image.getWidth();
 int h = image.getHeight();
 int pixel;
 temp=new int[256][256][256];;
// System.out.println("width, height: " + w + ", " + h);

for (int i = 0; i < h; i++) {
for (int j = 0; j < w; j++) {
//System.out.println("x,y: " + j + ", " + i);
pixel = image.getRGB(j, i);
// alpha = (pixel >> 24) & 0xff;
 red = (pixel >> 16) & 0xff;
 green = (pixel >> 8) & 0xff;
 blue = (pixel) & 0xff;
 temp[red][blue][green]=temp[red][blue][green]+1;
printPixelARGB(pixel);
 //System.out.println("value of K:"+k+ " value of  pixel: " + pixel[j][i]);
}
}

//System.out.println("loop is completed");

}
int[][][] value(){
	return temp;
	
}
public picture(String s) throws IOException {
//this is an image of a white spot on a black background.
//with the smoothing in the image it's of course not all black
//and white
BufferedImage image = ImageIO.read(new File(s));

final int w =(int) (image.getWidth()*0.05);
final int h = (int) (image.getHeight()*0.05);
BufferedImage scaledImage = new BufferedImage((w ),(h ), BufferedImage.TYPE_INT_ARGB);
final AffineTransform at = AffineTransform.getScaleInstance(0.05,0.05);
final AffineTransformOp ato = new AffineTransformOp(at, AffineTransformOp.TYPE_BICUBIC);
scaledImage = ato.filter(image, scaledImage);

marchThroughImage(image);
}
}