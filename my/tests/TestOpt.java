package my.tests;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.atan;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.tan;
import static my.tests.U.AS;
import static my.tests.U.AT;
import static my.tests.U.C;
import static my.tests.U.S;
import static my.tests.U.T;
import static my.tests.U.ac;
import static my.tests.U.as;
import static my.tests.U.defun;
import static my.tests.U.t;
import static my.tests.U.toa;
import static my.tests.U.trm;
import static my.tests.V.add;
import static my.tests.V.mR;
import static my.tests.V.mXY;
import static my.tests.V.mult;
import static my.tests.V.multD;
import static my.tests.V.sub;

import java.awt.Canvas;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Point;
import java.awt.geom.Area;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.lang.reflect.Method;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.swing.JFrame;

import org.omg.Messaging.SyncScopeHelper;

class Params{
	double z, Z, h, r;
	double R;

	double A;
	double n = 1;
	
	V centre;
	V obj;
	

	public Params(double z, double Z, double h, double r, double n, double A){
		this.z = z;
		this.Z = Z;
		this.h = h;
		this.r = r;
		this.A = A;
		this.n = n;
		calR();
		centre = new V(z*T(A), 0, z);
		obj = new V((z+Z)*T(A), 0, z+Z);
	}
	
	public Params init(double z, double Z, double h, double r, double n, double A){
		this.z = z;
		this.Z = Z;
		this.h = h;
		this.r = r;
		this.A = A;
		this.n = n;
		calR();
		return this;
	}
	
	void calR(){
		R = h*t(as(1/n))+r;
	}
	
	public V tr(V v){
		double mXY = mXY(v);
		double xy = h * tan(as(mXY/mR(v)/n)) ;
		return new V(v.x * xy/mXY, v.y * xy/mXY, h);
	}

	public V atr(V v){
		double mXY = mXY(v);
		if(mXY>(R-r)){
			return null;
		}
		double xy = Z * tan(as(mXY/mR(v)*n)) ;
		return new V(v.x * xy/mXY, v.y * xy/mXY, Z);
	}

	public V border_0(double B){
		V p = new V(C(B)*r, S(B)*r, 0);
		V li = add(centre, p);
		V im = add(li, tr(li));
		V li2 = sub(im, centre);
		V obAb = atr(li2);		
		if(obAb == null){
			return null;
		}
		V ob = add(centre, obAb);

		return mult(ob, z/(z+Z));
	}
}

class RDM extends Random{
	public static Random R = new Random();
	public static int nextI(int i, int j){
		if(i > j){
			i = i+j;
			j = i-j;
			i = i-j;
		}
		return R.nextInt(j-i) +i;
	}
	
	public static double nextIC(int c, int r){
		return R.nextInt(2*r)-r+c;
	}

	public static double nextD(double i, double j){
		return nextDC((i+j)/2, (i-j)/2);
	}

	public static double nextDC(double c, double r){
		return R.nextDouble()*2*r-r+c;
	}

};


class Drawer{
	static Color cl1 = new Color(240,248,248);
	static Color cl2 = new Color(255,240,240);
	static Color cl3 = new Color(252,250,250);
	static Color cl4 = new Color(240,240,240);
	
	static Color cla = new Color(250,250,255);
	static Color clb = new Color(250,200,200);

	int w, h;
	Graphics2D pt;
	
	public int getX(double x){
		return new Double(x+w/2).intValue();
	}
	
	public int getY(double y){
		return new Double(y+h/2).intValue();
	}
	
	public int inversY(int y){
		return h-y;
	}
	
	public int getI(double i){
		return new Double(i).intValue();
	}
	
	public Drawer(int w, int h, Graphics g){
		this.w = w;
		this.h = h;
		this.pt = (Graphics2D)g;
		g.setColor(cl1);

	}
	
	public void drPoint(double x, double y){
		int px = getX(x);
		int py = getY(y);
		pt.drawLine(px, inversY(py), px, inversY(py));
	}

	public void drPointXY(V v){
		int px = getX(v.x);
		int py = getY(v.y);
		pt.drawLine(px, inversY(py), px, inversY(py));
	}

	public void drPointYZ(V v){
		int px = getX(v.y);
		int py = getY(v.z);
		pt.drawLine(px, inversY(py), px, inversY(py));
	}

	public void drPointXZ(V v){
		int px = getX(v.x);
		int py = getY(v.z);
		pt.drawLine(px, inversY(py), px, inversY(py));
	}


	public void drCircle(double ox, double oy, double r){
		int oX = getX(ox);
		int oY = getY(oy);
		int R = getI(r);
		pt.drawOval(oX-R, inversY(oY+R), R*2, R*2);
	}
	
	public void drOval(double ox, double oy, double a, double b){
		int oX = getX(ox);
		int oY = getY(oy);
		int A = getI(a);
		int B = getI(b);
		pt.drawOval(oX-A, inversY(oY+B), A*2, B*2);
	}
	
	public void drRect(double ox, double oy, double width, double height){
		int oX = getX(ox);
		int oY = getY(oy);
		int wI = getI(width);
		int hI = getI(height);
		pt.drawRect(oX-wI/2, inversY(oY-hI/2), wI, hI);
	}
	
	public void fiRect(double ox, double oy, double width, double height){
		int oX = getX(ox);
		int oY = getY(oy);
		int wI = getI(width);
		int hI = getI(height);
		pt.fillRect(oX-wI/2, inversY(oY-hI/2), wI, hI);
	}
	
	public void flRect(double ox, double oy, double width, double height){
		int oX = getX(ox);
		int oY = getY(oy);
		int wI = getI(width);
		int hI = getI(height);
		pt.fillRect(oX-wI/2, inversY(oY-hI/2), wI, hI);
	}
	
	public void drLine(double x1, double y1, double x2, double y2){
		pt.drawLine(getX(x1), inversY(getY(y1)), getX(x2), inversY(getY(y2)));
	}

	public void drLXY(L l, double length){
		V v = add(l.o, mult(l.dir, length));
		drLineXY(l.o, v);
	}
	
	public void drLXZ(L l, double length){
		V v = add(l.o, mult(l.dir, length));
		drLineXZ(l.o, v);
	}
	
	public void drLYZ(L l, double length){
		V v = add(l.o, mult(l.dir, length));
		drLineYZ(l.o, v);
	}
	
	public void drSegXY(LSeg sg) {
		drLine(sg.p1.x, sg.p1.y, sg.p2.x, sg.p2.y);
//		drLXY(new L(sg.c, sg.dirN), 20);
	}

	public void drPL(List<LSeg> pl) {
		for(LSeg sg : pl){
			drSegXY(sg);
		}
		
	}
	
	public void drStr(Object val, double x, double y){
		pt.drawString(val.toString(), getX(x), inversY(getY(y)));
	}

	public void drLineXY(V v1, V v2){
		drLine(v1.x, v1.y, v2.x, v2.y);
	}

	public void drLineXZ(V v1, V v2){
		drLine(v1.x, v1.z, v2.x, v2.z);
	}
	
	public void drLineYZ(V v1, V v2){
		drLine(v1.y, v1.z, v2.y, v2.z);
	}

	public void drX(double x){
		drLine(x, -100*h, x, 100*h);
	}

	public void drY(double y){
		drLine(-100*w, y, 100*w, y);
	}

	// pp
	public V getVPP(V v){
//		return new V(v.x+v.y/4*Math.sqrt(2), 0, v.z+v.y/4*Math.sqrt(2));
		return new V(v.x, 0, v.z+v.y/4*Math.sqrt(2));
	}
	
	public void ppPoint(V v){
		drPointXZ(getVPP(v));
	}
	
	public void ppLine(V v1, V v2){
		drLineXZ(getVPP(v1), getVPP(v2));
	}

	public void ppLine(L l, double length){
		ppLine(l.o, add(l.o, mult(l.dir, length)));
	}
	
	public void ppCicleZ(L ax, V r){
		for(double a = 0; a<360; a++){
			V vi = V.ROLL(ax, r, a);
			ppPoint(vi);
		}
	}
	
	
	
	
	public void drFunXY(FunReal fun, V pO, V scale, double minX, double maxX, double dx){
		pt.setColor(new Color(0.5f,0.5f,0.5f,0.5f));
		drX(pO.x);
		drY(pO.y);

		drX(minX*scale.x+pO.x);
		drX(maxX*scale.x+pO.x);

		Double y0 = null;
		Double x0 = null;
		
		double x = minX;
		for(; x<maxX; x+=dx){
			double[] ys = fun.fun(new double[]{x});
			Double y = ys==null ? null : ys[0];
			
			if(y!=null){
				pt.setColor(Color.gray);
				if(x0!=null && y0!=null){
					drLine(x0*scale.x+pO.x, y0*scale.y+pO.y, x*scale.x+pO.x, y*scale.y+pO.y);
				}else{
					drPoint(x*scale.x+pO.x, y*scale.y+pO.y);
				}
			}else{
				pt.setColor(Color.pink);
				drPoint(x*scale.x+pO.x, pO.y);
			}
			
			x0 = x;
			y0 = y;
		}

		
		x = maxX;
		double[] ys = fun.fun(new double[]{x});
		Double y = ys==null ? null : ys[0];
		
		if(y!=null){
			pt.setColor(Color.gray);
			if(x0!=null && y0!=null){
				drLine(x0*scale.x+pO.x, y0*scale.y+pO.y, x*scale.x+pO.x, y*scale.y+pO.y);
			}else{
				drPoint(x*scale.x+pO.x, y*scale.y+pO.y);
			}
		}else{
			pt.setColor(Color.pink);
			drPoint(x*scale.x+pO.x, pO.y);
		}

		
	}
	
	public void drFun(FunReal fun, V pO, V scale){
		pt.setColor(new Color(0.5f,0.5f,0.5f,0.5f));
		drX(pO.x);
		drY(pO.y);

		fun.show(this, pO, scale);
	}
	
	public Color infoColor;
	public String infoString;
	
	public void show(){
//		Gra gra = new Gra(11, 20, 1);		
//		pt.setClip(gra.getArea(w/2, h/2, 2, 1, new V(0,0,2)));
		
//		Area a = new Area(new Rectangle2D.Float(4.5f,10,1,100));
//		a.add(new Area(new Rectangle2D.Float(4f,120,1,100)));
//		pt.setClip(a);
		
//		testClip();
		
//		testGroup();
		
		// tmp show oval
//		int a = 100; int b = 80;
//		double x0 = 0; 
//		double y0 = 0; 
//		for(double x=0; x<=a; x++){
//			double y = b*Math.pow(1-x*x/a/a, 0.5);
//			drLine(x0, y0, x, y);
//			if(x%2==0){
//				drLine(x, y, (1+b*b)*x, (1+a*a)*y);
//			}
//			x0 = x;
//			y0 = y;
//		}
		
		
//		test4();
//		test6();
		
//		test7();
//		test7V2();


		test7V3();
//		test7V3D();

//		test7V3DInverse();
		
//		testSegSurfaceCross();
//		testSegSurfaceCrossZ();
//		testSegSurfaceCrossZInverse();
//		testSegSurfaceCrossZPP();
		
//		testSegSurface();
//		statisticA();

		
//		testIPDirRec();
//		testSfcCllMtx();
		
		
		
//		showFunction();
//		showTrm();
		
//		showHSQX();
		
//		testPoly();
		
//		testADiscrete();
		

		
		
		
		
		
//		testSun();
//		testImg();
//		testFindR();
		
		
		if(infoColor!=null){
			pt.setColor(infoColor);
		}
		if(infoString!=null){
			pt.drawString(infoString, 10, 10);
		}

		
	}

	static Color[] objColors = {Color.darkGray, Color.red, Color.green, Color.blue, Color.pink, Color.lightGray, Color.magenta, Color.orange};
	
	
	int resolution = 200;
	double pixelRealSize = 0.10/1080;		
//	double pixelRealSize = 0.035/1080;		
	double cellLength = resolution*pixelRealSize ;
	double objRealL = cellLength * 20;
	double objReall = cellLength * 10;
	
	double objRealH = -15;
	double eyeRealH = 20;

	int objL = mToPixel(objRealL);
	int objl = mToPixel(objReall);

//	int objL = mToCellCount(objRealL);
//	int objl = mToCellCount(objReall);

	
	static Date crtTime = new Date();
//	static String imgPath = "C:/Users/lchen/tests/Equirectangular_projection_SW.jpg";
//	static String imgPath = "C:/Users/lchen/tests/soho2.jpg";
	static String imgPath = "C:/Users/lchen/tests/face.jpg";
//	static String imgPath = "C:/Users/lchen/tests/tree.jpg";
//	static String imgPath = "C:/Users/lchen/tests/round.jpg";

	static int edgeBorder = 5;
	
	static float DH = 0.05f; 
	static float DS = 0.4f;
	static float DB = 0.8f;

	static int objectMinR = 10;
	
	static double[] target = null;
	Thread thd = new Thread(){
		File file= new File(imgPath);
		BufferedImage image;
		Color cF = new Color(238,196,170);
		public void run(){
			while(true){
				doDetect();
				try {
					sleep(50);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		public void doDetect() {
			try {
				image = ImageIO.read(file);
			} catch (IOException e) {
				return;
			}		
			long t1 = System.currentTimeMillis();
			target = findTarget(cF, image);
//			double[] target = {0, 0, 50};
			long t2 = System.currentTimeMillis();
			System.out.println(t2-t1);
		}
	};

	
	public void testFindR(){
//		Color cF = Color.black;
//		Color cF = Color.red;
//		Color cF = new Color(35,177,77);
		Color cF = new Color(238,196,170);
//		Color cF = Color.yellow;
//		Color cF = new Color(254,230,0);
//		Color cF = new Color(54,33,30);

		
		if(!thd.isAlive()){
			thd.start();
		}
		
		pt.setColor(Color.white);
		pt.fillRect(0, 0, 500,500);
		if(target!=null){
			pt.setColor(Color.red);
			drCircle(target[0], target[1], 1);
			drCircle(target[0], target[1], 2);
			drCircle(target[0], target[1], 3);
			drCircle(target[0], target[1], target[2]);
		}
		
		
		
		
		
		
		
		
		
//		long t1 = System.currentTimeMillis();
//		
//		pt.setColor(Color.red);
//
//		List<int[]> edges = new ArrayList<int[]>(); 
//		for(int i=edgeBorder; i<image.getWidth()-edgeBorder; i++){
//			for(int j=edgeBorder; j<image.getHeight()-edgeBorder; j++){
//				if(isEdge(image, i, j, cF)){
//					edges.add(new int[]{i, j});
////					drPoint(i-image.getWidth()/2, j-image.getHeight()/2);
//					drPoint(i-image.getWidth()/2, j-image.getHeight()/2);
//				}
//			}
//		}
//		
//		long t2 = System.currentTimeMillis();
//
//		int[][] R = findR(edges);
//
//		long t3 = System.currentTimeMillis();
//		System.out.println(t2-t1);
//		System.out.println(t3-t2);
//		System.out.println(edges.size());
//		
//		drLineXY(new V(R[0][0]-image.getWidth()/2, R[0][1]-image.getHeight()/2, 0), new V(R[1][0]-image.getWidth()/2, R[1][1]-image.getHeight()/2, 0));
	}

	public  double[] findTarget(Color cF, BufferedImage image) {
		
		List<int[]> sameColors = new ArrayList<int[]>(); 
		List<Double> dists = new ArrayList<Double>(); 
		int[] distValues = new int[Math.max(image.getWidth(), image.getHeight())/5];
		Double centerX = 0d;
		Double centerY = 0d;
		Double everage = 0d;
		int maxAvailableDist = objectMinR;

		// center
		for(int i=0; i<image.getWidth(); i++){
			for(int j=0; j<image.getHeight(); j++){
				if(isSameColor(cF, new Color(image.getRGB(i, j)))){
					sameColors.add(new int[]{i, j});
					centerX += i;
					centerY += j;

					//test
//					pt.setColor(new Color(image.getRGB(i, j)));
////					pt.setColor(Color.blue);
//					drPoint(i-image.getWidth()/2, j-image.getHeight()/2);
				}else{
//					Color ci = new Color(image.getRGB(i, j));
//					float[] fi = Color.RGBtoHSB(ci.getRed(), ci.getGreen(), ci.getBlue(), new float[3]);
//
//					pt.setColor(new Color(Color.HSBtoRGB(fi[0], fi[1], fi[2]/1.5f)));
//					drPoint(i-image.getWidth()/2, j-image.getHeight()/2);
				}
			}
		}
		if(sameColors.size()==0){
			System.out.println("no same color");
			return null;
		}
		centerX /= sameColors.size();
		centerY /= sameColors.size();
		
		// dists, everage
		for(int i=0; i<sameColors.size(); i++){
			int[] p = sameColors.get(i);
			Double dist = Math.sqrt(Math.pow(p[0]-centerX, 2) + Math.pow(p[1]-centerY, 2));
			if(dist<=distValues.length){
				distValues[dist.intValue()] = 1;
			}
			dists.add(dist);
			everage += dist;
		}
		everage /= sameColors.size();

		// available max dist
		for(int i=objectMinR; i<distValues.length-3 && i<everage*2; i++){
			maxAvailableDist = i;
			if(distValues[i+3]==0 && distValues[i+2]==0 && distValues[i+1]==0){
				break;
			}
		}
		
		// center again
		double newCenterX = 0d;
		double newCenterY = 0d;
		for(int i=0; i<sameColors.size(); ){
			int[] p = sameColors.get(i);
			Double dist = Math.sqrt(Math.pow(p[0]-centerX, 2) + Math.pow(p[1]-centerY, 2));
			if(dist<=maxAvailableDist+1){
				newCenterX += p[0];
				newCenterY += p[1];
				
				//test
//				pt.setColor(new Color(image.getRGB(p[0], p[1])));
//				pt.setColor(cF);
//				drPoint(p[0]-image.getWidth()/2, p[1]-image.getHeight()/2);
				
				i++;
			}else{
				sameColors.remove(i);
			}
		}
		if(sameColors.size()==0){
			return null;
		}
		centerX = newCenterX/sameColors.size();
		centerY = newCenterY/sameColors.size();

		// max dists
		Double maxDist = 0d;
		for(int i=0; i<sameColors.size(); i++){
			int[] p = sameColors.get(i);
			Double dist = Math.sqrt(Math.pow(p[0]-centerX, 2) + Math.pow(p[1]-centerY, 2));
			maxDist = Math.max(maxDist, dist);
		}

		if(maxDist<0.5)return null;
		
		
		System.out.println(centerX);
		System.out.println(centerY);
		System.out.println(maxDist);
		
//		return new double[]{centerX-image.getWidth()/2, centerY-image.getHeight()/2, maxDist};
		return toTargetPosition(centerX-image.getWidth()/2, centerY-image.getHeight()/2, maxDist);
	}

	public static double camH = 500;
	public static double objR = 100;
	public double[] toTargetPosition(double x, double y, Double objr) {
		double r = Math.sqrt(x*x+y*y);
		double Z = camH * objR/objr;
		double R = r * objR/objr;
		double X = x/r*R;
		double Y = y/r*R;
		
		return new double[]{X, Y, Z};
//		return new double[]{-X, Y-screenY/2, -Z};
	}

	public boolean isEdge(BufferedImage image, int i, int j, Color cF, Color cB) {
		Color co = new Color(image.getRGB(i, j));
		return (isSameColor(co, cF) 
			&& (
				(isSameColor(co, new Color(image.getRGB(i+1, j))) != isSameColor(co, new Color(image.getRGB(i-1, j))))
			 || (isSameColor(co, new Color(image.getRGB(i, j+1))) != isSameColor(co, new Color(image.getRGB(i, j-1))))
			)
		);
	}

	public boolean isEdge(BufferedImage image, int i, int j, Color cF) {
		Color co = new Color(image.getRGB(i, j));
		return (isSameColor(co, cF) 
			&& (
				(isSameColor(co, new Color(image.getRGB(i+edgeBorder, j))) != isSameColor(co, new Color(image.getRGB(i-edgeBorder, j))))
			 || (isSameColor(co, new Color(image.getRGB(i, j+edgeBorder))) != isSameColor(co, new Color(image.getRGB(i, j-edgeBorder))))
			)
		);
	}

	public int[][] findR(List<int[]> edges) {
		int[][] R = null;
		double dist = 0;
		
		for(int i=1; i<edges.size(); i+=10){
			for(int j=0; j<i; j+=10){
				int[] p1 = edges.get(i);
				int[] p2 = edges.get(j);
				double distI = Math.pow(p1[0]-p2[0], 2) + Math.pow(p1[1]-p2[1], 2); 
				if(distI>dist){
					dist = distI;
					R = new int[][]{p1, p2};
				}
				
			}
		}
		
		return R;
	}

	public int countSameColor(Color[] cs, Color c) {
		int cnt = 0;
		for(Color ci : cs){
			if(isSameColor(ci,  c)){
				cnt++;
			}
		}
		return cnt;
	}

//	public boolean isSameColor(Color c1, Color c2){
//		return Math.abs(c1.getRed()-c2.getRed())<20 
//				&& Math.abs(c1.getGreen()-c2.getGreen())<20
//				&& Math.abs(c1.getBlue()-c2.getBlue())<20;
//	}
	
	public boolean isSameColor(Color c1, Color c2){
		float[] fc1 = Color.RGBtoHSB(c1.getRed(), c1.getGreen(), c1.getBlue(), new float[3]);
		float[] fc2 = Color.RGBtoHSB(c2.getRed(), c2.getGreen(), c2.getBlue(), new float[3]);
		
		float dH = Math.abs(fc1[0]-fc2[0]);
		dH = Math.min(dH, 1-dH);
		return dH<DH 
				&& Math.abs(fc1[1]-fc2[1])<DS
				&& Math.abs(fc1[2]-fc2[2])<DB;
	}
	
	

	public void testImg(){
		Drawer drr = new Drawer(w, h, pt);
		
		File file= new File(imgPath);
		BufferedImage image;
		try {
			image = ImageIO.read(file);
		} catch (IOException e) {
			return;
		}		
		
		pt.setColor(Color.black);
		fiRect(0, 0, w, h);
		int step = 10;
		int size = 1;
		for(int i=-w/2; i<=w/2; i++){
			for(int j=-h/2; j<=h/2; j++){
				if(i%step!=0 || j%step!=0){
					continue;
				}
				try{
					
					int cR = 0;
					int cG = 0;
					int cB = 0;

					int ct = 0;
					for(int si = 0; si<step; si++){
						for(int sj = 0; sj<step; sj++){
							Color c = new Color(image.getRGB(i+image.getWidth()/2+si, j+image.getHeight()/2)+sj);
							if(c != null){
								cR += c.getRed();
								cG += c.getGreen();
								cB += c.getBlue();
								ct++;
							}
						}
						
					}
					
					pt.setColor(new Color(cR/ct, cG/ct, cB/ct));
//					pt.setColor(new Color(cR, cG, cB));
//					fiRect(i*size/step, j*size/step, size, size);
					fiRect(i, j, size, size);
				}catch(Exception e){}
			}			
		}
		
	}

	public void testSun(){
		Sun sun = new Sun();
		
		Drawer drr = new Drawer(w, h, pt);
		
		int r = Math.max(w, h)/4-10;
		
		drr.drImg(r*2, imgPath, crtTime);
		

		pt.setColor(Color.black);
//		drr.drCircle(0, 0, r);
		drr.drLine(-r*2, 0, r*2, 0);
		
		double suna = sun.getSuna(crtTime);
		pt.setColor(Color.green);
		drr.drCircle(0, 0, 1);
		drr.drCircle(0, 0, r*2);
		pt.setColor(Color.red);
		V last = null;
		for(int i=0; i<=36000; i+=5){
			double b = U.toa(i*0.01);
			V si = sun.spToPl(sun.getSunP(suna, b));
			if(si!=null){
//				drr.drPoint(si.x*r, si.y*r);
				if(last != null){
					drr.drLine(si.x*r, si.y*r, last.x*r, last.y*r);
				}
				last = si;
			}
		}
		
		Calendar cal = Calendar.getInstance();
		cal.setTime(crtTime);
		cal.add(Calendar.HOUR, 48);
		crtTime = cal.getTime();

	}

	
	private void drImg(int r, String imgPath, Date tm) {
		System.out.println(r);
		File file= new File(imgPath);
		BufferedImage image;
		try {
			image = ImageIO.read(file);
		} catch (IOException e) {
			return;
		}		
		
		double rd = new Double(r);
		for(int i=-r; i<=r; i++){
			for(int j=-r; j<=r; j++){
				V p = Sun.plToTude(i/rd, j/rd, tm);
				if(p==null){
					continue;
				}
				
				int x = new Double(image.getWidth()*p.x).intValue();
				int y = new Double(image.getHeight()*p.y).intValue();
				y = Math.min(image.getHeight()-1, y);
				
				try{
					Color c = new Color(image.getRGB(x, y));
					pt.setColor(c);
					drPoint(i, j);
				}catch(Exception e){}
			}			
		}
	}

	public void testADiscrete(){
		// dddddddddddddddd

		double objH = mToPixel(objRealH);
		double eyeH = mToPixel(-eyeRealH);
		
		
		Drawer drr = new Drawer(w, h, pt);
//		pt.setColor(Color.black);		
		pt.setColor(Color.white);		
		drr.flRect(0d, 0d, w, h);
		
		int step = 1;
		PlView plR = new PlEmpty(new V(0,0,0), new V(0,0,1), step);
		PlView plI = new PlDisc(new V(0,0,0), new V(0,0,1), step, resolution);
		V o = new V(0,0,eyeH);
		
		System.out.println("w:"+new Double(cellCountToM(w)*10+0.5).intValue()/10d);
		System.out.println("cellLength:"+new Double(cellLength*1000+0.5).intValue()/1000d);
		System.out.println("obj:"+new Double(pixelToM(objH)*10+0.5).intValue()/10d +"("+objH +")");
		System.out.println("eye:"+new Double(pixelToM(-eyeH)*10+0.5).intValue()/10d);
		System.out.println("objl:"+new Double(pixelToM(objl)*1000+0.5).intValue()/1000d);
		System.out.println("objl in eye:"+AT(objl/(objH-eyeH))+"\b in img:"+pixelToM(objl/(eyeH-objH)*eyeH));
		
		
//		pt.setColor(Color.red);		
//		drr.drLine(-w, pl.o.z, w, pl.o.z);
//		pt.setColor(Color.lightGray);
		
		// draw plan xz:
//		PlDisc pl = plI;
//		for(int i=-w/2; i<=w/2; i+=plI.step){
//			V c = new V(i, pl.o.y, pl.o.z);
//			L in = new L(o, sub(c, o));
//			L out = pl.out(in, c, 0);
//			drr.drLineXZ(o, c);
//			drr.drLXZ(out, 1000);
//		}
		
		// draw obj : 
//		for(int i=-w/2; i<=w/2; i++){
//			for(int j=-h; j<=h; j++){
//				pt.setColor(discreteObjColor(i, j));
//				drr.drPoint(i, j);
//			}
//		}

		int mode = 0; // 1:theoric 2: cell as center
		
		// draw result
		for(int i=-w/2; i<=w/2; i+=plR.step){
//			for(int j=-h; j<=h; j+=plR.step){
			int j = 0;
				V c = new V(i*resolution, j*resolution, plR.o.z);

				if(true || j>h/6 || mode==1){// real
					int[] sumX = {0,0,0};
					int[] sumY = {0,0,0};
					V pTL = sub(c, new V(resolution/2-0.5, resolution/2-0.5, 0));
					for(int a = 0; a<resolution; a++){
						V caa = add(pTL, new V(a, a, 0));
						L in = new L(o, sub(caa, o));							
						L out = plR.out(in, caa, 0);
						Double objX = out.o.x+objH*out.dir.x/out.dir.z;
						Double objY = out.o.y+objH*out.dir.y/out.dir.z;
						Color clx = discreteObjColorX(objX.intValue());					
						Color cly = discreteObjColorY(objY.intValue());
						sumX[0] = sumX[0] + clx.getRed();
						sumX[1] = sumX[1] + clx.getGreen();
						sumX[2] = sumX[2] + clx.getBlue();
						sumY[0] = sumY[0] + cly.getRed();
						sumY[1] = sumY[1] + cly.getGreen();
						sumY[2] = sumY[2] + cly.getBlue();
					}
					pt.setColor(new Color(
							averageColor(sumX[0]/resolution, sumY[0]/resolution),
							averageColor(sumX[1]/resolution, sumY[1]/resolution),
							averageColor(sumX[2]/resolution, sumY[2]/resolution)
							));
					drr.drLine(i, -h/2, i, -h/6);

				}
				
				if(true || abs(j)<=h/6){
					// discrete theoric
					int[] sumX = {0,0,0};
					int[] sumY = {0,0,0};
					V pTL = sub(c, new V(resolution/2-0.5, resolution/2-0.5, 0));
					for(int a = 0; a<resolution; a++){
						V caa = add(pTL, new V(a, a, 0));
						L in = new L(o, sub(caa, o));							
						L out = plI.out(in, caa, 0);
						Double objX = out.o.x+objH*out.dir.x/out.dir.z;
						Double objY = out.o.y+objH*out.dir.y/out.dir.z;
						Color clx = discreteObjColorX(objX.intValue());					
						Color cly = discreteObjColorY(objY.intValue());
						sumX[0] = sumX[0] + clx.getRed();
						sumX[1] = sumX[1] + clx.getGreen();
						sumX[2] = sumX[2] + clx.getBlue();
						sumY[0] = sumY[0] + cly.getRed();
						sumY[1] = sumY[1] + cly.getGreen();
						sumY[2] = sumY[2] + cly.getBlue();
					}
					pt.setColor(new Color(
							averageColor(sumX[0]/resolution, sumY[0]/resolution),
							averageColor(sumX[1]/resolution, sumY[1]/resolution),
							averageColor(sumX[2]/resolution, sumY[2]/resolution)
							));
					drr.drLine(i, -h/6, i, h/6);
				}

				if(true || j<-h/6 || mode==2){
					// discrete cell as center
					int[] sumX = {0,0,0};
					int[] sumY = {0,0,0};
					V pTL = sub(o, new V(resolution/2-0.5, resolution/2-0.5, 0));
					for(int a = 0; a<resolution; a++){
						V oaa = add(pTL, new V(a, a, 0));
						L in = new L(oaa, sub(c, oaa));							
						L out = plI.out(in, c, 0);
						Double objX = out.o.x+objH*out.dir.x/out.dir.z;
						Double objY = out.o.y+objH*out.dir.y/out.dir.z;
						Color clx = discreteObjColorX(objX.intValue());					
						Color cly = discreteObjColorY(objY.intValue());
						sumX[0] = sumX[0] + clx.getRed();
						sumX[1] = sumX[1] + clx.getGreen();
						sumX[2] = sumX[2] + clx.getBlue();
						sumY[0] = sumY[0] + cly.getRed();
						sumY[1] = sumY[1] + cly.getGreen();
						sumY[2] = sumY[2] + cly.getBlue();
					}
					pt.setColor(new Color(
							averageColor(sumX[0]/resolution, sumY[0]/resolution),
							averageColor(sumX[1]/resolution, sumY[1]/resolution),
							averageColor(sumX[2]/resolution, sumY[2]/resolution)
							));
					drr.drLine(i, h/6, i, h/2);
				}
	
//			}
		}
		

		
		// center line
//		pt.setColor(Color.black);		
//		drr.drLine(-w, pl.o.y, w, pl.o.y);
//		drr.drLine(pl.o.x, -h, pl.o.x, h);
	}

	public int averageColor(int cx, int cy){
//		int bgValue = objBgColor.getRed();
//		
//		cx = bgValue-cx;
//		cy = bgValue-cy;
//		
//		return  bgValue - (cx+cy - cx*cy/255);

		return  cx;// (cx+cy) - cx*cy/255; //(cx+cy)/2;// 
	}
	
	public double cellCountToM(double count){
		return count*cellLength;
	}
	
	public int mToCellCount(double m){
		return new Double(m/cellLength).intValue();
	}
	
	public int mToPixel(double m){
		return new Double(m/pixelRealSize).intValue();
	}
	
	public double pixelToM(double pixel){
		return pixel*pixelRealSize;
	}
	
	public static Color objBgColor = Color.white;
	public Color discreteObjColorX(int x){
		x = abs(x) + objl/2;
		int mx = x%objL;
				
		if(mx<objl){
			return objColors[x/objL%objColors.length];
		}else{
			return objBgColor;
		}
	}
	
	public Color discreteObjColorY(int y){
		y = abs(y) + objl/2;
		int my = y%objL;
				
		if(my<objl){
			return objColors[y/objL%objColors.length];
		}else{
			return objBgColor;
		}
	}
	
	public Color discreteObjColor(int x, int y){
		x = abs(x) + objl/2;
		y = abs(y) + objl/2;
		int mx = x%objL;
		int my = y%objL;
				
//		if(mx<objl && my<objl){
//			Color cx = objColors[x/objL%objColors.length];;
//			Color cy = objColors[y/objL%objColors.length];;
//			return new Color((cx.getRed()+cy.getRed())/2, (cx.getGreen()+cy.getGreen())/2, (cx.getBlue()+cy.getBlue())/2);
//		}else if(mx<objl){
//			return objColors[x/objL%objColors.length];
//		}else if(my<objl){
//			return objColors[y/objL%objColors.length];					
//		}else{
//			return objBgColor;
//		}
		if(mx<objl){
			return objColors[x/objL%objColors.length];
		}else{
			return objBgColor;
		}
	}
	
	private void testGroup() {
		Group grp = new Group(3, 10, 21, 10, 7, 3, 1);
//		Group grp = new Group(3, 4, 9, 4, 3, 2, 1);
//		Group grp = new Group(3, 3, 3, 1, 3, 3, 1);
		grp.initImage(w, h, 100, 100, Color.red, Color.red);
		double ze = grp.hl*(3*300)+grp.H;
		V e = new V(grp.d*grp.r*10, 0, ze);
//		V e = new V(0, 0, ze);
		pt.setClip(grp.getArea(w/2, h/2, 1500, 3, e));
//		pt.setClip(null);
		pt.drawImage(grp.img, 0, 0, null);
		
		double d = e.z/(e.z-grp.obj.z);
		double tr = grp.objR * d;
		V tCenter = add(e, mult(sub(grp.obj, e), d));
		pt.setClip(null);
		new Drawer(w, h, pt).drRect(tCenter.x, tCenter.y, tr*2, tr*2);
	}

	public void testClip(){
		int r = 3;
		int m = 7; // x dir
		int d = 1; // y dir
		int l = 10;
		int H = l*(d*m+500);
		Area area = getArea(r, m, l, H);
		area.intersect(getArea(r, m, l*(d*m+1), H));
//		area.add(getArea(r, m, l*d, H));
		
		pt.setClip(area);
//		pt.setColor(cl1);
		pt.setColor(Color.lightGray);
		pt.fillRect(0,0,w,h);
		System.out.println("fin");
		
	}
	
	public Area getArea(int r, int m, int l, int H){
		int cx = w/2;
		int cy = h/2;

		int lx = w/(r*m);
		int ly = h/(r*m);
		
		System.out.println(lx+", "+ly);

		
		float c = H/new Float(H-l);
		float rf = r*c;
		float rm = r*m*c;
		float demiRf = (r/2)*c;
		
//		Area area = new Area();
//		for(int i=-lx; i<=lx; i++){
//			for(int j=-ly; j<=ly; j++){
//				area.add(new Area(new Rectangle2D.Float(cx+i*rm-demiRf, cy+j*rm-demiRf, rf, rf)));
//			}
//		}

		Area area = new Area();
		for(int i=-lx; i<=lx; i++){
			area.add(new Area(new Rectangle2D.Float(cx+i*rm-demiRf, 0, rf, h)));
		}

		Area area2 = new Area();
		for(int j=-ly; j<=ly; j++){
			area2.add(new Area(new Rectangle2D.Float(0, cy+j*rm-demiRf, w, rf)));
		}
		
		area.intersect(area2);
		return area;
	}
	
	static PlGenerator generator;
	static boolean threadOpen = false;
	
	public void testPoly(){
		pt.setColor(cl4);
		
		if(generator==null){
			generator = new PlGenerator(15, 15, 20, 100, 20, 60);
		}
		
		if(threadOpen){
			List<LSeg> pl = generator.getNextSymPl();
			if(pl==null){
				generator = null;
			}else{
				drPL(pl);
			}
		}else{
			int i=0;
			List<LSeg> pl = null;
			for(; i<RDM.R.nextInt(222401); i++){
				pl = generator.getNextSymPl();
				if(pl==null){
					System.out.println(i);
					break;
				}
//				drPL(pl);
			}
			System.out.println(i);
			pt.setColor(Color.pink);
			drPL(pl);
			for(LSeg s : pl){
				for(int a=15; a<=165; a+=30){
					L l = new L(s.c, new V(C(a), S(a), 0));
					if(l.crossBefore(pl)){
						pt.setColor(cl4);
					}else{
						pt.setColor(cl3);
					}
					drLXY(l, -100);
				}
				
			}

			
		}
		

	}

	public void test7(){
//		pt.setColor(Color.GRAY);
//		pt.fillRect(0, 0, w, h);
		double r = 0.5;
		double l = r*2;
		double h = r*1.9;
		
		double n=1.5;
		double am = 90;
		
		int dv = 800;

		double dY = -300;
		
		double[][] xyd = new double[dv+1][3];
		
		xyd[0][0] = 0;
		xyd[0][1] = 0;
		xyd[0][2] = 0;
		
		double dx = r/dv;
		for(int i=1; i<=dv; i++){
			double x = xyd[i-1][0]+dx;
			double y = xyd[i-1][1]+dx*xyd[i-1][2];
			
			double a1 = AT(x*(1+l/r)/(y+h));
			double a2 = am*x/r;
			
			xyd[i][0] = x;
			xyd[i][1] = y;
			xyd[i][2] = -(n*S(a1)-S(a2))/(n*C(a1)-C(a2));
		}
		
	
		float trans = 0.01f;
//		double xo = 4*l/4;
		Double hue = 0d;
//		for(double xo=l; xo>=-l; xo -= l/50)
		int lScale = 400;
		int rScale = 200;
		double dl = l/lScale;
		int stepl = 0;
//		for(double xo=-l; xo<=l; xo += dl)
		int a=-75;for(double xo=-l*a/50; xo<=-l*(a-1)/50; xo += dl)
//		for(double xo=-2*l; xo<=-1*l; xo += l/100) // out of range
		{
			V po = new V(xo, -h+dY, 0);
			
//			Double hue = 0.7d;
			Color color = Color.getHSBColor(hue.floatValue(), 1, 1);
			color = new Color(((float)color.getRed())/255, ((float)color.getGreen())/255, ((float)color.getBlue())/255, trans);
			pt.setColor(color);
			stepl++;
			if(stepl%(lScale/50)==0){
				hue+=0.333333333333333;
			}
			
//			int i = new Double(-xo/l*dv).intValue();
			for(int i=0; i<=dv; i+=dv/rScale)
			{
				for(int dr = -1; dr<=1; dr+=2){
					double x = xyd[i][0]*dr;
					double y = xyd[i][1]+dY;
					double dxy = xyd[i][2]*dr;
					V p1 = new V(x, y, 0);

					L lf = new L(p1, new V(-dxy, 1, 0));				
					
					V dir = sub(lf.o, po);
					L li = new L(po, dir);

					L lo = null;
					V ltV = trm(lf.dir, dir, 1/n);
					if(ltV != null){
						lo = new L(lf.o, ltV);					
					}

					
					drLineXY(po, p1);
					if(lo != null && lo.dir.y>=0){
						drLXY(lo, 900);
//						drLXY(lo, -900);
					}
					
				}
			}			
		}
		
		
		pt.setColor(new Color(0.2f,0.2f,0.2f,0.1f));
		drLine(-w, -h+dY, w, -h+dY);
		drLine(0, -900, 0, 900);
		for(int i=0; i<dv; i++){
			drLine(xyd[i][0], xyd[i][1]+dY, xyd[i+1][0], xyd[i+1][1]+dY);
			drLine(-xyd[i][0], xyd[i][1]+dY, -xyd[i+1][0], xyd[i+1][1]+dY);
		}
		
	}
	
	public void test7V2(){
		double r = 100;
//		SegCurve curve = new SegCurve(r, r*2, r*1.9, 1.5);
		curve = new SegCurve(r, r*2.4, r*2.7, 1.5);
		curve.setRange(90, 50);
		curve.setCenter(new V(0, -200, 0));
		curve.setHue(10, 50);
		curve.initPArray();
	
//		float trans = 0.002f;
////		double xo = 4*l/4;
////		for(double xo=l; xo>=-l; xo -= l/50)
//		int lScale = 400;
//		int rScale = 200;
//		double dl = curve.l/lScale;
//		for(double xo=-curve.l; xo<=curve.l; xo += dl)
////		int a=-75;for(double xo=-curve.l*a/50; xo<=-curve.l*(a-1)/50; xo += dl)
////		for(double xo=-2*curve.l; xo<=-1*curve.l; xo += curve.l/100) // out of range
//		{
//			V po = add(new V(xo, -curve.h, 0), curve.center);
//			
////			Double hue = 0.7d;
////			Color color = Color.getHSBColor(hue.floatValue(), 1, 1);
////			color = new Color(((float)color.getRed())/255, ((float)color.getGreen())/255, ((float)color.getBlue())/255, trans);
//			
//			Color color = curve.getColor(xo-curve.center.x, trans);
//			pt.setColor(color);
//			
////			int i = new Double(-xo/curve.l*curve.dv).intValue();
//			for(int i=0; i<=curve.dv; i+=curve.dv/rScale)
//			{
//				for(int dr = -1; dr<=1; dr+=2){
//					L seg = dr==-1?curve.getNSeg(abs(i)) : curve.getPSeg(abs(i));
//					
//					L lf = new L(seg.o, new V(-seg.dir.y, seg.dir.x, seg.dir.z));				
//					
//					V dir = sub(lf.o, po);
//					L li = new L(po, dir);
//
//					L lo = null;
//					V ltV = trm(lf.dir, dir, 1/curve.n);
//					if(ltV != null){
//						lo = new L(lf.o, ltV);					
//					}
//
//					
//					drLineXY(po, seg.o);
//					if(lo != null && lo.dir.y>=0){
//						drLXY(lo, 900);
////						drLXY(lo, -900);
//					}
//					
//				}
//			}			
//		}
		
		
		pt.setColor(new Color(0.2f,0.2f,0.2f,0.1f));
		drLine(-w, -curve.h+curve.center.y, w, -curve.h+curve.center.y);
		drLine(0, -900, 0, 900);
		for(int i=0; i<curve.dv; i++){
			drLineXY(curve.getPSeg(i).o, curve.getPSeg(i+1).o);
			drLineXY(curve.getNSeg(i).o, curve.getNSeg(i+1).o);
		}
		
	}
	
	public void test7V3(){
		double r = 1;
		if(curve==null){
			// prepare for  avegage :  find H by function
			double hFind = SegCurve.findAverageH(100, 1.5);
//			System.out.println("hFind="+hFind);
			
//			curve = new SegCurve(r, r*2, r*2*hFind, 1.5); // average
//			curve = new SegCurve(r, r*1.9, r*2*hFind*1.3, 1.5); // stat with h find
			curve = new SegCurve(r, r*2.4, r*2.7, 1.5); // stat

			
//			curve = new SegCurve(r, r*1.5, r*5, 1.5); // stat
			
			
			
			curve.setRange(50, 200);
			SegCurve.modeMiddle=true;
			
			
//			curve = new SegCurve(r, r*2.4, r*2.7, 1.5);
////		curve = new SegCurve(r, r*2, r*2, 1.5);
//			curve.setRange(90, 50);
			
//			curve = new SegCurve(r, r*2.4, r*2.7, 1/1.5);
//			curve.setRange(40, 200);
			
			curve.setCenter(new V(0, -200, 0));
			curve.setHue(2, 100);
			
			int dH = 0;
//			curve.setHue(7, 100, 0+dH, 10+dH, 20+dH, 30+dH, 40+dH, 50+dH, 60+dH, 70+dH, 80+dH, 90+dH, 100+dH, 110+dH, 120+dH, 130+dH, 140+dH, 150+dH);
			
//			curve.initPArray_simple();
			curve.initPArray();
//			curve.initPArray_average();
		}
		
		Color gray = new Color(0.2f,0.2f,0.2f,0.05f);
		pt.setColor(gray);

		drLXY(new L(curve.center, new V(1,1,0)), 800);
		drLXY(new L(curve.center, new V(-1,1,0)), 800);
		
//		double x=1; 
//		double y=curve.center.y+200;
		for(double x=-w/2; x<w/2; x+=1)
		{
			for(double y = curve.center.y; y<h/2; y+=1)
//			for(double y = curve.center.y+300; y<curve.center.y+310; y+=1)
			{
				V pXY = new V(x, y, 0);			
				
				// r condition
				double absR = sub(pXY, curve.center).abs();
//				if(absR<=h*0.1 || absR>=h*0.5+200){
//				if(absR<=h*0.2-5 || absR>=h*0.2+5){
				if(absR<=h*0.7 || absR>=h*0.7+20){
					continue;
				}
				
				V colorSumRGB = new V(0, 0, 0);
				double areaSum = 0d;
				
				for(int i=0; i<curve.dv; i+=1){
					// part of
//					if(i>curve.dv*0.7)continue;
					
					for(int dr = -1; dr<=1; dr+=2){
						L seg = dr==-1?curve.getNSeg(i) : curve.getPSeg(i);
						L pst = new L(dr==-1?curve.getPstN(i) : curve.getPstP(i), seg.dir);
//						L pst = new L(add(seg.o, mult(seg.dir, dr*0.5*curve.r/curve.dv/seg.dir.x)), seg.dir);
						V dir = sub(pst.o, pXY);
						
						if(curve.n<1){
							dir = trm(new V(0,-1,0), dir, 1/curve.n);
						}
						
//						Double tgIn = dir.tgXY();
//						Double tgSeg = seg.dir.tgXY();
//						if(tgIn!=null && tgSeg!=null && tgIn*tgSeg>0 && abs(tgSeg)>abs(tgIn)){
//							continue;
//						}

						V f = new V(-pst.dir.y, pst.dir.x, pst.dir.z);
						if(f.y<0){
							f = new V(-f.x, -f.y, f.z );
						}

						if(multD(dir, f)>=0){
							continue;
						}
					
						
						L lo = new L(pst.o, trm(f, dir, curve.n));
						if(lo!=null && lo.dir != null && lo.dir.y<0){
							// get cross point:
							double distY = curve.center.y-curve.h - lo.o.y;
							double xImg = lo.o.x + distY/lo.dir.y * lo.dir.x;
							V pImg = new V(xImg, curve.center.y-curve.h, 0);
							Color c = curve.getColor(xImg-curve.center.x, 1f);
							
							V pDest = add(lo.o, mult(lo.dir, (-curve.h-pst.o.y)/lo.dir.y));

//							double ASeg = tgSeg==null ? 90 : AT(tgSeg);
//							double AIn = tgIn==null ? 90 : AT(tgIn);
//							double l = abs(curve.r/curve.dv / C(ASeg) * S(ASeg-AIn));
							
							double l = abs(V.multDU(f, dir));

							colorSumRGB.x += c.getRed()*l;
							colorSumRGB.y += c.getGreen()*l;
							colorSumRGB.z += c.getBlue()*l;
							areaSum += l;
							
//							pt.setColor(gray);
//							drLXY(lo, 900);
//							drLineXY(pXY, seg.o);
							
//							pt.setColor(c);
//							drPoint(pImg.x, pImg.y);
//							drLXY(new L(seg.o, new V(-dir.y, dir.x, 0)), l);
						}
					}
				}
				
				double areaTotal = areaSum;
//				double areaTotal = 2*r;
				
				Color cXY = new Color(
						new Double(colorSumRGB.x/areaTotal).intValue(),
						new Double(colorSumRGB.y/areaTotal).intValue(),
						new Double(colorSumRGB.z/areaTotal).intValue()
						);

				pt.setColor(cXY);
				drPoint(x, y);
				
			}
		}
		
		
		pt.setColor(gray);
		drLine(-w, -curve.h+curve.center.y, w, -curve.h+curve.center.y);
		drLine(0, -900, 0, 900);
		for(int i=0; i<curve.dv; i++){
			drLineXY(add(mult(sub(curve.getPSeg(i).o, curve.center), 100/r), curve.center),
					 add(mult(sub(curve.getPSeg(i+1).o, curve.center), 100/r), curve.center));
			drLineXY(add(mult(sub(curve.getNSeg(i).o, curve.center), 100/r), curve.center),
					 add(mult(sub(curve.getNSeg(i+1).o, curve.center), 100/r), curve.center));
		}
		
	}
	
	public void test7V3D(){
		double r = 1;
		boolean limitZ = false;
		
		if(curve==null){
//			curve = new SegCurve(r, r*2, r*2, 1.5);
//			curve.setRange(90, 200);

			curve = new SegCurve(r, r*2.4, r*2.7, 1.5);
			curve.setRange(50, 50);
//			curve = new SegCurve(r, r*2.4, r*2.7, 1/1.5);
//			curve.setRange(40, 50);
			curve.setCenter(new V(0, 0, 0)); // do not change
//			curve.setHue(7, 100);
			
			int dH = 0;
			curve.setHue(7, 100, 0+dH, 10+dH, 20+dH, 30+dH, 40+dH, 50+dH, 60+dH, 70+dH, 80+dH, 90+dH, 100+dH);
//			curve.setHue(7, 100, 
//					0+dH, 2+dH, 4+dH, 6+dH, 8+dH,
//					10+dH, 12+dH, 14+dH, 16+dH, 18+dH,
//					20+dH, 22+dH, 24+dH, 26+dH, 28+dH,
//					30+dH, 32+dH, 34+dH, 36+dH, 38+dH,
//					40+dH, 42+dH, 44+dH, 46+dH, 48+dH,
//					50+dH, 52+dH, 54+dH, 56+dH, 58+dH,
//					60+dH, 62+dH, 64+dH, 66+dH, 68+dH,
//					70+dH, 72+dH, 74+dH, 76+dH, 78+dH,
//					80+dH, 82+dH, 84+dH, 86+dH, 88+dH,
//					90+dH, 92+dH, 94+dH, 96+dH, 98+dH,
//					100
//			);
			
//			int step = 4;
//			Integer[] hues = new Integer[100/step+1];
//			for(int i=0; i<hues.length; i++){
//				hues[i] = step * i;
//			}
//			curve.setHue(7, 100, hues);
			

//			curve.initPArray_simple();
			curve.initPArray();
		}
		
		Color gray = new Color(0.2f,0.2f,0.2f,0.05f);
		pt.setColor(gray);

//		double x=1; 
//		double y=curve.center.y+200;

		double rB = h/5; // r to show
		for(double xD=-rB; xD<=rB; xD+=1)
		{
			for(double zD = -rB; zD<=rB; zD+=1){
//				if(abs(xD)>7 && abs(zD)>7)continue;
//				if(abs(zD)>0)continue;
//				if(abs(xD)>0 || abs(zD)>0)continue;

				if(xD*xD+zD*zD>rB*rB)continue;
				
				V pXY = U.platToSpherePosition(1, rB, xD, 0, zD);			
								
				V colorSumRGB = new V(0, 0, 0);
				double areaSum = 0d;
				
				// i, j, k -> x, y, z
				for(double i=-curve.dv; i<=curve.dv; i+=4){
					for(double k=-curve.dv; k<=curve.dv; k+=4){
						
						Double rXY = new Double(pow(i*i+k*k, 0.5));
						if(rXY>curve.dv)continue;

						int iXY = rXY.intValue();
						L seg = curve.getPSeg(iXY);
						if(seg==null)continue;
						L pst = new L(curve.getPstP(iXY), seg.dir);

						
						V rollTo = new V(i, 0, k).unit();

						V f = new V(-pst.dir.y, pst.dir.x, pst.dir.z);
						if(f.y<0){
							f = new V(-f.x, -f.y, f.z );
						}
						
						f = new V(f.x*rollTo.x, f.y, f.x*rollTo.z);
						
						pst =new L(new V(pst.o.x*rollTo.x, pst.o.y, pst.o.x*rollTo.z), new V(pst.dir.x*rollTo.x, pst.o.y, pst.dir.x*rollTo.z));

						V dir = sub(pst.o, pXY);

						if(curve.n<1){
							dir = trm(new V(0,-1,0), dir, 1/curve.n);
						}

						
						if(multD(dir, f)>=0){
							continue;
						}
						
						L lo = new L(pst.o, trm(f, dir, curve.n));
						
						if(lo.dir!=null && lo.dir.y<0){
							V pDest = add(lo.o, mult(lo.dir, (-curve.h-pst.o.y)/lo.dir.y));
							// get cross point:
							double distY = curve.center.y-curve.h - lo.o.y;
							double xImg = lo.o.x + distY/lo.dir.y * lo.dir.x;
							V pImg = new V(xImg, curve.center.y-curve.h, 0);

							Color c = Color.white;
							if(limitZ){
								// plan 1
//								if(pDest.absXZ()== 0){
//									c = curve.getColor(pDest.absXZ(), 1f);									
//								}else{
//									Double arcL = pDest.absXZ() *100 * U.ac(pDest.x/pDest.absXZ())/curve.l;
//									if(arcL.intValue()%4==0){
//										c = curve.getColor(pDest.absXZ(), 1f);									
//									}
//								}
								
//								---------------------------------------
								
								// plan 2
//								if(new Double(pDest.x*100/curve.l+100).intValue()%4 == 0 || new Double(pDest.z*100/curve.l+100).intValue()%4 == 0){
////									c = curve.getColor(pDest.absXZ(), 1f);
//									c=Color.red;
//								}
								
//								---------------------------------------

								// plan 3
								if(new Double(100+pDest.x/curve.l*100).intValue()%4==0 ){
									c = Color.red;
								}else if(new Double(100+pDest.z/curve.l*100).intValue()%4==0 ){
									c = Color.blue;
								}
								
							}else{
								c = curve.getColor(pDest.absXZ(), 1f);
							}

							double l = abs(V.multDU(f, dir));
							colorSumRGB.x += c.getRed()*l;
							colorSumRGB.y += c.getGreen()*l;
							colorSumRGB.z += c.getBlue()*l;
							areaSum += l;							
							
							
//							pt.setColor(gray);
//							drLXY(lo, 900);
//							drLineXY(pXY, seg.o);
							
//							pt.setColor(c);
//							drPoint(pImg.x, pImg.y);
//							drLXY(new L(seg.o, new V(-dir.y, dir.x, 0)), l);
						}
					}
				}
				
				double areaTotal = areaSum;
//				double areaTotal = 2*r;
				
				Color cXY = new Color(
						new Double(colorSumRGB.x/areaTotal).intValue(),
						new Double(colorSumRGB.y/areaTotal).intValue(),
						new Double(colorSumRGB.z/areaTotal).intValue()
						);

				pt.setColor(cXY);
				drPoint(xD, zD);
//				drLine(xD, zD-100, xD, zD+100);

			}
		}
		
		pt.setColor(new Color(0.5f, 0.5f, 0.5f, 0.2f));
		for(int i=1; i<=6; i++){
			drOval(0, 0, rB*i/6, rB*i/6);
		}
		
		
//		pt.setColor(gray);
//		drLine(-w, -curve.h+curve.center.y, w, -curve.h+curve.center.y);
//		drLine(0, -900, 0, 900);
//		for(int i=0; i<curve.dv; i++){
//			drLineXY(add(mult(sub(curve.getPSeg(i).o, curve.center), 100/r), curve.center),
//					 add(mult(sub(curve.getPSeg(i+1).o, curve.center), 100/r), curve.center));
//			drLineXY(add(mult(sub(curve.getNSeg(i).o, curve.center), 100/r), curve.center),
//					 add(mult(sub(curve.getNSeg(i+1).o, curve.center), 100/r), curve.center));
//		}
		
	}
	
	
	public void test7V3DInverse(){
		pt.setColor(Color.black);
		pt.drawString("3D Inverse", 10, 10);

		boolean limitZ = false;
		SegCurve curve = getTestCurve(false);
		SegSurface surface = new SegSurface(curve, new V(0, 0, 0), V.AXIS_Y);
		
		double viewSphereR = 500;
		
		double rangeViewP = 70;
		int lineDimr = 25;
		
		
		int cellNumR = 50;
		
		V[][] colors = new V[U.toI(viewSphereR*2)][U.toI(viewSphereR*2)];
		
		
		double colorsCountTotal = 0; 

		int cellr = 3;
		int cellR = cellr * 2+1;
		int cellBR = cellR * cellNumR;
		
		double scale = curve.l/cellBR;
		double scaleR = curve.r/lineDimr;
		V imgCenter = curve.getImgCenter();
		
		// test color for cells
		Color[] cellColors = {Color.red, Color.green, Color.blue, Color.pink, Color.CYAN, Color.gray, Color.CYAN, Color.orange, Color.white};
//		float hue = 0.005f;
//		for(int i=0; i<cellColors.length; i++){
//			cellColors[i] = new Color(cellColors[i].getRed()/255f, cellColors[i].getGreen()/255f, cellColors[i].getBlue()/255f, hue);
//		}
		
		

		int a = 15;
		int b = 15;
		int jump = 0;
		int dab = 1;
		
		
		int countLine = 0;
		
		for(int i=0; i<=a; i+=dab){
			System.out.println("i = " + i);
			for(int j=0; j<=b; j+=dab){
				
//		for(int i=0; i<=cellNumR; i+=dab){
//			System.out.println("i = " + i);
//			for(int j=0; j<=cellNumR; j+=dab){
				// a cell
				V centerIJ = new V(i*cellR, 0, j*cellR);
				if(centerIJ.abs()>cellBR || i<j || i%(jump+1)!=0 || j%(jump+1)!=0){
					continue;
				}
				
				Color c = (cellColors[RDM.nextI(0, cellColors.length)]); // random
//				Color c = curve.getColor(centerIJ.abs()*scale, 1); // hue of curve
				
				int border = 0;
				for(int ic = -cellr+border; ic<=cellr-border; ic++){
					for(int jc = -cellr+border; jc<=cellr-border; jc++){
						// a pixel
						V pImg = add(mult(add(centerIJ, new V(ic, 0, jc)), scale), imgCenter);

//////////////////////////////////////////////////////////////////////////////////////////////
						// use lineDimr
//						for(double cvi=-lineDimr; cvi<=lineDimr; cvi+=1){
//							for(double cvj=-lineDimr; cvj<=lineDimr; cvj+=1){
//								// get f and seg
//								Double rXY = new Double(pow(cvi*cvi+cvj*cvj, 0.5));
//								if(rXY>lineDimr)continue;
//
//								countLine++;
//								
//								// get in dir
//								V in = sub(add(mult(new V(cvi, 0, cvj), scaleR*1), curve.getEndCenter()), pImg).unit();
//								L outLine = surface.out(new L(pImg, in), 1/curve.n);
//								V out = outLine==null? null : outLine.dir;		
//////////////////////////////////////////////////////////////////////////////////////////////
								
								
								
//////////////////////////////////////////////////////////////////////////////////////////////
						// use curve dv
						for(double cvi=-curve.dv; cvi<=curve.dv; cvi+=8){
							for(double cvj=-curve.dv; cvj<=curve.dv; cvj+=8){
								// get f and seg
								Double rXY = new Double(pow(cvi*cvi+cvj*cvj, 0.5));
								if(rXY>curve.dv)continue;
								countLine++;

								int iXY = rXY.intValue();
								L seg = curve.getPSeg(iXY);
								if(seg==null)continue;
								L pst = new L(curve.getPstP(iXY), seg.dir);

								V f = new V(-pst.dir.y, pst.dir.x, pst.dir.z);
								if(f.y<0){
									f = new V(-f.x, -f.y, f.z );
								}
								
								if(cvj!=0 || cvi!=0){
									V rollTo = new V(cvi, 0, cvj).unit();
									f = new V(f.x*rollTo.x, f.y, f.x*rollTo.z);
									pst =new L(new V(pst.o.x*rollTo.x, pst.o.y, pst.o.x*rollTo.z), new V(pst.dir.x*rollTo.x, pst.o.y, pst.dir.x*rollTo.z));
								}
								
								
								// get in dir
								V in = sub(pst.o, pImg);

//								old condition
								// of no cross // to change by using directely points in plat y0
								V end = curve.getEndPoint();
								double endR = sub(end, curve.center).absXZ();
								
								double yPp = (end.y-pImg.y)/(pst.o.y-pImg.y);
								V crossY0 = add(pImg, mult(in, yPp));
								if(sub(crossY0, curve.center).absXZ()>endR){
									continue;
								}
								//drPointXZ(crossY0);
								
								
								
								//////////////////////////////////////////////////////////////////////
								// get out line
								V out = trm(f, in, 1/curve.n);
								infoString = "old";
								//////////////////////////////////////////////////////////////////////


								//////////////////////////////////////////////////////////////////////
//								L outLine = surface.out(new L(pImg, in), 1/curve.n);
//								V out = outLine==null? null : outLine.dir;		
//								infoString = "mixed";
								//////////////////////////////////////////////////////////////////////

								
								
//////////////////////////////////////////////////////////////////////////////////////////////
								
								

								
								
								if(out==null || out.y<0){
									continue;
								}

//								if(U.is0(cvi)){
//									drLineXY(pImg, seg.o);
//									drLXY(new L(seg.o, out), 200);
//									drPointXY(seg.o);
//								}
								
								
								// get p on shpere
//								V pShp = add(curve.center, mult(out, viewSphereR)); 
//								V pXZ = U.spherePositionToPlat(1, viewSphereR, pShp);			

								// get h f dim=viewSphereR
								double h = viewSphereR/T(rangeViewP);
								V pXZ = add(curve.center, mult(out, h/out.y));

								
								int xInt = max(0, min(U.toI(viewSphereR*2)-1,new Double(pXZ.x+U.toI(viewSphereR)).intValue()));
								int zInt = max(0, min(U.toI(viewSphereR*2)-1,new Double(pXZ.z+U.toI(viewSphereR)).intValue()));
								if(colors[xInt][zInt] == null){
									colors[xInt][zInt] = new V(0,0,0);
								}
								colors[xInt][zInt] = add(colors[xInt][zInt], new V(c.getRed()/255d, c.getGreen()/255d, c.getBlue()/255d));
								colorsCountTotal++;
							}
						}
					}
				}
				
				
			}
		}
		
		
		
		double maxC=0;
		double minC=Double.MAX_VALUE;
		for(int i=0; i<U.toI(viewSphereR)*2; i++){
			for(int j=0; j<U.toI(viewSphereR)*2; j++){
				V cij = colors[i][j];
				if( cij!= null){
					maxC = max(maxC, max(cij.x, max(cij.y, cij.z)));
					minC = min(minC, min(cij.x, min(cij.y, cij.z)));
				}
			}
		}

		System.out.println("maxC = "+maxC);
		System.out.println("countLine = " + countLine+", countTotale = " + colorsCountTotal);
		for(int i=0; i<U.toI(viewSphereR)*2; i++){
			for(int j=0; j<U.toI(viewSphereR)*2; j++){
				V cij = colors[i][j];
				if( cij!= null){
					pt.setColor(new Color( U.toI(cij.x/maxC*255), U.toI(cij.y/maxC*255), U.toI(cij.z/maxC*255)));
					drPoint(i-U.toI(viewSphereR), j-U.toI(viewSphereR));
				}else{
					pt.setColor(Color.black);
					drPoint(i-U.toI(viewSphereR), j-U.toI(viewSphereR));
				}
			}
		}

		
		pt.setColor(Color.gray);
//		drCircle(curve.center.x, curve.center.y, curve.r);
		drPointXZ(curve.center);

		drX(0);
		drY(0);
		pt.setColor(new Color(0.5f, 0.5f, 0.5f, 0.2f));
		for(int i=1; i<=6; i++){
			drCircle(0, 0, viewSphereR*i/6);
		}
		

		
	}
		
	public static V[][] colors = null;
	public static Integer[][] DIPs = null;
	
	// for test
	public static Drawer crtDrawer;
	
	public SegCurve getInitedCurve(SegCurve curve, boolean forceCreation){
		String objName = curve.toString();
		boolean needCreation = forceCreation;
		if(needCreation==false){
			SegCurve curveInfile = (SegCurve)ObjSaver.readObj(SegCurve.class, objName);
			if(curveInfile==null){
				needCreation = true;				
			}else{
				curve=curveInfile;
			}
		}
		
		if(needCreation){
			curve.init();
			ObjSaver.saveObj(curve, curve.toString());
		}
		
		return curve;
	}
	
	public SegCurve getTestCurve(boolean forceCreation){
		String objName = "test Curve 200 50";
		boolean needCreation = forceCreation;
		if(needCreation==false){
			curve = (SegCurve)ObjSaver.readObjWhenNeed(SegCurve.class, curve, objName);
			needCreation = curve==null;
		}
		
		if(needCreation){
			double r = 1;
			
			double hFind = SegCurve.findAverageH(100, 1.5);
			System.out.println("hFind="+hFind);
			
//			curve = new SegCurve(r, r*2, r*2*hFind, 1.5); // average
//			curve = new SegCurve(r*0.9, r*2, r*2*hFind, 1.5); // stat with h find
			curve = new SegCurve(r, r*2.4, r*2.7, 1.5);


			curve.setRange(50, 200);
			infoString += "\n range="+curve.am;
			curve.setCenter(new V(0, 0, 0)); // do not change

			curve.setHue(7, 100);			
			int dH = 0;
//			curve.setHue(7, 100, 0+dH, 10+dH, 20+dH, 30+dH, 40+dH, 50+dH, 60+dH, 70+dH, 80+dH, 90+dH, 100+dH);
//			curve.setHue(7, 100, 
//					0+dH, 2+dH, 4+dH, 6+dH, 8+dH,
//					10+dH, 12+dH, 14+dH, 16+dH, 18+dH,
//					20+dH, 22+dH, 24+dH, 26+dH, 28+dH,
//					30+dH, 32+dH, 34+dH, 36+dH, 38+dH,
//					40+dH, 42+dH, 44+dH, 46+dH, 48+dH,
//					50+dH, 52+dH, 54+dH, 56+dH, 58+dH,
//					60+dH, 62+dH, 64+dH, 66+dH, 68+dH,
//					70+dH, 72+dH, 74+dH, 76+dH, 78+dH,
//					80+dH, 82+dH, 84+dH, 86+dH, 88+dH,
//					90+dH, 92+dH, 94+dH, 96+dH, 98+dH,
//					100
//			);
			
//			int step = 4;
//			Integer[] hues = new Integer[100/step+1];
//			for(int i=0; i<hues.length; i++){
//				hues[i] = step * i;
//			}
//			curve.setHue(7, 100, hues);
			

//			curve.initPArray_simple();
			curve.initPArray();
//			curve.initPArray_average();
			
			ObjSaver.saveObj(curve, objName);
			
		}
		
		return curve;
		
	}
	
	public void testSegSurfaceCross(){
		System.out.println("begin");
		
		crtDrawer = this;
		double r = 200;
		SegCurve curve = new SegCurve(r, r*2.4, r*2.7, 1.5);
		curve.setRange(90, 100);
		curve.setCenter(new V(0, 0, 0));
		curve = getInitedCurve(curve, false);

		
		
		Color[] segColors = {Color.black, Color.red};
		L segLeftP = curve.getNSeg(0);
		L segRightP = curve.getPSeg(0);
		for(int i=1; i<=curve.dv; i++){
			pt.setColor(segColors[i%2]);
			L segLeft = curve.getNSeg(i);
			L segRight = curve.getPSeg(i);
			
			drLineXY(segLeft.o, segLeftP.o);
			drLineXY(segRight.o, segRightP.o);
			
			segLeftP = segLeft;
			segRightP = segRight;
		}
		
		SegSurface surface = new SegSurface(curve, new V(0, 0, 0),  V.AXIS_Y);
		

		V EC = curve.getEndCenter();
		V EP = curve.getEndPoint();
		Double hLO = EC.y*2;
		pt.setColor(Color.lightGray);
		drX(0);
		drY(hLO);

		
		for(int c = 0; c<1; c++){
			
			segLeftP = null;
			segRightP = null;		
//			for(int i=0; i<1; i++){
//			for(int i=0; i<curve.dv; i++){
			for(int i=0; i<=curve.dv; i++){
				L segLeft = curve.getNSeg(i);
				L segRight = curve.getPSeg(i);
	
				// get random xy on seg
//				double dx = RDM.nextD(0, curve.r/curve.dv);

				// get int xy on seg
				double dx = 0;
				
				V vLeft = add(segLeft.o, mult(segLeft.dir, -dx/segLeft.dir.x));
				V vRight = add(segRight.o, mult(segRight.dir, dx/segRight.dir.x));
				
				// get random pO
				V poL = new V(RDM.nextD(-EP.x*2, EP.x*2), hLO, 0);
				V poR = new V(RDM.nextD(-EP.x*2, EP.x*2), hLO, 0);

				// roll random 
				vLeft = V.ROLL(V.AXIS_Y, vLeft, RDM.nextD(0, 360));
				vRight = V.ROLL(V.AXIS_Y, vRight, RDM.nextD(0, 360));
				poL = V.ROLL(V.AXIS_Y, poL, RDM.nextD(0, 360));
				poR = V.ROLL(V.AXIS_Y, poR, RDM.nextD(0, 360));
				

				if(c==0){
					pt.setColor(Color.green);
					drLineXY(vLeft, poL);
					drLineXY(vRight, poR);
					
				}
				
				L inL = new L(poL, sub(vLeft, poL));
				L inR = new L(poR, sub(vRight, poR));

//				inR = new L(new V(509.78010136766096,-254.45202682944546,-293.6146715834944), new V(-0.9482645240134459,0.3133431546521026,0.051092660238252124));
//				poR = inR.o;

				
				V crossL = surface.findCvCross(inL);
				V crossR = surface.findCvCross(inR);
//				V crossL = null;
				
//				System.out.println(crossL!=null && crossR!=null);
				
				pt.setColor(Color.red);
				if(crossL!=null){
					if(c==0){
						drLineXY(poL, crossL);						
					}
					double rInEndPL = sub(add(inL.o, mult(inL.dir, (EC.y-inL.o.y)/inL.dir.y)), EC).absXZ();
					if(U.aAbsSmallb(rInEndPL, EP.x) && !U.is0(V.dist(crossL, vLeft))){
						System.out.println(i+": not same p as vLeft\n" + crossL+"\n"+vLeft);
					}
				}else{
					double rInEndPL = sub(add(inL.o, mult(inL.dir, (EC.y-inL.o.y)/inL.dir.y)), EC).absXZ();
					if(U.aAbsSmallb(rInEndPL, EP.x)){
						pt.setColor(Color.blue);
						drLineXY(vLeft, poL);
						System.out.println(i+" l not outside line");
					}

					System.out.println(i+": not found for inL " + inL);
				}
				
				if(crossR!=null){
					if(c==0){
						drLineXY(poR, crossR);
					}
					
					double rInEndPL = sub(add(inR.o, mult(inR.dir, (EC.y-inR.o.y)/inR.dir.y)), EC).absXZ();
					if(U.aAbsSmallb(rInEndPL, EP.x) && !U.is0(V.dist(crossR, vRight))){
						System.out.println(i+": not same p as vRight " + crossR+", "+vRight);
					}
				}else{
					double rInEndPL = sub(add(inR.o, mult(inR.dir, (EC.y-inR.o.y)/inR.dir.y)), EC).absXZ();
					if(U.aAbsSmallb(rInEndPL, EP.x)){
						pt.setColor(Color.blue);
						drLineXY(vRight, poR);			
						System.out.println(i+" r not outside line");
					}
					System.out.println(i+": not found for inR " + inR);
				}
				
				
				segLeftP = segLeft;
				segRightP = segRight;
				
			}
		}
		
		
		
		
		
		
		
		
		
		
		
		System.out.println("end");
		
		
		
	}
	
	public void testSegSurfaceCrossZ(){
		System.out.println("begin");
		
		crtDrawer = this;
		double r = 200;
		SegCurve curve = new SegCurve(r, r*2.4, r*2.7, 1.5);
		curve.setRange(90, 100);
		curve.setCenter(new V(0, 0, 0));
		curve = getInitedCurve(curve, false);

		SegSurface surface = new SegSurface(curve, new V(0, 100, 0),  V.AXIS_Z);
		
//		V drawAxis = V.AXIS_X;
		V drawAxis = V.AXIS_Y;
//		V drawAxis = V.AXIS_Z;
		
		Color[] segColors = {Color.black, Color.red};
		
		L segLeftP = surface.toSurfaceAxisPos(curve.getNSeg(0));
		L segRightP = surface.toSurfaceAxisPos(curve.getPSeg(0));
		
		for(int i=1; i<=curve.dv; i++){
			pt.setColor(segColors[i%2]);
			L segLeft = surface.toSurfaceAxisPos(curve.getNSeg(i));
			L segRight = surface.toSurfaceAxisPos(curve.getPSeg(i));
			
			if(drawAxis==V.AXIS_X){
				drLineYZ(segLeft.o, segLeftP.o);
				drLineYZ(segRight.o, segRightP.o);
			}else if(drawAxis==V.AXIS_Y){
				drLineXZ(segLeft.o, segLeftP.o);
				drLineXZ(segRight.o, segRightP.o);
			}else if(drawAxis==V.AXIS_Z){
				drLineXY(segLeft.o, segLeftP.o);
				drLineXY(segRight.o, segRightP.o);
			}
			
			segLeftP = segLeft;
			segRightP = segRight;
		}
		
		V EC = surface.toSurfaceAxisPos(curve.getEndCenter());
		V EP = surface.toSurfaceAxisPos(curve.getEndPoint());

		Double hLO = EC.z*2;

		pt.setColor(Color.lightGray);
		drX(surface.center.x);
		drY(surface.center.y);
		drCircle(surface.center.x, surface.center.y, curve.r);
		drY(hLO);

		for(int c = 0; c<1; c++){
			
			segLeftP = null;
			segRightP = null;		
//			for(int i=0; i<3; i++){
//			for(int i=0; i<curve.dv; i++){
			for(int i=0; i<curve.dv; i+=10){
				L segLeft = surface.toSurfaceAxisPos(curve.getNSeg(i));
				L segRight = surface.toSurfaceAxisPos(curve.getPSeg(i));
	
				// get random xy on seg
				double dx = RDM.nextD(0, curve.r/curve.dv);

				// get int xy on seg
//				double dx = 0;
				
				V vLeft = add(segLeft.o, mult(segLeft.dir, -dx/segLeft.dir.x));
				V vRight = add(segRight.o, mult(segRight.dir, dx/segRight.dir.x));
				
				// get random pO
//				V poL = add(surface.center, new V(RDM.nextD(-EP.x*2, EP.x*2), 0, hLO));
//				V poR = add(surface.center, new V(RDM.nextD(-EP.x*2, EP.x*2), 0, hLO));

				// get random pO from top
				V poL = add(surface.center, new V(RDM.nextD(-EP.x*2, EP.x*2), 0, -hLO));
				V poR = add(surface.center, new V(RDM.nextD(-EP.x*2, EP.x*2), 0, -hLO));

				
				// roll random 
				vLeft = V.ROLL(surface.sfCenterAxis, vLeft, RDM.nextD(0, 360));
				vRight = V.ROLL(surface.sfCenterAxis, vRight, RDM.nextD(0, 360));
				poL = V.ROLL(surface.sfCenterAxis, poL, RDM.nextD(0, 360));
				poR = V.ROLL(surface.sfCenterAxis, poR, RDM.nextD(0, 360));
				

				if(c==0){
					pt.setColor(Color.green);
					if(drawAxis==V.AXIS_X){
						drLineYZ(vLeft, poL);
						drLineYZ(vRight, poR);
					}else if(drawAxis==V.AXIS_Y){
						drLineXZ(vLeft, poL);
						drLineXZ(vRight, poR);
					}else if(drawAxis==V.AXIS_Z){
						drLineXY(vLeft, poL);
						drLineXY(vRight, poR);
					}
					
				}
				
//				if(true)continue;
				L inL = new L(poL, sub(vLeft, poL));
				L inR = new L(poR, sub(vRight, poR));
				
//				inR = new L(new V(509.78010136766096,-254.45202682944546,-293.6146715834944), new V(-0.9482645240134459,0.3133431546521026,0.051092660238252124));
//				poR = inR.o;

				
				V crossL = surface.findCvCrossReal(inL);
				V crossR = surface.findCvCrossReal(inR);

				// from bottom
//				L outL = surface.out(inL, 1/curve.n);
//				L outR = surface.out(inR, 1/curve.n);
				
				// from top
				L outL = surface.out(inL, curve.n);
				L outR = surface.out(inR, curve.n);
				
				//				V crossL = null;
				
//				System.out.println(crossL!=null && crossR!=null);
				
				
				pt.setColor(Color.black);
				if(outL!=null && outL.dir!=null){
					if(drawAxis==V.AXIS_X){
						drLYZ(outL, 300);
					}else if(drawAxis==V.AXIS_Y){
						drLXZ(outL, 300);
					}else if(drawAxis==V.AXIS_Z){
						drLXY(outL, 300);
					}
				}
				if(outR!=null && outR.dir!=null){
					if(drawAxis==V.AXIS_X){
						drLYZ(outR, 300);
					}else if(drawAxis==V.AXIS_Y){
						drLXZ(outR, 300);
					}else if(drawAxis==V.AXIS_Z){
						drLXY(outR, 300);
					}
				}
				
				if(crossL!=null){
					if(c==0){
						pt.setColor(Color.red);
						if(drawAxis==V.AXIS_X){
							drLineYZ(poL, crossL);						
						}else if(drawAxis==V.AXIS_Y){
							drLineXZ(poL, crossL);						
						}else if(drawAxis==V.AXIS_Z){
							drLineXY(poL, crossL);						
						}
					}
					double rInEndPL = sub(add(inL.o, mult(inL.dir, (EC.z-inL.o.z)/inL.dir.z)), EC).absXZ();
					if(U.aAbsSmallb(rInEndPL, EP.x) && !U.is0(V.dist(crossL, vLeft))){
//						System.out.println(i+": not same p as vLeft\n" + crossL+"\n"+vLeft);
					}
				}else{
					double rInEndPL = sub(add(inL.o, mult(inL.dir, (EC.z-inL.o.z)/inL.dir.z)), EC).absXZ();
					if(U.aAbsSmallb(rInEndPL, EP.x)){
						pt.setColor(Color.blue);
						if(drawAxis==V.AXIS_X){
							drLineYZ(vLeft, poL);
						}else if(drawAxis==V.AXIS_Y){
							drLineXZ(vLeft, poL);
						}else if(drawAxis==V.AXIS_Z){
							drLineXY(vLeft, poL);
						}
						System.out.println(i+" l not outside line");
					}

					System.out.println(i+": not found for inL " + inL);
				}
				
				if(crossR!=null){
					if(c==0){
						pt.setColor(Color.red);
						if(drawAxis==V.AXIS_X){
							drLineYZ(poR, crossR);
						}else if(drawAxis==V.AXIS_Y){
							drLineXZ(poR, crossR);
						}else if(drawAxis==V.AXIS_Z){
							drLineXY(poR, crossR);
						}
					}
					
					double rInEndPL = sub(add(inR.o, mult(inR.dir, (EC.z-inR.o.z)/inR.dir.z)), EC).absXZ();
					if(U.aAbsSmallb(rInEndPL, EP.x) && !U.is0(V.dist(crossR, vRight))){
//						System.out.println(i+": not same p as vRight \n" + crossR+"\n"+vRight);
					}
				}else{
					double rInEndPL = sub(add(inR.o, mult(inR.dir, (EC.z-inR.o.z)/inR.dir.z)), EC).absXZ();
					if(U.aAbsSmallb(rInEndPL, EP.x)){
						pt.setColor(Color.blue);
						if(drawAxis==V.AXIS_X){
							drLineYZ(vRight, poR);			
						}else if(drawAxis==V.AXIS_Y){
							drLineXZ(vRight, poR);			
						}else if(drawAxis==V.AXIS_Z){
							drLineXY(vRight, poR);			
						}
						System.out.println(i+" r not outside line");
					}
					System.out.println(i+": not found for inR " + inR);
				}
				
				
				segLeftP = segLeft;
				segRightP = segRight;
				
			}
		}
		
		
		System.out.println("end");
		
		
		
	}
	
	public void testSegSurfaceCrossZInverse(){
		System.out.println("begin");
		
		crtDrawer = this;
		double r = 200;
		SegCurve curve = new SegCurve(r, r*2.4, r*2.7, 1.5);
		curve.setRange(90, 100);
		curve.setCenter(new V(0, 0, 0));
		curve = getInitedCurve(curve, false);

		SegSurface surface = new SegSurface(curve, new V(0, 100, 0),  V.AXIS_Z);
		
//		V drawAxis = V.AXIS_X;
		V drawAxis = V.AXIS_Y;
//		V drawAxis = V.AXIS_Z;
		
		Color[] segColors = {Color.black, Color.red};
		
		L segLeftP = surface.toSurfaceAxisPos(curve.getNSeg(0));
		L segRightP = surface.toSurfaceAxisPos(curve.getPSeg(0));
		
		for(int i=1; i<=curve.dv; i++){
			pt.setColor(segColors[i%2]);
			L segLeft = surface.toSurfaceAxisPos(curve.getNSeg(i));
			L segRight = surface.toSurfaceAxisPos(curve.getPSeg(i));
			
			if(drawAxis==V.AXIS_X){
				drLineYZ(segLeft.o, segLeftP.o);
				drLineYZ(segRight.o, segRightP.o);
			}else if(drawAxis==V.AXIS_Y){
				drLineXZ(segLeft.o, segLeftP.o);
				drLineXZ(segRight.o, segRightP.o);
			}else if(drawAxis==V.AXIS_Z){
				drLineXY(segLeft.o, segLeftP.o);
				drLineXY(segRight.o, segRightP.o);
			}
			
			segLeftP = segLeft;
			segRightP = segRight;
		}
		
		V EC = surface.toSurfaceAxisPos(curve.getEndCenter());
		V EP = surface.toSurfaceAxisPos(curve.getEndPoint());

		Double hLO = EC.z*2;

		pt.setColor(Color.lightGray);
		drX(surface.center.x);
		drY(surface.center.y);
		drCircle(surface.center.x, surface.center.y, curve.r);
		drY(hLO);

		for(int c = 0; c<1; c++){
			
			segLeftP = null;
			segRightP = null;		
			for(int i=0; i<10; i++){
				// get random pO
				V pTop = add(surface.center, new V(RDM.nextD(-EP.x*2, EP.x*2), 0, hLO));
				V pBottom = add(surface.center, new V(RDM.nextD(-EP.x*2, EP.x*2), 0, -hLO));

				// roll random 
				pTop = V.ROLL(surface.sfCenterAxis, pTop, RDM.nextD(0, 360));
				pBottom = V.ROLL(surface.sfCenterAxis, pBottom, RDM.nextD(0, 360));
				

				if(c==0){
					pt.setColor(Color.green);
					if(drawAxis==V.AXIS_X){
						drLineYZ(pTop, pBottom);
					}else if(drawAxis==V.AXIS_Y){
						drLineXZ(pTop, pBottom);
					}else if(drawAxis==V.AXIS_Z){
						drLineXY(pTop, pBottom);
					}
					
				}
				
//				if(true)continue;
				L inL = new L(pTop, sub(pBottom, pTop));
				
				
				V crossL = surface.findCvCrossReal(inL);

				
				L outL = surface.out(inL, curve.n);

				pt.setColor(Color.black);
				if(outL!=null && outL.dir!=null){
					if(drawAxis==V.AXIS_X){
						drLYZ(outL, 300);
					}else if(drawAxis==V.AXIS_Y){
						drLXZ(outL, 300);
					}else if(drawAxis==V.AXIS_Z){
						drLXY(outL, 300);
					}
				}

				if(crossL!=null){
					if(c==0){
						pt.setColor(Color.red);
						if(drawAxis==V.AXIS_X){
							drLineYZ(pTop, crossL);						
						}else if(drawAxis==V.AXIS_Y){
							drLineXZ(pTop, crossL);						
						}else if(drawAxis==V.AXIS_Z){
							drLineXY(pTop, crossL);						
						}
					}
				}
				
			}
		}
		
		
		System.out.println("end");
		
		
		
	}
	
	public void testSegSurfaceCrossZPP(){
		System.out.println("begin");
		
		
		crtDrawer = this;
		double r = 300;
		SegCurve curve = new SegCurve(r, r*2.4, r*2.7, 1.5);
		curve.setRange(90, 10);
		curve.setCenter(new V(0, 0, 0));
		curve = getInitedCurve(curve, false);

		SegSurface surface = new SegSurface(curve, new V(0, 0, 0),  V.AXIS_Z);
		
		Color[] segColors = {Color.black, Color.red};
		
		for(int i=1; i<=curve.dv; i++){
			pt.setColor(segColors[i%2]);
			L segRight = surface.toSurfaceAxisPos(curve.getPSeg(i));
			ppCicleZ(surface.sfCenterAxis, segRight.o);
		}

		V EC = surface.toSurfaceAxisPos(curve.getEndCenter());
		V EP = surface.toSurfaceAxisPos(curve.getEndPoint());

		Double hLO = EC.z*2;

		pt.setColor(Color.lightGray);
		drX(getVPP(surface.center).x);
		drY(getVPP(surface.center).z);
		drY(-hLO);
		drY(hLO);

		
		
		
		for(int c = 0; c<1; c++){
			
			V segLeftP = null;
			V segRightP = null;		
			for(int i=0; i<10; i++){
				// get random pO
				V pTop = add(surface.center, new V(RDM.nextD(-EP.x*2, EP.x*2), 0, -hLO));
				V pBottom = add(surface.center, new V(RDM.nextD(-EP.x*2, EP.x*2), 0, hLO));

				// roll random 
				pTop = V.ROLL(surface.sfCenterAxis, pTop, RDM.nextD(0, 360));
				pBottom = V.ROLL(surface.sfCenterAxis, pBottom, RDM.nextD(0, 360));
				

				if(c==0){
					pt.setColor(Color.green);
					ppLine(pTop, pBottom);
					
				}
				
//				if(true)continue;
				L inL = new L(pTop, sub(pBottom, pTop));
				
				
				V crossP = surface.findCvCrossReal(inL);

				
//				L outL = surface.out(inL, curve.n);
				
//				pt.setColor(Color.black);
//				if(outL!=null && outL.dir!=null){
//						ppLine(outL, 300);
//				}

				if(crossP!=null){
					if(c==0){
						pt.setColor(Color.black);
//						pt.setColor(Color.red);
						ppLine(pTop, crossP);						
					}
				}
				
			}
		}
		
		
		System.out.println("end");
		
		
		
	}
	
	public void testSegSurface(){
		infoString = "seg surface";
		infoColor = Color.WHITE;
		
		crtDrawer = this;
		
		pt.setColor(Color.lightGray);
		fiRect(0, 0, w, h); 

		SegCurve curve = getTestCurve(false);

//		pt.setColor(Color.gray);
//		for(int i=0; i<=curve.dv; i++){
//			drPointXY(curve.getPSeg(i).o);
//		}
		
		SegSurface surface = new SegSurface(curve, new V(0, 0, 0), V.AXIS_Y);
		
		double viewSphereR = 500;
		
		double rangeViewP = 70;
		
		int cellNumR = 50;

		int cellr = 3;
		int cellBorder = 0;

		int cellR = cellr * 2+1;
		int cellBR = cellR * cellNumR;
		
		V imgCenter = curve.getImgCenter();
		
		// test color for cells
		Color[] cellColors = {Color.red, Color.green, Color.blue, Color.pink, Color.CYAN, Color.gray, Color.CYAN, Color.orange, Color.white};
//		float hue = 0.005f;
//		for(int i=0; i<cellColors.length; i++){
//			cellColors[i] = new Color(cellColors[i].getRed()/255f, cellColors[i].getGreen()/255f, cellColors[i].getBlue()/255f, hue);
//		}
		
		
		int dab = 1;
		int jump = 0;
		int a = 20;
		int b = 20;

		int lineDimr = 10;
		
//		double scale = curve.l*1.5/cellBR; // 1.5 for cell voisin 
		double scale = curve.l/cellBR; // 1.5 for cell voisin 
		double scaleR = curve.r/lineDimr;

		int countLine = 0;
		int countValid = 0;
		
		boolean showColorByDensity = true;
		double colorPowerIndex = 1;

		DirCounts dirCounts = new DirCounts(U.toI(viewSphereR*1.2), cellNumR);
		
		
//		[cellNumR*2+1][cellNumR*2+1]		
		if(colors == null){
			colors = new V[U.toI(viewSphereR*2+1)][U.toI(viewSphereR*2+1)];

			
			for(int i=-a; i<=a; i+=dab){
				System.out.println("i = " + i);
				for(int j=-b; j<=b; j+=dab){
					
//			for(int i=0; i<=cellNumR; i+=dab){
//				System.out.println("i = " + i);
//				for(int j=0; j<=cellNumR; j+=dab){
					// a cell
					if(j%10==0){
						System.out.println("\tj = "+j);
					}
					
					V centerIJ = new V(i*cellR, 0, j*cellR);
					if(centerIJ.abs()>cellBR){
						continue;
					}

					// sec 0.5
					if(i<0 || j<0 || i<j){
						continue;
					}
//					
//					// sec 1
//					if(i<0 || j<0){
//						continue;
//					}
//
//					// jump
//					if(i%(jump+1)!=0 || j%(jump+1)!=0){
//						continue;
//					}
					
					
//					Color c = Color.green; // fix
					Color c = (cellColors[RDM.nextI(0, cellColors.length)]); // random
//					Color c = curve.getColor(centerIJ.abs()*scale, 1); // hue of curve
					
					int border = 0;
					int maxIJ = cellr-border;
					for(int ic = -maxIJ; ic<=maxIJ; ic++){
						for(int jc = -maxIJ; jc<=maxIJ; jc++){

							// center
//							if(ic != 0 || jc !=0)continue;
							
							// i=j
//							if(ic != jc)continue;
							
							// end
//							if(abs(ic)!=maxIJ || abs(jc)!=maxIJ)continue;
							
							// border
//							if(abs(ic)<maxIJ-1 && abs(jc)<maxIJ-1)continue;

							// a pixel
							V pImg = add(mult(add(centerIJ, new V(ic, 0, jc)), scale), imgCenter);
							
							for(double cvi=-lineDimr; cvi<=lineDimr; cvi+=1){
//								for(double cvj=0; cvj<=0; cvj+=1){
								for(double cvj=-lineDimr; cvj<=lineDimr; cvj+=1){
									
									// get f and seg
									Double rXY = new Double(pow(cvi*cvi+cvj*cvj, 0.5));
									if(rXY>lineDimr)continue;

									countLine++;
									
									// get in dir
									V in = sub(add(mult(new V(cvi, 0, cvj), scaleR*1), curve.getEndCenter()), pImg).unit();
									L outLine = surface.out(new L(pImg, in), 1/curve.n);
//									L outLine = surface.outFor0(new L(pImg, in), 1/curve.n);
									V out = outLine==null? null : outLine.dir;		

//									pt.setColor(Color.lightGray);
//									drLXY(new L(pImg, in), 100);
//									if(U.is0(cvj) && outLine!=null){
//										pt.setColor(Color.gray);
//										drLineXY(pImg, outLine.o);
//										if(out!=null){
//											drLXY(new L(outLine.o, out), 500);
//										}
//										pt.setColor(Color.red);
//										drPointXY(outLine.o);
//									}
//									pt.setColor(Color.green);
//									drPointXY(add(mult(new V(cvi, 0, cvj), scaleR), curve.getEndCenter()));

									if(out==null){
										continue;
									}


									
									// get p on shpere
//									V pShp = add(curve.center, mult(out, viewSphereR));
//									V pXZ = U.spherePositionToPlat(1, viewSphereR, pShp);			

									// get h f dim=viewSphereR
									double h = viewSphereR/T(rangeViewP);
									V pXZ = add(curve.center, mult(out, h/out.y));
									
									int xInt = max(0, min(U.toI(viewSphereR*2),U.toI45(pXZ.x)+U.toI(viewSphereR)));

									int zInt = max(0, min(U.toI(viewSphereR*2),U.toI45(pXZ.z)+U.toI(viewSphereR)));
									if(colors[xInt][zInt] == null){
										colors[xInt][zInt] = new V(0,0,0);
									}
									
//									dirCounts.addValue(U.toI(pXZ.x), U.toI(pXZ.z), i, j, 1);

									dirCounts.addValue(U.toI45(pXZ.x), U.toI45(pXZ.z), i, j, multD(in, V.AXIS_Y));
//									dirCounts.addValueOct(U.toI45(pXZ.x), U.toI45(pXZ.z), i, j, multD(in, V.AXIS_Y));
									
									countValid++;
									colors[xInt][zInt] = add(colors[xInt][zInt], new V(c.getRed()/255d, c.getGreen()/255d, c.getBlue()/255d));
//									colors[xInt][zInt] = add(colors[xInt][zInt], mult(new V(c.getRed()/255d, c.getGreen()/255d, c.getBlue()/255d), multD(in, V.AXIS_Y)));

								}
							}
						}
					}
					
					
				}
			}

		}


		
		
		
		
		
		double maxC=0;
		double minC=Double.MAX_VALUE;
		for(int i=0; i<=U.toI(viewSphereR)*2; i++){
			for(int j=0; j<=U.toI(viewSphereR)*2; j++){
				V cij = colors[i][j];
				if( cij!= null){
					maxC = max(maxC, max(cij.x, max(cij.y, cij.z)));
					minC = min(minC, min(cij.x, min(cij.y, cij.z)));
				}
			}
		}

		System.out.println("maxC = " + maxC);
		System.out.println("countLine = " + countLine+", countTotale = " + countValid);

		
		pt.setColor(Color.black);
		pt.fillRect(0, 0, w, h);

	

		
		// show area shp
//		if(DIPs==null){
//			DIPs = new Integer[dirCounts.dirR*2+1][dirCounts.dirR*2+1];
//			for(int i= -dirCounts.dirR; i<= dirCounts.dirR; i++){
//				for(int j= -dirCounts.dirR; j<= dirCounts.dirR; j++){
//					Integer maxV = null;
//					double maxCt = 1;
//					double totalCount = 0;
//					double advance = 1;
//					int[][] dirOcts = Oct.getOctAll(i, j);
//					for(int octI=0; octI<8; octI++){
//						IPCounts ipCounts = dirCounts.getValue(dirOcts[octI][0], dirOcts[octI][1]);
//						
//						if(ipCounts != null && ipCounts.ipMap != null){
//							
//							for(Map.Entry<Integer, Double> entry : ipCounts.ipMap.entrySet()){
//								Integer ip = entry.getKey();
//								int[] ipOct = Oct.toIDis(ipCounts.ijk.getI(ip), ipCounts.ijk.getJ(ip),-octI);
//								if(ipOct != null){
//									Double d = entry.getValue();
//									totalCount += d;
//									if(d>maxCt){
//										advance = d/maxCt;
//										maxCt = d;
//										maxV = ipCounts.ijk.getK(ipOct[0], ipOct[1]);										
//									}else{
//										advance = min(advance, maxCt/d);
//									}
//								}
//							}
//
//						}
//					}
//
////					maxCt*3>=totalCount ||
//					if(maxV != null && ( advance>1.2 )){
//						DIPs[i+dirCounts.dirR][j+dirCounts.dirR] = maxV;				
//					}else{
//						DIPs[i+dirCounts.dirR][j+dirCounts.dirR] = null;				
//					}
//				}
//			}
//		}
//		
//		Map<Integer, Color> clForArea = new HashMap<Integer, Color>();
//		for(int i= -dirCounts.dirR; i<= dirCounts.dirR; i++){
//			for(int j= -dirCounts.dirR; j<= dirCounts.dirR; j++){
//				Integer maxV = DIPs[i+dirCounts.dirR][j+dirCounts.dirR];
//				if(maxV != null){
//					Color c = clForArea.get(maxV);
//					if(c==null){
//						int mij = (300000+dirCounts.ijk.getI(maxV)%2*2+dirCounts.ijk.getJ(maxV))%4;
//						c =  mij==0?Color.red : mij==1? Color.green : mij==2? Color.white : Color.black; // fix
////						c = (cellColors[RDM.nextI(0, cellColors.length)]); // random
//						clForArea.put(maxV, c);
//					}
//					
//					pt.setColor(c);
//					drPoint(i, j);
//				}	
//					
//			}
//		}
		
		
		// show light sph
		for(int i=0; i<=U.toI(viewSphereR)*2; i++){
			for(int j=0; j<=U.toI(viewSphereR)*2; j++){
				V cij = colors[i][j];
				if( cij!= null){
					pt.setColor(new Color( U.toI(cij.x/maxC*255), U.toI(cij.y/maxC*255), U.toI(cij.z/maxC*255)));
					if(showColorByDensity){
						pt.setColor(vToColor(cij, maxC, colorPowerIndex));
					}else{
						pt.setColor(vToColor2(cij, maxC));
					}
					drPoint(i-U.toI(viewSphereR), j-U.toI(viewSphereR));
				}else{
				}
			}
		}
		

		// positions
		pt.setColor(new Color(1f, 0.3f, 0.3f, 0.8f));
//		drCircle(curve.center.x, curve.center.y, curve.r);
		drPointXZ(curve.center);

		drX(0);
		drY(0);
		pt.setColor(new Color(0.8f, 0.2f, 0.2f, 0.8f));
		for(int i=1; i<=rangeViewP/10; i++){
			drCircle(0, 0, viewSphereR/T(rangeViewP)*T(i*10));
		}
		
		
		
		
		
		
		// count analyse
//
//		int maxDirCount = 0;
//		
//		int[] counts = new int[10];
//		int countN = 0;
//		int count0 = 0;
//		for(int i= -dirCounts.dirR; i<= dirCounts.dirR; i++){
//			for(int j= -dirCounts.dirR; j<= dirCounts.dirR; j++){
//				IPCounts ipCounts = dirCounts.getValue(i, j);
//				if(ipCounts != null && ipCounts.ipMap != null){
//					int s = dirCounts.getValue(i, j).ipMap.size();
//					if(s==0){
//						count0++;
//					}else{
//						counts[U.toI(Math.log10(s))] ++;				
//					}
//					maxDirCount = max(maxDirCount, s);
//				}else{
//					countN++;
//				}
//			}
//		}
//		
//		System.out.println("maxDirCount = "+maxDirCount);
//		System.out.println("N : "+countN);
//		System.out.println("0 : "+count0);
//		for(int i=0; i<counts.length; i++){
//			System.out.println(">"+pow(10, i)+" : "+counts[i]);
//		}

	}
		
	
	public void testIPDirRec(){
		SegCurve curve = getTestCurve(false);

		IPDirRec dataRec = new IPDirRec(500, 50, 5, 20, 60);
		dataRec.setCurve(curve);
		
		VList[][] ipDirs = dataRec.getIpDirs();
		dataRec.getDirIPs();
		
		dataRec.dealFarDirs();
		dataRec.dealNullDirs();
		
		Color[][] colorPave = {
				{Color.red, Color.white}, {Color.green, Color.blue}
		};

		pt.setColor(Color.gray);
		fiRect(0, 0, dataRec.rDir, dataRec.rDir);
		
		for(int i=0; i<=dataRec.rIp*2; i++){
			for(int j=0; j<=dataRec.rIp*2; j++){
				
//				if(i!=0 || j!=20)continue;
				
				pt.setColor(colorPave[i%2][j%2]);
				if(ipDirs[i][j]!=null){
					for(V dir : ipDirs[i][j].dirs){
						drPointXY(dir);
					}
				}
			}
		}

		pt.setColor(Color.black);
		for(int A=0; A<= dataRec.rangeA; A+=5){
			pt.setColor(A%10==0 ? Color.black : Color.yellow);
			drCircle(0, 0, dataRec.getHDir()*T(A));
		}
	
	
		System.out.println("fin");
	}

	
	public void testSfcCllMtx(){
		SegCurve curve = getTestCurve(false);
		IPDirRec ipDirRec = new IPDirRec(500, 50, 5, 20, 60);
		ipDirRec.setCurve(curve);
		
		SfcCllMtx mtx = new SfcCllMtx(curve, ipDirRec, null);
		mtx.initObs(1, 300, 45);
		mtx.initClls();
		
		mtx.calMtxScnColors(1);
		
		double maxSeg = 0;
		for(V[] vs: mtx.scnColors){
			if(vs==null)continue;
			for(V v : vs){
				if(v==null)continue;
				maxSeg = max( maxSeg, max(v.x, max(v.y, v.z)));
			}
		}
		
		for(int i=-mtx.pxScnR; i<mtx.pxScnR; i++){
			for(int j=-mtx.pxScnR; j<mtx.pxScnR; j++){
				V vij = mtx.scnColors[i+mtx.pxScnR][j+mtx.pxScnR];
				Color cij = vij==null? mtx.mtxBgColor : vToColor(vij, maxSeg, 1);
				pt.setColor(cij);
				drPoint(i, j);
			}
		}
		
	}
	
	
	public Color vToColor(V vColor, double max, double colorPowerIndex){
//		max = max(vColor.x, max(vColor.y, vColor.z));
		
		if(colorPowerIndex == 1){
			int R = U.toI(vColor.x/max *255);
			int G = U.toI(vColor.y/max *255);
			int B = U.toI(vColor.z/max *255);
			return new Color(R, G, B);
		}else{
			int R = U.toI(pow(vColor.x/max, colorPowerIndex) *255);
			int G = U.toI(pow(vColor.y/max, colorPowerIndex) *255);
			int B = U.toI(pow(vColor.z/max, colorPowerIndex) *255);
			return new Color(R, G, B);
		}
	}
	
	public Color vToColor2(V vColor, double max){
		max = max(vColor.x, max(vColor.y, vColor.z));
		
		int R = U.toI(vColor.x/max *255);
		int G = U.toI(vColor.y/max *255);
		int B = U.toI(vColor.z/max *255);
		return new Color(R, G, B);
	}
	
	public void showTrm(){
		double r = 200;
		double h = 1.9*r;
		double l = 1.69*r;
		double dr = 20;
		double n = 1.5;
		double r0 = 1; double b = -30;

		V O = new V(r0*r,0,0);
		V f = V.ROLL(new V(0, 0, 1), new V(0, 1, 0), b);
		pt.setColor(Color.black);
		drLXY(new L(O, f), 200);
					
		double[][] douts = new double[new Double(dr).intValue()][2];
		double lastA = 0;
		for(int dx=0; dx<=dr; dx++){
			double x = dx*l/dr;
			V pO = new V(-x, -h, 0);
			
			V out = trm(f, sub(O, pO), 1/n);
			
			V out2 = null;

			double S = n*S(AT(x/h+r0*r/h)+b);
			if(!(abs(S)>1)){
				double y = AS(S)-b;
				y = 90-y;
				out2 = new V(C(y), S(y), 0);
			}
			
			pt.setColor(Color.black);
			drLineXY(pO, O);
			if(out!=null){
				pt.setColor(Color.gray);
				drLXY(new L(O, out), 800);
				if(dx>0){
					System.out.println(out.AXY()-lastA);
					douts[dx-1] = new double[]{out.AXY()-lastA, 1};
					
				}
				lastA = out.AXY();
			}else{
				System.out.println("out of range");
			}
			
			if(out2!=null){
				pt.setColor(Color.red);
				drLXY(new L(O, out2), 200);
				System.out.println("  "+out2);
			}
			
		}
		
		System.out.println(U.SD(douts));
	}
	
	public void statisticA(){
		double n = 1.5;
		
		int accu = 200;
		int accu2 = 50;
		
//		int vPerAccu2 = 100;
		
		double maxOut = 90;
		
		double r=200;
		double h = r*2;
		double R = r*2;
		
		
		// SegCurve
		if(curve==null){
			double hFind = SegCurve.findAverageH(100, 1.5);
			System.out.println(hFind);
			
//			curve = new SegCurve(r*0.95, r*2, r*2*hFind*1.3, 1.5); // stat with h find
			curve = new SegCurve(r, r*2.4, r*2.7, 1.5);
			curve.setRange(80, accu);
			curve.setCenter(new V(0, 0, 0));
//			curve.setHue(7, 100);
			
			h = curve.h;
			R = curve.l;

			curve = getInitedCurve(curve, false);
			
			int dH = 0;
			curve.setHue(7, 100, 0+dH, 20+dH, 40+dH, 60+dH, 80+dH, 100+dH, 120+dH, 140+dH, 150+dH);
			
		}
		

		V[] Ps = new V[accu+1]; // x, y, tg
		Ps[0] = new V(0, 0, 0);

		for(int i=1; i<=accu; i++){
			L seg = curve.getPSeg(i);
			Ps[i] = new V(seg.o.x, seg.o.y, seg.dir.tgXY());
		}
		
		pt.setColor(Color.black);
		// curve
//		for(int i=1; i<=accu; i++){
//			V pO = new V(-i*R/accu, -h, 0);
//			V pXY = new V(i*r/accu, Ps[i-1].y+r/accu*Ps[i-1].z, 0);
//			
//			V dirIn = sub(pXY, pO).unit();
//
//			double in = AT(dirIn.x/dirIn.y) ;
//			
//			double out1 = pow(i*1d/accu, 1 );
////			double out1 = i*1d/accu;
////			double out = (i/10)*10 *maxOutd/accu;
//			
//			double out = out1 *  maxOut;
////			double out = new Double(out1*100d).intValue()/100d  *  maxOutd;
////			double out = new Double(out1 * maxOutd).intValue();
//
//					
//			double Tc = (n*S(in) -S(out)) / (n*C(in) -C(out));
//			
//			Ps[i] = new V(pXY.x, pXY.y, -Tc);			
//		}

		pt.setColor(Color.pink);
		for(int i=0; i<accu; i++){
			drLine(Ps[i].x, Ps[i].y, Ps[i+1].x, Ps[i+1].y);
		}
		
		// oval
//		double rx = h*0.5;
//		double ry = rx;		
//		double rx = h*0.60;
//		double ry = h*0.85;		
//		drOval(0, -ry, rx, ry);
//		for(int i=1; i<=accu; i++){
//			V pXY = new V(i*r/accu, Ps[i-1].y+r/accu*Ps[i-1].z, 0);
//			double Tc = abs(pXY.x/r)<0.49 ? (pXY.x*ry*ry/(pXY.y+ry)/rx/rx) : 0.5d;
//			Ps[i] = new V(pXY.x, pXY.y, -Tc);						
//		}

		
		boolean showSingle = false;
//		showSingle = true;
		
		double scaleX = 1.5*600d/accu;
		double scaleY = 2*pow(accu2/100d, 0.3333);

		pt.setColor(new Color(230,230,230));
		for(double x = -accu2; x<=accu2; x++){
//			double xscale = (x) * 10;
			drLine(scaleX*x*accu/accu2, -h, scaleX*x*accu/accu2, h);
		}
		pt.setColor(new Color(220,220,220));
		drLine(0, -h, 0, h);

		pt.setColor(Color.black);
		for(int i=0; i<accu; i++){
			drLine(Ps[i].x, Ps[i].y, Ps[i+1].x, Ps[i+1].y);
		}

//		if(true)return;

		
		
		float trs = 1;
		Color[] colors = {
				new Color(1f, 0f, 0f, trs),
				new Color(0f, 1f, 0f, trs),
				new Color(0f, 0f, 1f, trs),
				new Color(0f, 1f, 1f, trs),
				new Color(0f, 0f, 0f, trs),
		};

		int colorIndex = 0;

		XYStatA statI = null;
		int dirY = 1;
		for(int i=0; i<=accu; i+=1){
			if(new Double(i/(accu/accu2)).intValue()%4!=0){
//				continue;
			}
			
			V pO = new V(-i*R/accu, -h, 0);
			
			boolean changeAccu2 = i%(accu/accu2) == 0;
			if(changeAccu2){
				if(statI != null){
					if(showSingle){
						pt.setColor(Color.white);
						pt.fillRect(0, 0, w, 2000);
					}

					pt.setColor(colors[(colorIndex++)%colors.length]);
					statI.show(this, scaleX, scaleY*accu2/accu*dirY);
//					dirY *= -1;
				}
//				statI = new XYStatA(-maxOut, maxOut, accu2*2*vPerAccu2);
				statI = new XYStatA(-maxOut, maxOut, accu*2);
			}
			for(int j=-accu+1; j<accu;j++){
				double tg = j<0? -Ps[-j].z : Ps[j].z;
				V pXY = j<0? new V(-(Ps[-j].x+Ps[-j+1].x)/2, (Ps[-j].y+Ps[-j+1].y)/2, 0) : new V((Ps[j].x+Ps[j+1].x)/2, (Ps[j].y+Ps[j+1].y)/2, 0);
				V dirIn = sub(pXY, pO).unit();
				V f = new V(-tg, 1, 0);
				V dirOut = trm(f, dirIn, 1/n);
				if(dirOut!=null){
					double a = AT(dirOut.x/dirOut.y);
//					statI.addValue(a);

					double b = AT(tg);
					double lj = 1/C(b);
					statI.addValue(a, lj * U.s(U.ac(V.multDU(dirOut, new V(1, tg, 0)))));
				}
			}
			
		}
		if(statI != null){
			pt.setColor(colors[(colorIndex++)%colors.length]);
			statI.show(this, scaleX, scaleY*accu2/accu);
		}

		
	}
	
	public void showFunction(){
			P1P2AOutToADir fun = new P1P2AOutToADir(1/1.5);
			drFun(fun, new V(0,-320,0), new V(2,2,0));
		
	}
	
	public static class XYStatA{
		double vMin;
		double vMax;
		int accu;
		
		double[] count;
		
		public XYStatA(double vMin, double vMax, int accu){
			this.accu = accu;
			this.vMax = vMax;
			this.vMin = vMin;
			count = new double[accu+1];
		}
		
		public void addValue(double v){
			if(v<=vMax && v>=vMin){
				count[new Double((v-vMin)/(vMax-vMin)*accu).intValue()]++;
			}
		}
		
		public void addValue(double v, double c){
			if(v<=vMax && v>=vMin){
				count[new Double((v-vMin)/(vMax-vMin)*accu).intValue()]+=c;
			}
		}
		
		public void show(Drawer drawer, double scalex, double scaley){
			for(int i=1; i<count.length; i++){
				drawer.drLine((i-1-count.length/2) * scalex, count[i-1] * scaley, (i-count.length/2) * scalex, count[i] * scaley);
			}
			
		}
		
	}
	


	
	static SegCurve curve = null;
	public void test7V4(){
		double r = 100;
		if(curve==null){
			curve = new SegCurve(r, r*2, r*1.9, 1.5);
			curve.setRange(90, 5);
			curve.setCenter(new V(0, -300, 0));
//			curve.setHue(7, 100);
			
			int dH = 0;
			curve.setHue(7, 100, 0+dH, 20+dH, 40+dH, 60+dH, 80+dH, 100+dH, 120+dH, 140+dH, 150+dH);
			
			curve.initPArray();
		}

		Color gray = new Color(0.2f,0.2f,0.2f,0.05f);
		pt.setColor(gray);

		if(false)
		{
//		double x=1; 
//		double y=curve.center.y+200;
		for(double x=-w/2; x<w/2; x+=1)
		{
			for(double y = curve.center.y; y<h/2; y+=1)
//			for(double y = curve.center.y+300; y<curve.center.y+310; y+=1)
			{
				V pXY = new V(x, y, 0);			
				
				// r condition
				double absR = sub(pXY, curve.center).abs();
				if(absR<=h/2 || absR>=h/2+10){
					continue;
				}
				
				V colorSumRGB = new V(0, 0, 0);
				double areaSum = 0d;
				
				for(int i=0; i<=curve.dv; i+=1){
					for(int dr = -1; dr<=1; dr+=2){
						L seg = dr==-1?curve.getNSeg(abs(i)) : curve.getPSeg(abs(i));
						
						V dir = sub(seg.o, pXY);
						
						Double tgIn = dir.tgXY();
						Double tgSeg = seg.dir.tgXY();
						if(tgIn!=null && tgSeg!=null && tgIn*tgSeg>0 && abs(tgSeg)>abs(tgIn)){
							continue;
						}

						L lo = new L(seg.o, trm(new V(-seg.dir.y, seg.dir.x, seg.dir.z), dir, curve.n));
						if(lo!=null && lo.dir.y<0){
							// get cross point:
							double distY = curve.center.y-curve.h - lo.o.y;
							double xImg = lo.o.x + distY/lo.dir.y * lo.dir.x;
							V pImg = new V(xImg, curve.center.y-curve.h, 0);
							Color c = curve.getColor(xImg-curve.center.x, 1f);
							
							double ASeg = tgSeg==null ? 90 : AT(tgSeg);
							double AIn = tgIn==null ? 90 : AT(tgIn);
							double l = abs(curve.r/curve.dv / C(ASeg) * S(ASeg-AIn));
							
							colorSumRGB.x += c.getRed()*l;
							colorSumRGB.y += c.getGreen()*l;
							colorSumRGB.z += c.getBlue()*l;
							areaSum += l;
							
//							pt.setColor(gray);
//							drLXY(lo, 900);
//							drLineXY(pXY, seg.o);
							
//							pt.setColor(c);
//							drPoint(pImg.x, pImg.y);
//							drLXY(new L(seg.o, new V(-dir.y, dir.x, 0)), l);
						}
					}
				}
				
				double areaTotal = areaSum;
//				double areaTotal = 2*r;
				
				Color cXY = new Color(
						new Double(colorSumRGB.x/areaTotal).intValue(),
						new Double(colorSumRGB.y/areaTotal).intValue(),
						new Double(colorSumRGB.z/areaTotal).intValue()
						);

				pt.setColor(cXY);
				drPoint(x, y);
				
			}
		}
		}
		
		
		pt.setColor(gray);
		drLine(-w, -curve.h+curve.center.y, w, -curve.h+curve.center.y);
		drLine(0, -900, 0, 900);
		for(int i=0; i<curve.dv; i++){
			drLineXY(curve.getPSeg(i).o, curve.getPSeg(i+1).o);
			drLineXY(curve.getNSeg(i).o, curve.getNSeg(i+1).o);
		}
		
	}

	public void showHSQX(){
		double n = 1.5;
		double r0 = -0.5;
		double h = 200;
		double b = 0;

		
		pt.setColor(Color.lightGray);
		for(double i=-100; i<=100; i+=10){
			drLine(0, i, w, i);
		}
		pt.setColor(Color.gray);
		drLine(0, 0, w, 0);
		drLine(0, 0, w, 0);
		
		pt.setColor(Color.black);

		for(; b<40; b+=3){
//		for(; h>150; h-=2){
//		for(r0=-0.5; r0<0.5; r0+=0.05, h-=U.unit(r0)*4){
			Double y0 = null;
			for(double x=0; x<h*(1-r0)*h; x++){
				double y = (r0>=0) ? AS(n*S(AT(x/h+r0)+b))-b  :  AS(n*S(AT(x/h+r0)-b))+b;
				
				if(y0 != null){
					drLine(x-1, y0, x, y);
				}
				y0 = y;
				System.out.println(y);
			}
		}
//		}
		
	}
	
	public void test6(){
//		n²x²  + 2n²T(a0)hx  =  y²[ S²(a1) - n²]     +   y {nS(a1)*[(r+T(a0)h)² + n²(h)²]° - n²h}  +    n²(r+T(a0)h)² -n²T²(a0)h² 
		
		
		double n = 1.5d; 
		double r = 200;
		double h = r * 2;
		double br = r/1;

		double sf = 1;
		
//		double drn = n>1? 1 : -1;
		double drn = 1;
		
		double filtreA = 90;
		
		double n2 = 1.5;
		double yP = 0;
		double filtreP = AS(1/n2);
		
//		double odAMax = 90;
		
		double dY = -0;
		
		double scr = 1;
		double scr1 = 1;

		double A0M = AS(1/n);
//		double A0 = 0;		
//		for(double odAx = 0; odAx<=2*r*scr; odAx+= 2*r*scr/10){
			
//		for(double odA = 45; odA<=odAMax; odA=odAMax+1){
		for(double A0 = 0; A0<=A0M; A0+=A0M/20)
//		for(double odA = odAMax; odA>=0; odA-=odAMax/20){
		{
			double A1 = 90-AS(S(A0)*n);
			double odAx =h*T(A0);
			
//			Double hue = (odA*4/odAMax) % 1;
//			Double hue = (odAx*2/(2*r*scr)) % 1;
//			Double hue = 0.7d;
			Double hue = (A0*3/A0M) % 1;
			Color color = Color.getHSBColor(hue.floatValue(), 1, 1);
			color = new Color(((float)color.getRed())/255, ((float)color.getGreen())/255, ((float)color.getBlue())/255, 1f);

			
//			double od = h * T(odA);
			double od = odAx;
			double m = 1/n; 
//			double s = pow(br*br+h*h, 0.5);
			
			
			double x0 = r+T(A0)*h;
			double L = n*pow(pow(x0, 2) + h*h, 0.5) - r*C(A1);
			
			double xa = n*n - pow(C(A1), 2);
			double xb = 2*(n*n*T(A0)*h - L*C(A1));

			double a = S(A1)*S(A1) - n*n;
			double yb = 2* (L*S(A1) - n*n*h);
			double xya = 2*S(A1)*C(A1);
			double c0 = n*n*( pow(r+T(A0)*h, 2) -pow(T(A0)*h, 2) );
			
			c0 = L*L - n*n *h*h*( pow(T(A0),2)+1);
//			double a = (m*m-1) ;
//			double b = 2*(drn*s*m - h);
//			double c0 = (s*s-h*h);

			
			
			V o = new V(-od, -h*sf + dY, 0);

			boolean pr = true;
			
			for(double x=-r*scr1; x<=r*scr1; x+= r*scr1/40){
				
				double b = yb+xya*x;
				double c = c0-x*x*xa - x*xb;
				
				double D = b*b-4*a*c;
				if(D<0){
					System.out.println("D<0");
					continue;
				}
				double d = pow(D, 0.5);
				
//				double y = (-b-drn*d)/2/a ;
//				double y = (-b-d)/2/a ;
				
				
				double y1 = (-b-d)/2/a ;
				double y2 = (-b+d)/2/a ;
				double y;
				if(n>1){
					y = y1>=0? y1 : y2;
				}else{
					y = y2<=0? y2 : y1;
				}

				
				
//				System.out.println(x+","+y);
				
				double dxy = (2*x*(n*n - pow(C(A1), 2))+2*(n*n*T(A0)*h - L*C(A1))-y*2*S(A1)*C(A1)) /(  2*y*(pow(S(A1),2)-n*n)+2*(L*S(A1)-n*n*h)+x*2*S(A1)*C(A1) );

				L lf = new L(new V(x, y+ dY,0), new V(-dxy, 1, 0));				
				
				V dir = sub(lf.o, o);
				L li = new L(o, dir);

				L lo = null;
				V ltV = trm(lf.dir, dir, m);
				if(ltV != null){
					lo = new L(lf.o, ltV);					
				}
				
				pt.setColor(color);

//				drCircle(0, dY, 25);

//				drLXY(lf, 100);

/** /				
				if(lo!=null && AT(lo.dir.x/lo.dir.y)<filtreA){
					drLXY(lo, 900);
					drLXY(lo, -900);
					drLXY(li, dir.abs());
//					drLXY(lf, 200);
					if(pr){
//						System.out.println("od="+od+"\nodA="+odA);
						pr = false;
					}
				}
				
				/* /				
				if(lo!=null && AT(lo.dir.x/lo.dir.y)<filtreP){

					double dyP = yP+dY-lo.o.y;
					double dxP = lo.dir.x/lo.dir.y*dyP;
					V oP = add(lo.o, new V(dxP, dyP, 0));
					
					V outP = trm(new V(0,1,0), lo.dir, 1/n2);
					if(outP != null && abs(AT(outP.x/outP.y))<filtreA){
						L lop = new L(oP, outP);
						
						drLXY(lo, new V(dxP, dyP, 0).abs());
						
						drLXY(lop, 900);
						drLXY(lop, -900);

						drLXY(li, dir.abs());

					}
				}
				/**/				

//				pt.setColor(new Color(Color.gray.getRed()/255f, Color.gray.getRed()/255f, Color.gray.getRed()/255f, 0.5f));
				drPoint(x, y+ dY);
				
				pt.setColor(new Color(0.2f,0.2f,0.2f,0.01f));
				drLine(-x, dY, x, dY);

			}
			
		}
		
		
		
	}

	
	public void test4(){
//		n²(x+d1)² + n²(y+h)² =  ( L + y*S(a1) + x*C(a1) )²
//		(x+d1)² + (y+h)² =  m²( n²*d1/C(a1)  + y*S(a1) + x*C(a1))²
//
//		L = n²*d1/C(a1)  or  C=90°, d1 = 0, L = ns

//		(x)² 	=  ( L + y)² / n²  -  (y+h)²
//		dy/dx   =  x/( y(m²-1) + sm - h)
//		(x)² 	=  ( s + my)²   -  (y+h)²
		
//		if(this.equals(null))
		
		double n = 1.5; 
		double r = 5;
		double h = 10;
		double filtreA = 90;
		
		double odAMax = 90;
		
		double dX = -440;
		double dY = -430;
		
		float trans = 0.01f;
		double dr = 200;
		double dA = 0.25;

		for(double odA = 0; odA<=odAMax; odA+=dA){
//		for(double odA = odAMax; odA>=0; odA-=dA){
			
			Double hue = (odA/10) % 1;
			Color color = Color.getHSBColor(hue.floatValue(), 1, 1);
			color = new Color(((float)color.getRed())/255, ((float)color.getGreen())/255, ((float)color.getBlue())/255, trans);

			
			double od = h * t(toa(odA));
			double m = 1/n; 
			double s = pow(r*r+h*h, 0.5);
			
			double a = (m*m-1) ;
			double b = 2*(s*m - h);
			double c0 = (s*s-h*h);
			
			V o = new V(-od + dX, -h + dY, 0);

			boolean pr = true;
//			for(double x=r; x>=-r; x-= r/dr){
			for(double x=-r; x<=r; x+= r/dr){
				
				double c = c0-x*x;
				
				double D = b*b-4*a*c;
				if(D<0){
					continue;
				}
				double d = pow(D, 0.5);
				
				double y = (-b-d)/2/a ;

				
				
//				System.out.println(x+","+y);
				
				double dxy = x/( y*(m*m-1) + s*m - h);
				L lf = new L(new V(x + dX, y+ dY,0), new V(-dxy, 1, 0));

				V dir = sub(lf.o, o);
				L li = new L(o, dir);

				L lo = null;
				V ltV = trm(lf.dir, dir, m);
				if(ltV != null){
					lo = new L(lf.o, ltV);					
				}
				
				pt.setColor(color);
				
				if(lo!=null && AT(lo.dir.x/lo.dir.y)<filtreA && lo.dir.y>0){
					drLXY(lo, 900);
					drLXY(li, dir.abs());
//					drLXY(lf, 200);
					if(pr){
						System.out.println("od="+od+"\nodA="+odA);
						pr = false;
					}
				}

				pt.setColor(Color.gray);
				drPoint(x + dX, y+ dY);

			}
			
		}
		

		pt.setColor(Color.black);
		Font font = pt.getFont();
		pt.setFont(new Font(font.getName(), font.getStyle(), 8));
		for(double ruleA = 0; ruleA<=90; ruleA+=5){
			V o = new V(dX, dY, 0);
			V dir = new V(S(ruleA), C(ruleA), 0);
			V o2 = add(o, mult(dir, 400));
			V o3 = add(o, mult(dir, 390));
			L l = new L(o2, dir);
			drLXY(l, ruleA%10==0?20:10);
			
			drStr(new Double(ruleA).intValue(), o3.x-4, o3.y);
		}

		
		
	}
	
	public void test2(){
		Sph s0 = new Sph(new V(0,-100000,0), new V(0, 1, 0), 100110, toa(20));
//		Sph s1 = new Sph(new V(0,-100000,0), new V(0, 1, 0), 100105, toa(30));
//		Sph s1 = new Sph(new V(0,-100,0), new V(0, 1, 0), 205, toa(30));
		Sph s1 = new Sph(new V(0, 307,0), new V(0, -1, 0), 205, toa(50));
		Sph s2 = new Sph(new V(0, 30, 0), new V(0, 1, 0), 70, toa(50));

		Sph sa = new Sph(new V(0, 55, 0), new V(0, 1, 0), 75, toa(50));
//		Sph sa = new Sph(new V(0, 180, 0), new V(0, -1, 0), 75, toa(50));
		Ovl ovl = new Ovl(new V(0, 50, 0), new V(0, 1, 0), 50, 100, toa(100));
		
		int d = 20;
		double n = 1.5;

//		Surface[] surfaces = {s0, s1, s2};
//		double[] ns = {n, 1/n, n};

		Surface[] surfaces = {ovl};
		double[] ns = {n};

		for(int arc = 20; arc<=170; arc+=10){
			
			if((arc<90) != ((arc/10)%2==0)){
				continue;
			}
			
			
			V dir = new V(C(arc), -S(arc), 0);
			V o = sub(new V(0,110,0), mult(dir, 300));
			
			
			List<L> ls = new ArrayList<L>();
			ls.add(new L(o, dir));
			for(int i=1; i<2; i++){
				ls.add(new L(add(o, new V(d*i, 0, 0)), dir));
				ls.add(new L(add(o, new V(-d*i, 0, 0)), dir));
			}
			

			pt.setColor(cl4);
			
			for(L l : ls){
				int i=0;
				L la = l;
				L lb = null;
				surfaces[i].show(this);
				do{
					lb = surfaces[i].out(la, ns[i]);
					if(lb!=null){
						drLineXY(la.o, lb.o);
					}
					la = lb;
					i++;
				}while(lb!=null && lb.dir!=null && i<surfaces.length);

				if(lb!=null && lb.dir!=null){
					drLXY(lb, 500);
				}
			}
		}
	}
	
	public void test1(){
		double r = 1;
		
		System.out.println(AT(40d*4/(1000 * 5))*2);
		
//		for(Double a : new Double[]{0d, 9d, 18d, 27d, 36d, 45d}){
		for(Double a : new Double[]{0d, 9d, 18d, 27d, 36d, 45d, 54d, 63d}){
			pt.setColor(cl1);
//			Params prm = new Params(r*1000 * 10, -r*1000 * 2.5, r*40, r/2, 1.4, a);
			Params prm = new Params(r*1000 * 10, -r*1000 * 2, r*40, r, 1.3, a);
			System.out.println(prm.R);

			V v = null;
			for(double b=0; b<=360; b++){
				V border = prm.border_0(b); 
				V vb = border==null ? null : sub(border, prm.centre);
				
				if(v != null && vb != null){
					drLineXY(v, vb);
				}
				v = vb;
			}
			pt.setColor(cl2);
			drCircle(0, 0, prm.r);
//			drCircle(0, 0, prm.h);
			drCircle(0, 0, prm.R);
		}
	}
	
}

public class TestOpt {
//	static int w = 1900;
//	static int h = 300;

//	static int w = 1900;
//	static int h = 1000;
	static int w = 1260;
	static int h = 670;

	static Canvas cvs;
	public static void main(String[] args) {
		ObjSaver.setFolderPath(args[0]);
		
		//T ijk
//		int ri=20000, rj = 20000;
//		IJtoK ijk = new IJtoK(ri, rj);
//		for(int i=-ri; i<=ri; i++){
//			for(int j=-rj; j<=rj; j++){
//				int k = ijk.getK(i, j);
//				
//				if(i!=ijk.getI(k) || j!=ijk.getJ(k)){
//					System.out.println(i+","+j);
//					System.out.println(ijk.getI(k)+","+ijk.getJ(k));
//				}
//			}
//		}
//		System.out.println("fin");
//		if(true)return;

//		System.out.println(AS(1/1.33));
//		System.out.println(AS(1/1.4));
//		System.out.println(AS(1/1.5));
//
//		System.out.println(AS(S(72)/1.5));
//		if(true)return;
		
		JFrame jf = new JFrame();
		jf.setSize(w+15, h+40);
		jf.setLayout(null);
		jf.setLocation(0, 0);
		
		cvs = new Canvas(){
			public void paint(Graphics g){
//				g.setColor(Color.white);
//				g.fillRect(0, 0, w, h);
				Drawer dr = new Drawer(w, h, g);
				dr.show();
			}
		};
		cvs.setSize(w, h);
		cvs.setLocation(0, 0);
		jf.add(cvs);
		jf.setVisible(true);
		jf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		Thread th = new Thread(){
			public void run(){
				for(int i=0; i<10000000; i++){
					long t = System.currentTimeMillis();
					cvs.repaint();
					try {
						t = System.currentTimeMillis() - t;
						sleep(Math.max(100-t, 1));
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		};
		if(Drawer.threadOpen){
			th.start();
		}
		
	}


}
		
class U{
	public static double toa(double arc){
		return PI*arc/180d;
	}

	public static V spherePositionToPlat(int axIndex, double rpl, V vSph) {
		double r = rpl*2/PI;
		
		V pole = null;
		if(axIndex==0){
			pole = new V(1, 0, 0);
		}else if(axIndex==1){
			pole = new V(0, 1, 0);
		}else if(axIndex==2){
			pole = new V(0, 0, 1);
		}else{
			return null;
		}
		
		double ca = V.multDU(pole, vSph);
		double r12 = r*ac(ca);
		
		if(axIndex==0){
			return mult(new V(0, vSph.y, vSph.z).unit(), r12);
		}else if(axIndex==1){
			return mult(new V(vSph.x, 0, vSph.z).unit(), r12);
		}else if(axIndex==2){
			return mult(new V(vSph.x, vSph.y, 0).unit(), r12);
		}else{
			return null;
		}		
	}

	public static V platToSpherePosition(int axIndex, double rpl, V vpl) {
		return platToSpherePosition(axIndex, rpl, vpl.x, vpl.y, vpl.z);
	}
	
	public static V platToSpherePosition(int axIndex, double rpl, double ... vpl) {
		double v1 = 0;
		double v2 = 0;
		
		if(axIndex==0){
			v1 = vpl[1]+0;
			v2 = vpl[2]+0;
		}else if(axIndex==1){
			v1 = vpl[0]+0;
			v2 = vpl[2]+0;
		}else if(axIndex==2){
			v1 = vpl[0]+0;
			v2 = vpl[1]+0;
		}else{
			return null;
		}
		
		double r = rpl*2/PI;
		double r12 = pow(v1*v1+v2*v2, 0.5);
		double a = r12/r; 
		
		double v3 = r*c(a);
		double v12 = r*sin(a);
		v1 = U.is0(v1)? 0 : v1/r12 * v12;
		v2 = U.is0(v2)? 0 : v2/r12 * v12;
		
		if(axIndex==0){
			return new V(v3, v1, v2);
		}else if(axIndex==1){
			return new V(v1, v3, v2);
		}else if(axIndex==2){
			return new V(v1, v2, v3);
		}else{
			return null;
		}
		
	}

	public static double pow2(double d){
		return d*d;
	}
	
	public static int toI(double d){
		return new Double(d).intValue();
	}

	public static int toI45(double d){
		return new Double(d>0? d+0.5 : (d<0 ? d-0.5 : 0)).intValue();
	}

	public static double toA(double arc){
		return arc*180d/PI;
	}

	public static double c(double arc){
		return cos(arc);
	}
	
	public static double s(double arc){
		return sin(arc);
	}
	
	public static double t(double arc){
		return tan(arc);
	}
	
	public static double ac(double v){
		return acos(v);
	}
	
	public static double as(double v){
		return asin(v);
	}
	
	public static double at(double v){
		return atan(v);
	}

	
	public static double C(double arc){
		return c(toa(arc));
	}
	
	public static double S(double arc){
		return s(toa(arc));
	}
	
	public static double T(double arc){
		return t(toa(arc));
	}
	
	public static double AC(double v){
		return toA(acos(v));
	}
	
	public static double AS(double v){
		return toA(asin(v));
	}
	
	public static double AT(double v){
		return toA(atan(v));
	}
	
	public static V trm(V f, V in, double n){
		V uf = f.unit();
		V uIn = in.unit();
		double c = multD(uf, uIn);
		double a = ac(c);
		double s2 = s(a)/n;
		if(s2 >1){
			return null;
		}
		
		double c2 = abs(c(as(s2)))*(c>0 ? 1 : -1);
		V x = is0(s2)? new V(0,0,0) : sub(uIn, mult(uf, c)).unit();

		return add(mult(uf, c2), mult(x, s2)).unit();
	}
	
	public static boolean is0(double d){
		return d==0 || abs(d)<0.00000001;
	}

	public static boolean absEqual(double a, double b){
		return equal(abs(a), abs(b));
	}
	
	public static boolean aAbsSmallEqualb(double a, double b){
		return aSmallEqualb(abs(a), abs(b));
	}
	
	public static boolean aAbsBigEqualb(double a, double b){
		return aBigEqualb(abs(a), abs(b));
	}

	public static boolean aAbsSmallb(double a, double b){
		return !aAbsBigEqualb(a, b);
	}
	
	public static boolean aAbsBigb(double a, double b){
		return !aAbsSmallEqualb(a, b);
	}

	public static boolean equal(double a, double b){
		return a==b
				|| (is0(a) && b==0)
				|| (is0(b) && a==0)
				|| is0((a-b)/max(abs(a),abs(b)))
		;
	}
	
	public static boolean aSmallb(double a, double b){
		return !aBigEqualb(a, b);
	}
	
	public static boolean aBigb(double a, double b){
		return !aSmallEqualb(a, b);
	}

	public static boolean aSmallEqualb(double a, double b){
		return a<b || equal(a, b);
	}
	
	public static boolean aBigEqualb(double a, double b){
		return a>b || equal(a, b);
	}

	public static boolean aInBC(double a, double b, double c){
		double small = min(b, c);
		double big = max(b, c);
		return aBigEqualb(a, small) && aSmallEqualb(a, big);
	}
	public static boolean aInBCStrict(double a, double b, double c){
		double small = min(b, c);
		double big = max(b, c);
		return !aSmallEqualb(a, small) && !aBigEqualb(a, big);
	}
	
	public static boolean aAbsInBC(double a, double b, double c){
		return aInBC(abs(a), abs(b), abs(c));
	}
	public static boolean aAbsInBCStrict(double a, double b, double c){
		return aInBCStrict(abs(a), abs(b), abs(c));
	}

	public static double[] defun(double a, double b, double c){
		if(is0(a) && is0(b)){
			return null;
		}else if(is0(a)){
			return new double[]{-c/b};
		}
		double d = b*b-4*a*c;
		if(is0(d)){
			return new double[]{-b/2/a};
		}else if(d<0){
			return null;
		}else {
			double rd = pow(d, 0.5);
			double r1 = (rd-b)/2/a;
			double r2 = (-rd-b)/2/a;
			return new double[]{min(r1, r2), max(r1, r2)};
		}
	}
	
	public static double unit(Double x){
		if(x==null || is0(x)){
			return 0;
		}
		
		return x/abs(x);
	}
	
	// data : {{value, count}, {value, count}}
	public static Double SD(double[][] datas){
		double N = 0;
		double sum = 0;
		double psum = 0;
		for(double[] data : datas){
			N += data[1];
			sum += data[0]*data[1];
			psum += data[0]*data[0]*data[1];
		}
		
		if(is0(N)){
			return null;
		}
		
		return pow(psum/N - pow(sum/N, 2), 0.5);
	}

	public static Double SD(double[][] datas, double v){
		double N = 0;
		double psum = 0;
		for(double[] data : datas){
			N += data[1];
			psum += pow(data[0]-v, 2)*data[1];
		}
		
		if(is0(N)){
			return null;
		}
		
		return pow(psum/N, 0.5);
	}

	public static Double SD_B(double[][] datas, double v){
		double N = 0;
		double psum = 0;
		int i=0;
		for(double[] data : datas){
			if(abs(data[0])<abs(v)){
				N += data[1];
			}
			psum += data[1]/(1+abs(data[0]/abs(v)));
//			psum += (i<datas.length/3?1 : 1) /(1+abs(data[0]/abs(v)));
//			psum += data[1];
			i++;
		}
		
		if(is0(N)){
			return null;
		}
		
		return psum;
	}


}

class V implements Serializable{
	private static final long serialVersionUID = "V".hashCode();

	static V AXIS_X = new V(1, 0, 0);
	static V AXIS_Y = new V(0, 1, 0);
	static V AXIS_Z = new V(0, 0, 1);

	
	double x;
	double y;
	double z;
	public V(double x, double y, double z){
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public Double AXY() {
		Double a = aXY();
		return a==null ? null :U.toA(a);
	}

	public Double aXY() {
		if(U.is0(x) && U.is0(y)){
			return null;
		}
		
		if(U.is0(x)){
			return y>0? PI/2 :-PI/2;
		}
		
		if(U.is0(y)){
			return x>0? 0 :PI;
		}

		double a = U.at(y/x);
		
		return x<0? a+PI : a;
	}

	public V unit(){
		double r = mR(this);
		if(r==0){
			return null;
		}
		return new V(x/r, y/r, z/r);
	}
	
	public V(V v){
		this.x = v.x;
		this.y = v.y;
		this.z = v.z;
	}
	
	public boolean equals(Object o){
		if(o instanceof V){
			V v = (V)o;
			return U.is0(v.x-x) && U.is0(v.y-y) && U.is0(v.z-z);
		}
		
		return false;
	}
	
	public int hashCode(){
		return toString().hashCode();
	}
	
	public Double tgXY(){
		if(x==0){
			return null;
		}
		
		return y/x;
	}
	
	public V projXY(){
		return new V(x, y, 0);
	}
	
	public V projYZ(){
		return new V(0, y, z);
	}
	
	public V projXZ(){
		return new V(x, 0, z);
	}
	
	public String toString(){
		return "("+x+","+y+","+z+")";
	}
	
	public static V N(V v){
		return new V(-v.x, -v.y, -v.z);
	}

	public static V add(V v1, V ... vs){
		V res = new V(v1);
		
		for(V v:vs){
			res.x += v.x;
			res.y += v.y;
			res.z += v.z;
		}
		
		return res;
	}

	public static V sub(V v1, V v2){
		return add(v1, N(v2));
	}
	
	public static V mult(V v, double d){
		return new V(v.x * d, v.y * d, v.z * d);
	}

	public static V mult(V v1, double d1, double d2, double d3){
		return new V(v1.x * d1, v1.y * d2, v1.z * d3);
	}

	public static V mult(V v1, V v2){
		return mult(v1, v2.x, v2.y, v2.z);
	}

	public static V multV(V v1, V v2){
		return new V(v1.y*v2.z-v2.y*v1.z, v1.z*v2.x-v2.z*v1.x, v1.x*v2.y-v2.x*v1.y);
	}

	public static V multVU(V v1, V v2){
		return multV(v1, v2).unit();
	}
	
	public static double multD(V v1, V v2){
		return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	}

	public static double multD2(V v){
		return multD(v, v);
	}

	public static double multDU(V v1, V v2){
		return multD(v1.unit(), v2.unit());
	}

	public static V prodV(V v1, V v2){
		return new V(
				v1.y*v2.z - v1.z*v2.y,
				v1.z*v2.x - v1.x*v2.z,
				v1.x*v2.y - v1.y*v2.x
		);
	}
	
	public static V prodVU(V v1, V v2){
		return prodV(v1.unit(), v2.unit());
	}
	
	public static double mR(V v){
		return pow(v.x*v.x + v.y*v.y + v.z*v.z, 0.5);
	}

	public static double mXY(V v){
		return pow(v.x*v.x + v.y*v.y, 0.5);
	}

	public static double moXZ(V v){
		return pow(v.x*v.x + v.z*v.z, 0.5);
	}

	public static double moYZ(V v){
		return pow(v.y*v.y + v.z*v.z, 0.5);
	}
	
	public static double a(V v1, V v2){
		return ac(multDU(v1, v2));
	}

	public double abs() {
		return pow(multD(this, this), 0.5);
	}

	public double absXY() {
		return pow(this.x*this.x+this.y*this.y, 0.5);
	}

	public double absXZ() {
		return pow(this.x*this.x+this.z*this.z, 0.5);
	}

	public double absYZ() {
		return pow(this.z*this.z+this.y*this.y, 0.5);
	}

	public static double dist(V p1, V p2){
		return sub(p1, p2).abs();
	}

	public static double dist(L l, V p){
		return prodV(sub(p, l.o), l.dir.unit()).abs();
	}
	
	public static double distFoot(L l, V p){
		return multD(sub(p, l.o), l.dir.unit());
	}
	
	public static double dist(L l1, L l2){
		return Math.abs(multD(sub(l1.o, l2.o), prodV(l1.dir, l2.dir).unit()));
	}

	public static V ROLL(L axis, V v, double A){
		return roll(axis, v, U.toa(A));
	}

	public static V roll(L axis, V v, double a){
		return add(axis.o, roll(axis.dir, sub(v, axis.o), a));
	}

	public static V ROLL(V axis, V v, double A){
		return roll(axis, v, U.toa(A));
	}

	public static V roll(V axis, V v, double a){
		if(v==null){
			return null;
		}
		
		if(axis==null){
			return v;
		}
		
		axis = axis.unit();
		double s = U.s(a);
		double c = U.c(a);
		double x = axis.x;
		double y = axis.y;
		double z = axis.z;
		double xx = axis.x * axis.x;
		double yy = axis.y * axis.y;
		double zz = axis.z * axis.z;
		double xy = axis.x * axis.y;
		double yz = axis.y * axis.z;
		double zx = axis.z * axis.x;
		
		V[] matrix = {
				new V(c+(1-c)*xx, (1-c)*xy-s*z, (1-c)*zx+s*y),
				new V((1-c)*xy+s*z, c+(1-c)*yy, (1-c)*yz-s*x),
				new V((1-c)*zx-s*y, (1-c)*yz+s*x, c+(1-c)*zz)
		};
		
		return new V(multD(matrix[0], v), multD(matrix[1], v), multD(matrix[2], v));		
	}
	
	public static L ROLL(V axis, L l, double A){
		return roll(axis, l, U.toa(A));
	}

	public static L roll(V axis, L l, double a) {
		return new L(roll(axis, l.o, a), roll(axis, l.dir, a));
	}

	public static L ROLL(L axis, L l, double A){
		return roll(axis, l, U.toa(A));
	}

	public static L roll(L axis, L l, double toa) {
		return new L(roll(axis, l.o, toa), roll(axis.dir, l.dir, toa));
	}



}

abstract class FunReal{
	public abstract double[] fun(double[] xs);

	public abstract void show(Drawer drawer, V pO, V scale);
}

abstract class FunRealXY extends FunReal{
	public abstract Double fun(double x);	
	public double[] fun(double[] xs){
		Double y = fun(xs[0]);
		return y==null? null : new double[]{y};
	} 
}

class OutToADir extends FunRealXY{
	double n;
	
	public OutToADir(double n){
		this.n = n;
	}
	
	@Override
	public Double fun(double x) {
		double y = AT(S(x)/(C(x)-1/n));

		if(abs(x-y)>90 || abs(y)>90){
			return null;
		}
		
		if(abs(S(y)/n)>1){
			return null;
		}
		
		if(!U.is0(S(y)-S(y-x)*n)){
			return null;
		}
		return y;
	}

	@Override
	public void show(Drawer drawer, V pO, V scale) {
		// TODO Auto-generated method stub
		
	}
}

class P1P2AOutToADir extends FunReal{
	OutToADir fOutToADir;
	Double filtreOutA = null;
	double n;
	
	public P1P2AOutToADir(double n) {
		this.n = n;
		this.fOutToADir = new OutToADir(n);
	}
	
	public P1P2AOutToADir(double n, Double filtreOutA) {
		this.n = n;
		this.fOutToADir = new OutToADir(n);
		this.filtreOutA = filtreOutA;
	}
	
	public Double fun(V pO, V p, double outA){
		V in = sub(p, pO);
		if(U.is0(in.absXY())){
			return null;
		}
		
		double AIn = in.AXY();
		Double fOut = fOutToADir.fun(outA-AIn);
		return fOut==null ? null : fOut+AIn-90;
		
	}
	
	@Override
	public double[] fun(double[] xs) {
		Double ADir = fun(new V(xs[0], xs[1], xs[2]), new V(xs[3], xs[4], xs[5]), xs[6]);
		if(ADir==null){
			return null;
		}else{
			return new double[]{ADir};
		}
	}
	
	@Override
	public void show(Drawer drawer, V pC, V scale) {
		double rangeXY = 200;
		double dXY = 5;
		
//		for(double i=-rangeXY; i<=0; i+=dXY){
		for(double i=-rangeXY; i<=rangeXY; i+=dXY){
			filtreOutA = 180d * (i+rangeXY)/rangeXY/2;
			V pOI = add(pC, V.mult(new V(i, 0, 0), scale));
			showSingle(drawer, pOI, pC, scale);
		}
		System.out.println("finished");
	}
	
	public void showSingle(Drawer drawer, V pO, V pC, V scale) {
		double rangeXY = 320;
		double dXY = 10;
		double dA = 4.5/2;
		
		for(double outA=0; outA<=180; outA += dA){
			if(filtreOutA != null && !U.is0(outA-filtreOutA)){
				continue;
			}
			float hue = new Double(outA/200).floatValue();
			Color cl = Color.getHSBColor(hue, 1, 1);
//			cl = new Color(cl.getRed()/255f, cl.getGreen()/255f, cl.getBlue()/255f, 0.02f);
//			cl = new Color(0f,0f,0f, 0.02f);

			for(double x = -rangeXY; x<=rangeXY; x+=dXY){
				for(double y = 0; y<=rangeXY; y+=dXY){
					V p = add(mult(scale, new V(x, y, 0)), pC);
					Double dirA = fun(pO, p, outA);
					
					if(dirA!=null){						
						drawer.pt.setColor(cl);
					
					// draw null aussi
//					}else{
//						drawer.pt.setColor(Color.black);
//					}
//					{
						
						// show direction
						L l = new L(p, new V(C(dirA), S(dirA), 0));
						drawer.drLXY(l, dXY*scale.x/3);
						drawer.drLXY(l, -dXY*scale.y/3);

						
						
						// show function line
//						double xA = p.x+ (90-outA)/90 * dXY*0.4*scale.x;
//						double yA = dirA != null
//								? p.y+(dirA)/90 * dXY*0.4*scale.y
//								: p.y;
//
//						if(new Double(abs(outA)).intValue()%45==0){
//							drawer.drLine(xA, yA+2, xA, yA-2);
//						}
//						drawer.drPoint(xA, yA);
//						
//						drawer.pt.setColor(Color.black);
//						drawer.drPoint(p.x, yA);
						
					}

					drawer.pt.setColor(Color.black);
					drawer.drPoint(p.x, p.y);		
				}
				
			}
		}
		
	}

}



interface Surface{
	public L out(L in, double n);
	public void show(Drawer drawer);
}

class Sph implements Surface{
	V o;
	V df;
	double r;
	double range;
	
	Sph(V o, V df, double r, double a){
		this.o = o;
		this.df = df.unit();
		this.r = r;
		this.range = a;
	}

	@Override
	public L out(L in, double n) {
		V dis = sub(in.o, o);
		double a = multD(in.dir, in.dir);
		double b = multD(in.dir, dis)*2;
		double c = multD(dis, dis) - pow(r, 2);
		
		double[] rt = defun(a, b, c);
		if(rt == null){
			return null;
		}

		for(double l : rt){
			if(l>0){
				V pcr = add(mult(in.dir, l), in.o);
				V dfl = sub(pcr, o).unit();
				double ar = ac(multD(dfl, df));
				if(ar<=range){
					V dir = trm(dfl, in.dir, n);
					return new L(pcr, dir);
				}
			}
		}
		
		return null;
	}

	@Override
	public void show(Drawer drawer) {
		drawer.drCircle(o.x, o.y, r);
	}
}

class SegSurface implements Surface{
	SegCurve curve;
	V center;
	V axis;

	V rollAxis;
	L cvCenterAxis;
	L sfCenterAxis;
	double rolla;
	
	int crossMode = CMODE_OPT;
	
	static final int CMODE_OPT = 0;
	static final int CMODE_GEO = 1;
	
	static final V curveAxis = V.AXIS_Y;
	
	public SegSurface(SegCurve curve, V center, V axis) {
		this.curve = curve;
		this.center = center;
		this.axis = axis;
		
		rollAxis = V.prodVU(axis, curveAxis).unit();
		rolla = ac(V.multDU(axis, curveAxis));
		cvCenterAxis = new L(curve.center, curveAxis);
		sfCenterAxis = new L(center, axis);
		
	}

	//use curve's axis and center
	public L outCv(L in, double n) {
		V crossP = findCvCross(in);
		if(crossP==null){
			return null;
		}

		V f = findCvF(crossP);
		if(f==null){
			return new L(crossP, null);
		}
		
		if( (in.dir.y<0 && V.a(in.dir, f)<PI/2) || (in.dir.y>=0 && V.a(in.dir, f)>PI/2) ){
			return new L(crossP, null);
		}
		
		V out = trm(f, in.dir, n);

		if(out!=null && out.y*in.dir.y<=0){
			return new L(crossP, null);
		}
		
		if(out!=null && out.y<0){
			double distPlanEnd = sub(add(crossP, mult(out, (curve.getEndCenter().y-crossP.y)/out.y)), curve.getEndCenter()).absXZ();
			if(distPlanEnd > abs(curve.getEndPoint().x)){
				return new L(crossP, null);
			}
		}
		
		return new L(crossP, out);
	}

	public V findCvCrossReal(L in){
		return toSurfaceAxisPos(findCvCross(toCurveAxisPos(in)));
	}

	private V findSegCross(L in, Integer sg) {
		L segI0 = curve.getPSeg(sg);
		L segI1 = curve.getPSeg(sg+1);

		if(sg==0){
			V p = add(in.o, mult(in.dir, ((segI1.o.y-in.o.y)/in.dir.y)));
			if(U.aAbsSmallEqualb(p.absXZ(), (segI1.o.x))){
				return p;
			}else{
				return null;
			}
		}
		
		double tg = segI0.dir.y/segI0.dir.x;
		double yc = segI0.o.y+(curve.center.x-segI0.o.x)*tg;
		V c = new V(curve.center.x, yc, curve.center.z);
		
		V cToO = sub(in.o, c);
		
		double A = U.pow2(tg) * V.multD2((in.dir.projXZ()))-U.pow2(in.dir.y);
		double B = 2*( U.pow2(tg)*multD(cToO.projXZ(), in.dir.projXZ()) -in.dir.y*cToO.y );
		double C = U.pow2(tg)* V.multD2(cToO.projXZ()) - U.pow2(cToO.y);
		
		double[] ls = U.defun(A, B, C);
		if(ls==null){
			return null;
		}
		for(double l : ls){
			if(l<0){
				continue;
			}

			V p = add(in.o, mult(in.dir, l));
			if(segI1==null){
				System.out.println(sg);
			}
			if( U.aInBC(p.y, segI0.o.y, segI1.o.y)){
				double xz = sub(p, curve.center).absXZ();
				double r0 = abs(segI0.o.x);
				double r1 = abs(segI1.o.x);
				if( U.aInBC(xz, r0, r1 )){
					return p;
				}else{
//					System.out.println("..");
				}
			}else{
//				Drawer.crtDrawer.pt.setColor(Color.red);
//				System.out.println(p.y+",   ("+segI0.o.y+", "+segI1.o.y+")");
//				Drawer.crtDrawer.fiRect(p.x, p.y, 4, 4);
			}
		}
		
		return null;
	}

	public void setCrossMode(int mode){
		this.crossMode = mode;
	}
	
	public Integer[] findCvDrossBetween(L in, int sgT, int sgB){
		boolean fromTop = in.dir.y<0;
		double distTop = distSegO(in, sgT);
		double distBottom = distSegO(in, sgB);		
		boolean inTop = U.aSmallEqualb(distTop, 0);
		boolean inBottom = U.aSmallEqualb(distBottom, 0);

		if(!fromTop && crossMode==CMODE_OPT && !inBottom){
			return null;
		}
		
		
		if(inTop && inBottom){
			return new Integer[]{sgT-1};
		}else if(inTop && !inBottom){
			if(fromTop){
				return new Integer[]{sgT-1};
			}else{
				if(sgT+1==sgB){
					return new Integer[]{sgT};
				}
				
				int sgM = (sgT+sgB)/2;
				double distMi = distSegO(in, sgM);		
				boolean inMi = U.aSmallEqualb(distMi, 0);
				if(inMi){
					return findCvDrossBetween(in, sgM, sgB);
				}else{
					return findCvDrossBetween(in, sgT, sgM);
				}
			}
		}else if(!inTop && inBottom){
			if(sgT+1==sgB){
				return new Integer[]{sgT};
			}
			
			int sgM = (sgT+sgB)/2;
			double distMi = distSegO(in, sgM);		
			boolean inMi = U.aSmallEqualb(distMi, 0);
			if(inMi){
				return findCvDrossBetween(in, sgT, sgM);
			}else{
				return findCvDrossBetween(in, sgM, sgB);
			}
		}else if(!inTop && !inBottom){
			if(sgT+1==sgB){
				return new Integer[]{sgT};
			}
			
			int sgM = (sgT+sgB)/2;
			double distMi = distSegO(in, sgM);		
			boolean inMi = U.aSmallEqualb(distMi, 0);
			if(inMi){
				return fromTop? findCvDrossBetween(in, sgT, sgM) : findCvDrossBetween(in, sgM, sgB);
			}
			
			if(sgT+2==sgB){
				return new Integer[]{sgM, sgM+1};
			}
			
			int sgM2 = sgM+1;
			double distMi2 = distSegO(in, sgM2);		
			boolean inMi2 = U.aSmallEqualb(distMi2, 0);

			if(inMi2){
				return fromTop? new Integer[]{sgM} : findCvDrossBetween(in, sgM2, sgB);
			}else{
				if(U.aSmallb(distMi, distMi2)){
					return findCvDrossBetween(in, sgT, sgM2);
				}else if(U.aBigb(distMi, distMi2)){
					return findCvDrossBetween(in, sgM, sgB);
				}else {
					return new Integer[]{sgM};
				}
				
			}

		}
		
		return null;
	}
	
	// for segI>0, in.dir.y!=0
	public double distSegO(L in, int sgI){
		L segI = curve.getPSeg(sgI);
		V pL = add(in.o, mult(in.dir, ((segI.o.y-in.o.y)/in.dir.y)));
		return pL.absXZ() - abs(segI.o.x);
	}
	
	// new
	public V findCvCross(L in){
		Integer[] possibleSg = findCvDrossBetween(in, 1, curve.dv);
		if(possibleSg==null){
			return null;
		}else if(possibleSg.length==1){
			return findSegCross(in, possibleSg[0]);
		}else{
			if(in.dir.y<0){
				for(int i=0; i<possibleSg.length; i++){
					if(possibleSg[i]==null){
						System.out.println("length="+possibleSg.length);
					}
					V res = possibleSg[i]<curve.dv? findSegCross(in, possibleSg[i]) : null;
					if(res != null){
						return res;
					}
				}
			}else if(in.dir.y>0){
				for(int i=possibleSg.length-1; i>=0; i--){
					V res = possibleSg[i]<curve.dv? findSegCross(in, possibleSg[i]) : null;
					if(res != null){
						return res;
					}
				}
			}
		}
		
		return null;
	}
	

	public V findCvCrossB(L in){
		int i = 0;
		int iBorder = curve.dv;
		int di = 1;
		if(in.dir.y>0){
			i = curve.dv-1;
			iBorder = -1;
			di = -1;
		}
		
		// find seg
		Integer findI = null;
		L segI0 = null;
		L segI1 = null;
		for(; di*i<di*iBorder; i+=di){
			segI0 = curve.getPSeg(i);
			segI1 = curve.getPSeg(i+1);
			
			// in.dir.y == 0
			if(U.is0(in.dir.y)){
				if(U.aInBC(in.o.y, segI0.o.y, segI1.o.y)){
					findI = i;
					break;
				}else{
					findI = null;
					segI0 = null;
					segI1 = null;
					continue;
				}
			}
			
			double dist0 = V.dist(cvCenterAxis, add(in.o, mult(in.dir, (segI0.o.y-in.o.y)/in.dir.y)));
			double dist1 = V.dist(cvCenterAxis, add(in.o, mult(in.dir, (segI1.o.y-in.o.y)/in.dir.y)));
			double segX0 = abs(segI0.o.x -curve.center.x);
			double segX1 = abs(segI1.o.x -curve.center.x);

//			if((dist0<segX0) != (dist1<=segX1)){
//			if(U.aInBC(0, dist0-segX0, dist1-segX1)){
			if((U.aBigEqualb(dist0, segX0) && U.aSmallEqualb(dist1,segX1)) || (U.aSmallEqualb(dist0, segX0) && U.aBigEqualb(dist1,segX1))){
				findI = i;
//				System.out.println((dist0-segX0)+", "+(dist1-segX1));
				break;
			}else if(dist0<segX0 && dist1<segX1){
				findI = null;
				segI0 = null;
				segI1 = null;
				continue;
			}else {				
				if(segI0.dir.y!=0){
					V segI00 = add(segI0.o, mult(segI0.dir, -curve.r/curve.dv/segI0.dir.x));
					V segI11 = add(segI0.o, mult(segI0.dir, 2*curve.r/curve.dv/segI0.dir.x));
					double dist00 = V.dist(cvCenterAxis, add(in.o, mult(in.dir, (segI00.y-in.o.y)/in.dir.y)));
					double dist11 = V.dist(cvCenterAxis, add(in.o, mult(in.dir, (segI11.y-in.o.y)/in.dir.y)));

					if(dist00-segI00.x>dist0-segX0 && dist11-segI11.x>dist1-segX1){
						findI = i;
						break;
					}
					
				}
				
			}
		}
		
		if(findI==null){
//			System.out.println("i not found for " + in);
			return null;
		}else{
//			System.out.println("i fount = "+findI);
		}
		
		double tg = segI0.dir.y/segI0.dir.x;
		double yc = segI0.o.y+(curve.center.x-segI0.o.x)*tg;
		V c = new V(curve.center.x, yc, curve.center.z);
		
		V cToO = sub(in.o, c);
		
		double A = U.pow2(tg) * V.multD2((in.dir.projXZ()))-U.pow2(in.dir.y);
		double B = 2*( U.pow2(tg)*multD(cToO.projXZ(), in.dir.projXZ()) -in.dir.y*cToO.y );
		double C = U.pow2(tg)* V.multD2(cToO.projXZ()) - U.pow2(cToO.y);
		
		double[] ls = U.defun(A, B, C);
		if(ls==null){
//			System.out.println("ls null for" + in);
//			System.out.println("A="+A+",B="+B+",C="+C);
//			System.out.println("D="+(U.pow2(B)-4*A*C));
//			System.out.println("ls==null");
			return null;
		}
//		System.out.println("\n"+findI);
		for(double l : ls){
			if(l<0){
//				System.out.println("l:"+l);
				continue;
			}
//			System.out.println(".");
			V p = add(in.o, mult(in.dir, l));
//			System.out.println(p.y+"\n"+segI0.o.y+"\n"+segI1.o.y);
//			if( (p.y<=segI0.o.y && p.y>=segI1.o.y) || ( U.is0(segI0.dir.y) && U.is0(p.y-segI0.o.y) ) ){
			if( U.aInBC(p.y, segI0.o.y, segI1.o.y)){
				double xz = sub(p, curve.center).absXZ();
				double r0 = abs(segI0.o.x);
				double r1 = abs(segI1.o.x);
//				if( (xz>=r0 && xz<=r1 )                  || U.is0(xz-r1) ){
				if( U.aInBC(xz, r0, r1 )){
					return p;
				}else{
//					System.out.println("..");
				}
			}else{
//				Drawer.crtDrawer.pt.setColor(Color.red);
//				System.out.println(p.y+",   ("+segI0.o.y+", "+segI1.o.y+")");
//				Drawer.crtDrawer.fiRect(p.x, p.y, 4, 4);
			}
		}
//		System.out.println("...");
		return null;
	}
	
	public V findCvF(V p){
		V cP = sub(p, curve.center);
//		int iXY = U.toI45(cP.absXZ()/curve.r*curve.dv);
		int iXY = U.toI(cP.absXZ()/curve.r*curve.dv);
		L seg = curve.getPSeg(iXY);
		if(seg==null){
			return null;
		}

		V f = new V(-seg.dir.y, seg.dir.x, seg.dir.z);
		if(f.y<0){
			f = V.N(f);
		}
		
		if(!U.is0(cP.x) || !U.is0(cP.z)){
			V rollTo = new V(cP.x, 0, cP.z).unit();
			f = new V(f.x*rollTo.x, f.y, f.x*rollTo.z);
		}

		return f;
		
	}

	public V toCurveAxis(V v){
		if(v==null){
			return null;
		}
		return V.roll(rollAxis, v, rolla);
	}
	
	public V toSurfaceAxis(V v){
		if(v==null){
			return null;
		}
		return V.roll(rollAxis, v, -rolla);
	}
	
	public L toCurveAxis(L l){
		if(l==null){
			return null;
		}
		return V.roll(rollAxis, l, rolla);
	}
	
	public L toSurfaceAxis(L l){
		if(l==null){
			return null;
		}
		return V.roll(rollAxis, l, -rolla);
	}

	public V toSurfaceAxisPos(V v){
		if(v==null){
			return null;
		}
		return add(sfCenterAxis.o, toSurfaceAxis(v));
	}
	
	public V toCurveAxisPos(V v){
		if(v==null){
			return null;
		}
		return V.roll(rollAxis, sub(v, sfCenterAxis.o), rolla);
	}

	public L toSurfaceAxisPos(L l){
		if(l==null){
			return null;
		}
		L l2 = V.roll(rollAxis, l, -rolla);
		l2.o = add(l2.o, sfCenterAxis.o);
		return l2;
	}

	public L toCurveAxisPos(L l){
		if(l==null){
			return null;
		}
		l = new L(sub(l.o, sfCenterAxis.o), l.dir);
		return V.roll(rollAxis, l, rolla);
	}
	

	
	@Override
	public L out(L in, double n) {
		in = toCurveAxisPos(in);
		
		L out = outFor0(in, n);
		if(out==null){
			return null;
		}
		
		return toSurfaceAxisPos(out);
	}
	
	public L outFor0(L in, double n) {
		if(V.dist(in, curve.getEndCenter())>curve.getMaxDistToEC()){
			return null;
		}

		L out = outCv(in, n);
		if(out==null){
			return null;
		}
		
		return out;
	}

	@Override
	public void show(Drawer drawer) {
		// TODO Auto-generated method stub
		
	}
}

class SfcCllMtx{
	SegCurve curve;
	IPDirRec ipDirRec;	
	int pxCllR;

	int pxCllW;		// cal
	double dPx;		// cal
	double dCllW;	// cal
	
	int clObsDist;
	int pxScnR;
	double rangeA;

	int clRMtx; // cal
	double pxScnH; // cal (to obs)
	
	V obsP; //cal
	
	int clTgtDist;	
	double clTgtRso;
	
	SegSurface[][] clls;
	V[][] scnColors;
	
	Color mtxBgColor = Color.black;
	Color cllBgColor = Color.gray;
	
	public SfcCllMtx(SegCurve curve, IPDirRec ipDirRec, Integer pxCllR){
		this.curve = curve;
		this.ipDirRec = ipDirRec;
		ipDirRec.setCurve(curve);
		
		dPx = curve.l/ipDirRec.rIp;
		
		pxCllR = pxCllR!=null ? pxCllR : ipDirRec.rIp;
		pxCllW = pxCllR*2+1;
		dCllW = pxCllW*dPx;
	}

	public void initObs(int clObsDist, int pxScnR, double rangeA){
		this.clObsDist = clObsDist;
		this.pxScnR = pxScnR;
		this.rangeA = rangeA;
		this.clRMtx = U.toI45(clObsDist*T(rangeA));
		this.pxScnH = pxScnR/T(rangeA);
		this.obsP = new V(0, 0, clObsDist*dCllW);
		
		initClls();
	}
	
	public void initClls() {
		clls = new SegSurface[clRMtx*2+1][clRMtx*2+1];
		for(int x=-clRMtx; x<=clRMtx; x++){
			for(int y=-clRMtx; y<=clRMtx; y++){
				SegSurface surface = new SegSurface(curve, new V(x*dCllW, y*dCllW, 0), V.AXIS_Z);
				setSurface(x, y, surface);
			}
		}
		
	}
	
	public SegSurface getSurface(int i, int j){
		return clls[i+clRMtx][j+clRMtx];
	}
	
	public void setSurface(int i, int j, SegSurface surface){
		clls[i+clRMtx][j+clRMtx] = surface;
	}
	
	public void calMtxScnColors(int pxStep){
		scnColors = new V[pxScnR*2+1][pxScnR*2+1];
		for(int x=-clRMtx; x<=clRMtx; x++){
			for(int y=-clRMtx; y<=clRMtx; y++){
				calCllScnColors(x, y, pxStep);
			}

		}
	}

	public void calCllScnColors2(int sfX, int sfY, int pxStep) {
		for(int y=0; true; y+= pxStep){
			int countValidPx = 0;
			// area 1
			for(int x=0; true; x+= pxStep){
				boolean res = calPxColor(sfX, sfY, x, y);
				if(res){
					countValidPx++;
				}else{
					break;
				}
			}
			
			// area 4
			for(int x=-1; true; x-= pxStep){
				boolean res = calPxColor(sfX, sfY, x, y);
				if(res){
					countValidPx++;
				}else{
					break;
				}
			}

			if(countValidPx==0){
				break;
			}
		}

		for(int y=-1; true; y-= pxStep){
			int countValidPx = 0;
			// area 2
			for(int x=0; true; x+= pxStep){
				boolean res = calPxColor(sfX, sfY, x, y);
				if(res){
					countValidPx++;
				}else{
					break;
				}
			}
			
			// area 4
			for(int x=-1; true; x-= pxStep){
				boolean res = calPxColor(sfX, sfY, x, y);
				if(res){
					countValidPx++;
				}else{
					break;
				}
			}

			if(countValidPx==0){
				break;
			}
		}

	}

	public void calCllScnColors(int sfX, int sfY, int pxStep) {
		for(int x=0; true; x+= pxStep){
			int countValidPx = 0;
			// area 1
			for(int y=0; true; y+= pxStep){
				boolean res = calPxColor(sfX, sfY, x, y);
				if(res){
					countValidPx++;
				}else{
					break;
				}
			}
			
			// area 4
			for(int y=-pxStep; true; y-= pxStep){
				boolean res = calPxColor(sfX, sfY, x, y);
				if(res){
					countValidPx++;
				}else{
					break;
				}
			}

			if(countValidPx==0){
				break;
			}
		}

		for(int x=-pxStep; true; x-= pxStep){
			int countValidPx = 0;
			// area 2
			for(int y=0; true; y+= pxStep){
				boolean res = calPxColor(sfX, sfY, x, y);
				if(res){
					countValidPx++;
				}else{
					break;
				}
			}
			
			// area 4
			for(int y=-pxStep; true; y-= pxStep){
				boolean res = calPxColor(sfX, sfY, x, y);
				if(res){
					countValidPx++;
				}else{
					break;
				}
			}

			if(countValidPx==0){
				break;
			}
		}

	}

	public boolean calPxColor(int sfX, int sfY, int pxX, int pxY) {
		SegSurface surface = getSurface(sfX, sfY);
		double endH = abs(curve.getEndPoint().y);
		V p0 = add(surface.center, mult(new V(pxX, pxY, 0), dPx));
		L in = new L(obsP, sub(p0, obsP));
		V pE = add(obsP, mult(in.dir, (endH+obsP.z)/in.dir.z));

		// out of rangeA
		V scnPx = mult(in.dir, -pxScnH/in.dir.z);
		int scnX = U.toI45(scnPx.x);
		int scnY = U.toI45(scnPx.y);
		
		if(abs(scnX)>pxScnR || abs(scnY)>pxScnR){
			return false;
		}
		

		
		int[][] range = getIndexRangeByP1P2(p0, pE);
		if(hasInterrupt(sfX, sfY, in, range)){
			return false;
		}
		
		L out = surface.out(in, surface.curve.n);
		if(out==null){
			return false;
		}
		
		Color c = null;
		if(out.dir != null){
			V imP = add(out.dir, mult(out.dir, (-surface.curve.h-out.o.z)/out.dir.z));
			c = getImColor(imP);
		}else{
			c = cllBgColor;
		}
		
		addColorCount(scnX, scnY, c);
		
		return true;
	}

	public void addColorCount(int x, int y, Color cl){
		x += pxScnR;
		y += pxScnR;
		
		V colorXY = scnColors[x][y];
		if(colorXY==null){
			colorXY = new V(0, 0, 0);
			scnColors[x][y] = colorXY;
		}
		
		scnColors[x][y] = add(colorXY, new V(cl.getRed(), cl.getGreen(), cl.getBlue()));
	}
	
	public Color getImColor(V imP) {
		int x = U.toI45(imP.x/dPx);
		int y = U.toI45(imP.y/dPx);
		
		return Color.green;
	}

	private boolean hasInterrupt(int sfX, int sfY, L in, int[][] range) {
		for(int i=range[0][0]; abs(i)<=abs(range[0][1]); i+=range[0][2]){
			for(int j=range[1][0]; abs(j)<=abs(range[1][1]); j+=range[1][2]){
				if(i==sfX && j==sfY){
					return false;
				}
					
				V cross = getSurface(i, j).findCvCrossReal(in);
				if(cross != null && abs(i)<=abs(sfX) && abs(j)<=abs(sfY)){
					return true;
				}
			}			
		}
		
		return false;
	}
	
	public int[] getIndexByP(V p){
		return new int[]{
				U.toI45(p.x/dCllW),
				U.toI45(p.y/dCllW),
		};
	}
	
	public int[][] getIndexRangeByP1P2(V p1, V p2){
		int xP1 = U.toI45(p1.x/dCllW);
		int yP1 = U.toI45(p1.y/dCllW);

		int xP2 = U.toI45(p2.x/dCllW);
		int yP2 = U.toI45(p2.y/dCllW);

		return new int[][]{
			getMinMaxD(xP1, xP2),
			getMinMaxD(yP1, yP2)
		};
		
	}

	private int[] getMinMaxD(int i1, int i2) {
		boolean order12 = abs(i1)<=abs(i2);
		int iMin = order12? i1 : i2;
		int iMax = order12? i2 : i1;
		int d = iMax>=iMin?1 : -1;
		return new int[]{iMin, iMax, d};
	}
	

	
	
	
}


// z=0;
abstract class PlView implements Surface{
	V o;
	V df;
	double step;

	@Override
	public L out(L in, double n) {
		return null;
	}

	@Override
	public void show(Drawer drawer) {
	}

	abstract public L out(L in, V c, double n);
}

class PlEmpty extends PlView{
	public PlEmpty(V o, V df, double step){
		this.o = o;
		this.df = df.unit();
		this.step = step;
	}

	public L out(L in, V c, double n) {
		L out = new L(c, new V((in.dir.x/in.dir.z)*in.dir.z, (in.dir.y/in.dir.z)*in.dir.z, in.dir.z));
		return out;
	}
}

class PlDisc extends PlView{
	int resolution;
	int ua; 
	
	static int resoUnit = 1000;
	static double resoUnitD = new Double(resoUnit);

	public PlDisc(V o, V df, double step, int resolution){
		this.o = o;
		this.df = df.unit();
		this.step = step;
		this.resolution = resolution;
		this.ua = 180*resoUnit / resolution;
	}

	public double getDiscretA(double A){
		return (new Double( A*resoUnit + (A>0 ? 0.5 : -0.5)*ua).intValue())/ua*ua/resoUnitD;
	}
	
	public L out(L in, V c, double n) {
		double AXOut = getDiscretA(AT(in.dir.x/in.dir.z));
		double AYOut = getDiscretA(AT(in.dir.y/in.dir.z));
		
		L out = new L(c, new V(T(AXOut)*in.dir.z, T(AYOut)*in.dir.z, in.dir.z));
		return out;
	}
}

class PlDiscN extends PlView{
	int resolution;
	int ua; 
	
	static int resoUnit = 1000;
	static double resoUnitD = new Double(resoUnit);

	public PlDiscN(V o, V df, double step, int resolution){
		this.o = o;
		this.df = df.unit();
		this.step = step;
		this.resolution = resolution;
		this.ua = 180*resoUnit / resolution;
	}

	public double getDiscretA(double A){
		return (new Double( A*resoUnit + (A>0 ? 0.5 : -0.5)*ua).intValue())/ua*ua/resoUnitD;
	}
	
	public L out(L in, V c, double n) {
		double AXOut = getDiscretA(AT(in.dir.x/in.dir.z));
		double AYOut = getDiscretA(AT(in.dir.y/in.dir.z));
		
		L out = new L(c, new V(T(AXOut)*in.dir.z, T(AYOut)*in.dir.z, in.dir.z));
		return out;
	}
}

class Ovl implements Surface{
//	x�/a� + y�/b� = 1
//  y=b(1-x�/a�)�
//  y'=-b�x / a�y
//  faxian xielv = a�y/b�x
	 
	V o;
	V df;
	double a;
	double b;
	double range;
	
	Ovl(V o, V df, double a, double b, double range){
		this.o = o;
		this.df = df.unit();
		this.a = a;
		this.b = b;
		this.range = range;
	}

	@Override
	public L out(L in, double n) {
		V dis = sub(in.o, o);
		V ratio = new V(1/a/a, 1/b/b, 1/a/a);
		
		double da = multD(in.dir, mult(in.dir, ratio));
		double db = multD(dis, mult(in.dir,ratio))*2;
		double dc = multD(dis, mult(dis,ratio)) - 1;
		
//		e�/a�  +  f�/b�
//		2((p-r)e/a�+(q-s)f/b�)
//		(p-r)�/a�  +  (q-s)�/b� -1
		
		double[] rt = defun(da, db, dc);
		if(rt == null){
			return null;
		}

		for(double l : rt){
			if(l>0){
				V pcr = add(mult(in.dir, l), in.o);
				V dfp = sub(pcr, o).unit();
				V xyz = sub(pcr, o);
			//  faxian xielv = a�y/b�x
				V fp = new V(xyz.x*b*b, xyz.y*a*a, 0);
				
				double ar = ac(multD(df, dfp));
				if(ar<=range){
					V dir = trm(fp, in.dir, n);
					return new L(pcr, dir);
				}
			}
		}
		
		return null;
	}

	@Override
	public void show(Drawer drawer) {
		drawer.drOval(o.x, o.y, a, b);
	}


}

class SegCurve implements Serializable{
	private static final long serialVersionUID = "SegCurve".hashCode();

	double r = 0.5;
	double l = r*2;
	double h = r*1.9;	

	double n=1.5;

	double am = 0;
	int dv = 800;

	String initMethod;
	
	int hueNum;
	int hueAccu;
	Integer[] validLs = null;
	
	V center;

	Double maxDistToEC = null;
	
	L[] pArray;

	static boolean modeMiddle = true;

	
	public String toString(){
		return r+"_"+l+"_"+h+"_"+n+"_"+am+"_"+dv+"__"+center.x+"_"+center.y+"_"+center.z+"_"+getInitMethod();
	}

	public String getStdName(){
		return (l/r)+"_"+(h/r)+"_"+n+"_"+am+"_"+dv+"__"+getInitMethod();
		
	}

	public void init() {
		boolean inited = false;
		if(initMethod==null){
			setInitMethod("initPArray");
		}
		
		try {
			Method m = this.getClass().getMethod(initMethod);
			if(m !=null){
				m.invoke(this);
				inited = true;
			}
		} catch (Exception e) {}
		
		if(!inited){
			initPArray();
		}
	}

	public SegCurve(double r, double l, double h, double n){
		this.r = r;
		this.l = l;
		this.h = h;
		this.n = n;
	}
	
	public void setHue(int hueNum, int hueAccu){
		setHue(hueNum, hueAccu, (Integer[])null);
	}
	
	public void setHue(int hueNum, int hueAccu, Integer ... validLs){
		this.hueAccu = hueAccu;
		this.hueNum = hueNum;
		this.validLs = validLs;
	}

	public void setInitMethod(String method){
		this.initMethod = method;
	}
	
	public String getInitMethod(){
		return initMethod==null ? "initPArray" : initMethod;
	}

	
	public Color getColor(double dl, float trans){
		int dlNb = new Double((100*l + dl)*hueAccu/l).intValue();
//		int dlRealNb = dlNb-100*hueAccu;
		int dlRealNb = dlNb%hueAccu;
		if(validLs!=null){
			boolean isValid = false;
			for(int validL : validLs){
				if(validL==dlRealNb){
					isValid = true;
					break;
				}
			}
			
			if(!isValid){
				return Color.white;
			}
		}
		
		float hue = (abs(dlNb)%hueNum)/(float)hueNum;
		
		Color color = Color.getHSBColor(hue, 1, 1);
		return new Color(((float)color.getRed())/255, ((float)color.getGreen())/255, ((float)color.getBlue())/255, trans);
	}
	
	public void setRange(double am, int dv){
		this.am = am;
		this.dv = dv;
	}
	
	public void setCenter(V center){
		this.center = center;
	}
	
	public V getImgCenter(){
		return add(center, new V(0, -h, 0));
	}
	
	public V getEndPoint(){
		return getPSeg(dv).o;
	}
	
	public V getEndCenter(){
		return add(center, new V(0, getEndPoint().y ,0));
	}
	
	public double getMaxDistToEC(){
		if(maxDistToEC == null){
			V p = getEndCenter();
			double maxD = 0;
			
			for(int i=0; i<=dv; i++){
				maxD = Math.max(V.dist(getPSeg(i).o, p), maxD);
			}
			maxDistToEC = maxD;
		}
		
		return maxDistToEC;
	}
	
	public void initPArray_simple(){
		pArray = new L[dv+1];
		
		pArray[0] = new L(new V(0, 0, 0), new V(1, 0, 0));
		
		double dx = r/dv;
		for(int i=1; i<=dv; i++){
			double tg0 = pArray[i-1].dir.y / pArray[i-1].dir.x;
			double x = pArray[i-1].o.x+dx;
			double y = pArray[i-1].o.y+dx*tg0;
			
			double a1 = AT(x*(1+l/r)/(y+h));
			double a2 = am*x/r;
			
			V p = new V(x, y, pArray[i-1].o.z);
			double tg1 = -(n*S(a1)-S(a2))/(n*C(a1)-C(a2));
			pArray[i] = new L(p, new V(1, tg1, 0));
		}

	}
	
	public static double findAverageH(int c, double n){
		P1P2AOutToADir fun = new P1P2AOutToADir(1/n);

		Double selectedH = null;
		Double selectedAvgSumA = null;
				
		for(double y=0; y<5*c; y+=0.1){
			double sumA = 0;
			double count = 0;
			V p1 = new V(0, y, 0);
			for(double i=0; i<=c; i++){
				V pO = new V(-i, 0, 0);
				Double dirA = fun.fun(pO, p1, 90-90*i/c);
				if(dirA != null){
					sumA += abs(dirA);
					count++;
				}
			}
			sumA /= count;
			if(selectedAvgSumA==null || selectedAvgSumA>sumA){
				selectedH = y;
				selectedAvgSumA = sumA;
			}
		}
		
		return selectedH/c;
	}
	
	public void initPArray_average(){
		pArray = new L[dv+1];

		P1P2AOutToADir fun = new P1P2AOutToADir(1/this.n);
		
		pArray[0] = new L(new V(0, 0, 0), new V(1, 0, 0));
		
		double dx = r/dv;
		for(int i=1; i<=dv; i++){
			double tg0 = pArray[i-1].dir.y / pArray[i-1].dir.x;
			double x = pArray[i-1].o.x+dx;
			double y = pArray[i-1].o.y+dx*tg0;
			
			V p = new V(x, y, pArray[i-1].o.z);
			
			double tg1 = averageForI(fun, i, p);
			pArray[i] = new L(p, new V(1, tg1, 0));

		}

	}

	private double averageForI(P1P2AOutToADir fun, int i, V p) {
		double avgA = 0;
		double count = 0;
		for(double d=-dv; d<=dv; d++){
			V pO = new V(d*l/dv, -h, 0);
			
			if(i==1)System.out.println(d+"::"+(90+90*d/dv));
			
			Double dirA = fun.fun(pO, p, 90+90*d/dv);
			if(dirA != null){
				avgA += dirA;
				count++;
			}
		}
		avgA /= count;
		double tg1 = T(avgA);
		return tg1;
	}
	
	public void initPArray(){
		
		pArray = new L[dv+1];
		
		pArray[0] = new L(new V(0, 0, 0), new V(1, 0, 0));

		// failed attempt : combine init and init average A
//		P1P2AOutToADir fun = new P1P2AOutToADir(1/this.n);

		double dx = r/dv;
//		double swichP = 1; 		// failed attempt : combine init and init average B
		for(int i=1; i<=dv; i++){			
			double tg0 = pArray[i-1].dir.y / pArray[i-1].dir.x;
			double x = pArray[i-1].o.x+dx;
			double y = pArray[i-1].o.y+dx*tg0;
			V p = new V(x, y, pArray[i-1].o.z);
			
			// failed attempt : combine init and init average C fin
//			if(i>dv*swichP){
//				double tg1 = averageForI(fun, i, p);
//				pArray[i] = new L(p, new V(1, tg1, 0));
//				continue;
//			}
			
			V pLeft = new V(-p.x, p.y, p.z);
			
			// get tg1
//			double a1 = AT(x*(1+l/r)/(y+h));
//			double a2 = am*x/r;
//			double tg1 = -(n*S(a1)-S(a2))/(n*C(a1)-C(a2));

//////////////////////////////////
			
			// get tg1 by stat
			double tg1 = 0;
			
			Double aSelected = null;
			Double bestSD = null;
			
			double[][] SDDatas = new double[dv*2+2][2];
			for(double a=-85; a<=85; a+=0.25){
				double tga = T(a);
				L segL = new L(pLeft,  new V(-1, tga, 0));
				L setR = new L(p,  new V(1, tga, 0));
				
				V fL = new V(tga, 1, 0);
				V fR = new V(tga, -1, 0);

				//				getSD for a
				int countValieOut = 0;
				for(int d=0; d<=dv; d++){
//					double aOut = 90-d*am/dv; // average
//					double aOut = 90-d*90/dv; // average 90
					double aOut = 90 - AS(min(1, abs(n*S(AT(d*l/dv/h))))); // trm value in O
					
					
					
					V pO = new V(-d*l/dv, -h, 0);
					
					V outL = null;
					V outR = null;
					if(modeMiddle){
						outL = trm(fL, sub(add(pLeft, mult(segL.dir, -0.5*dx/segL.dir.x)), pO), 1/n);
						outR = trm(fR, sub(add(p, mult(setR.dir, 0.5*dx/setR.dir.x)), pO), 1/n);
					}else{
						outL = trm(fL, sub(pLeft, pO), 1/n);
						outR = trm(fR, sub(p, pO), 1/n);
					}
					

					
//					System.out.println("i="+i+", a="+a+", d="+d+", aOut="+aOut+", left=("+sub(pLeft, pO).AXY()+","+(outL==null?"null":outL.AXY())+"), right=("+sub(p, pO).AXY()+","+(outR==null?"null":outR.AXY())+")");
					
					double di = 2-(abs(i-d)*1d/dv);
					double indiceDi = 0;
					if(outL!=null && outL.y>=0){
						SDDatas[d*2][0] = outL.AXY()-aOut;
						SDDatas[d*2][1] = abs(V.multDU(fL, outL)) * pow(di, indiceDi);
						countValieOut++;
					}else{
						SDDatas[d*2][0] = 0;
						SDDatas[d*2][1] = 0;
					}
					
					if(outR!=null && outR.y>=0){
						SDDatas[d*2+1][0] = outR.AXY()-aOut;
						SDDatas[d*2+1][1] = abs(V.multDU(fR, outR)) * pow(di, indiceDi);
						countValieOut++;
					}else{
						SDDatas[d*2+1][0] = 0;
						SDDatas[d*2+1][1] = 0;
					}
				}
				
				
				
				
				// selection strategie 1
//				Double sda = U.SD(SDDatas, 0);
//				if(countValieOut<(dv+1)*2*0.2)continue;				
////				System.out.println("i="+i+",a="+a+", countValieOut="+countValieOut+", sda="+sda);
//				if((bestSD==null || bestSD>sda)){
//					bestSD = sda;
//					aSelected = a;
//				}
				

				
				
				// selection strategie 2
//				Double sda = U.SD_B(SDDatas, am*0.5/dv);
				Double sda = U.SD_B(SDDatas, am*1/dv);
//				Double sda = U.SD_B(SDDatas, am*2/dv);
//				if(countValieOut<(dv+1)*2*0.2)continue;				
//				System.out.println("i="+i+",a="+a+", countValieOut="+countValieOut+", sda="+sda);
				if((bestSD==null || (sda!=null && bestSD<sda))){
					bestSD = sda;
					aSelected = a;
				}

			}

			System.out.println("i = " +i +", aSelected = " +aSelected +", bestSD = " +bestSD);
			tg1 = T(aSelected);

			
			
//////////////////////////////////
			
			
			pArray[i] = new L(p, new V(1, tg1, 0));
		}

	}

	public L getPSeg(int i){
		if(i>dv || i<0){
			return null;
		}
		
		if(i>=0){
			return new L(add(pArray[i].o, center), pArray[i].dir);
		}else{
			return null;
		}
	}

	public L getNSeg(int i){
		if(i>dv || i<0){
			return null;
		}
		
		if(i>=0){
			L l = new L(add(pArray[i].o, center), pArray[i].dir);
			l.o.x = -l.o.x;
			l.dir.x = -l.dir.x;
			return l;
		}else{
			return null;
		}
	}

	public V getPstP(int i){
		if(i>dv || i<0){
			return null;
		}
		
		if(modeMiddle){
			return add(center, mult(add(pArray[i].o, pArray[i+1].o), 0.5));  // new 1
//			return mult(add(getPSeg(i).o, getPSeg(i+1).o), 0.5);  //new 2
		}else{
			return getPSeg(i).o; // old
		}
	}
	
	public V getPstN(int i){
		if(i>dv || i<0){
			return null;
		}
		
		return mult(getPstP(i), new V(-1, 1, 1));
	}
	
}

class L implements Serializable{
	private static final long serialVersionUID = "L".hashCode();

	V o;
	V dir;
	
	L(V o, V dir){
		this.o = o;
		if(dir!=null){
			this.dir = dir.unit();
		}else{
			this.dir = null;
		}
	}
	
	public static int CROSS_NO = 0;
	public static int CROSS_BEFORE = 1;
	public static int CROSS_AFTER = 2;
	
	public String toString(){
		return "L:o"+o +", dir:"+dir;
	}
	
	public L projXY(){
		return new L(o.projXY(), dir.projXY());
	}
	
	public L projYZ(){
		return new L(o.projYZ(), dir.projYZ());
	}
	
	public L projXZ(){
		return new L(o.projXZ(), dir.projXZ());
	}

	
	public int cross(LSeg s){
		if(s.c.equals(o)){
			return CROSS_NO;
		}
		
		double a = dir.x*(s.p1.y-o.y) - dir.y*(s.p1.x-o.x);
		double b = dir.x*(s.p2.y-o.y) - dir.y*(s.p2.x-o.x);
		if(a*b>=0){
			return CROSS_NO;
		}
		
		double R = (o.y-s.p1.y) * s.lx - (o.x-s.p1.x)*s.ly;
		return (s.ly*dir.x - s.lx*dir.y)*R>0 ? CROSS_AFTER : CROSS_BEFORE;
	}
	
	public boolean crossBefore(List<LSeg> pl){
		for(LSeg s : pl){
			if(CROSS_BEFORE == cross(s)){
				return true;
			}
		}
		return false;
	}

	public boolean crossAfter(List<LSeg> pl){
		for(LSeg s : pl){
			if(CROSS_AFTER == cross(s)){
				return true;
			}
		}
		return false;
	}

}

class LSeg{
	V p1; 
	V p2;
	V dir;	
	double len;
	double lx;
	double ly;
	
	V c;
	V dirN;
	
	LSeg(V p1, V p2){
		this.p1 = new V(p1);
		this.p2 = new V(p2);
		dir = sub(p2, p1).unit();

		len = sub(p2, p1).abs();
		lx = this.p2.x - this.p1.x;
		ly = this.p2.y - this.p1.y;
		
		calC();
		calNXY();
	}
	
	LSeg(V p1, V dir, double len){
		this.p1 = new V(p1);
		this.dir = dir.unit();
		this.len = len;
		this.p2 = add(p1, mult(this.dir, this.len));

		lx = this.p2.x - this.p1.x;
		ly = this.p2.y - this.p1.y;

		calC();
		calNXY();
	}
	
	V calC(){
		c = mult(add(p1, p2), 0.5);
		return c;
	}
	
	V calNXY(){
		dirN = new V(-dir.y, dir.x, dir.z);
		return dirN;
	}
}

class PlGenerator{
	int segNum;
	double segLen;
	double aStep;
	double limitX;
	double limitY;
	double limitA;
	
	Double[] crtAs;	
	LSeg lsBegin;

	public PlGenerator(int segNum, double segLen, double aStep, double limitX, double limitY, double limitA){
		this.segNum = segNum;
		this.segLen = segLen;
		this.aStep = aStep;
		this.limitX = limitX;
		this.limitY = limitY;
		this.limitA = limitA;
		
		crtAs = new Double[segNum];
		lsBegin = new LSeg(new V(-segLen/2, 0, 0), new V(1, 0, 0), segLen);
	}
	
	public List<LSeg> getNextPl(){
		List<LSeg> pl = new ArrayList<LSeg>();
		pl.add(lsBegin);
		
		int iLast = 0;
		for(; iLast<segNum; iLast++){
			if(crtAs[iLast]==null){
				break;
			}else{
				pl.add(new LSeg(getLast(pl).p2, getDir(crtAs[iLast]), segLen));
			}
		}
		iLast = max(iLast -1, 0);
		
		int tour = 0;
		while(tour++<segNum*3){
			double preA = iLast==0? 0 : crtAs[iLast-1];
			
			// new
			if(crtAs[iLast] == null){
				crtAs[iLast] = max(-limitA, preA-60);
//				crtAs[iLast] = -limitA;
				pl.add(new LSeg(getLast(pl).p2, getDir(crtAs[iLast]), segLen));
			}else{

				double a = crtAs[iLast]+aStep;
				if(a>limitA || abs(preA-a)>60){
					// <--
					crtAs[iLast] = null;
					removeLast(pl);
					if(iLast == 0){
						return null;
					}
					iLast--;
					continue;
				}else{
					setLast(pl, new LSeg(getLast(pl).p1, getDir(a), segLen));
					crtAs[iLast] = a;
				}
			}

			if(isEndable(pl)){
				return pl;
			}else{
				iLast++;
			}
		}
	
		System.out.println("error");
		return null;
	}

	public List<LSeg> getNextSymPl(){
		List<LSeg> pl = getNextPl();
		if(pl == null){
			return null;
		}
		int l = pl.size();
		V sym = new V(-1, 1, 1);
		for(int i=1; i<l; i++){
			LSeg si = pl.get(i);
			pl.add(new LSeg(mult(si.p1, sym), mult(si.p2, sym)));
		}
		return pl;
	}

	
	public V getDir(double a){
		return new V(C(a), S(a), 0).unit();
	}
	
	public boolean isEndable(List<LSeg> poly){
		if(poly.size()>segNum){
			return true;
		}else{
			LSeg last = getLast(poly);
			return poly.size()>= segNum || abs(last.p2.x)>limitX  || abs(last.p2.y)>limitY ;
		}
	}
	
	public LSeg getLast(List<LSeg> poly){
		return poly.get(poly.size()-1);
	}

	public void removeLast(List<LSeg> poly){
		poly.remove(poly.size()-1);
	}

	public void setLast(List<LSeg> poly, LSeg lSeg){
		poly.set(poly.size()-1, lSeg);
	}
}

class Group{
	int r; // dir 1
	float z1; // dir 3

	int d; // dir 1
	int md; // dir 1

	int l;
	int n;
	int h;

	float zmd;
	int t;
	double a;
	float H;
	float hl;
	Gra g0;
	Gra gmd;
	Gra gl;
	Gra gn;
	
	public Image img;
	int objR;
	V obj;

	public Group(int r, float z1, int d, int md, int l, int n, int h){
		this.r = r;
		this.z1 = z1;
		this.d = d;
		this.md = md;
		this.l = l;
		this.n = n;
		this.h = h;

		zmd = z1*md/d;
		t = l*n+h;
		a = AT((d-1)/2/z1);

		g0 = new Gra(r, d*r, (l*n+h)*z1);
		gmd = new Gra(md*r, d*r, (l*n+h)*z1-zmd);
		gl = new Gra(r, l*r, (l*n-l+h)*z1);
		gn = new Gra(r, l*r, h*z1);
		
		H = (l*n+h)*z1;
		hl = l*z1;
	}
	
	public void initImage(int width, int height, int objDis, int oR, Color cO, Color cB) {
		img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		Graphics g = img.getGraphics();
		Drawer dr = new Drawer(width, height, g);
		
		V objC = new V(0,0,objDis * H);

		objR = oR * r;
		obj = new V(0,0,H-objC.z);

		
		int cx = width/r/d; 
		int cy = height/r/d; 
		int cd = d/2;
//		List<Integer> ss = new ArrayList<Integer>();
		for(int i=-cx; i<=cx; i++){
			for(int j=-cy; j<=cy; j++){
				V ctij = new V(i*r*d, j*r*d, 0);
//				if(abs(i)>0 || abs(j)>0)continue;
//				System.out.println(objC.x+","+objC.y+","+objC.z+","+objR);
				for(int p=-cd; p<=cd; p++){
					for(int q=-cd; q<=cd; q++){
						V dirijpq = new V(p*t*r,q*t*r,H);
						V imijpq = add(ctij, dirijpq);
//						System.out.println(p+","+q+"  "+dirijpq.x+", "+dirijpq.y);
						
						V rej = sub(add(ctij, mult(dirijpq, objC.z/imijpq.z)), objC);
						if(abs(rej.x)<objR && abs(rej.y)<objR){
							g.setColor(cO);
							dr.flRect(imijpq.x, imijpq.y, r, r);
						}else{
							g.setColor(cB);
							dr.flRect(imijpq.x, imijpq.y, r, r);
						}
						
//						int sss = (i*d+p*t)*100000+(j*d+q*t);
////							ss.add(sss);
//							ss.add(new Double(imijpq.x).intValue()*10000000+new Double(imijpq.y).intValue());
					}
				}
			}
		}
//		System.out.println(ss.size());
//		Collections.sort(ss);
//		for(int i=0; i<ss.size()-1; i++){
////			if(ss.get(i)==ss.get(i+1)){
//				System.out.println(ss.get(i));
////			}
//		}
		
	}

	public Area getArea(int cx, int cy, int lx, int ly, V e){
		System.out.println("H="+H+",a="+a+",R="+r*d);
		System.out.println("E="+e.x+", "+e.z+", S="+obj.z+", a="+AT(e.x/e.z));
		System.out.println("S="+obj.z);
		System.out.println("e="+e.z/r/d+", s="+obj.z/r/d);

		Area a = g0.getArea(cx, cy, lx, ly, e);
//		a.intersect(gmd.getArea(cx, cy, lx, ly, e));
		
		int b = new Float(d*(e.z-gl.h)/l/(e.z-g0.h)).intValue();
		a.intersect(gl.getArea(cx, cy, lx*b, ly*b, e));
		
		b = new Float(d*(e.z-gn.h)/l/(e.z-g0.h)).intValue();
		a.intersect(gn.getArea(cx, cy, lx*b, ly*b, e));
		return a;
	}
}

class Gra{
	int r;
	int d;
	float h;
	
	public Gra(	int r, int d, float h){
		this.r = r;
		this.d = d;
		this.h = h;
	}

	public Area getArea(int cx, int cy, int lx, int ly, V e){
		
		float c = new Float(e.z/(e.z-h));
		float c2 = new Float(h/(e.z-h));
		Float rf = r*c;
		int ri = rf.intValue();
		float demiRf = (r-1)*c/2;
//		float demiRf = (r)*c/2;
		float rm = d*c;
		float ox = cx + new Float(-e.x*c2);
		float oy = cy + new Float(-e.y*c2);
		
		float x0 = ox-rm*(lx+1);
		float wx = 2*rm*(lx+1);
		
		float y0 = oy-rm*(ly+1);
		float hy = 2*rm*(ly+1);
		
		Area area = new Area();
		for(int i=-lx; i<=lx; i++){
//			area.add(new Area(new Rectangle2D.Float(ox+i*rm-demiRf, y0, rf, hy)));
			area.add(new Area(new Rectangle2D.Float(new Float(ox+i*rm-demiRf+0.5).intValue(), y0, ri, hy)));
		}

		Area area2 = new Area();
		for(int j=-ly; j<=ly; j++){
//			area2.add(new Area(new Rectangle2D.Float(x0, oy+j*rm-demiRf, wx, rf)));
			area2.add(new Area(new Rectangle2D.Float(x0, new Float(oy+j*rm-demiRf+0.5).intValue(), wx, ri)));
		}
		
		area.intersect(area2);
		return area;
	}	
}

class IPDirRec{
	public int rDir;
	public int rIp;
	public int rMnIp;
	public int rCvDiv;
	public double rangeA;
	
	private double hDir;
	private String objName;
	SegCurve curve;

	VList[][] ipDirs;
	DirCounts dirCounts;
	Integer[][] dirIPs;

	public IPDirRec(int rDir, int rIp, int rMnIp, int rCvDiv, double rangeA) {
		super();
		this.rDir = rDir; 		// 500
		this.rIp = rIp;			// 50
		this.rMnIp = rMnIp;
		this.rCvDiv = rCvDiv;	// 10
		this.rangeA = rangeA;	// 60
		
		hDir = rDir/T(rangeA);
	}

	public double getHDir(){
		return hDir;
	}

	public void setCurve(SegCurve curve){
		this.curve = curve;
		objName = calObjName();

	}
	
	public String calObjName(){
		return rDir+"_"+rIp+"_"+rMnIp+"_"+rCvDiv+"_"+rangeA+"__"+curve.getStdName();
	}

	public String ipDirsObjName(){
		objName = calObjName();
		return("IpDirs_"+objName);
	}

	public String dirIpsObjName(){
		objName = calObjName();
		return("DirIPs_"+objName);
	}

	public String dirCountsObjName(){
		objName = calObjName();
		return("DirCounts_"+objName);
	}

	public void rec(){
		calDatas();
		objName = calObjName();
		
		ObjSaver.saveObj(dirCounts, dirCountsObjName());
		ObjSaver.saveObj(ipDirs, ipDirsObjName());
		ObjSaver.saveObj(dirIPs, dirIpsObjName());
	}

	public VList[][] getIpDirs() {
		ipDirs = (VList[][])ObjSaver.readObjWhenNeed(VList[][].class, ipDirs, ipDirsObjName());
		if(ipDirs == null){
			rec();
		}
		
		return ipDirs;
	}

	public DirCounts getDirCounts() {
		dirCounts = (DirCounts)ObjSaver.readObjWhenNeed(DirCounts.class, dirCounts, dirCountsObjName());
		if(dirCounts == null){
			rec();
		}
		
		return dirCounts;
	}

	public Integer[][] getDirIPs() {
		dirIPs = (Integer[][])ObjSaver.readObjWhenNeed(Integer[][].class, dirIPs, dirIpsObjName());
		if(dirIPs == null){
			rec();
		}

		return dirIPs;
	}

	public void setDirIPs(Integer[][] dirIPs) {
		this.dirIPs = dirIPs;
	}

	public void calDatas(){
		SegSurface surface = new SegSurface(curve, new V(0, 0, 0), V.AXIS_Y);
		
		int cellR = rMnIp * 2+1;
		int cellBR = cellR * rIp;
		
		V imgCenter = curve.getImgCenter();
		
//		double scale = curve.l*1.5/cellBR; // 1.5 for cell voisin 
		double scale = curve.l/cellBR;  
		double scaleCvr = curve.r/rCvDiv;

		dirCounts = (DirCounts)ObjSaver.readObj(DirCounts.class, dirCountsObjName());
		if(dirCounts==null){
			System.out.println("dirCounts not found");
			dirCounts = new DirCounts(rDir, rIp);
			
			for(int i=0; i<=rIp; i++){
				System.out.println("rec : i = "+i);
				for(int j=0; j<=rIp; j++){
					// a IP
					V centerIJ = new V(i*cellR, 0, j*cellR);
					
					// sec 0.5 in cicle
					if(centerIJ.abs()>cellBR || i<j){
						continue;
					}

					if(j%5==0){
						System.out.println("rec : \tj = "+j);
					}
					
					for(int ic = -rMnIp; ic<=rMnIp; ic++){
						for(int jc = -rMnIp; jc<=rMnIp; jc++){
							// a mn ip
							V pImg = add(mult(add(centerIJ, new V(ic, 0, jc)), scale), imgCenter);
							
							for(double cvi=-rCvDiv; cvi<=rCvDiv; cvi++){
								for(double cvj=-rCvDiv; cvj<=rCvDiv; cvj++){							
									Double rXY = new Double(pow(cvi*cvi+cvj*cvj, 0.5));
									if(rXY>rCvDiv)continue;

									// get in dir
									V in = sub(add(mult(new V(cvi, 0, cvj), scaleCvr*1), curve.getEndCenter()), pImg).unit();
//									L outLine = surface.out(new L(pImg, in), 1/curve.n);
									L outLine = surface.outFor0(new L(pImg, in), 1/curve.n);
									V out = outLine==null? null : outLine.dir;		

									if(out==null){
										continue;
									}

									V pXZ = add(curve.center, mult(out, hDir/out.y));
									
//									dirCounts.addValue(U.toI(pXZ.x), U.toI(pXZ.z), i, j, 1);
									dirCounts.addValue(U.toI45(pXZ.x), U.toI45(pXZ.z), i, j, multD(in, V.AXIS_Y));
//									dirCounts.addValueOct(U.toI45(pXZ.x), U.toI45(pXZ.z), i, j, multD(in, V.AXIS_Y));
									
								}
							}
						}
					}
					
					
				}
			}
		}
		
		dirIPs = new Integer[dirCounts.dirR*2+1][dirCounts.dirR*2+1];
		ipDirs = new VList[rIp*2+1][rIp*2+1];

		IJtoK ipIJK = new IJtoK(rIp);
		for(int i= -dirCounts.dirR; i<= dirCounts.dirR; i++){
			for(int j= -dirCounts.dirR; j<= dirCounts.dirR; j++){
				
				Integer maxV = null;
				double maxCt = 1;
//				double totalCount = 0;
				double advance = 1;
				int[][] dirOcts = Oct.getOctAll(i, j);
				for(int octI=0; octI<8; octI++){
					IPCounts ipCounts = dirCounts.getValue(dirOcts[octI][0], dirOcts[octI][1]);
					
					if(ipCounts != null && ipCounts.ipMap != null){
						
						for(Map.Entry<Integer, Double> entry : ipCounts.ipMap.entrySet()){
							Integer ip = entry.getKey();
							int[] ipOct = Oct.toIDis(ipCounts.ijk.getI(ip), ipCounts.ijk.getJ(ip),-octI);
							if(ipOct != null){
								Double d = entry.getValue();
//								totalCount += d;
								if(d>maxCt){
									advance = d/maxCt;
									maxCt = d;
									maxV = ipCounts.ijk.getK(ipOct[0], ipOct[1]);										
								}else{
									advance = min(advance, maxCt/d);
								}
							}
						}

					}
				}
				

//					maxCt*3>=totalCount ||
//				if(maxV != null && ( advance>1.2 )){
				if(maxV != null){
					dirIPs[i+dirCounts.dirR][j+dirCounts.dirR] = maxV;	
					
					int ipI = ipIJK.getI(maxV);
					int ipJ = ipIJK.getJ(maxV);
					
					VList ipDirList = ipDirs[ipI+rIp][ipJ+rIp];
					if(ipDirList==null){
						ipDirList = new VList(ipI, ipJ);
						ipDirs[ipI+rIp][ipJ+rIp] = ipDirList;
					}
					ipDirList.dirs.add(new V(i, j, getHDir()));
					
				}else{
					dirIPs[i+dirCounts.dirR][j+dirCounts.dirR] = null;				
				}
			}
		}
		
	}
	
	public void dealFarDirs() {
		for(VList[] ipdirArr : ipDirs){
			for(VList vList : ipdirArr){
				if(vList != null){
					vList.deleteFarValues(this);
				}
			}
		}
	}
	
	public void dealNullDirs() {
		Integer[][] dirIpsNew = new Integer[rDir*2+1][rDir*2+1];
		IJtoK ipIJK = new IJtoK(rIp);

		for(int i=0; i<= rDir*2; i++){
			for(int j=0; j<= rDir*2; j++){
				if(true){

					Map<Integer, Integer> ipCounts = new HashMap<Integer, Integer>();
					int nullCount = 0;
					
					int paveR = 5;
					for(int p=max(0, i-paveR); p<=min(rDir*2, i+paveR); p++){
						for(int q=max(0, j-paveR); q<=min(rDir*2, j+paveR); q++){
							Integer vpq = dirIPs[p][q];
							if(vpq==null){
								nullCount++;
							}else{
								if(ipCounts.containsKey(vpq)){
									ipCounts.put(vpq, ipCounts.get(vpq)+1);
								}else{
									ipCounts.put(vpq, 1);
								}
							}							
						}
					}
					int maxCount = 0;
					Integer maxIp = null;
					for(Map.Entry<Integer, Integer> entry : ipCounts.entrySet()){
						if(entry.getValue()>maxCount){
							maxCount = entry.getValue();
							maxIp = entry.getKey();
						}
					}
					
					if(maxCount<nullCount){
						maxIp = null;
					}
					
					dirIpsNew[i][j] = maxIp;

					if((maxIp!=null && !maxIp.equals(dirIPs[i][j])) || maxIp!=dirIPs[i][j]){
						if(dirIPs[i][j] != null){
							int ipI = ipIJK.getI(dirIPs[i][j]);
							int ipJ = ipIJK.getJ(dirIPs[i][j]);
							VList ipDirList = ipDirs[ipI+rIp][ipJ+rIp];
							if(ipDirList != null){
								ipDirList.dirs.remove(new V(i-rDir, j-rDir, getHDir()));
							}
						}
						
						if(maxIp!=null){
							int ipI = ipIJK.getI(maxIp);
							int ipJ = ipIJK.getJ(maxIp);
							
							VList ipDirList = ipDirs[ipI+rIp][ipJ+rIp];
							if(ipDirList==null){
								ipDirList = new VList(ipI, ipJ);
								ipDirs[ipI+rIp][ipJ+rIp] = ipDirList;
							}
							ipDirList.dirs.add(new V(i-rDir, j-rDir, getHDir()));
						}
	
					}
					
					
				}else{
					dirIpsNew[i][j] = dirIPs[i][j];
				}
			}
		}
		
		dirIPs = dirIpsNew;
	}

}

class DirCounts implements Serializable{
	private static final long serialVersionUID = 3006051378467273177L;

	IPCounts[][] dirsIPs;
	int dirR;
	int ipR;

	IJtoK ijk;

	public DirCounts(int dirR, int ipR){
		this.dirR = dirR;
		this.ipR = ipR;
		dirsIPs = new IPCounts[dirR*2+1][dirR*2+1];
		ijk = new IJtoK(ipR);
	}

	public IPCounts getValue(int i, int j) {
		return dirsIPs[i+dirR][j+dirR];
	}

	public IPCounts getValue(double x, double y) {
		return getValue(U.toI45(x), U.toI45(y));
	}

	public Double getValue(int i, int j, int x, int y) {
		IPCounts ipc = dirsIPs[i+dirR][j+dirR];
		if(ipc==null || ipc.ipMap==null){
			return null;
		}
		return ipc.getValue(x, y);
	}

	public Double getValue(double i, double j, double x, double y) {
		return getValue(U.toI45(i), U.toI45(j), U.toI45(x), U.toI45(y));
	}

	public void addValue(double dirX, double dirY, double ipX, double ipY, double value){
		addValue(U.toI45(dirX), U.toI45(dirY), U.toI45(ipX), U.toI45(ipY), value);
	}

	public void addValue(int dirX, int dirY, int ipX, int ipY, double value){
		int dirI = dirX + dirR;
		int dirJ = dirY + dirR;
		
		IPCounts currentIPCounts = null;
		try{
			currentIPCounts = dirsIPs[dirI][dirJ];
		}catch(ArrayIndexOutOfBoundsException e){
			return;
		}
		
		if(currentIPCounts==null){
			currentIPCounts = new IPCounts(ijk);
			dirsIPs[dirI][dirJ] = currentIPCounts;
		}
		
		currentIPCounts.addValue(ipX, ipY, value);
	}

	public void addValueOct(int dirX, int dirY, int ipX, int ipY, double value) {
		int[][] dirOcts = Oct.getOctAll(dirX, dirY);
		int[][] ipOcts = Oct.getOctDis(ipX, ipY);
		
		for(int i=0; i<8; i++){
			if(ipOcts[i]!=null){
				addValue(dirOcts[i][0], dirOcts[i][1], ipOcts[i][0], ipOcts[i][1], value);
				
			}
		}
	}


}

class VList implements Serializable{
	private static final long serialVersionUID = 1095558864259654260L;

	int i;
	int j;
	List<V> dirs;
	
	V avgUnit = null;
	Integer countV = null; 
	Double aAvg = null;
	
	public VList(int i, int j){
		this.i = i;
		this.j = j;
		dirs = new ArrayList<V>();
	}
	
	public void calAvgUnit(){
		if(dirs==null && dirs.isEmpty()){
			avgUnit = null;
			countV = null;
		}
		
		V avg = new V(0, 0, 0);
		int count = 0;
		for(V dir : dirs){
			if(dir!=null){
				avg = add(avg, dir.unit());
				count++;
			}
		}
		countV = count;
		avgUnit =  count==0? null : V.mult(avg, 1d/count).unit();
	}
	
	public void calaAvg(){
		double aDist = 0;
		if(countV!=null && countV!=0 && avgUnit != null){
			for(V dir : dirs){
				if(dir!=null){
					aDist += abs(ac(V.multD(dir.unit(), avgUnit)));
				}
			}
			aAvg = aDist/countV;
		}
	}
	
	public void deleteFarValues(IPDirRec rec){
		int countDel = 0;
		calAvgUnit();
		calaAvg();
		do{
			countDel = 0;
			if(countV!=null && countV>0){
				for(int i=0; i<dirs.size();){
					V dir = dirs.get(i);
					
					if(dir!=null && abs(ac(V.multD(dir.unit(), avgUnit)))>aAvg*2){
						dirs.remove(i);
						if(rec != null){
							rec.dirIPs[U.toI45(dir.x)+rec.rDir][U.toI45(dir.y)+rec.rDir] = null;
						}
						countDel++;
					}else{
						i++;
					}
				}
			}
			calAvgUnit();
			calaAvg();
		}while(countDel != 0);
	}
}

class IPCounts implements Serializable{
	private static final long serialVersionUID = -2619409562435070912L;

	Map<Integer, Double> ipMap;
	IJtoK ijk;
	public IPCounts(IJtoK ijk){
		ipMap = new HashMap<Integer, Double>();
		this.ijk = ijk;
	}

	public Double getValue(int x, int y) {
		if(ipMap==null){
			return null;
		}
		Integer ipK = ijk.getK(x, y);
		return ipMap.get(ipK);
	}

	public void addValue(int x, int y, double value) {
		Integer ipK = ijk.getK(x, y);
		if(ipMap.containsKey(ipK)){
			ipMap.put(ipK, ipMap.get(ipK)+value);
		}else{
			ipMap.put(ipK, value);
		}
	}
}

class Oct{
	public static int[][] getOctAll(int x, int y){
		return getOct(x, y, true);
	}
	
	public static int[][] getOctDis(int x, int y){
		return getOct(x, y, false);
	}

	public static int[][] getOct(int x, int y, boolean canRep) {
		int[][] octs = new int[8][];
		octs[0] = new int[]{x, y};
		if(canRep || (x!=0 || y!=0)){
			octs[2] = new int[]{-y, x};
			octs[4] = new int[]{-x, -y};
			octs[6] = new int[]{y, -x};
		}

		if(canRep || (abs(x)!=abs(y) && y!=0)){
			octs[1] = new int[]{y, x};
			octs[3] = new int[]{-x, y};
			octs[5] = new int[]{-y, -x};
			octs[7] = new int[]{x, -y};
		}
		
		return octs;
	}
	
	public static int[] toI(int x, int y, int index, boolean canRep){
		if(index==0){
			return new int[]{x, y}; 
		}

		if(canRep || (x!=0 || y!=0)){
			if(index==2 || index==-6){
				return new int[]{-y, x};
			}else if(index==4 || index==-4){
				return new int[]{-x, -y};
			}else if(index==6 || index==-2){
				return new int[]{y, -x};
			}
		}

		if(canRep || (abs(x)!=abs(y) && y!=0)){
			if(index==1 || index==-1){
				return new int[]{y, x};
			}else if(index==3 || index==-3){
				return new int[]{-x, y};
			}else if(index==5 || index==-5){
				return new int[]{-y, -x};
			}else if(index==7 || index==-7){
				return new int[]{x, -y};
			}
		}
		
		return null;
	}

	public static int[] toIDis(int x, int y, int index){
		return toI(x, y, index, false);
	}
	
	public static int[] toIAll(int x, int y, int index){
		return toI(x, y, index, true);
	}

}

class MapArrayAccess<T>{
	Map<Integer, T>[][] datas;
	IJtoK ijToK;

	public MapArrayAccess(Map<Integer, T>[][] datas, int ri, int rj){
		this.datas = datas;
		ijToK = new IJtoK(ri, rj);
	}

	public MapArrayAccess(Map<Integer, T>[][] datas, int r){
		this.datas = datas;
		ijToK = new IJtoK(r);
	}
	
	public void setData(int arrX, int arrY, int ri, int rj, T val){
		if(datas[arrX][arrY]==null){
			datas[arrX][arrY] = new HashMap<Integer, T>();
		}
		datas[arrX][arrY].put(ijToK.getK(ri, rj), val);
	}
}

class IJtoK implements Serializable{
	private static final long serialVersionUID = 8086376179961220082L;

	int ri;
	int Ri;
	
	int rj;
	int Rj;
	
	int dist;
	
	public IJtoK(int r){
		this(r, r);
	}

	public IJtoK(int ri, int rj){
		this.ri = ri;
		this.rj = rj;

		this.Ri = ri*2+1;
		this.Rj = rj*2+1;

		this.dist = ri*Rj + rj;
	}
	
	public int getK(int i, int j){
		return i*Rj +j;
	}
	
	public int getI(int k){
		return (k+dist)/Rj - ri;
	}

	public int getJ(int k){
		return (k+dist)%Rj - rj;
	}

}

class Sun{
	SimpleDateFormat DF = new SimpleDateFormat("yyyyMMdd-HHmmss");
	double obliquity = U.toa(23+4d/60);
	double daysOfYear = 365d;
	public Date getEquinoxN(){
		try {
//			return DF.parse("20180320-161500");
//			return DF.parse("20180105-000000");
			return DF.parse("20180405-175800");
		} catch (ParseException e) {
			return null;
		}
	}
	
	public static V plToTude(double i, double j, Date tm) {
		double l = Math.sqrt(i*i+j*j);
		
		if(l>2){
			return null;
		}

		double a = U.ac(-j/l);
		if(i<0){
			a=2*PI-a;
		}
		
		if(useDD){
			double d = dd;
			l = U.at(l*d/dm)*2/d;		
			if(l>range){
				return null;
			}
			
		}
		
		return new V((a/2/PI+0.0)%1, 1-l, 0);
	}

	public double getSuna(Date date){
		double diffa = PI*2*(getEquinoxN().getTime()-date.getTime())/(3600000*24)/daysOfYear;
		V axEq = new V(sin(obliquity), 0, cos(obliquity));
		V axAct = new V(axEq.x*cos(diffa), axEq.x*sin(diffa), axEq.z);
		
		return V.a(axAct, new V(0,1,0));
	}
	
	public V getSunP(double sana, double b){
		V x = new V(1, 0, 0);
		V y = new V(0, sin(sana), -cos(sana));
		V z = new V(0, cos(sana), sin(sana));
		V vb = new V(cos(b), 0, sin(b));
		return new V(multD(x, vb), multD(y, vb), multD(z, vb)).unit();
	}

//	public static double dd = 2.435;
//	public static double dd = 3;
	public static double dd = 1.5;
	public static double range = 0.970;
	public static double dm = dd/(U.t(range*dd/2));
	public static boolean useDD = true;

	public V spToPl(V v){
		double l = V.a(v, new V(0,0,1))/PI*2;
		double mXY = mXY(v);
		if(!U.is0(mXY)){
			
			if(useDD){
				if(l>range*4){
					return null;
				}
				double d = dd;
				l = U.t(l*d/4d)*2/d * dm;
			}

			return new V(l/mXY*v.x, l/mXY*v.y, 0);
		}else{
			return null;
		}
	}

}
