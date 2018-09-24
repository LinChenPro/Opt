package my.tests;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.awt.Canvas;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.JFrame;

public class Phy extends Drawer {
	static int dt = 5;
	static double dts = dt/1000d;
	
	public Phy(int w, int h, Graphics g) {
		super(w, h, g);
	}
	static int w = 1260;
	static int h = 670;

	static Canvas cvs;

	
	double scaleM = 1000;

	double yFlr = -300/scaleM;
	double hShl = 400/scaleM;
	
	double rShl = hShl/2;
	double yShl = yFlr + rShl;
	double yCore = yFlr + rShl;
	double minYShl = yShl;
	double minYCore = yFlr;
	
	double vShl = 0;
	double vCore = -3000/scaleM;

	double kEls = 0.5;
	double g=-9.8;
	
	double mShl =  0.01;
	double mCore = 0.1;
	
	
	public void show(){
		// cal y;
		yCore = max(minYCore, yCore + vCore * dts);
		yShl = min(yCore+rShl, (max(yCore-rShl, max(minYShl, yShl+ vShl * dts))));
//		yShl = max(minYShl, yShl+ vShl * dts);
		
		pt.setColor(Color.gray);
		drY(yFlr * scaleM);
		drCircle(0, yCore * scaleM, 10);
		drCircle(0, yShl * scaleM, rShl*scaleM);
		
		// cal f;
		double fElsT = (rShl + yShl-yCore)*kEls;
		double fElsB = (rShl + yCore-yShl)*kEls;
		double fCore = fElsT - fElsB;
		double fShl = -fCore;
		
		// cal v
		vCore = vCore + (fCore/mCore + g)*dts;
		if(yCore<=yFlr){
//			vCore = max(0, vCore);
			vCore = -vCore;
		}
		
		vShl = vShl + (fShl/mShl + g)*dts;
		if(fShl>0){
			System.out.println(fShl/mShl);
		}
		if(yShl<=minYShl){
			vShl = max(0, vShl);
		}
		if(yShl<=yCore-rShl){
			vShl = max(vCore, vShl);
		}
		if(yShl>=yCore+rShl){
			vShl = min(vCore, vShl);
		}
		
	}
	
	static Phy phy = null;
	public static void main(String[] args) {
		JFrame jf = new JFrame();
		jf.setSize(w+15, h+40);
		jf.setLayout(null);
		jf.setLocation(0, 0);
		Drawer.threadOpen = true;
		cvs = new Canvas(){
			public void paint(Graphics g){
//				g.setColor(Color.white);
//				g.fillRect(0, 0, w, h);
				if(phy==null){
					phy = new Phy(w, h, g);
				}
				phy.pt = (Graphics2D)g;
				phy.show();
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
						sleep(Math.max(Phy.dt-t, 1));
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
