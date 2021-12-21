/*
 * Brandon Li 
 * Microwave simulation for 2210
 * 11/14/2021
 */

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.RenderingHints;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class MicrowaveObstacleSimulation extends JFrame implements Runnable, MouseListener, MouseMotionListener, KeyListener {
	
	//Units: cm, ns
	
	double[][] display;
	double[][] u;
	double[][] up;
	double[][] upp;
	double[][] e;
	double[][] bx;
	double[][] by;
	double[][] eavg;
	double[][] phix;
	double[][] phiy;
	int[][] obstacles;
	int detectornum = 41;
	ArrayList<Obstacle> obstaclelist;
	MeasurementApparatus[] detectors;
	public static double xmin = -8;
	public static double xmax = 122;
	public static double ymin = 0;
	public static double ymax = 50;
	public static double width = xmax-xmin;
	public static double height = ymax-ymin;
	public static double detectorx = 122-6;
	public static double detectordeltay = 5.0;
	public static double ds = 0.1;
	public static double dt = 0.002;
	public static int nx = (int)(width/ds);
	public static int ny = (int)(height/ds);
	double c = 30.0;
	double sigma = c*(dt/ds);
	double sigmasq = sigma*sigma;
	double time = 0.0;
	int measurementsteps = 0;
	int maxmeasurementsteps = 100;
	double avgintensity = 0;
	double intensityfactor = 10.0;
	int iterationmultiplier = 1;

	double obstconfinementwidth = 61;
	double obstconfinementheight = 29;
	double obstconfinementwidth_random = 61*0.8;
	double obstconfinementheight_random = 29*0.8;
	double obstacleradius = 1.25;
	
	boolean paused = true;
	boolean frame = false;
	boolean clear = false;
	boolean reset = false;
	boolean measure = false;
	int brightness = 0;
	int view = 1;
	RenderCanvas r;
	
	boolean experimentended = false;

	public static void main(String[] args) {
		MicrowaveObstacleSimulation w = new MicrowaveObstacleSimulation();
		w.run();
	}
	
	public MicrowaveObstacleSimulation() {
		this.setSize(800,800);
		display = new double[nx][ny];
		u = new double[nx][ny];
		up = new double[nx][ny];
		upp = new double[nx][ny];
		e = new double[nx][ny];
		bx = new double[nx][ny];
		by = new double[nx][ny];
		eavg = new double[nx][ny];
		phix = new double[nx][ny];
		phiy = new double[nx][ny];
		obstacles = new int[nx][ny];
		obstaclelist = new ArrayList<Obstacle>();
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				display[i][j]=0.0;
				u[i][j]=0.0;
				up[i][j]=0.0;
				upp[i][j]=0.0;
				e[i][j]=0.0;
				bx[i][j]=0.0;
				by[i][j]=0.0;
				eavg[i][j]=0.0;
				phix[i][j]=0.0;
				phiy[i][j]=0.0;
				obstacles[i][j] = 1;
			}
		}
		
		detectors = new MeasurementApparatus[detectornum];
		for (int i = 0; i < detectornum; i++) {
			double f = i/((double)detectornum-1.0);
			double y = 5.0*(1-f) + (ymax - 5.0)*f;
			detectors[i] = new MeasurementApparatus((int)((y-detectordeltay/2.0)/ds),(int)((y+detectordeltay/2.0)/ds), i+5);
		}

/** Configurations **/
//		goldenRatio2();
//		SquareGridWithRandomOffset();
		HexGrid();
//		SquareGrid();
//		Vshape();

		this.setVisible(true);
		r= new RenderCanvas(this);
		this.add(r);
		r.addMouseListener(this);
		r.addMouseMotionListener(this);
		this.addKeyListener(this);
		this.setVisible(true);
		this.setTitle("Microwave diffraction");
		r.setPreferredSize(new Dimension(r.imgwidth-1, r.imgheight));
		this.pack();
		for (double x = -10; x < 2; x += 0.2) {
			double y = (2.0-x)/10.0*1.7 +(x+8.0)/10.0*5.0;
			placeobstacle(x, ymax/2.0+y, 0.3, 0, false);
			placeobstacle(x, ymax/2.0-y, 0.3, 0, false);
		}
	}
	
	
	public void SquareGrid() {
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				placeobstacle(xmax/2.0+obstconfinementwidth_random*(i-1.5)/3.0, ymax/2.0+obstconfinementheight_random*(j-1.0)/2.0, obstacleradius, 4.0, true);
			}
		}
	}

	//seeds: 50023 (Random configuration #1)
	// 1522 (Random config #2)
	// 314159 (Random config #3)
	public void SquareGridWithRandomOffset() {
		Random rand = new Random(314159);
		
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				placeobstacle(xmax/2.0+obstconfinementwidth_random*((i-1.5 + ((rand.nextDouble()-0.5)*0.5))/3.0), ymax/2.0+obstconfinementheight_random*((j-1.0+((rand.nextDouble()-0.5)*0.5))/2.0), obstacleradius, 4.0, true);
			}
		}
	}
	
	public void arrangement2() {
		for (int i = 0; i < 15;)
		{
			if (placeobstacle(xmax*(Math.random()*0.6+0.2), ymax*(Math.random()*0.6+0.2), obstacleradius, 4.0, true))
				i++;
		}
	}
	
	public void HexGrid() {
		placegrid(0,0);
		placegrid(0,1);
		placegrid(1,-1);
		placegrid(1,0);
		placegrid(1,1);
		placegrid(2,-1);
		placegrid(2,0);
		placegrid(3,-2);
		placegrid(3,-1);
		placegrid(3,0);
		placegrid(4,-2);
		placegrid(4,-1);
	}
	
	public void goldenRatio() {
		for (int n = 0; n < 3; n++) {
			for (int m = 0; m < 4; m++) {
				double theta = n;
				double r = 10.0*Math.exp(0.30635*theta);
				theta *= 0.4;
				theta += m*Math.PI/2;
				double x = r*Math.cos(theta);
				double y = r*Math.sin(theta);
				placeobstacle(x+xmax/2, y+ymax/2, obstacleradius, 0.0, true);
			}
		}
	}
	

	public void goldenRatio2() {
		for (int n = 0; n < 12; n++) {
				double theta = n;
				double r = 5.0*Math.exp(0.15*theta);
				theta *= 1.1;
				double x = r*Math.cos(theta);
				double y = r*Math.sin(theta);
				placeobstacle(x+xmax/2, y+ymax/2, obstacleradius, 0.0, true);
		}
	}

	public void Vshape() {
		double dx = 0.5;
		double dy = 0.4;
		double dx2 = 1.0;
		placeobstaclenormalized(-1.0, 0.0, obstacleradius, 4.0);
		placeobstaclenormalized(-1.0+dx, dy, obstacleradius, 4.0);
		placeobstaclenormalized(-1.0+2*dx, 2*dy, obstacleradius, 4.0);
		placeobstaclenormalized(-1.0+3*dx, 3*dy, obstacleradius, 4.0);
		placeobstaclenormalized(-1.0, 0.0, obstacleradius, 4.0);
		placeobstaclenormalized(-1.0+dx, -dy, obstacleradius, 4.0);
		placeobstaclenormalized(-1.0+2*dx, -2*dy, obstacleradius, 4.0);
		placeobstaclenormalized(-1.0+3*dx, -3*dy, obstacleradius, 4.0);

		placeobstaclenormalized(dx2-1.0, 0.0, obstacleradius, 4.0);
		placeobstaclenormalized(dx2-1.0+dx, dy, obstacleradius, 4.0);
		placeobstaclenormalized(dx2-1.0+2*dx, 2*dy, obstacleradius, 4.0);
		placeobstaclenormalized(dx2-1.0, 0.0, obstacleradius, 4.0);
		placeobstaclenormalized(dx2-1.0+dx, -dy, obstacleradius, 4.0);
		placeobstaclenormalized(dx2-1.0+2*dx, -2*dy, obstacleradius, 4.0);
	}
	
	public void placegrid(int n1, int n2) {
		double d = 17.5;
		double b1x = d*Math.cos(Math.PI/6.0);
		double b1y = d*Math.sin(Math.PI/6.0);
		double b2x = 0;
		double b2y = d;
		double v0x = 32.5;
		double v0y = 16.25;
		placeobstacle((int) (b1x*n1 + b2x*n2 + v0x), (int) (b1y*n1 + b2y*n2 + v0y), obstacleradius, 0, true);
	}

	public boolean placeobstaclenormalized(double x, double y, double r, double mindist) {
		return placeobstacle(x*obstconfinementwidth/2.0 + xmax/2.0, y*obstconfinementheight/2.0+ymax/2.0, r, mindist, true);
	}
	
	public boolean placeobstacle(double x, double y, double r, double mindist, boolean addtolist) {
		for (Obstacle o: obstaclelist) {
			if ((o.x - x)*(o.x - x) + (o.y - y)*(o.y - y) < mindist*mindist) {
				return false;
			}
		}
		if (addtolist) {
			obstaclelist.add(new Obstacle(x, y, r));
		}
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				double dist = Math.sqrt((i*ds+xmin-x)*(i*ds+xmin-x) + (j*ds+ymin-y)*(j*ds+ymin-y));
				if (dist <= r)
					obstacles[i][j] = 0;
			}
		}
		if (addtolist) {
			//System.out.println("Obstacle #" + obstacle + ": x = " + String.format("%2.2f", x) + ", y = " + String.format("%2.2f", y));
		}
		return true;
	}
	int obstacle = 1;

	public void Physics() {
		
		if (clear || reset) {
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					display[i][j]=0.0;
					u[i][j]=0.0;
					up[i][j]=0.0;
					upp[i][j]=0.0;
					e[i][j]=0.0;
					eavg[i][j]=0.0;
					phix[i][j]=0.0;
					phiy[i][j]=0.0;
					
					if(reset)
						obstacles[i][j]=1;
					//obstacles[i][j] = 1;
				}
			}
			time = 0.0;
			clear = false;
			reset = false;
		}
		

		if (mouseIsPressed) {

			int ip = (int)(mouseX/RenderCanvas.scalefactor);
			int jp = ny - (int)(mouseY/RenderCanvas.scalefactor);
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					double dx = ds*Math.sqrt((i-ip)*(i-ip) + (j-jp)*(j-jp));
					if (dx <= obstacleradius) {
						if (mouseButton == MouseEvent.BUTTON1)
							obstacles[i][j] = 0;
						else if (mouseButton == MouseEvent.BUTTON3)
							obstacles[i][j] = 1;
					}
				}
			}
			//u[mouseX/RenderCanvas.INTERVAL][mouseY/RenderCanvas.INTERVAL] = 10.0;//10.0*Math.sin(time*20.0);
			//10.0*Math.sin(time*20.0);
		}
		
		if (!paused || frame) {
			
			if (!experimentended && time > 15.0) {
				measure = true;
				paused = true;
				experimentended = true;
			}
			
			/*Interior*/
			for (int i = 1; i < nx-1; i++)
			{
				for (int j = 1; j < ny-1; j++)
				{
					u[i][j] = sigmasq*(up[i+1][j] + up[i-1][j] + up[i][j+1] + up[i][j-1] - 4*up[i][j]) + 2*up[i][j] - upp[i][j]; // - 0.1*dt*(up[i][j] - upp[i][j]);
					//if (!newenergy) {
					//phix[i][j] = -(up[i][j] - upp[i][j])*(up[i+1][j]-up[i-1][j])/(2.0*ds*dt);
					//phiy[i][j] = -(up[i][j] - upp[i][j])*(up[i][j+1]-up[i][j-1])/(2.0*ds*dt);
					//e[i][j] = (c*c*((up[i+1][j] - up[i-1][j])*(up[i+1][j] - up[i-1][j]) + (up[i][j+1] - up[i][j-1])*(up[i][j+1] - up[i][j-1]))/(4*ds*ds) + (up[i][j]-upp[i][j])*(up[i][j]-upp[i][j])/(dt*dt));
					//}
					bx[i][j] += -dt/(2*ds)*(up[i][j+1] - up[i][j-1]);
					by[i][j] += dt/(2*ds)*(up[i+1][j] - up[i-1][j]);
					phix[i][j] = -u[i][j]*by[i][j];
					phiy[i][j] = u[i][j]*bx[i][j];
					e[i][j] = c*c*u[i][j]*u[i][j] + (by[i][j]*by[i][j] + by[i][j]*by[i][j]);
					eavg[i][j] = 0.99*eavg[i][j] + 0.01*e[i][j];
				}
			}

			/*Boundary*/
			for (int j = 1; j < ny-1; j++) {
				u[0][j] = 2.0*up[0][j]-upp[0][j]+sigma*(up[0+1][j]-up[0][j]-upp[0+1][j]+upp[0][j]) + (sigmasq/2.0)*(up[0][j+1]-2*up[0][j]+up[0][j-1]);
			}

			for (int j = 1; j < ny-1; j++) {
				u[nx-1][j] = 2.0*up[nx-1][j]-upp[nx-1][j]-sigma*(up[nx-1][j]-up[nx-2][j]-upp[nx-1][j]+upp[nx-2][j]) + (sigmasq/2.0)*(up[nx-1][j+1]-2*up[nx-1][j]+up[nx-1][j-1]);
			}

			int ctr = ny/2;
			int halfwidth = (int)(0.66/ds);
			//halfwidth = (int)(24.9/ds);
			for (int j = ctr-halfwidth; j <= ctr+halfwidth; j++) {
				u[2][j] += dt*3.0*1500*(Math.sin(2*Math.PI*10.5*time))*Math.atan(4*time*time);
				//u[2][j] = Math.atan(4*time*time)*200*(Math.sin(2*Math.PI*1.5*time));
				//u[2][j] = 20.0*(Math.random()-0.5);
			}

			measurementsteps++;
			if (measurementsteps == maxmeasurementsteps) {
				avgintensity = 0;
				for (int i = 0; i < detectornum; i++) {
					detectors[i].measure(this, true);
					avgintensity += detectors[i].avgintensity;
				}
				avgintensity /= detectornum;

				measurementsteps = 0;
			} else {
				for (int i = 0; i < detectornum; i++) {
					detectors[i].measure(this, false);
				}
			}
			if (measure) {
				for (int i = 0; i < detectornum; i++) {
					System.out.println(detectors[i].number + "\t" + (detectors[i].avgintensity*intensityfactor));
				}
				measure=false;
			}

			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					u[i][j] *= obstacles[i][j];
					upp[i][j] = up[i][j];
					up[i][j] = u[i][j];
					display[i][j] = u[i][j];
				}
			}

			time += dt;
			frame = false;
		}
	}
	
	@Override
	public void run() {
		while (true) {
			try {
				Thread.sleep(17);
			} catch (InterruptedException e) {}
			for (int i = 0; i < iterationmultiplier ; i++)
				this.Physics();
			r.repaint();
		}
	}

	int mouseX = 0;
	int mouseY = 250;
	int mouseButton = 0;

	@Override
	public void mouseClicked(MouseEvent arg0) {
		//updateMousePosition = !updateMousePosition;
	}
	boolean mouseIsPressed = false;
	PointerInfo p = MouseInfo.getPointerInfo();
	@Override
	public void mouseEntered(MouseEvent arg0) {}
	@Override
	public void mouseExited(MouseEvent arg0) {}
	@Override
	public void mousePressed(MouseEvent arg0) {mouseIsPressed = true;
	mouseButton = arg0.getButton();
	}
	@Override
	public void mouseReleased(MouseEvent arg0) {mouseIsPressed = false;}

	@Override
	public void mouseDragged(MouseEvent arg0) {
			mouseX = arg0.getX();
			mouseY = arg0.getY();
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
			mouseX = arg0.getX();
			mouseY = arg0.getY();
	}

	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void keyPressed(KeyEvent e) {
		if (e.getKeyCode() == KeyEvent.VK_P) {
			paused = !paused;
		} else if (e.getKeyCode() == KeyEvent.VK_F) {
			frame = true;
		} else if (e.getKeyCode() == KeyEvent.VK_1) {
			view = 1;
		} else if (e.getKeyCode() == KeyEvent.VK_2) {
			view = 2;
		} else if (e.getKeyCode() == KeyEvent.VK_3) {
			view = 3;
		} else if (e.getKeyCode() == KeyEvent.VK_4) {
			view = 4;
		} else if (e.getKeyCode() == KeyEvent.VK_C) {
			clear = true;
		} else if (e.getKeyCode() == KeyEvent.VK_R) {
			reset = true;
		} else if (e.getKeyCode() == KeyEvent.VK_M) {
			measure = true;
		} else if (e.getKeyCode() == KeyEvent.VK_EQUALS) {
			brightness++;
		} else if (e.getKeyCode() == KeyEvent.VK_MINUS) {
			brightness--;
		}// else if (e.getKeyCode() == KeyEvent.VK_E) {
		//	newenergy = !newenergy;
		//}
	}

	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	
}

class Obstacle {
	double x;
	double y;
	double r;
	
	Obstacle(double x, double y, double r) {
		this.x = x;
		this.y = y;
		this.r = r;
	}
}

class MeasurementApparatus {
	public int yi = 0;
	public int yf = 0;
	public int x = 0;
	public float totintensity = 0;
	public float avgintensity = 0;
	public int steps = 0;
	public int number = 0;
	
	public MeasurementApparatus(int yi, int yf, int number) {
		x = (int)((MicrowaveObstacleSimulation.detectorx-MicrowaveObstacleSimulation.xmin)/MicrowaveObstacleSimulation.ds);
		
		this.yi = yi;
		if (this.yi < 0)
			this.yi = 0;
		
		this.yf = yf;
		if (yf >= MicrowaveObstacleSimulation.ny)
			this.yf = MicrowaveObstacleSimulation.ny-1;
		
		this.number = number;
	}
	
	void measure(MicrowaveObstacleSimulation s, boolean calculateavg) {
		double intensity = 0.0;
		for (int y = yi; y <= yf; y++) {
			intensity += s.phix[x][y];
		}
		totintensity += intensity/(yf-yi);

		steps++;
		if (calculateavg) {
			avgintensity = totintensity/steps;
			totintensity = 0;
			steps = 0;
		}
	}
}

class RenderCanvas extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = -226046103513143999L;
	Image screen;
	MicrowaveObstacleSimulation parent;
	public static final double scalefactor = MicrowaveObstacleSimulation.ds*12.0;
	@Override
	public void paintComponent(Graphics real) {
		Graphics g = screen.getGraphics();
		g.setColor(new Color(50, 50, 50));
		g.fillRect(0, 0, imgwidth, imgheight);
		for (int i = 0; i < MicrowaveObstacleSimulation.nx; i++) {
			for (int j = 0; j < MicrowaveObstacleSimulation.ny; j++) {
				double scalingconstant = Math.pow(2.0, parent.brightness/2.0);
				int energy = (int) Math.round(parent.e[i][j]*0.0005*scalingconstant);
				int energyavg = (int) Math.round(parent.eavg[i][j]*0.0015*scalingconstant);
				int displacement = (int) Math.round(parent.u[i][j]*10*scalingconstant);
				int fluxx = (int) Math.round(parent.phix[i][j]*20*scalingconstant);
				int fluxy = (int) Math.round(parent.phiy[i][j]*20*scalingconstant);
				int bx = (int) Math.round(parent.bx[i][j]*50*scalingconstant);
				int by = (int) Math.round(parent.by[i][j]*50*scalingconstant);
				Color c = Color.WHITE;
				if (parent.view == 1) {
					c = new Color(clamp(Math.abs(fluxx)+128, 0, 255), clamp(127+energy, 0, 255), clamp(Math.abs(fluxy) + 128, 0, 255));
				} else if (parent.view == 2) {
					c = new Color(clamp(energyavg, 0, 255), 0, 0);
				} else if (parent.view == 3) {
					c = new Color(clamp(displacement, 0, 255), 0, clamp(-displacement, 0, 255));
				} else if (parent.view == 4) {
					c = new Color(clamp(Math.abs(bx)+128, 0, 255), 127, clamp(Math.abs(by) + 128, 0, 255));
				}
				if (parent.obstacles[i][j] == 1)
					g.setColor(c);
				else
					g.setColor(Color.WHITE);
				g.fillRect((int)(i*scalefactor), (int)((MicrowaveObstacleSimulation.ny-j-1)*scalefactor), (int)Math.ceil(scalefactor), (int)Math.ceil(scalefactor));
				//g.drawLine(i*7, (int)(parent.pos[i]*10.0+250), (i+1)*7, (int)(parent.pos[i+1]*10.0+250));
				//g.setColor(Color.BLUE);
				//g.fillRect(i*7, (int)(parent.pos[i]*10.0+250), 2, 2);
			}
		}
		
		g.setColor(new Color(1.0f, 0.4f, 1.0f, 0.5f));
		for (int i = 0; i < parent.detectornum; i++) {
			int x = (int)(parent.detectors[i].x*scalefactor);
			int y = (int)((parent.detectors[i].yi+parent.detectors[i].yf)*0.5*scalefactor);
			g.fillOval(x-3, y-3, 7, 7);
		}
		//g.setColor(.)
		
		double Imax = 0.1;

		for (int i = 0; i < parent.detectornum; i++) {
			if (parent.detectors[i].avgintensity > Imax)
				Imax = parent.detectors[i].avgintensity;
		}
		((Graphics2D)g).setRenderingHint(
		        RenderingHints.KEY_TEXT_ANTIALIASING,
		        RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		g.setFont(new Font("Arial Bold", Font.BOLD, 12));
		g.setColor(new Color(255, 50, 50));
		g.drawString("Average intensity: " + parent.avgintensity*parent.intensityfactor, 20, imgheight-20);
		g.drawString("Time: " + String.format("%2.2f", parent.time), 20, imgheight-40);
		
		g.setColor(Color.BLACK);
		g.fillRect(graphcx - 15*(parent.detectornum/2), imgheight - graphheight, 15*(parent.detectornum-1), graphheight);
		g.setColor(new Color(255, 50, 50));
		g.drawString("Expected intensity distribution", graphcx - 100, imgheight-graphheight+15);
		for (int i = 0; i < parent.detectornum - 1; i++) {
			g.setColor(new Color(155, 0, 0));
			g.drawLine(graphcx + 15*(i - parent.detectornum/2), imgheight-2-(int)(225.0*parent.detectors[i].avgintensity/Imax), graphcx + 15*(1 + i - parent.detectornum/2), imgheight-2-(int)(225.0*parent.detectors[i+1].avgintensity/Imax));
			g.setColor(new Color(255, 50, 50));
			g.fillRect(graphcx + 15*(i - parent.detectornum/2) - 2, imgheight-2-(int)(225.0*parent.detectors[i].avgintensity/Imax) - 2, 4, 4);
			if (i == parent.detectornum - 2)
				g.fillRect(graphcx + 15*(i + 1 - parent.detectornum/2) - 2, imgheight-2-(int)(225.0*parent.detectors[i+1].avgintensity/Imax) - 2, 4, 4);
		}
		real.drawImage(screen, 0, 0, parent);
	}
	
	private int clamp(int val, int min, int max) {
		if (val < min) return min;
		if (val > max) return max;
		return val;
	}
	int graphcx = 0;
	int graphcy = 0;
	int graphheight = 250;
	int imgwidth = 0;
	int imgheight = 0;
	public RenderCanvas(MicrowaveObstacleSimulation w) {
		parent = w;
		imgwidth = (int)Math.ceil(scalefactor*MicrowaveObstacleSimulation.nx);
		imgheight = (int)Math.ceil(scalefactor*MicrowaveObstacleSimulation.ny) + graphheight;
		screen = parent.createImage(imgwidth, imgheight);
		graphcy = (int)Math.ceil(scalefactor*MicrowaveObstacleSimulation.ny) + graphheight/2;
		graphcx = (int)Math.ceil(scalefactor*MicrowaveObstacleSimulation.nx) / 2;
	}
}
