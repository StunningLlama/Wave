import java.awt.Color;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class SutureSimulation extends JFrame implements Runnable {

	double[] y;
	double[] vy;
	double[] s;
	double[] x;
	double[] deq;
	double[] ay;
	double[] alpha;
	double[] phi;
	double[] tensions;

	double mp = 1.0; /*Pulley mass*/
	double T0 = 20; /*Tension*/
	double mu = 0.3;
	//double beta = Math.PI;
	//double c = Math.exp(mu*beta);
	double k = 0.3;
	double ks = 0.0;
	double s_even = 0.0;
	double s_odd = 0.0;
	double timestep = 0.04;
	double time = 0.0;
	double wound_length = 25.0;
	double wound_width = 4.0;
	double frac = 0.2;
	double min_spacing = 0.0;
	int pulleys = 25;
	
	
	RenderCanvas r;

	public static void main(String[] args) {
		SutureSimulation w = new SutureSimulation();
		w.run();
	}
	
	public SutureSimulation() {
		this.setSize(800,800);
		resetSuture();
		this.setVisible(true);
		r= new RenderCanvas(this);
		this.add(r);
	}
	
	public void resetSuture() {
		ay = new double[pulleys];
		y = new double[pulleys];
		x = new double[pulleys];
		vy = new double[pulleys];
		deq = new double[pulleys];
		alpha = new double[pulleys+1];
		s = new double[pulleys+1];
		phi = new double[pulleys];
		tensions = new double[pulleys+1];
		
		double stot = 2.0*wound_length/pulleys;
		s_even = stot*frac;
		s_odd = stot*(1.0-frac);
		k = 0.3*25.0/pulleys;
		
		s[0] = 1.0;
		x[0] = 0.0;
		for (int i = 0; i < pulleys; i++)
		{
			double r = ((double)i)/pulleys;
			ay[i]=0.0;
			vy[i]=0.0;
			if (i % 2 == 0) {
				//deq[i]=5.0;
				//y[i]=r*(1.0-r)*5.0;
				y[i] = Math.sqrt(1-(2.0*r-1.0)*(2.0*r-1.0))*(wound_width/2.0);
				deq[i] = y[i];
				s[i+1]=s_odd;
				//y[i] = 5.0;
			} else {
				//deq[i]=-5.0;
				//y[i]=-r*(1.0-r)*5.0;
				y[i] = -Math.sqrt(1-(2.0*r-1.0)*(2.0*r-1.0))*(wound_width/2.0);
				deq[i] = y[i];
				s[i+1]=s_even;
			}
		}

		for (int i = 1; i < pulleys; i++)
		{
			x[i] = x[i-1]-s[i];
		}
	}
	
	public void Physics() {
		alpha[0] = Math.PI/2;
		if (pulleys%2 == 1) {
			alpha[pulleys] = Math.atan((y[pulleys-1])/s[pulleys]);
		} else {
			alpha[pulleys] = -Math.atan((y[pulleys-1])/s[pulleys]);
		}
		for (int i = 1; i < pulleys; i++)
		{
			if (i % 2 == 0) {
				alpha[i] = -Math.atan((y[i-1]-y[i])/s[i]);
			} else {
				alpha[i] = Math.atan((y[i-1]-y[i])/s[i]);
			}
		}
		for (int i = 0; i < pulleys; i++) {
			phi[i] = alpha[i] + alpha[i+1];
		}
		tensions[0] = T0;
		for (int i = 1; i <= pulleys; i++) {
			tensions[i] = tensions[i-1]*Math.exp(-phi[i-1]*mu);
		}
		/*for (int i = 0; i < POINTS; i++)
		{
			if (i % 2 == 0) {
				ay[i] = -(T0*(1.0/Math.pow(c, (double)(i+1)) + 1.0/Math.pow(c, (double)(i)))) + (deq[i] - y[i])*k;
			} else {
				ay[i] = (T0*(1.0/Math.pow(c, (double)(i+1)) + 1.0/Math.pow(c, (double)(i)))) - (y[i] - deq[i])*k;
			}
			ay[i] /= mp;
		}*/
		for (int i = 0; i < pulleys; i++)
		{
			if (i % 2 == 0) {
				ay[i] = -(tensions[i]*Math.sin(alpha[i]) + tensions[i+1]*Math.sin(alpha[i+1])) + (deq[i] - y[i])*k;
			} else {
				ay[i] = (tensions[i]*Math.sin(alpha[i]) + tensions[i+1]*Math.sin(alpha[i+1])) - (y[i] - deq[i])*k;
			}
		}
		for (int i = 0; i < pulleys; i++)
		{
			if (i+2 < pulleys) {
				ay[i] += ks*(y[i+2]-y[i]);
			}
			if (i-2 >= 0) {
				ay[i] += ks*(y[i-2]-y[i]);
			}
		}
		for (int i = 0; i < pulleys; i++)
		{
			ay[i] /= mp;
		}
		for (int i = 0; i < pulleys; i++)
		{
			vy[i] += ay[i]*timestep;
		}
		for (int i = pulleys-1; i >= 0; i--)
		{
			double vend = 0.0;
			for (int j = pulleys - 1; j > pulleys-1; j--) {
				vend = vy[j]*(Math.sin(alpha[i])+Math.sin(alpha[i+1])) - vend;
			}
			
			if (i%2 == 0) {
				if (vy[i]*Math.sin(alpha[i+1]) > vend) {
					vy[i] = vend/Math.sin(alpha[i+1]);
				}
			} else {
				if (vy[i]*Math.sin(alpha[i+1]) < vend) {
					vy[i] = vend/Math.sin(alpha[i+1]);
				}
			}
			
			

			if (i % 2 == 0) {
				if (y[i] < min_spacing) {
					y[i] = min_spacing;
					vy[i] = 0.0;
				}
			} else {
				if (y[i] > -min_spacing) {
					y[i] = -min_spacing;
					vy[i] = 0.0;
				}
			}
		}
		

		for (int i = 0; i < pulleys; i++)
		{
			y[i] += vy[i]*timestep;
		}
		
		time += timestep;
	}
	
	public void compute() {


		FileWriter out = null;

		try {
			out = new FileWriter("C:\\Users\\Brandon\\Desktop\\output3.txt");
			for (pulleys = 8; pulleys <= 30; pulleys++) {
				for (frac = 0.05; frac <=0.95; frac += 0.05) {
					resetSuture();
					double alpha_even = Math.atan(wound_width/s_even);
					double alpha_odd = Math.atan(wound_width/s_odd);
					out.write(alpha_even + " " + alpha_odd);
					for (mu = 0.01; mu < 1.0; mu += 0.01) {
						resetSuture();
						while (true) {
							Physics();
							//try {
							//	Thread.sleep(50);
							//} catch (InterruptedException e) {}
							//r.repaint();
							if (Math.abs(y[pulleys-1]) < min_spacing+0.01 || time > 100.0) {
								out.write(" " + time );
								time = 0.0;
								break;
							}
						}
					}
					out.write("\n");
				}
			}
			if (out != null) {
				out.close();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	int i = 0;
	@Override
	public void run() {
		compute();
		/*while (true) {
			try {
				Thread.sleep(50);
			} catch (InterruptedException e) {}
			this.Physics();
			r.repaint();
		}*/
		System.out.println("Done.");
	}

	
	public static double absmin(double x, double y) {
		if (Math.abs(x) < Math.abs(y)) return x;
		return y;
	}
}

class RenderCanvas extends JPanel {
	Image screen;
	SutureSimulation parent;
	public static final int INTERVAL = 4;
	@Override
	public void paintComponent(Graphics real) {
		Graphics g = screen.getGraphics();
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, 800, 800);
		for (int i = 0; i < parent.pulleys; i++) {
			g.setColor(Color.RED);
			int radius = 2;
			g.drawOval(transformX(parent.x[i]) - radius, transformY(parent.y[i]) - radius, 2*radius, 2*radius);
			g.drawLine(transformX(parent.x[i]), transformY(parent.deq[i]), transformX(parent.x[i]), transformY(parent.y[i]));
			g.setColor(Color.BLUE);
			if (i < parent.pulleys - 1)
				g.drawLine(transformX(parent.x[i]), transformY(parent.y[i]), transformX(parent.x[i+1]), transformY(parent.y[i+1]));
		}
		real.drawImage(screen, 0, 0, parent);
	}
	
	int transformX(double x) {
		return (int)(700+20*x);
	}
	
	int transformY(double y) {
		return (int)(400-20*y);
	}
	
	private int clamp(int val, int min, int max) {
		if (val < min) return min;
		if (val > max) return max;
		return val;
	}
	public RenderCanvas(SutureSimulation w) {
		parent = w;
		screen = parent.createImage(800,800);
	}
}
