import processing.core.*;
import processing.data.*;
import processing.event.*;
import processing.opengl.*;

import java.nio.*;
import processing.core.PMatrix3D;

import java.util.HashMap;

import java.util.ArrayList;
import java.io.File;
import java.io.BufferedReader;
import java.io.PrintWriter;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.IOException;

public class swirl extends PApplet {

	class correspondence
	{
		public int a;
		public int b;
		public correspondence(int x,int y){a=x;b=y;}
		float r_inflation;
	};
	
	// Skate dancer on moving terrain
	float dz = 0; // distance to camera. Manipulated with wheel or when
	// float rx=-0.06*TWO_PI, ry=-0.04*TWO_PI; // view angles manipulated when
	// space pressed but not mouse
	float rx = 0, ry = 0; // view angles manipulated when space pressed but not
							// mouse
	Boolean twistFree = false, animating = false, tracking = false, center = true, gouraud = true,
			showControlPolygon = false, showNormals = false;
	int show_inflation = 0;
	float t = 0f, s = 0;
	
	
	boolean viewpoint = false;
	pt Viewer = P();
	pt F = P(50, -100, -100); // focus point: the camera is looking at it (moved when
						// 'f or 'F' are pressed
	pt Of = P(100, 100, 0), Ob = P(110, 110, 0); // red point controlled by the
													// user via mouseDrag : used
													// for inserting vertices
													// ...
	pt Vf = P(0, 0, 0), Vb = P(0, 0, 0);

	pt[] control_point_A = new pt[4];
	pt[] control_point_B = new pt[4];
	pt[] sample_A = new pt[1001];
	pt[] sample_B = new pt[1001];
	
	pt[] medial_axis = new pt[1001];
	correspondence[] corres =new correspondence[1001];
	pt[] pt_inflation=new pt[1001];
	pt[] intermediate=new pt[1001];
	
	pt[] ABintersect = new pt[1001];
	vec orientation;
	
	vec[] normal_A=new vec[1000];
	vec[] normal_B=new vec[1000];
	vec[] normal_intermediate=new vec[1000];
	vec[] normal_inflation=new vec[1000];
	
	int medial_axis_size=1001;

	public void setup() {
		myFace = loadImage("data/Shao.jpg"); // load image from file pic.jpg in
											// folder data *** replace that file
											// with your pic of your own face
		myFace2=loadImage("data/wang.jpg");
		
		
		textureMode(NORMAL);
		// p3D means that we will do 3D graphics
		P.declare();
		Q.declare();
		PtQ.declare(); // P is a polyloop in 3D: declared in pts
		// P.resetOnCircle(12,100); // used to get started if no model exists on
		// file
		P.loadPts("data/A.pts"); // loads saved model from file
		//Q.loadPts("data/B.pts"); // loads saved model from file

		control_point_A[0] = P.G[0];
		control_point_A[3] = P.G[3];
		control_point_A[1] = P.G[1];
		control_point_A[2] = P.G[2];

		control_point_B[0] = P.G[0];
		control_point_B[3] = P.G[3];
		control_point_B[1] = P.G[4];
		control_point_B[2] = P.G[5];

		calculate_sample();

		// test ABintersect
		ABintersect[0] = control_point_A[0];
		ABintersect[1000] = control_point_A[3];

		medial_axis[0] = control_point_A[0];
		medial_axis[1000] = control_point_A[3];
		
		
		intermediate[0] = control_point_A[0];
		intermediate[1000] = control_point_A[3];
		
		pt_inflation[0] = control_point_A[0];
		pt_inflation[1000] = control_point_A[3];
		
		corres[0]=new correspondence(0,0);
		corres[1000]=new correspondence(1000,1000);

		compute_medial_axis();
		find_normal(sample_A,normal_A,1000);
		find_normal(sample_B,normal_B,1000);
		
	}

	public void draw() {

		background(255);
		pushMatrix(); // to ensure that we can restore the standard view before
						// writing on the canvas

		float fov = PI / 3.0f;
		float cameraZ = (height / 2.0f) / tan(fov / 2.0f);
		camera(0, 0, cameraZ, 0, 0, 0, 0, 1, 0); // sets a standard perspective
		perspective(fov, 1.0f, 0.1f, 10000);

		translate(0, 0, dz); // puts origin of model at screen center and moves
								// forward/away by dz
		lights(); // turns on view-dependent lighting
		rotateX(rx);
		rotateY(ry); // rotates the model around the new origin (center of
						// screen)
		rotateX(PI / 2); // rotates frame around X to make X and Y basis vectors
							// parallel to the floor
		if (center)
			translate(-F.x, -F.y, -F.z);
		noStroke(); // if you use stroke, the weight (width) of it will be
					// scaled with you scaleing factor
		showFrame(50); // X-red, Y-green, Z-blue arrows
		// fill(yellow); pushMatrix(); translate(0,0,-1.5f); box(400,400,1);
		// popMatrix(); // draws floor as thin plate
		// fill(magenta); show(F,4); // magenta focus point (stays at center of
		// screen)
		// fill(magenta,100); showShadow(F,5); // magenta translucent shadow of
		// focus point (after moving it up with 'F'

		computeProjectedVectors(); // computes screen projections I, J, K of
									// basis vectors (see bottom of pv3D): used
									// for dragging in viewer's frame

		stroke(255, 0, 0);
		noFill();
		showSphere(control_point_A[0], 15);
		showSphere(control_point_A[1], 15);
		showSphere(control_point_A[2], 15);
		showSphere(control_point_A[3], 15);
		//stroke(0, 0, 255);
		strokeWeight(3);
		noFill();
		bezier(control_point_A);

		stroke(0, 255, 0);
		showSphere(control_point_B[0], 15);
		showSphere(control_point_B[1], 15);
		showSphere(control_point_B[2], 15);
		showSphere(control_point_B[3], 15);
		//stroke(0, 255, 0);
		strokeWeight(3);
		bezier(control_point_B);

		stroke(0, 0, 255);
		strokeWeight(7);
		for (int i = 0; i < medial_axis_size; i++) {
			showSphere(medial_axis[i], 4);
		}	

		draw_curve_quad(8,6);
		draw_net();		
		
		if(show_inflation != 0){
			find_normal_inflation(pt_inflation,normal_inflation,medial_axis_size-1);
			draw_quad_inflation(8);		
		}
		
	
		for(int i=1;i<medial_axis_size;i++)
		{
			intermediate[i]=L(L(sample_A[corres[i].a],t,medial_axis[i]),t,L(medial_axis[i],t,sample_B[corres[i].b]));
		}		
		
		
		if(animating) 
		{
			find_normal(intermediate,normal_intermediate,medial_axis_size-1);
			draw_quad(8,8);			
			
			t+=0.01f;
			if(t>1)
				t=0;
		}

		pp = P.idOfVertexWithClosestScreenProjectionTo(Mouse()); 
		
		if (mousePressed && !keyPressed) {
			Of = pick(mouseX, mouseY);
			for (int i = 0; i < 10; i++)
				attractFront(PtQ.G, 0.001f);
		} // modifies bu & bv
		else {
			for (int i = 0; i < 10; i++)
				slide(PtQ.G, 0.0002f);
		}
		for (int i = 0; i < 10; i++)
			attractBack(PtQ.G, 50, 0.001f);
		/*
		 * Vf.setTo(coons(PtQ.G,fu,fv)); vec Nf =
		 * CoonsNormal(PtQ.G,fu,fv,0.01f); Vb.setTo(coons(PtQ.G,bu,bv)); vec Nb
		 * = CoonsNormal(PtQ.G,bu,bv,0.01f);
		 * 
		 * fill(green); show(Of,3); // fill(red,100); showShadow(Of,5); } //
		 * show ret tool point and its shadow showMan(PtQ.G,50);
		 * 
		 * // fill(red); arrow(Vc,V(50,Nc),5); fill(yellow,100);
		 * show(P(Vc,15,Nc),15); if(tracking) F.setTo(P(F,0.01f,Vf));
		 * 
		 * if(showNormals){ strokeWeight(2); stroke(magenta);
		 * showNormals(PtQ.G,0.1f,0.01f); stroke(orange);
		 * showNormals(PtQ.G,0.1f,0.2f); } noFill(); stroke(blue);
		 * strokeWeight(2);
		 */
		popMatrix(); // done with 3D drawing. Restore front view for writing
						// text on canvas

		// for demos: shows the mouse and the key pressed (but it may be hidden
		// by the 3D model)
		// if(keyPressed) {stroke(red); fill(white);
		// ellipse(mouseX,mouseY,26,26); fill(red);
		// text(key,mouseX-5,mouseY+4);}
		if (scribeText) {
			fill(black);
			displayHeader();
		} // dispalys header on canvas, including my face
		if (scribeText && !filming)
			displayFooter(); // shows menu at bottom, only if not filming
		//if (animating) {
		//	t += PI / 180 / 2;
		//	if (t >= TWO_PI)
		//		t = 0;
		//	s = (cos(t) + 1.f) / 2;
		//} // periodic change of time
		if (filming && (animating || change))
			saveFrame("FRAMES/F" + nf(frameCounter++, 4) + ".tif"); // save next
																	// frame to
																	// make a
																	// movie
		change = false; // to avoid capturing frames when nothing happens
						// (change is set uppn action)
		uvShow(); // **<01
		// scribe("F1 x: "+global_I.x+" F1 y: "+global_I.y+" F1 z:
		// "+global_I.z,20,100);
		// scribe("F0 x: "+F0.I.x+"F0 y: "+F0.I.y+"F0 z: "+F0.I.z,20,150);
	}

	public void calculate_sample() {

		for (int i = 0; i <= 1000; i++) {
			sample_A[i] = bezierPoint(control_point_A, ((float) i) / 1000);
			sample_B[i] = bezierPoint(control_point_B, ((float) i) / 1000);
		}

	}
/*
	public void calculate_MA() {

		for (int i = 0; i < 100; i++) {

			// println("i = " + i);
			float min_temp = Integer.MAX_VALUE, min_final = Integer.MAX_VALUE;
			for (int i_1 = 0; i_1 <=1000; i_1++) {

				// println("i_1 = " + i_1);
				min_temp = d(sample_A[i], sample_B[i_1]);
				if (min_final >= min_temp) {
					// println("i = " + i);
					min_final = min_temp;
					medial_axis[i] = P(sample_A[i], sample_B[i_1]);
				}
			}
		}
	}
	*/
	public void find_normal(pt[] sample,vec[] sample_normal,int size)
	{
		float parameter;
		vec temp1=U(V(sample[0],sample[1]));//tangent
		vec X1=new vec(1,0,0);//find a random vector to form a plane
		if( abs(dot(X1,temp1))>0.999 )//check parallel
			X1=V(0,1,0);
		vec J1=U( N(X1,temp1) );//get normal vector to tangent
		sample_normal[0]=J1;
	
		for (int i = 1; i < size; i++) {
			vec AC=U(V( sample[i-1] , sample[i+1] ));//U
			parameter=dot(sample_normal[i-1],AC)/dot(AC,AC);//x=(Na.U)/(U.U)
			sample_normal[i]=U(   M(sample_normal[i-1], V(parameter,AC) )  );//Nb=Na-xU
		}

	}
	
	
	public void find_normal_inflation(pt[] sample,vec[] sample_normal,int size)
	{
		float parameter;
		vec temp1=U(V(sample[0],sample[1]));//tangent
		vec X1=V(sample_A[corres[1].a],sample_B[corres[1].b]);//find a random vector to form a plane
		if( abs(dot(X1,temp1))>0.999 )//check parallel
			X1=V(0,1,0);
		vec J1=U( N(X1,temp1) );//get normal vector to tangent
		sample_normal[0]=J1;
	
		for (int i = 1; i < size; i++) {
			vec AC=U(V( sample[i-1] , sample[i+1] ));//U
			parameter=dot(sample_normal[i-1],AC)/dot(AC,AC);//x=(Na.U)/(U.U)
			sample_normal[i]=U(   M(sample_normal[i-1], V(parameter,AC) )  );//Nb=Na-xU
		}

	}

	
	
	public void draw_tube(pt[] sample,vec[] normal_sample,int r,int s,int size)
	{
		
		for(int i=0;i<size;i++)
		{
			if(i%100<=50)
				{stroke(0, 255, 255);
				fill(0,255,255);}
			else
				{stroke(255,255,0);
				fill(255,255,0);
				}

			vec I=U(V(sample[i],sample[i+1]));
			vec J=normal_sample[i];
			vec K=U(N(I,J));
			beginShape(QUAD_STRIP);
			for (float t=0; t<TWO_PI; t+=TWO_PI/s) 
			{
				v(P(P(sample[i],r*cos(t),J),r*sin(t),K)); v(P(P(sample[i+1],r*cos(t),J),r*sin(t),K));
			}
			endShape();
			
		}
	}
	
	//display input curves as tubes using sample points of curves and parallel transport norms
	public void draw_curve_quad(int r,int s)
	{
		int j=0;
		for(int i=0;i<medial_axis_size-1;i+=1)
		{
			vec I=U(V(sample_A[corres[i].a],sample_A[corres[i+1].a]));
			vec J=normal_A[i];
			vec K=U(N(I,J));
			
			vec I1=U(V(sample_B[corres[i].b],sample_B[corres[i+1].b]));
			vec J1=normal_B[i];
			vec K1=U(N(I1,J1));
			
			
			for (float t=0; t<TWO_PI; t+=TWO_PI/s,j++) 
			{
				if(j%2==0)
				{					
						//stroke(255, 0, 0);
					noStroke();
					fill(255, 0, 0);
				}
				else
					{
					noStroke();
					//stroke(0,0,0);
					fill(0, 0, 0);
					}
				beginShape(QUAD);
				v(P(P(sample_A[corres[i].a],r*cos(t),J),r*sin(t),K)); 
				v(P(P(sample_A[corres[i+1].a],r*cos(t),J),r*sin(t),K));
				
				v(P(P(sample_A[corres[i+1].a],r*cos(t-TWO_PI/s),J),r*sin(t-TWO_PI/s),K));
				v(P(P(sample_A[corres[i].a],r*cos(t-TWO_PI/s),J),r*sin(t-TWO_PI/s),K)); 

				endShape();
				
				if(j%2==0)
				{				noStroke();	
						//stroke(0, 255, 0);
						fill(0,255,0);
				}
				else{
					noStroke();
					//stroke(0,0,0);
					fill(0,0,0);
					}
				beginShape(QUAD);
				v(P(P(sample_B[corres[i].b],r*cos(t),J1),r*sin(t),K1)); 
				v(P(P(sample_B[corres[i+1].b],r*cos(t),J1),r*sin(t),K1));
				
				v(P(P(sample_B[corres[i+1].b],r*cos(t-TWO_PI/s),J1),r*sin(t-TWO_PI/s),K1));
				v(P(P(sample_B[corres[i].b],r*cos(t-TWO_PI/s),J1),r*sin(t-TWO_PI/s),K1)); 

				endShape();
	
			}
			if(i%((medial_axis_size-1)/10)!=((medial_axis_size-1)/20))
				j++;

		}
	}	
	
	// draw ball morph between two curves
	public void draw_quad(int r,int s)
	{
		int j=0;
		for(int i=0;i<medial_axis_size-1;i+=1)
		{
			vec I=U(V(intermediate[i],intermediate[i+1]));
			vec J=normal_intermediate[i];
			vec K=U(N(I,J));
			
			for (float t=0; t<TWO_PI; t+=TWO_PI/s,j++) 
			{
				if(j%2==0)
				{
					noStroke();
					fill(0,0,255);
				}
				else
				{
					noStroke();
					fill(0,0,0);
				}

				beginShape(QUAD);
				v(P(P(intermediate[i],r*cos(t),J),r*sin(t),K)); 
				v(P(P(intermediate[i+1],r*cos(t),J),r*sin(t),K));
				
				v(P(P(intermediate[i+1],r*cos(t-TWO_PI/s),J),r*sin(t-TWO_PI/s),K));
				v(P(P(intermediate[i],r*cos(t-TWO_PI/s),J),r*sin(t-TWO_PI/s),K)); 

				endShape();
			}
			if(i%((medial_axis_size-1)/10)!=((medial_axis_size-1)/20))
				j++;
		}
	}

	
	public void draw_quad_inflation(int s)
	{
		int j=0;
		for(int i=0;i<medial_axis_size-1;i+=1)
		{
			vec I=U(V(pt_inflation[i],pt_inflation[i+1]));
			vec J=normal_inflation[i];
			vec K=U(N(I,J));
			
			vec line_correspts = U( V(sample_A[corres[i].a],sample_B[corres[i].b]) ); 
			//calculate dot of line between correspongding points and normal, in order to figure out rotation angle
			float my_temp = asin(dot(line_correspts,J));
			
			for (float t= -TWO_PI/s - my_temp; t <TWO_PI*show_inflation/2 - TWO_PI/s - my_temp; t+=TWO_PI/s,j++)
				//show_inflation control what part to show
			{
				if(j%2==0)
				{	
					strokeWeight(2);
					stroke(0,0,0);
					noFill();

				}
				else
				{
					strokeWeight(2);
					noFill();
					stroke(180,180,180);
				}
				beginShape(QUAD);
				v(P(P(pt_inflation[i],corres[i].r_inflation*cos(t),J),corres[i].r_inflation*sin(t),K)); 
				v(P(P(pt_inflation[i+1],corres[i+1].r_inflation*cos(t),J),corres[i+1].r_inflation*sin(t),K));
				
				v(P(P(pt_inflation[i+1],corres[i+1].r_inflation*cos(t-TWO_PI/s),J),corres[i+1].r_inflation*sin(t-TWO_PI/s),K));
				v(P(P(pt_inflation[i],corres[i].r_inflation*cos(t-TWO_PI/s),J),corres[i].r_inflation*sin(t-TWO_PI/s),K)); 

				endShape();
			}
			if(i%((medial_axis_size-1)/10)!=((medial_axis_size-1)/20))
				j++;
		}
	}
	
	// use quadratic Bezier curve to simulate the transversal arcs
	public void draw_net()
	{
		for(int i=0;i<medial_axis_size;i++)
		{
			noFill();
			strokeWeight(2);
			stroke(255,255,0);
			for(float t=0;t<1f;t+=0.01f)
			{
				pt a=L(L(sample_A[corres[i].a],t,medial_axis[i]),t,L(medial_axis[i],t,sample_B[corres[i].b]));
				pt b=L(L(sample_A[corres[i].a],t+0.01f,medial_axis[i]),t+0.01f,L(medial_axis[i],t+0.01f,sample_B[corres[i].b]));
				line(a.x,a.y,a.z,b.x,b.y,b.z);
			}
		}
	}
	
	
	/*
	public void draw_medial_tube()
	{	
		
		for(int i=0;i<100;i++)
		{
			float radius_1=d( sample_A[corres[i].a],sample_B[corres[i].b]    )/2;
			pt center_1=P(sample_A[corres[i].a],sample_B[corres[i].b]);
			
			float radius_2=d( sample_A[corres[i+1].a],sample_B[corres[i+1].b]    )/2;
			pt center_2=P(sample_A[corres[i+1].a],sample_B[corres[i+1].b]);
			
			
			vec I=U(V(sample[i],sample[i+1]));
			vec J=normal_sample[i];
			vec K=U(N(I,J));
			
			for (float t=0; t<TWO_PI; t+=TWO_PI/s,j++) 
			{
				if(j%2==0)
				{
					if(color==1)
						stroke(255, 0, 0);
					else if(color==2)
						stroke(0, 0, 255);
					else
						stroke(0, 255, 0);
				}
				else
					stroke(0,0,0);
				beginShape(QUAD);
				v(P(P(sample[i],r*cos(t),J),r*sin(t),K)); 
				v(P(P(sample[i+1],r*cos(t),J),r*sin(t),K));
				
				v(P(P(sample[i],r*cos(t-TWO_PI/s),J),r*sin(t-TWO_PI/s),K)); 
				v(P(P(sample[i+1],r*cos(t-TWO_PI/s),J),r*sin(t-TWO_PI/s),K));

				endShape();
			}		
		}
	}*/
	
	
	// find closest point
	public int find_closest_projection(pt p, pt[] sample) {
		float distance = d(sample[0], p), temp;
		int index = 0;
		for (int i = 1; i <=1000; i++) {
			temp = d(sample[i], p);
			if (temp < distance) {
				distance = temp;
				index = i;
			}
		}
		return index;
	}

	public void compute_medial_axis() {

		float temp_ptinplane;
		pt AinBplane, BinAplane, median_AB;
		vec line_medianAB;
		vec t1 = V(), t2 = V();

		orientation = U(A(U(bezierTangent(control_point_A, 0)), U(bezierTangent(control_point_B, 0))));

		for (int i = 1; i < 1000; i++) 
		{
			medial_axis[i] = P(medial_axis[i - 1], 7f, orientation);
			
			for (int ii = 0; ii < 3; ii++) // three rounds of guess to obtain better location of next point
			{
				int temp_A = find_closest_projection(medial_axis[i], sample_A);
				int temp_B = find_closest_projection(medial_axis[i], sample_B);

				
				//see whether is parallel
				t1 = U(bezierTangent(control_point_A, ((float) temp_A) / 1000));
				t2 = U(bezierTangent(control_point_B, ((float) temp_B) / 1000));

				if(parallel(t1,t2)){
					medial_axis[i] = P(sample_A[temp_A], sample_B[temp_B]);
					ABintersect[i] = P(sample_A[temp_A], sample_B[temp_B]);
					break;
				}
				
				
				// get modified guess
				vec normA = V(medial_axis[i], sample_A[temp_A]);
				vec normB = V(medial_axis[i], sample_B[temp_B]);

				// tangent A intersect planeB
				temp_ptinplane = ((sample_B[temp_B].x - sample_A[temp_A].x) * normB.x
						+ (sample_B[temp_B].y - sample_A[temp_A].y) * normB.y
						+ (sample_B[temp_B].z - sample_A[temp_A].z) * normB.z)
						/ (normB.x * t1.x + normB.y * t1.y + normB.z * t1.z);

				AinBplane = P(sample_A[temp_A].x + t1.x * temp_ptinplane, sample_A[temp_A].y + t1.y * temp_ptinplane,
						sample_A[temp_A].z + t1.z * temp_ptinplane);

				// tangent B intersect planeA
				temp_ptinplane = ((-sample_B[temp_B].x + sample_A[temp_A].x) * normA.x
						+ (-sample_B[temp_B].y + sample_A[temp_A].y) * normA.y
						+ (-sample_B[temp_B].z + sample_A[temp_A].z) * normA.z)
						/ (normA.x * t2.x + normA.y * t2.y + normA.z * t2.z);

				BinAplane = P(sample_B[temp_B].x + t2.x * temp_ptinplane, sample_B[temp_B].y + t2.y * temp_ptinplane,
						sample_B[temp_B].z + t2.z * temp_ptinplane);
				
				// average of two intersects + average of two lines between it and tangent points
				
				ABintersect[i] = P(AinBplane, BinAplane);
				
				vec new_guessLine = U(V(V(ABintersect[i], sample_A[temp_A]), V(ABintersect[i], sample_B[temp_B])));

				medial_axis[i] = P(ABintersect[i],
						V(dot(V(ABintersect[i], medial_axis[i]), new_guessLine), new_guessLine));
				
				
				corres[i] = new correspondence(temp_A,temp_B);
				pt_inflation[i] = P(sample_A[temp_A],sample_B[temp_B]);
				corres[i].r_inflation = Math.max((d(sample_A[temp_A],sample_B[temp_B]) / 2) -4 , 0);
			}			
						
			// just media axis, not enough
			// medial_axis[i] = P(sample_A[temp_A], sample_B[temp_B]);

			orientation = U(A(t1, t2));			
			
			// judge whether MA reaches the end point of curve
			// points exceed end point can cause unpredictable errors
			if(  Float.isNaN(d(medial_axis[i],medial_axis[1000]))    )
			{
				medial_axis_size=i;
				break;
			}

		}
	}

	public void penFill(int c, float w) {
		stroke(c);
		strokeWeight(w);
		fill(c);
	}

	public void keyPressed() {
		if (key == '`')
			picking = true;
		if (key == '?')
			scribeText = !scribeText;
		if (key == '!')
			snapPicture();
		if (key == '~')
			filming = !filming;
		if (key == ']')
			showControlPolygon = !showControlPolygon;
		if (key == '|')
			showNormals = !showNormals;
		if (key == 'G')
			gouraud = !gouraud;
		if (key == 'q')
			Q.copyFrom(P);
		if (key == 'p')
			P.copyFrom(Q);
		if (key == 'e') {
			PtQ.copyFrom(Q);
			Q.copyFrom(P);
			P.copyFrom(PtQ);
		}
		if (key == '=') {
			bu = fu;
			bv = fv;
		}
		// if(key=='.') F=P.Picked(); // snaps focus F to the selected vertex of
		// P (easier to rotate and zoom while keeping it in center)
		if (key == 'c')
			center = !center; // snaps focus F to the selected vertex of P
								// (easier to rotate and zoom while keeping it
								// in center)
		if (key == 't')
			tracking = !tracking; // snaps focus F to the selected vertex of P
									// (easier to rotate and zoom while keeping
									// it in center)
		if (key == 'x' || key == 'z' || key == 'd')
			P.setPickedTo(pp); // picks the vertex of P that has closest
								// projeciton to mouse
		if (key == 'd')
			P.deletePicked();
		if (key == 'i')
			P.insertClosestProjection(Of); // Inserts new vertex in P that is
											// the closeset projection of O
		if (key == 'W') {
			P.savePts("data/pts");
			Q.savePts("data/pts2");
		} // save vertices to pts2
		if (key == 'L') {
			P.loadPts("data/pts");
			Q.loadPts("data/pts2");
		} // loads saved model
		if (key == 'w')
			P.savePts("data/pts"); // save vertices to pts
		if (key == 'l')
			P.loadPts("data/pts");
		if (key == 'a')
			animating = !animating; // toggle animation
		if (key == 'A')
			show_inflation = (show_inflation + 1)%3; // show whole,half,none of inflation of average curve.
		if (key == ',')
			viewpoint = !viewpoint;
		if (key == '#')
			exit();
		change = true;
	}

	public void mouseWheel(MouseEvent event) {
		dz += event.getAmount() * 10;
		change = true;
	}

	public void mousePressed() {
		if (!keyPressed)
			picking = true;
	}

	public void mouseMoved() {
		if (keyPressed && key == ' ') {
			rx -= PI * (mouseY - pmouseY) / height;
			ry += PI * (mouseX - pmouseX) / width;
		}
		;
		if (keyPressed && key == 's')
			dz += (float) (mouseY - pmouseY); // approach view (same as wheel)
		if (keyPressed && key == 'v') { // **<01
			u += (float) (mouseX - pmouseX) / width;
			u = max(min(u, 1), 0);
			v += (float) (mouseY - pmouseY) / height;
			v = max(min(v, 1), 0);
		}
	}

	public void mouseDragged() {
		if (!keyPressed) {
			Of.add(ToIJ(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));
		}
		if (keyPressed && key == CODED && keyCode == SHIFT) {
			Of.add(ToK(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));
		}
		;
		if (keyPressed && key == 'x')
			P.movePicked(ToIJ(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));
		if (keyPressed && key == 'z')
			P.movePicked(ToK(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));

		// runs every time user moves control point of two curves
		// re-initialize 
		control_point_A[0] = P.G[0];
		control_point_A[3] = P.G[3];
		control_point_A[1] = P.G[1];
		control_point_A[2] = P.G[2];

		control_point_B[0] = P.G[0];
		control_point_B[3] = P.G[3];
		control_point_B[1] = P.G[4];
		control_point_B[2] = P.G[5];
		
		calculate_sample();
		
		ABintersect[0] = control_point_A[0];
		ABintersect[100] = control_point_A[3];

		medial_axis[0] = control_point_A[0];
		medial_axis[100] = control_point_A[3];

		intermediate[0] = control_point_A[0];
		intermediate[100] = control_point_A[3];		
		
		pt_inflation[0] = control_point_A[0];
		pt_inflation[100] = control_point_A[3];		
		
		corres[0]=new correspondence(0,0);
		corres[100]=new correspondence(1000,1000);

		compute_medial_axis();
		find_normal(sample_A,normal_A,1000);
		find_normal(sample_B,normal_B,1000);
		
		
		
		//if (keyPressed && key == 'X')
			//P.moveAll(ToIJ(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));
		//if (keyPressed && key == 'Z')
			//P.moveAll(ToK(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));
		//if (keyPressed && key == 'f') { // move focus point on plane
			//if (center)
				//F.sub(ToIJ(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));
			//else
				//F.add(ToIJ(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));
		//}
		/*if (keyPressed && key == 'F') { // move focus point vertically
			if (center)
				F.sub(ToK(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));
			else
				F.add(ToK(V((float) (mouseX - pmouseX), (float) (mouseY - pmouseY), 0)));
		}*/
	}

	// **** Header, footer, help text on canvas
	public void displayHeader() { // Displays title and authors face on screen
		scribeHeader(title, 0);
		scribeHeaderRight(name);
		fill(white);
		image(myFace2, width - myFace.width / 4, 25, myFace.width / 4, myFace.height / 4);
		image(myFace,width - myFace.width / 4, 25+myFace.height / 4, myFace.width / 4, myFace.height / 4);
	}

	public void displayFooter() { // Displays help text at the bottom
		scribeFooter(guide, 1);
		scribeFooter(menu, 0);
	}

	String title = "6491 P3 2015: Curve average in 3D", name = "   Shao         Wang",
			menu = "space:rotate, s/wheel:closer, a:anim, A:show (half) inflation, AA: show whole inflation #:quit",
			guide = "CURVES x/z:select&edit"; // user's
																									// guide
	/********
	 * Editor of an Animated Coons Patch
	 * 
	 * Implementation steps: <01 Manual control of (u,v) parameters. <02 Draw 4
	 * boundary curves CT(u), CB(u), SL(v), CR(v) using proportional Neville <03
	 * Compute and show Coons point C(u,v) <04 Display quads filed one-by-one
	 * for the animated Coons patch <05 Compute and show normal at C(u,v) and a
	 * ball ON the patch
	 * 
	 */
	// **<01: mouseMoved; 'v', draw: uvShow()
	float u = 0, v = 0;
	float bu = 0.5f, bv = 0.5f; // back foot
	float fu = 0.55f, fv = 0.55f; // front foot

	public void uvShow() {
		fill(red);
		if (keyPressed && key == 'v')
			text("u=" + u + ", v=" + v, 10, 30);
		noStroke();
		fill(blue);
		ellipse(u * width, v * height, 5, 5);
	}

	/*
	 * 0 1 2 3 11 4 10 5 9 8 7 6
	 */
	public pt coons(pt[] P, float s, float t) {
		pt Lst = L(L(P[0], s, P[3]), t, L(P[9], s, P[6]));
		pt Lt = L(N(0, P[0], 1.f / 3, P[1], 2.f / 3, P[2], 1, P[3], s), t,
				N(0, P[9], 1.f / 3, P[8], 2.f / 3, P[7], 1, P[6], s));
		pt Ls = L(N(0, P[0], 1.f / 3, P[11], 2.f / 3, P[10], 1, P[9], t), s,
				N(0, P[3], 1.f / 3, P[4], 2.f / 3, P[5], 1, P[6], t));
		return P(Ls, V(Lst, Lt));
	}

	public pt B(pt A, pt B, pt C, float s) {
		return L(L(A, s, B), s, L(B, s, C));
	}

	public pt B(pt A, pt B, pt C, pt D, float s) {
		return L(B(A, B, C, s), s, B(B, C, D, s));
	}

	public pt B(pt A, pt B, pt C, pt D, pt E, float s) {
		return L(B(A, B, C, D, s), s, B(B, C, D, E, s));
	}

	public pt N(float a, pt A, float b, pt B, float t) {
		return L(A, (t - a) / (b - a), B);
	}

	public pt N(float a, pt A, float b, pt B, float c, pt C, float t) {
		return N(a, N(a, A, b, B, t), c, N(b, B, c, C, t), t);
	}

	public pt N(float a, pt A, float b, pt B, float c, pt C, float d, pt D, float t) {
		return N(a, N(a, A, b, B, c, C, t), d, N(b, B, c, C, d, D, t), t);
	}

	public pt N(float a, pt A, float b, pt B, float c, pt C, float d, pt D, float e, pt E, float t) {
		return N(a, N(a, A, b, B, c, C, d, D, t), e, N(b, B, c, C, d, D, e, E, t), t);
	}

	public void drawBorders(pt[] P) {
		float e = 0.01f;
		beginShape();
		for (float t = 0; t < 1.001f; t += e)
			v(coons(P, 0, t));
		endShape();
		beginShape();
		for (float t = 0; t < 1.001f; t += e)
			v(coons(P, 1, t));
		endShape();
		beginShape();
		for (float t = 0; t < 1.001f; t += e)
			v(coons(P, t, 0));
		endShape();
		beginShape();
		for (float t = 0; t < 1.001f; t += e)
			v(coons(P, t, 1));
		endShape();
	}

	public vec CoonsNormal(pt[] P, float u, float v, float e) {
		vec Tu = V(coons(P, u - e, v), coons(P, u + e, v));
		vec Tv = V(coons(P, u, v - e), coons(P, u, v + e));
		return U(N(Tu, Tv));
	}

	public void shadeSurface(pt[] P, float e) {
		for (float s = 0; s < 1.001f - e; s += e)
			for (float t = 0; t < 1.001f - e; t += e) {
				beginShape();
				v(coons(P, s, t));
				v(coons(P, s + e, t));
				v(coons(P, s + e, t + e));
				v(coons(P, s, t + e));
				endShape(CLOSE);
			}
	}

	public void shadeSurfaceTextured(pt[] P, float e) {
		fill(white);
		for (float s = 0; s < 1.001f - e; s += e)
			for (float t = 0; t < 1.001f - e; t += e) {
				beginShape();
				texture(myFace);
				vTextured(coons(P, s, t), s, t);
				vTextured(coons(P, s + e, t), s + e, t);
				vTextured(coons(P, s + e, t + e), s + e, t + e);
				vTextured(coons(P, s, t + e), s, t + e);
				endShape(CLOSE);
			}
	}

	public void shadeSurfaceGouraud(pt[] P, float e, float ee) {
		Boolean col = true;
		for (float s = 0; s < 1.001f - e; s += e) {
			col = !col;
			for (float t = 0; t < 1.001f - e; t += e) {
				if (col)
					fill(cyan);
				else
					fill(magenta);
				col = !col;
				beginShape();
				nv(CoonsNormal(P, s, t, ee));
				v(coons(P, s, t));
				nv(CoonsNormal(P, s + e, t, ee));
				v(coons(P, s + e, t));
				nv(CoonsNormal(P, s + e, t + e, ee));
				v(coons(P, s + e, t + e));
				nv(CoonsNormal(P, s, t + e, ee));
				v(coons(P, s, t + e));
				endShape(CLOSE);
			}
		}
	}

	public void showNormals(pt[] P, float e, float ee) {
		for (float s = 0; s < 1.001f - e; s += e)
			for (float t = 0; t < 1.001f - e; t += e)
				show(coons(P, s, t), 50, CoonsNormal(P, s, t, ee));
	}

	public void slide(pt[] P, float e) {
		float nfu = fu, nfv = fv;
		float z = coons(P, fu, fv).z;
		for (float a = 0; a <= TWO_PI; a += PI / 20) {
			float nu = fu + e * cos(a), nv = fv + e * sin(a);
			float nz = coons(P, nu, nv).z;
			if (nz < z) {
				z = nz;
				nfu = nu;
				nfv = nv;
			}
		}
		fu = nfu;
		fv = nfv;
	}

	public void attractFront(pt[] P, float e) {
		float nfu = fu, nfv = fv;
		float z = d(coons(P, fu, fv), Of);
		for (float a = 0; a <= TWO_PI; a += PI / 20) {
			float nu = fu + e * cos(a), nv = fv + e * sin(a);
			float nz = d(coons(P, nu, nv), Of);
			if (nz < z) {
				z = nz;
				nfu = nu;
				nfv = nv;
			}
		}
		fu = nfu;
		fv = nfv;
	}

	public void attractBack(pt[] P, float r, float e) {
		float nbu = bu, nbv = bv;
		pt O = coons(PtQ.G, fu, fv);
		float z = abs(d(coons(PtQ.G, bu, bv), O) - r);
		for (float a = 0; a <= TWO_PI; a += PI / 20) {
			float nu = bu + e * cos(a), nv = bv + e * sin(a);
			float nz = abs(d(coons(P, nu, nv), O) - r);
			if (nz < z) {
				z = nz;
				nbu = nu;
				nbv = nv;
			}
		}
		bu = nbu;
		bv = nbv;
	}

	public void showMan(pt[] P, float h) {
		pt Of = coons(PtQ.G, fu, fv);
		vec Nf = CoonsNormal(PtQ.G, fu, fv, 0.01f);
		pt Ff = P(Of, 5, Nf);
		pt Kf = P(Of, h, Nf);
		pt Ob = coons(PtQ.G, bu, bv);
		vec Nb = CoonsNormal(PtQ.G, bu, bv, 0.01f);
		pt Fb = P(Ob, 5, Nb);
		pt Kb = P(Ob, h, Nb);
		float d = d(Kf, Kb) / 2, b = sqrt(sq(h) - sq(d));
		vec V = V(Nf, Nb);
		pt B = P(P(Kf, Kb), b, V);
		pt T = P(B, h, V);
		fill(red);
		show(Ff, 5);
		collar(Of, V(h, Nf), 1, 5);
		show(Kf, 5);
		collar(Kf, V(Kf, B), 5, 10);
		fill(blue);
		show(Fb, 5);
		collar(Ob, V(h, Nb), 1, 5);
		show(Kb, 5);
		collar(Kb, V(Kb, B), 5, 10);
		fill(orange);
		show(B, 10);
		show(T, 15);
		collar(B, V(B, T), 10, 15);
	}

	pt PP = P(); // picked point
	Boolean picking = false;

	public pt pick(int mX, int mY) { // returns point on visible surface at
										// pixel (mX,My)
		PGL pgl = beginPGL();
		FloatBuffer depthBuffer = ByteBuffer.allocateDirect(1 << 2).order(ByteOrder.nativeOrder()).asFloatBuffer();
		pgl.readPixels(mX, height - mY - 1, 1, 1, PGL.DEPTH_COMPONENT, PGL.FLOAT, depthBuffer);
		float depthValue = depthBuffer.get(0);
		depthBuffer.clear();
		endPGL();

		// get 3d matrices
		PGraphics3D p3d = (PGraphics3D) g;
		PMatrix3D proj = p3d.projection.get();
		PMatrix3D modelView = p3d.modelview.get();
		PMatrix3D modelViewProjInv = proj;
		modelViewProjInv.apply(modelView);
		modelViewProjInv.invert();

		float[] viewport = { 0, 0, p3d.width, p3d.height };
		float[] normalized = new float[4];
		normalized[0] = ((mX - viewport[0]) / viewport[2]) * 2.0f - 1.0f;
		normalized[1] = ((height - mY - viewport[1]) / viewport[3]) * 2.0f - 1.0f;
		normalized[2] = depthValue * 2.0f - 1.0f;
		normalized[3] = 1.0f;

		float[] unprojected = new float[4];

		modelViewProjInv.mult(normalized, unprojected);
		return P(unprojected[0] / unprojected[3], unprojected[1] / unprojected[3], unprojected[2] / unprojected[3]);
	}

	public pt pick(float mX, float mY, float mZ) {
		// get 3d matrices
		PGraphics3D p3d = (PGraphics3D) g;
		PMatrix3D proj = p3d.projection.get();
		PMatrix3D modelView = p3d.modelview.get();
		PMatrix3D modelViewProjInv = proj;
		modelViewProjInv.apply(modelView);
		modelViewProjInv.invert();
		float[] viewport = { 0, 0, p3d.width, p3d.height };
		float[] normalized = new float[4];
		normalized[0] = ((mX - viewport[0]) / viewport[2]) * 2.0f - 1.0f;
		normalized[1] = ((height - mY - viewport[1]) / viewport[3]) * 2.0f - 1.0f;
		normalized[2] = mZ * 2.0f - 1.0f;
		normalized[3] = 1.0f;
		float[] unprojected = new float[4];
		modelViewProjInv.mult(normalized, unprojected);
		return P(unprojected[0] / unprojected[3], unprojected[1] / unprojected[3], unprojected[2] / unprojected[3]);
	}

	public pt viewPoint() {
		return pick(0, 0, (height / 2) / tan(PI / 6));
	}

	/*
	 * in draw, before popMatrix, insert
	 * 
	 * if(picking) {PP = pick( mouseX, mouseY ); picking=false;} else
	 * {fill(yellow); show(PP,3);}
	 * 
	 * in keyPressed,
	 * 
	 * if(key=='`') picking=true;
	 * 
	 */
	int pp = 1; // index of picked vertex
	pts P = new pts(); // polyloop in 3D
	pts Q = new pts(); // second polyloop in 3D
	pts PtQ = new pts(); // inbetweening polyloop L(P,t,Q);

	class pts { // class for manipulaitng and sisplaying polyloops
		Boolean loop = true;
		int pv = 0, // picked vertex index,
				iv = 0, // insertion vertex index
				nv = 0; // number of vertices currently used in P
		int maxnv = 16000; // max number of vertices
		pt[] G = new pt[maxnv]; // geometry table (vertices)

		pts() {
		}

		public pts declare() {
			for (int i = 0; i < maxnv; i++)
				G[i] = P();
			return this;
		} // init all point objects

		public pts empty() {
			nv = 0;
			pv = 0;
			return this;
		} // resets P so that we can start adding points

		public pts addPt(pt P) {
			G[nv].setTo(P);
			pv = nv;
			nv++;
			return this;
		} // adds a point at the end

		public pts addPt(float x, float y) {
			G[nv].x = x;
			G[nv].y = y;
			pv = nv;
			nv++;
			return this;
		}

		public pts copyFrom(pts Q) {
			empty();
			nv = Q.nv;
			for (int v = 0; v < nv; v++)
				G[v] = P(Q.G[v]);
			return this;
		}

		public pts setToL(pts P, float t, pts Q) { // lerp (linear interpolation
													// betwen P and Q
			empty();
			nv = min(P.nv, Q.nv);
			for (int v = 0; v < nv; v++)
				G[v] = L(P.G[v], t, Q.G[v]);
			return this;
		}

		public pts resetOnCircle(int k, float r) { // makes new polyloo[p with k
													// points on a circle around
													// origin
			empty(); // resert P
			pt C = P(); // center of circle
			for (int i = 0; i < k; i++)
				addPt(R(P(C, V(0, -r, 0)), 2.f * PI * i / k, C)); // points on
																	// z=0 plane
			pv = 0; // picked vertex ID is set to 0
			return this;
		}

		public int idOfVertexWithClosestScreenProjectionTo(pt M) { // for
																	// picking a
																	// vertex
																	// with the
																	// mouse
			pp = 0;
			for (int i = 1; i < nv; i++)
				if (d(M, ToScreen(G[i])) <= d(M, ToScreen(G[pp])))
					pp = i;
			return pp;
		}

		public pt closestProjectionOf(pt M) { // for picking inserting O.
												// Returns projection but also
												// CHANGES iv !!!!
			pt C = P(G[0]);
			float d = d(M, C);
			for (int i = 1; i < nv; i++)
				if (d(M, G[i]) <= d) {
					iv = i;
					C = P(G[i]);
					d = d(M, C);
				}
			for (int i = nv - 1, j = 0; j < nv; i = j++) {
				pt A = G[i], B = G[j];
				if (projectsBetween(M, A, B) && disToLine(M, A, B) < d) {
					d = disToLine(M, A, B);
					iv = i;
					C = projectionOnLine(M, A, B);
				}
			}
			return C;
		}

		public pts insertPt(pt P) { // inserts new vertex after vertex with ID
									// iv
			for (int v = nv - 1; v > iv; v--)
				G[v + 1].setTo(G[v]);
			iv++;
			G[iv].setTo(P);
			nv++; // increments vertex count
			return this;
		}

		public pts insertClosestProjection(pt M) {
			pt P = closestProjectionOf(M); // also sets iv
			insertPt(P);
			return this;
		}

		public pts deletePicked() {
			for (int i = pv; i < nv; i++)
				G[i].setTo(G[i + 1]);
			pv = max(0, pv - 1);
			nv--;
			return this;
		}

		public pts setPt(pt P, int i) {
			G[i].setTo(P);
			return this;
		}

		public pts showPicked() {
			show(G[pv], 13);
			return this;
		}

		public pts drawBalls(float r) {
			for (int v = 0; v < nv; v++)
				show(G[v], r);
			return this;
		}

		public pts showPicked(float r) {
			show(G[pv], r);
			return this;
		}

		public pts drawClosedCurve(float r) {
			for (int v = 0; v < nv - 1; v++)
				stub(G[v], V(G[v], G[v + 1]), r, r / 2);
			stub(G[nv - 1], V(G[nv - 1], G[0]), r, r / 2);
			return this;
		}

		public pts setPickedTo(int pp) {
			pv = pp;
			return this;
		}

		public pts movePicked(vec V) {
			G[pv].add(V);
			return this;
		} // moves selected point (index p) by amount mouse moved recently

		public pts moveAll(vec V) {
			for (int i = 0; i < nv; i++)
				G[i].add(V);
			return this;
		};

		public pt Picked() {
			return G[pv];
		}

		public void savePts(String fn) {
			String[] inppts = new String[nv + 1];
			int s = 0;
			inppts[s++] = str(nv);
			for (int i = 0; i < nv; i++) {
				inppts[s++] = str(G[i].x) + "," + str(G[i].y) + "," + str(G[i].z);
			}
			saveStrings(fn, inppts);
		};

		public void loadPts(String fn) {
			println("loading: " + fn);
			String[] ss = loadStrings(fn);
			String subpts;
			int s = 0;
			int comma, comma1, comma2;
			float x, y;
			int a, b, c;
			nv = PApplet.parseInt(ss[s++]);
			print("nv=" + nv);
			for (int k = 0; k < nv; k++) {
				int i = k + s;
				float[] xy = PApplet.parseFloat(split(ss[i], ","));
				G[k].setTo(xy[0], xy[1], xy[2]);
			}
			pv = 0;
		};

	} // end of pts class
		// points, vectors, frames in 3D

	class vec {
		float x = 0, y = 0, z = 0;

		vec() {
		};

		vec(float px, float py, float pz) {
			x = px;
			y = py;
			z = pz;
		};

		vec(float px, float py) {
			x = px;
			y = py;
		};

		public vec set(float px, float py, float pz) {
			x = px;
			y = py;
			z = pz;
			return this;
		};

		public vec setTo(vec V) {
			x = V.x;
			y = V.y;
			z = V.z;
			return this;
		};

		public vec set(vec V) {
			x = V.x;
			y = V.y;
			z = V.z;
			return this;
		};

		public vec add(vec V) {
			x += V.x;
			y += V.y;
			z += V.z;
			return this;
		};

		public vec add(float s, vec V) {
			x += s * V.x;
			y += s * V.y;
			z += s * V.z;
			return this;
		};

		public vec sub(vec V) {
			x -= V.x;
			y -= V.y;
			z -= V.z;
			return this;
		};

		public vec mul(float f) {
			x *= f;
			y *= f;
			z *= f;
			return this;
		};

		public vec div(float f) {
			x /= f;
			y /= f;
			z /= f;
			return this;
		};

		public vec div(int f) {
			x /= f;
			y /= f;
			z /= f;
			return this;
		};

		public vec rev() {
			x = -x;
			y = -y;
			z = -z;
			return this;
		};

		public float norm() {
			return (sqrt(sq(x) + sq(y) + sq(z)));
		};

		public vec normalize() {
			float n = norm();
			if (n > 0.000001f) {
				div(n);
			}
			;
			return this;
		};

		public vec rotate(float a, vec I, vec J) { // Rotate this by angle a
													// parallel in plane (I,J)
													// Assumes I and J are
													// orthogonal
			float x = d(this, I), y = d(this, J); // dot products
			float c = cos(a), s = sin(a);
			add(x * c - x - y * s, I);
			add(x * s + y * c - y, J);
			return this;
		};
	} // end class vec

	class pt {
		float x = 0, y = 0, z = 0;

		pt() {
		};

		pt(float px, float py) {
			x = px;
			y = py;
		};

		pt(float px, float py, float pz) {
			x = px;
			y = py;
			z = pz;
		};

		public pt set(float px, float py, float pz) {
			x = px;
			y = py;
			z = pz;
			return this;
		};

		public pt set(pt P) {
			x = P.x;
			y = P.y;
			z = P.z;
			return this;
		};

		public pt setTo(pt P) {
			x = P.x;
			y = P.y;
			z = P.z;
			return this;
		};

		public pt setTo(float px, float py, float pz) {
			x = px;
			y = py;
			z = pz;
			return this;
		};

		public pt add(pt P) {
			x += P.x;
			y += P.y;
			z += P.z;
			return this;
		};

		public pt add(vec V) {
			x += V.x;
			y += V.y;
			z += V.z;
			return this;
		};

		public pt sub(vec V) {
			x -= V.x;
			y -= V.y;
			z -= V.z;
			return this;
		};

		public pt add(float s, vec V) {
			x += s * V.x;
			y += s * V.y;
			z += s * V.z;
			return this;
		};

		public pt sub(pt P) {
			x -= P.x;
			y -= P.y;
			z -= P.z;
			return this;
		};

		public pt mul(float f) {
			x *= f;
			y *= f;
			z *= f;
			return this;
		};

		public pt div(float f) {
			x /= f;
			y /= f;
			z /= f;
			return this;
		};

		public pt div(int f) {
			x /= f;
			y /= f;
			z /= f;
			return this;
		};
	}

	// ===== vector functions
	public vec V() {
		return new vec();
	}; // make vector (x,y,z)

	public vec V(float x, float y, float z) {
		return new vec(x, y, z);
	}; // make vector (x,y,z)

	public vec V(vec V) {
		return new vec(V.x, V.y, V.z);
	}; // make copy of vector V

	public vec A(vec A, vec B) {
		return new vec(A.x + B.x, A.y + B.y, A.z + B.z);
	}; // A+B

	public vec A(vec U, float s, vec V) {
		return V(U.x + s * V.x, U.y + s * V.y, U.z + s * V.z);
	}; // U+sV

	public vec M(vec U, vec V) {
		return V(U.x - V.x, U.y - V.y, U.z - V.z);
	}; // U-V

	public vec M(vec V) {
		return V(-V.x, -V.y, -V.z);
	}; // -V

	public vec V(vec A, vec B) {
		return new vec((A.x + B.x) / 2.0f, (A.y + B.y) / 2.0f, (A.z + B.z) / 2.0f);
	} // (A+B)/2

	public vec V(vec A, float s, vec B) {
		return new vec(A.x + s * (B.x - A.x), A.y + s * (B.y - A.y), A.z + s * (B.z - A.z));
	}; // (1-s)A+sB

	public vec V(vec A, vec B, vec C) {
		return new vec((A.x + B.x + C.x) / 3.0f, (A.y + B.y + C.y) / 3.0f, (A.z + B.z + C.z) / 3.0f);
	}; // (A+B+C)/3

	public vec V(vec A, vec B, vec C, vec D) {
		return V(V(A, B), V(C, D));
	}; // (A+B+C+D)/4

	public vec V(float s, vec A) {
		return new vec(s * A.x, s * A.y, s * A.z);
	}; // sA

	public vec V(float a, vec A, float b, vec B) {
		return A(V(a, A), V(b, B));
	} // aA+bB

	public vec V(float a, vec A, float b, vec B, float c, vec C) {
		return A(V(a, A, b, B), V(c, C));
	} // aA+bB+cC

	public vec V(pt P, pt Q) {
		return new vec(Q.x - P.x, Q.y - P.y, Q.z - P.z);
	}; // PQ

	public vec U(vec V) {
		float n = V.norm();
		vec temp = V(V);
		if (n < 0.0000001f)
			return V(0, 0, 0);
		else
			return temp.div(n);
	}; // V/||V||

	public vec U(pt P, pt Q) {
		return U(V(P, Q));
	}; // PQ/||PQ||

	public vec U(float x, float y, float z) {
		return U(V(x, y, z));
	}; // make vector (x,y,z)

	public vec N(vec U, vec V) {
		return V(U.y * V.z - U.z * V.y, U.z * V.x - U.x * V.z, U.x * V.y - U.y * V.x);
	}; // UxV cross product (normal to both)

	public vec N(pt A, pt B, pt C) {
		return N(V(A, B), V(A, C));
	}; // normal to triangle (A,B,C), not normalized (proportional to area)

	public vec B(vec U, vec V) {
		return U(N(N(U, V), U));
	}

	public vec Normal(vec V) {
		if (abs(V.z) <= min(abs(V.x), abs(V.y)))
			return V(-V.y, V.x, 0);
		if (abs(V.x) <= min(abs(V.z), abs(V.y)))
			return V(0, -V.z, V.y);
		return V(V.z, 0, -V.x);
	}

	// ===== point functions
	public pt P() {
		return new pt();
	}; // point (x,y,z)

	public pt P(float x, float y, float z) {
		return new pt(x, y, z);
	}; // point (x,y,z)

	public pt P(float x, float y) {
		return new pt(x, y);
	}; // make point (x,y)

	public pt P(pt A) {
		return new pt(A.x, A.y, A.z);
	}; // copy of point P

	public pt P(pt A, float s, pt B) {
		return new pt(A.x + s * (B.x - A.x), A.y + s * (B.y - A.y), A.z + s * (B.z - A.z));
	}; // A+sAB

	public pt L(pt A, float s, pt B) {
		return new pt(A.x + s * (B.x - A.x), A.y + s * (B.y - A.y), A.z + s * (B.z - A.z));
	}; // A+sAB

	public pt P(pt A, pt B) {
		return P((A.x + B.x) / 2.0f, (A.y + B.y) / 2.0f, (A.z + B.z) / 2.0f);
	} // (A+B)/2

	public pt P(pt A, pt B, pt C) {
		return new pt((A.x + B.x + C.x) / 3.0f, (A.y + B.y + C.y) / 3.0f, (A.z + B.z + C.z) / 3.0f);
	}; // (A+B+C)/3

	public pt P(pt A, pt B, pt C, pt D) {
		return P(P(A, B), P(C, D));
	}; // (A+B+C+D)/4

	public pt P(float s, pt A) {
		return new pt(s * A.x, s * A.y, s * A.z);
	}; // sA

	public pt A(pt A, pt B) {
		return new pt(A.x + B.x, A.y + B.y, A.z + B.z);
	}; // A+B

	public pt P(float a, pt A, float b, pt B) {
		return A(P(a, A), P(b, B));
	} // aA+bB

	public pt P(float a, pt A, float b, pt B, float c, pt C) {
		return A(P(a, A), P(b, B, c, C));
	} // aA+bB+cC

	public pt P(float a, pt A, float b, pt B, float c, pt C, float d, pt D) {
		return A(P(a, A, b, B), P(c, C, d, D));
	} // aA+bB+cC+dD

	public pt P(pt P, vec V) {
		return new pt(P.x + V.x, P.y + V.y, P.z + V.z);
	} // P+V

	public pt P(pt P, float s, vec V) {
		return new pt(P.x + s * V.x, P.y + s * V.y, P.z + s * V.z);
	} // P+sV

	public pt P(pt O, float x, vec I, float y, vec J) {
		return P(O.x + x * I.x + y * J.x, O.y + x * I.y + y * J.y, O.z + x * I.z + y * J.z);
	} // O+xI+yJ

	public pt P(pt O, float x, vec I, float y, vec J, float z, vec K) {
		return P(O.x + x * I.x + y * J.x + z * K.x, O.y + x * I.y + y * J.y + z * K.y,
				O.z + x * I.z + y * J.z + z * K.z);
	} // O+xI+yJ+kZ

	public void makePts(pt[] C) {
		for (int i = 0; i < C.length; i++)
			C[i] = P();
	}

	public pt ToScreen(pt P) {
		return P(screenX(P.x, P.y, P.z), screenY(P.x, P.y, P.z), 0);
	} // O+xI+yJ+kZ

	public pt ToModel(pt P) {
		return P(modelX(P.x, P.y, P.z), modelY(P.x, P.y, P.z), modelZ(P.x, P.y, P.z));
	} // O+xI+yJ+kZ

	public void showSphere(pt P, float r) {
		pushMatrix();
		translate(P.x, P.y, P.z);
		sphereDetail(3);
		sphere(r);
		popMatrix();
	}

	// ===== mouse
	public pt Mouse() {
		return P(mouseX, mouseY, 0);
	}; // current mouse location

	public pt Pmouse() {
		return P(pmouseX, pmouseY, 0);
	};

	public vec MouseDrag() {
		return V(mouseX - pmouseX, mouseY - pmouseY, 0);
	}; // vector representing recent mouse displacement

	public pt ScreenCenter() {
		return P(width / 2, height / 2);
	} // point in center of canvas

	// ===== measures
	public float d(vec U, vec V) {
		return U.x * V.x + U.y * V.y + U.z * V.z;
	}; // U*V dot product

	public float dot(vec U, vec V) {
		return U.x * V.x + U.y * V.y + U.z * V.z;
	}; // U*V dot product

	public float det2(vec U, vec V) {
		return -U.y * V.x + U.x * V.y;
	}; // U|V det product

	public float det3(vec U, vec V) {
		return sqrt(d(U, U) * d(V, V) - sq(d(U, V)));
	}; // U|V det product

	public float m(vec U, vec V, vec W) {
		return d(U, N(V, W));
	}; // (UxV)*W mixed product, determinant

	public float m(pt E, pt A, pt B, pt C) {
		return m(V(E, A), V(E, B), V(E, C));
	} // det (EA EB EC) is >0 when E sees (A,B,C) clockwise

	public float n2(vec V) {
		return sq(V.x) + sq(V.y) + sq(V.z);
	}; // V*V norm squared

	public float n(vec V) {
		return sqrt(n2(V));
	}; // ||V|| norm

	public float d(pt P, pt Q) {
		return sqrt(sq(Q.x - P.x) + sq(Q.y - P.y) + sq(Q.z - P.z));
	}; // ||AB|| distance

	public float area(pt A, pt B, pt C) {
		return n(N(A, B, C)) / 2;
	}; // area of triangle

	public float volume(pt A, pt B, pt C, pt D) {
		return m(V(A, B), V(A, C), V(A, D)) / 6;
	}; // volume of tet

	public boolean parallel(vec U, vec V) {
		return n(N(U, V)) < n(U) * n(V) * 0.00001f;
	} // true if U and V are almost parallel

	public float angle(vec U, vec V) {
		return acos(d(U, V) / n(V) / n(U));
	}; // angle(U,V)

	public boolean cw(vec U, vec V, vec W) {
		return m(U, V, W) > 0;
	}; // (UxV)*W>0 U,V,W are clockwise

	public boolean cw(pt A, pt B, pt C, pt D) {
		return volume(A, B, C, D) > 0;
	}; // tet is oriented so that A sees B, C, D clockwise

	public boolean projectsBetween(pt P, pt A, pt B) {
		return dot(V(A, P), V(A, B)) > 0 && dot(V(B, P), V(B, A)) > 0;
	};

	public float disToLine(pt P, pt A, pt B) {
		return det3(U(A, B), V(A, P));
	};

	public pt projectionOnLine(pt P, pt A, pt B) {
		return P(A, dot(V(A, B), V(A, P)) / dot(V(A, B), V(A, B)), V(A, B));
	}

	// ===== rotate
	public vec R(vec V) {
		return V(-V.y, V.x, V.z);
	} // rotated 90 degrees in XY plane

	public pt R(pt P, float a, vec I, vec J, pt G) {
		float x = d(V(G, P), I), y = d(V(G, P), J);
		float c = cos(a), s = sin(a);
		return P(P, x * c - x - y * s, I, x * s + y * c - y, J);
	}; // Rotated P by a around G in plane (I,J)

	public vec R(vec V, float a, vec I, vec J) {
		float x = d(V, I), y = d(V, J);
		float c = cos(a), s = sin(a);
		return A(V, V(x * c - x - y * s, I, x * s + y * c - y, J));
	}; // Rotated V by a parallel to plane (I,J)

	public pt R(pt Q, pt C, pt P, pt R) { // returns rotated version of Q by
											// angle(CP,CR) parallel to plane
											// (C,P,R)
		vec I0 = U(C, P), I1 = U(C, R), V = V(C, Q);
		float c = d(I0, I1), s = sqrt(1.f - sq(c));
		if (abs(s) < 0.00001f)
			return Q;
		vec J0 = V(1.f / s, I1, -c / s, I0);
		vec J1 = V(-s, I0, c, J0);
		float x = d(V, I0), y = d(V, J0);
		// stroke(red); show(C,400,I0); stroke(blue); show(C,400,I1);
		// stroke(orange); show(C,400,J0); stroke(magenta); show(C,400,J1);
		// noStroke();
		return P(Q, x, M(I1, I0), y, M(J1, J0));
	}

	public pt R(pt Q, float a) {
		float dx = Q.x, dy = Q.y, c = cos(a), s = sin(a);
		return P(c * dx + s * dy, -s * dx + c * dy, Q.z);
	}; // Q rotated by angle a around the origin

	public pt R(pt Q, float a, pt C) {
		float dx = Q.x - C.x, dy = Q.y - C.y, c = cos(a), s = sin(a);
		return P(C.x + c * dx - s * dy, C.y + s * dx + c * dy, Q.z);
	}; // Q rotated by angle a around point P

	// ===== render
	public void normal(vec V) {
		normal(V.x, V.y, V.z);
	}; // changes normal for smooth shading

	public void vertex(pt P) {
		vertex(P.x, P.y, P.z);
	}; // vertex for shading or drawing

	public void v(pt P) {
		vertex(P.x, P.y, P.z);
	}; // vertex for shading or drawing

	public void nv(vec N) {
		normal(N.x, N.y, N.z);
	}; // vertex for shading or drawing

	public void vTextured(pt P, float u, float v) {
		vertex(P.x, P.y, P.z, u, v);
	}; // vertex with texture coordinates

	public void show(pt P, pt Q) {
		line(Q.x, Q.y, Q.z, P.x, P.y, P.z);
	}; // draws edge (P,Q)

	public void show(pt P, vec V) {
		line(P.x, P.y, P.z, P.x + V.x, P.y + V.y, P.z + V.z);
	}; // shows edge from P to P+V

	public void show(pt P, float d, vec V) {
		fill(black);
		line(P.x, P.y, P.z, P.x + d * V.x, P.y + d * V.y, P.z + d * V.z);
	}; // shows edge from P to P+dV

	public void show(pt A, pt B, pt C) {
		beginShape();
		vertex(A);
		vertex(B);
		vertex(C);
		endShape(CLOSE);
	}; // volume of tet

	public void show(pt A, pt B, pt C, pt D) {
		beginShape();
		vertex(A);
		vertex(B);
		vertex(C);
		vertex(D);
		endShape(CLOSE);
	}; // volume of tet

	public void show(pt P, float r) {
		pushMatrix();
		translate(P.x, P.y, P.z);
		sphere(r);
		popMatrix();
	}; // render sphere of radius r and center P

	public void show(pt P, float s, vec I, vec J, vec K) {
		noStroke();
		fill(yellow);
		show(P, 5);
		stroke(red);
		show(P, s, I);
		stroke(green);
		show(P, s, J);
		stroke(blue);
		show(P, s, K);
	}; // render sphere of radius r and center P

	public void show(pt P, String s) {
		text(s, P.x, P.y, P.z);
	}; // prints string s in 3D at P

	public void show(pt P, String s, vec D) {
		text(s, P.x + D.x, P.y + D.y, P.z + D.z);
	}; // prints string s in 3D at P+D

	public void showShadow(pt P, float r) {
		pushMatrix();
		translate(P.x, P.y, 0);
		scale(1, 1, 0.01f);
		sphere(r);
		popMatrix();
	}

	public void show(FR F) {
		strokeWeight(5);
		fill(yellow);
		show(F.O, 3);
		show(F.O1, 3);
		stroke(red);
		show(F.O, F.I);
		stroke(green);
		show(F.O, F.J);
		stroke(blue);
		show(F.O, F.K);
	}// directly show frame

	public void showotherFrame(FR F) {
		strokeWeight(5);
		fill(yellow);
		show(F.O, 3);
		show(F.O1, 3);
		stroke(brown);
		show(F.O, F.I);
		stroke(grey);
		show(F.O, F.J);
		stroke(black);
		show(F.O, F.K);
	}

	public void drawball(FR F0) {
		pushMatrix();
		translate(F0.O.x, F0.O.y, F0.O.z);
		stroke(grey);
		noFill();
		sphere(n(F0.I) / 2);
		popMatrix();
	}

	class FR {
		pt O;
		pt O1;
		vec I;
		vec J;
		vec K;

		FR() {
			O = P();
			I = V(1, 0, 0);
			J = V(0, 0, 1);
			K = V(0, 1, 0);
		}

		FR(vec II, vec JJ, vec KK, pt OO) {
			I = V(II);
			J = V(JJ);
			K = V(KK);
			O = P(OO);
		}

		FR(pt A, pt B) {
			O = P(A);
			O1 = P(B);
			I = V(A, B);
			J = R(I);
			K = V(n(I), K = U(N(I, J)));
		}// make a point in A, create a vector from A to B. I turns right 90
			// degrees.

		FR(pt A, pt B, pt C, pt D) {
			O = P(A);
			I = V(A, B);
			J = V(A, C);
			K = V(A, D);
		}
	}

	public void show3Dframe(FR F) {
		show(F);
	}
	/*
	 * //normalized normal_N public vec calNormal(FR F1, FR F2) { return U(
	 * A(A(N( M(U(F2.I),U(F1.I)), M(U(F2.J),U(F1.J))),N( M(U(F2.J),U(F1.J)),
	 * M(U(F2.K),U(F1.K)))),N( M(U(F2.K),U(F1.K)), M(U(F2.I),U(F1.I)))) ); }
	 * 
	 * public pt spiralCenter(float a, float z, pt A, pt C) { float c=cos(a),
	 * s=sin(a); float D = sq(c*z-1)+sq(s*z); float ex = c*z*A.x - C.x -
	 * s*z*A.y; float ey = c*z*A.y - C.y + s*z*A.x; float x=(ex*(c*z-1) +
	 * ey*s*z) / D; float y=(ey*(c*z-1) - ex*s*z) / D; return P(x,y,0); }
	 * 
	 * 
	 * public pt spiral_3d_center(float a, float s, FR F0,FR F1) { vec
	 * FP0_P=V(dot(F0.I,normal_N)/dot(normal_N,normal_N),normal_N); vec
	 * FP0_V=M(F0.I, FP0_P);
	 * 
	 * pt local_F0_O=P(0,0);//local coordination origin vec
	 * i1=U(V(FP0_V));//local axis i vec j1=U(N(normal_N,i1));//local axis j vec
	 * k1=U(V(normal_N));//local axis k vec O0_O1=V(F0.O,F1.O); pt
	 * local_F1_O=P(dot(i1,O0_O1),dot(j1,O0_O1)); pt
	 * local_fixedpoint=spiralCenter( a, s, local_F0_O, local_F1_O); pt
	 * global_fixedpoint=P(F0.O, local_fixedpoint.x, i1, local_fixedpoint.y,
	 * j1);
	 * 
	 * float h=dot(O0_O1,k1); float z=h/(1-s); pt
	 * real_fixedpoint=P(global_fixedpoint, z, k1);
	 * 
	 * return real_fixedpoint;
	 * 
	 * }
	 * 
	 * 
	 * public float angle2d (vec U, vec V) {return atan2(det2(U,V),dot(U,V)); };
	 * public FR swirl_animatingFrame(FR F0,float t,FR F1) { float a =
	 * angle2d(F0.I,F1.I); float s = n(F1.I)/n(F0.I); System.out.println(
	 * "angle= "+a); float cos_v=cos(a*t), sin_v=sin(a*t);
	 * 
	 * fixedpoint=spiral_3d_center(a,s,F0,F1); show(fixedpoint,10);
	 * 
	 * pt new_O=calculate_new_point(F0.O,cos_v,sin_v,s,t); pt
	 * new_O1=calculate_new_point(F0.O1,cos_v,sin_v,s,t);
	 * 
	 * return new FR(new_O,new_O1);
	 * 
	 * } public pt calculate_new_point(pt p,float cos_v,float sin_v,float
	 * s,float t) { //rodrigues' rotation formula vec FP0=V(fixedpoint,p); vec
	 * part1=V(cos_v, FP0); vec part2=V(sin_v,N(normal_N, FP0)); vec
	 * part3=V(dot(normal_N,FP0)*(1-cos_v),normal_N); vec
	 * after_rotate=A(A(part1,part2),part3); vec new_vector=V(pow(s,t),
	 * after_rotate); pt new_O=P(fixedpoint, new_vector); return new_O;
	 * 
	 * }
	 * 
	 */

	public String toText(vec V) {
		return "(" + nf(V.x, 1, 5) + "," + nf(V.y, 1, 5) + "," + nf(V.z, 1, 5) + ")";
	}

	// ==== curve
	public void bezier(pt A, pt B, pt C, pt D) {
		bezier(A.x, A.y, A.z, B.x, B.y, B.z, C.x, C.y, C.z, D.x, D.y, D.z);
	} // draws a cubic Bezier curve with control points A, B, C, D

	public void bezier(pt[] C) {
		bezier(C[0], C[1], C[2], C[3]);
	} // draws a cubic Bezier curve with control points A, B, C, D

	public pt bezierPoint(pt[] C, float t) {
		return P(bezierPoint(C[0].x, C[1].x, C[2].x, C[3].x, t), bezierPoint(C[0].y, C[1].y, C[2].y, C[3].y, t),
				bezierPoint(C[0].z, C[1].z, C[2].z, C[3].z, t));
	}

	public vec bezierTangent(pt[] C, float t) {
		return V(bezierTangent(C[0].x, C[1].x, C[2].x, C[3].x, t), bezierTangent(C[0].y, C[1].y, C[2].y, C[3].y, t),
				bezierTangent(C[0].z, C[1].z, C[2].z, C[3].z, t));
	}

	public void PT(pt P0, vec T0, pt P1, vec T1) {
		float d = d(P0, P1) / 3;
		bezier(P0, P(P0, -d, U(T0)), P(P1, -d, U(T1)), P1);
	} // draws cubic Bezier interpolating (P0,T0) and (P1,T1)

	public void PTtoBezier(pt P0, vec T0, pt P1, vec T1, pt[] C) {
		float d = d(P0, P1) / 3;
		C[0].set(P0);
		C[1].set(P(P0, -d, U(T0)));
		C[2].set(P(P1, -d, U(T1)));
		C[3].set(P1);
	} // draws cubic Bezier interpolating (P0,T0) and (P1,T1)

	public vec vecToCubic(pt A, pt B, pt C, pt D, pt E) {
		return V((-A.x + 4 * B.x - 6 * C.x + 4 * D.x - E.x) / 6, (-A.y + 4 * B.y - 6 * C.y + 4 * D.y - E.y) / 6,
				(-A.z + 4 * B.z - 6 * C.z + 4 * D.z - E.z) / 6);
	}

	public vec vecToProp(pt B, pt C, pt D) {
		float cb = d(C, B);
		float cd = d(C, D);
		return V(C, P(B, cb / (cb + cd), D));
	};

	// ==== perspective
	public pt Pers(pt P, float d) {
		return P(d * P.x / (d + P.z), d * P.y / (d + P.z), d * P.z / (d + P.z));
	};

	public pt InverserPers(pt P, float d) {
		return P(d * P.x / (d - P.z), d * P.y / (d - P.z), d * P.z / (d - P.z));
	};

	// ==== intersection
	public boolean intersect(pt P, pt Q, pt A, pt B, pt C, pt X) {
		return intersect(P, V(P, Q), A, B, C, X);
	} // if (P,Q) intersects (A,B,C), return true and set X to the intersection
		// point

	public boolean intersect(pt E, vec T, pt A, pt B, pt C, pt X) { // if ray
																	// from E
																	// along T
																	// intersects
																	// triangle
																	// (A,B,C),
																	// return
																	// true and
																	// set X to
																	// the
																	// intersection
																	// point
		vec EA = V(E, A), EB = V(E, B), EC = V(E, C), AB = V(A, B), AC = V(A, C);
		boolean s = cw(EA, EB, EC), sA = cw(T, EB, EC), sB = cw(EA, T, EC), sC = cw(EA, EB, T);
		if ((s == sA) && (s == sB) && (s == sC))
			return false;
		float t = m(EA, AC, AB) / m(T, AC, AB);
		X.set(P(E, t, T));
		return true;
	}

	public boolean rayIntersectsTriangle(pt E, vec T, pt A, pt B, pt C) { // true
																			// if
																			// ray
																			// from
																			// E
																			// with
																			// direction
																			// T
																			// hits
																			// triangle
																			// (A,B,C)
		vec EA = V(E, A), EB = V(E, B), EC = V(E, C);
		boolean s = cw(EA, EB, EC), sA = cw(T, EB, EC), sB = cw(EA, T, EC), sC = cw(EA, EB, T);
		return (s == sA) && (s == sB) && (s == sC);
	};

	public boolean edgeIntersectsTriangle(pt P, pt Q, pt A, pt B, pt C) {
		vec PA = V(P, A), PQ = V(P, Q), PB = V(P, B), PC = V(P, C), QA = V(Q, A), QB = V(Q, B), QC = V(Q, C);
		boolean p = cw(PA, PB, PC), q = cw(QA, QB, QC), a = cw(PQ, PB, PC), b = cw(PA, PQ, PC), c = cw(PQ, PB, PQ);
		return (p != q) && (p == a) && (p == b) && (p == c);
	}

	public float rayParameterToIntersection(pt E, vec T, pt A, pt B, pt C) {
		vec AE = V(A, E), AB = V(A, B), AC = V(A, C);
		return -m(AE, AC, AB) / m(T, AC, AB);
	}

	public float angleDraggedAround(pt G) { // returns angle in 2D dragged by
											// the mouse around the screen
											// projection of G
		pt S = P(screenX(G.x, G.y, G.z), screenY(G.x, G.y, G.z), 0);
		vec T = V(S, Pmouse());
		vec U = V(S, Mouse());
		return atan2(d(R(U), T), d(U, T));
	}

	public float scaleDraggedFrom(pt G) {
		pt S = P(screenX(G.x, G.y, G.z), screenY(G.x, G.y, G.z), 0);
		return d(S, Mouse()) / d(S, Pmouse());
	}

	// FANS, CONES, AND ARROWS
	public void disk(pt P, vec V, float r) {
		vec I = U(Normal(V));
		vec J = U(N(I, V));
		disk(P, I, J, r);
	}

	public void disk(pt P, vec I, vec J, float r) {
		float da = TWO_PI / 36;
		beginShape(TRIANGLE_FAN);
		v(P);
		for (float a = 0; a <= TWO_PI + da; a += da)
			v(P(P, r * cos(a), I, r * sin(a), J));
		endShape();
	}

	public void fan(pt P, vec V, float r) {
		vec I = U(Normal(V));
		vec J = U(N(I, V));
		fan(P, V, I, J, r);
	}

	public void fan(pt P, vec V, vec I, vec J, float r) {
		float da = TWO_PI / 36;
		beginShape(TRIANGLE_FAN);
		v(P(P, V));
		for (float a = 0; a <= TWO_PI + da; a += da)
			v(P(P, r * cos(a), I, r * sin(a), J));
		endShape();
	}

	public void collar(pt P, vec V, float r, float rd) {
		vec I = U(Normal(V));
		vec J = U(N(I, V));
		collar(P, V, I, J, r, rd);
	}

	public void collar(pt P, vec V, vec I, vec J, float r, float rd) {
		float da = TWO_PI / 36;
		beginShape(QUAD_STRIP);
		for (float a = 0; a <= TWO_PI + da; a += da) {
			v(P(P, r * cos(a), I, r * sin(a), J, 0, V));
			v(P(P, rd * cos(a), I, rd * sin(a), J, 1, V));
		}
		endShape();
	}

	public void cone(pt P, vec V, float r) {
		fan(P, V, r);
		disk(P, V, r);
	}

	public void stub(pt P, vec V, float r, float rd) {
		collar(P, V, r, rd);
		disk(P, V, r);
		disk(P(P, V), V, rd);
	}

	public void arrow(pt P, vec V, float r) {
		stub(P, V(.8f, V), r * 2 / 3, r / 3);
		cone(P(P, V(.8f, V)), V(.2f, V), r);
	}

	// **************************** PRIMITIVE
	public void showFrame(float d) {
		noStroke();
		fill(metal);
		sphere(d / 10);
		fill(blue);
		showArrow(d, d / 10);
		fill(red);
		pushMatrix();
		rotateY(PI / 2);
		showArrow(d, d / 10);
		popMatrix();
		fill(green);
		pushMatrix();
		rotateX(-PI / 2);
		showArrow(d, d / 10);
		popMatrix();
	}

	public void showFan(float d, float r) {
		float da = TWO_PI / 36;
		beginShape(TRIANGLE_FAN);
		vertex(0, 0, d);
		for (float a = 0; a <= TWO_PI + da; a += da)
			vertex(r * cos(a), r * sin(a), 0);
		endShape();
	}

	public void showCollar(float d, float r, float rd) {
		float da = TWO_PI / 36;
		beginShape(QUAD_STRIP);
		for (float a = 0; a <= TWO_PI + da; a += da) {
			vertex(r * cos(a), r * sin(a), 0);
			vertex(rd * cos(a), rd * sin(a), d);
		}
		endShape();
	}

	public void showCone(float d, float r) {
		showFan(d, r);
		showFan(0, r);
	}

	public void showStub(float d, float r, float rd) {
		showCollar(d, r, rd);
		showFan(0, r);
		pushMatrix();
		translate(0, 0, d);
		showFan(0, rd);
		popMatrix();
	}

	public void showArrow() {
		showArrow(1, 0.08f);
	}

	public void showArrow(float d, float r) {
		float dd = d / 5;
		showStub(d - dd, r * 2 / 3, r / 3);
		pushMatrix();
		translate(0, 0, d - dd);
		showCone(dd, r);
		popMatrix();
	}

	public void showBlock(float w, float d, float h, float x, float y, float z, float a) {
		pushMatrix();
		translate(x, y, h / 2);
		rotateZ(TWO_PI * a);
		box(w, d, h);
		popMatrix();
	}

	// *********** PICK
	vec I = V(1, 0, 0), J = V(0, 1, 0), K = V(0, 0, 1); // screen projetions of
														// global model frame

	public void computeProjectedVectors() {
		pt O = ToScreen(P(0, 0, 0));
		pt A = ToScreen(P(1, 0, 0));
		pt B = ToScreen(P(0, 1, 0));
		pt C = ToScreen(P(0, 0, 1));
		I = V(O, A);
		J = V(O, B);
		K = V(O, C);
	}

	public vec ToIJ(vec V) {
		float x = det2(V, J) / det2(I, J);
		float y = det2(V, I) / det2(J, I);
		return V(x, y, 0);
	}

	public vec ToK(vec V) {
		float z = dot(V, K) / dot(K, K);
		return V(0, 0, z);
	}

	// ************************************ IMAGES & VIDEO
	int pictureCounter = 0, frameCounter = 0;
	Boolean filming = false, change = false;
	PImage myFace; // picture of author's face, should be: data/pic.jpg in
					// sketch folder
	PImage myFace2;
	
	public void snapPicture() {
		saveFrame("PICTURES/P" + nf(pictureCounter++, 3) + ".jpg");
	}

	// ******************************************COLORS
	int black = 0xff000000, white = 0xffFFFFFF, // set more colors using Menu >
												// Tools > Color Selector
			red = 0xffFF0000, green = 0xff00FF01, blue = 0xff0300FF, yellow = 0xffFEFF00, cyan = 0xff00FDFF,
			magenta = 0xffFF00FB, grey = 0xff818181, orange = 0xffFFA600, brown = 0xffB46005, metal = 0xffB5CCDE,
			dgreen = 0xff157901;

	public void pen(int c, float w) {
		stroke(c);
		strokeWeight(w);
	}

	// ******************************** TEXT , TITLE, and USER's GUIDE
	Boolean scribeText = true; // toggle for displaying of help text

	public void scribe(String S, float x, float y) {
		fill(0);
		text(S, x, y);
		noFill();
	} // writes on screen at (x,y) with current fill color

	public void scribeHeader(String S, int i) {
		fill(0);
		text(S, 10, 20 + i * 20);
		noFill();
	} // writes black at line i

	public void scribeHeaderRight(String S) {
		fill(0);
		text(S, width - 7.5f * S.length(), 20);
		noFill();
	} // writes black on screen top, right-aligned

	public void scribeFooter(String S, int i) {
		fill(0);
		text(S, 10, height - 10 - i * 20);
		noFill();
	} // writes black on screen at line i from bottom

	public void scribeAtMouse(String S) {
		fill(0);
		text(S, mouseX, mouseY);
		noFill();
	} // writes on screen near mouse

	public void scribeMouseCoordinates() {
		fill(black);
		text("(" + mouseX + "," + mouseY + ")", mouseX + 7, mouseY + 25);
		noFill();
	}

	// **************************** FILE SELECTION FOR SAVING AND LOADING MODELS
	// String fileName="data/points";

	// String path="data/pts";
	// void saveToFile(File selection) {
	// if (selection == null) println("Window was closed or the user hit
	// cancel.");
	// else path=selection.getAbsolutePath();
	// println(" save path = "+path);
	// }
	//
	// void readFromFile(File selection) {
	// if (selection == null) println("Window was closed or the user hit cancel
	// or file not found.");
	// else path=selection.getAbsolutePath();
	// println(" read path = "+path);
	// }
	//
	//
	// void fileSelected(File selection) {
	// if (selection == null) println("Window was closed or the user hit
	// cancel.");
	// else {
	// fileName = selection.getAbsolutePath();
	// println("User selected " + fileName);
	// }
	// }
	//
	public void settings() {
		size(900, 900, P3D);
		noSmooth();
	}

	static public void main(String[] passedArgs) {
		String[] appletArgs = new String[] { "swirl" };
		if (passedArgs != null) {
			PApplet.main(concat(appletArgs, passedArgs));
		} else {
			PApplet.main(appletArgs);
		}
	}
}
