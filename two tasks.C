//
//
//

//
// We need a routine to calculate derivatives.
// It will take as an argument the number of ordinary diff. eqs.
// the value at which we want to evaluate the derivatives.
// An array "y" to hold the nODE functions in our system of Diff. eqs.
// An array "f" to hold the nODE derivatives of the functions. 
//
#include <algorithm>
#include <TMarker.h>
#include <TLine.h>
#include <TArrow.h>
void derivative(double x, double y[], double f[]) {
  // Note, in general, the derivative can depend on the independent variable, x.
  // For Free-Fall, it doesn't (where time is the independent variable).
  // Therefore, the independent variable is not used in the calculation below.
  // The position of the object will be stored in y[0].
  // The velocity of the object will be stored in y[1].
  // The derivatives of each will be stored in the f arrays.
  
  double g=9.81;
  double mass=65;
  double b = 7;//pick b value 
  double k = 25; // 
 
  f[0]=y[1]; 
  if (y[0] < 94)//when spring force is added
  {
      f[1] = -(g)-b * y[1] / mass+k* (134-y[0]-40)/mass;// chose L=40
      f[2] = -b * f[1] / mass - k * f[0] / mass;// a'=(-ba-kv)/m

  }
  else {
      
      f[1] = -(g)-b * y[1] / mass;
      f[2] = -b * f[1] / mass;
  }
 
  return;
}

// Based on the algorithm of Klein-Godunov p. 75, with bugs fixed.
void rk4_step(double* y, double* ydiff, double* y_old, int nODE, double step, double x0) {
  // Definition of arguments:
  // y : array of the calculated fifth order integral for the step (output of routine).
  // ydiff : array of the difference between the fifth and sixth orders (output of routine).
  // y_old : array of starting values for this step (input to routine).
  // x0 : the current value of the independent variable for this step.
  
  // Setup of constans for the RK method.
  double const a2(1./5.), a3(3./10), a4(3./5.), a5(1.0), a6(7./8.);
  double const b21(1./5.);
  double const b31(3./40.),       b32(+9./40.);
  double const b41(3./10.),       b42(-9./10.),   b43(6./5.);
  double const b51(-11./54.),     b52(5./2.),     b53(-70./27),     b54(-35./27.);
  double const b61(1631./55296.), b62(175./512.), b63(575./13824.), b64(44275./110592.), b65(253./4096.);

  // coefficients for calculating y(nth-order) (i.e. c's in Eq.8.83)
  double const c[6] = {37./378., 0, 250./621., 125./594., 0, 512./1771.};
  // coefficients for calculating y(n-1-th order) (i.e. c*'s in Eq. 8.84)
  double const cs[6] = {2825./27648., 0, 18575./48384., 13525./55296., 277./14336., 1./4.};
  
  // Order of Differential equation = number of 1st order equations
  // array to store the values of the derivatives
  // Again, since we don't know how many Diff. Eqs. we will need,
  // best is to use a vector.  So we will use a 6 dimensional array
  // for the 6 different calculations, and each element of the array
  // will have a vector of nODE doubles, the first element of the vector
  // for the dependent variable, and the other elements for the
  // derivatives.
  // We could do it with a call to vector::resize(n), but 
  // root doesn't support it. Will have to do it
  // by 
  vector<double> f[6];
  for (int iFunc=0;iFunc<6;++iFunc) {
    for (int iN=0; iN<nODE; ++iN) {
      f[iFunc].push_back(0.0);
    }
    // at this point the size of each vector in f[iFunc] should be 2
    //cout << "Size of vector #" << iFunc << " is " << f[iFunc].size() << endl;
  }
  // 1st step
  // With x0, y_old, calculate derivatives f[0]
  derivative(x0,y_old,&f[0][0]);
  // With the f0 derivatives, calculate new y's.
  for (int i1=0;i1<nODE;++i1) {
    y[i1]=y_old[i1]+b21*step*f[0][i1];
  }
  //cout << "Step 1, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  double x = x0 + a2*step;
  
  // 2nd step
  // With the x, and y's just calculated, calculate derivatives f[1]
  derivative(x,y,&f[1][0]);
  // With the f0 and f1, derivatives, calculate new y's.
  for (int i2=0;i2<nODE;++i2) {
    y[i2]=y_old[i2]+b31*step*f[0][i2]+b32*step*f[1][i2];    
  }
  //cout << "Step 2, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  x = x0 + a3*step;
  
  //3d step
  // With the x, and y's just calculated, calculate derivatives f[2]
  derivative(x,y,&f[2][0]);
  // With the f0, f1, and f2 derivatives, calculate new y's.
  for (int i3=0;i3<nODE;++i3) {
    y[i3]=y_old[i3]+b41*step*f[0][i3]+b42*step*f[1][i3]+b43*step*f[2][i3];    
  }
  //cout << "Step 3, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  x = x0 + a4*step;

  //4th step  double b_wind = 7;

  // With the x, and y's just calculated, calculate derivatives f[3]
  derivative(x,y,&f[3][0]);
  // With the f0, f1, ,f2 and f3 derivatives, calculate new y's.
  for (int i4=0;i4<nODE;++i4) {
    y[i4]=y_old[i4]+b51*step*f[0][i4]+b52*step*f[1][i4]+b53*step*f[2][i4]+b54*step*f[3][i4];    
  }
  //cout << "Step 4, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  x = x0 + a5*step;

  //5th step
  // With the x, and y's just calculated, calculate derivatives f[4]
  derivative(x,y,&f[4][0]);
  // With the f0, f1, ,f2, f3 and f4 derivatives, calculate new y's.
  for (int i5=0;i5<nODE;++i5) {
    y[i5]=y_old[i5]+b61*step*f[0][i5]+b62*step*f[1][i5]+b63*step*f[2][i5]+b64*step*f[3][i5]+b65*step*f[4][i5];    
  }
  //cout << "Step 5, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  x = x0 + a6*step;

  //6th step
  // With the x, and y's just calculated, calculate derivatives f[5]
  derivative(x,y,&f[5][0]);
  // With the f0, f1, ,f2, f3, f4 and f5 derivatives, calculate final y's. 
  // There will be two calculations, one using the c's and one using the c*'s.
  
  // create a vector to hold the 2nd solution
  vector<double> yp(nODE);
  for (int i6=0;i6<nODE;++i6) {
    // first add the old values to both solutions.
    y[i6]=y_old[i6];
    yp[i6]=y_old[i6];
    // next add the sum of the c[i]*f[i] for all 6 products
    // do the same for the solution using c*'s.
    for (int j1=0; j1<6; ++j1) {
      y[i6]+=step*c[j1]*f[j1][i6];
      yp[i6]+=step*cs[j1]*f[j1][i6];
    }
    // then calculate the difference between the solutions.
    ydiff[i6]=yp[i6]-y[i6];
  }
  //cout << "Step 6, y and y': " << y[0] << ", " << y[1] << endl;

  return;
}

void rk4_stepper(double* y_old, const int nODE, double xstart, double xmax, double hstep, double eps, int& nxstep) {
  // Arguments:
  // y_old, starting values of the dependent variable, and the derivatives of the dependent
  //        variable.
  // nODE: degree of the diff. eq. For the SHO, nODE=2.
  // xtart: starting value of the independent variable.
  // xmax:  ending value of the independent variable
  // hstep: beginning step size
  // eps  : chosen precision
  // nxstep: number of steps taken
  cout << "nODE   : " << nODE <<  endl;
  cout << "xstart : " << xstart << endl;
  cout << "xmax   : " << xmax   << endl;
  cout << "hstep  : " << hstep << endl;
  cout << "eps    : " << eps   << endl;
  
  
  double heps=hstep*eps;
  double yerrmax=0.99; // max error in any given step
  double redPow = 1./5.; //power used in formula to reduce the step size, eq. 8.79
  double esmall = eps/100.; // lower limit of precision, if the difference is smaller, increase step size.
  // For the y and ydiff arrays, since the size of the array is given as a parameter, 
  // a static array is not the best solution.  It is better to use a vector<double>
  // which can have a size given at runtime.  We can still handle it just like an array.
  vector<double> ydiff(nODE); // store the difference between the two R-K methods.
  vector<double> y(nODE);  // store the result of the calculation of the new dep. variable at each step.
  double hnew;
  double x0 = xstart; // value of the independent variable at each step. Initialize to xstart.

  double step_lolim = hstep/1000.; // lower limit for step size
  double step_hilim = hstep*10.; // upper limit for step size
  
  vector<double> xVals,yVals[3]; // vector of values of the indep and dep variables.  These are the values that will go into drawing the graphs of y(x) and y'(x), etc.
  
  // First, store the starting values.
  xVals.push_back(x0);
  for (int i=0; i<nODE; ++i) {
    yVals[i].push_back(y_old[i]); // store into final vector of results
    y[i] = y_old[i]; // store into array for current step, maybe not needed?
  }
  // Next, we enter into a loop and call the helper function rk4_step
  // to find the next value.  If the difference is big, we reduce the 
  // step size and try again, until the difference is smaller than our expected
  // precision.
  while (x0<xmax) {
    yerrmax=0.99; //reset max error for each step in the loop
    while (yerrmax<=0.99) { // when error is of this size, we found the necessary step size.
      //cout << "Current Step size " << hstep << endl;
      //cout << "Current Values : " << endl;
      //cout << "x = " << x0 << ", y= " << y[0] << ", y'= " << y[1] << endl; 
      // calculate the new values using rk4 for the current step size.
      // obtain the new y's and the difference between the 5th and 6th order rk values
      rk4_step(&y[0], &ydiff[0], &y_old[0], nODE, hstep, x0);
      
      // compare the difference to the step*chosen error for the function and the
      // derivatives. If the difference is large, the ratio heps/ydiff will be small
      // so it will not be greater than yerrmax.
      for (int j=0;j<nODE;++j) {
	yerrmax=max(yerrmax,fabs(heps/ydiff[j]));
	//cout << "Error for Eq# " << j << " : " << fabs(heps/ydiff[j]) << endl;
      }
      // at this point, yerrmax should be 0.99 if ydiff is big (need to reduce step size)
      // or greater than 1 if ydiff is small (we're done)
      if (yerrmax==0.99) {
	// reduce the step size.
	hstep=0.5*hstep*pow(fabs(heps/ydiff[0]),redPow);
	//heps=hstep*eps;
	cout << "Reducing step size to  " << hstep << endl;
      }
      
      // We will repeat the calculation, unless the step size is too small.
      // Protect against it.
      if (hstep<step_lolim) {
	cout << "rk4_stepper: step size requested is smaller than lowest stepsize allowed" << endl;
	cout << "Fix by trying a lower step size of lowering the step_lolim" << endl;
	cout << "Exiting..." << endl;
	return;
      }
    }// loop until yerrmax>0.99
    
    // If in last step the difference is very small, we could also increase the step size.
    // If ydiff is too small, then yerrmax will be large, we can compare this to
    // 1/esmall, our lower limit of precision, which should tell us how large should
    // yerrmax is allowed to be with our given precision.  If ydiff is too small,
    // yerrmax will be much bigger than 1/esmall, so we can use this as a criterion
    // to increase the step size. Also check that we don't go beyond an upper limit.
    
    if (yerrmax>1./esmall) hstep*=2.;
    if (hstep>step_hilim) hstep=step_hilim;
    
    // Store new data into vector. Set the new values y[i] into the old values y_old
    // so that we can call the routine again.
    // increase the counter for the number of steps we've taken.
    
    for (int k=0;k<nODE; ++k) {
      y_old[k]=y[k];
      yVals[k].push_back(y_old[k]); // store into final vector of results for dep. variable.
    }
    x0+=hstep; // propagate the indep. varialbe by one step.
    xVals.push_back(x0); // store into vector of indep variable values.
    nxstep++;
    //cout << "Step " << nxstep << "------" << endl;
    //cout << "x = " << x0 << ", y= " << y_old[0] << ", diff= " << ydiff[0] << endl;
  }// loop from x0=xstart, until x0>xmax
  cout << "Total steps " << nxstep << endl;
  cout << "Size of x vector " << xVals.size() << endl;
  cout << "Size of y vector " << yVals[0].size() << endl;
 
  // Make TGraphs out of the output vectors.
  TGraph* yGraph = new TGraph(xVals.size(),&xVals[0],&(yVals[0][0])); //position vs. time
  TGraph* ypGraph = new TGraph(xVals.size(),&xVals[0],&(yVals[1][0])); //velocity vs. time
  TGraph* yppGraph = new TGraph(xVals.size(), &xVals[0], &(yVals[2][0]));
  // Make the graphs look better, add color
  yGraph->SetLineColor(2);
  ypGraph->SetLineColor(4);
  yGraph->SetMinimum(-100);
  yGraph->SetMaximum(+150);
  // Plot them in a canvas
  TCanvas* freeFallAdapStep = new TCanvas("freeFallAdapStep","Adaptive Step Size RK",500,500);
  yGraph->Draw("AL");
  ypGraph->Draw("L");
  yppGraph->Draw("L");
  yGraph->GetHistogram()->SetXTitle("t (sec)");
  ypGraph->GetHistogram()->SetYTitle("y (m), v (m/s)");
  double tmin = 0.0;
  double tmax = 45;
      double step = 0.010; // 10 ms time step, which matches the time delay of the gif file.
  int Nsteps = static_cast<int>((tmax - tmin) / step);
  const int maxPoints = 4501;
  double vVal[maxPoints];
  double xVal[maxPoints];
  double aVal[maxPoints];
  double tVal[maxPoints];
  tVal[0] = 0.0;
  for (int i = 0; i < yVals[0].size(); ++i) {
      tVal[i] = xVals[i];
      vVal[i] = yVals[1][i]; 
      xVal[i] = yVals[0][i]; 
      aVal[i] = yVals[2][i];
  }
  TCanvas* ygraph = new TCanvas("ygraph", "position Animation", 500, 500);
 
 
  TH1D* dummy = new TH1D("dummy", "ygraph", 1000, 0, 45);
  dummy->SetMinimum(-40);
  dummy->SetMaximum(150);
  gStyle->SetOptStat(0); //Don't need a stat box in the canvas
  // Create the marker we are going to animate
  TMarker* massMarker = new TMarker(0, 0, 21);
  massMarker->SetMarkerColor(kBlue);
  massMarker->SetMarkerSize(1); 
  TLine* lineAnim = new TLine(0, 134, 0, 134);
  lineAnim->SetLineWidth(2);
  TArrow* vgraph = new TArrow(0, 134, 0, 134);
      TArrow * agraph = new TArrow(0, 134, 0, 134);
  char timeMessage[50];
  TLatex* ltx = new TLatex();
  // Use TLatex to draw a message on the canvas stating the value of t for each frame
  

  // Loop over the values of the x-coordinate and draw the marker at that x-value
  // and draw the object in the canvas, frame by frame
  for (size_t iFrame = 0; iFrame < 1000; iFrame += 10) {
      // clear the canvas in each frame
      gPad->Clear();

      // Draw the coordinate system.
      dummy->Draw();

      // Draw the marker at the right location
      massMarker->DrawMarker(0, xVal[iFrame]);
      lineAnim->DrawLine(0, 134, 0, xVal[iFrame]);
      vgraph->DrawArrow(tVal[iFrame], vVal[iFrame], tVal[iFrame+1],vVal[iFrame + 1],0.1);
      agraph->DrawArrow(tVal[iFrame], aVal[iFrame], tVal[iFrame+1], aVal[iFrame + 1],0.05);
    
      ypGraph->Draw("L");
      yppGraph->Draw("L");
      sprintf(timeMessage, "t=%2.1f", tVal[iFrame]);
      ltx->DrawLatex(-1, 1, timeMessage);
      ygraph->Modified();
      ygraph->Update();
      ygraph->SaveAs("compundgraph.gif+1");
  }
  
  return;
} 

void freeFallAdaptiveStep() {
  
  // Set up parameters for running the Runge-Kutta stepper
  // routine.
  
  //initial position and initial velocity, hence 2 differential equations
  const int numberOfDiffEqs(3);

  double y_old[numberOfDiffEqs]; // set up the initial conditions
  // i.e. the starting position, y_old[0] should be set to 134, and the starting velocity y_old[1] should be set to 0, and acceleration is -9.8.
  y_old[0]=134; y_old[1]=0; 
  y_old[2] = -9.8;
  
  // Initial and final values of the indep. variable (it is really time, but in the text, the independent variable is labeled x, so that is what we use here)
  double xstart = 0;
  double xmax = 45;
  
  // Step size and desired precision eps (given as percent of step size).
  // i.e. 0.1=10% of hstep will be used to compare against the difference
  // in the two solutions given by the 5th and 6th order r-k.

  double hstep = 0.01;
  double eps = 0.1; 
  int NumberOfSteps;
  rk4_stepper(y_old, numberOfDiffEqs, xstart, xmax, hstep, eps, NumberOfSteps);
  
  return;
}
