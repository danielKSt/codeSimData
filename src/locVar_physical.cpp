#include <Rcpp.h>
#include <vector>

// [[Rcpp::plugins("cpp11")]]

class LocvarCalculatorPhysical {
private:
  int xDim, yDim, winDim, shellDim, divisor, divisorLeft;
  double xs, x2s, xsLeft, x2sLeft;
  std::vector<int> win;
  std::vector<int> shell;
  std::vector<double> sDiv;
  std::vector<double> locvar;
public:
  LocvarCalculatorPhysical(int inxDim, int inyDim, int inwinDim, int inShellDim, std::vector<int> inWin, std::vector<int> inShell, std::vector<double> insDiv);

  // The main shift functions
  void initialLocvarCalculation();
  void sideshift(int x, int y);
  void upshift(int x, int y);

  // Some sub-functions for shifting
  void sideshift_sub(int sx_curr, int sy_curr, int sy_next, int sy_prev, int x, int y);
  void upshift_sub(int sy_curr, int sx_curr, int sx_next, int sx_prev, int x, int y);

  // The export funciton
  Rcpp::NumericVector returnLocVar();

  // Some funcitons for debugging with r
  void rprintsdiv();
  void rprintlcov();
  void rptintwin();
  void rprintshell();
};


LocvarCalculatorPhysical::LocvarCalculatorPhysical(int inxDim, int inyDim, int inwinDim, int inShellDim, std::vector<int> inWin, std::vector<int> inShell, std::vector<double> insDiv){
  xDim = inxDim;
  yDim = inyDim;
  winDim = inwinDim;
  shellDim = inShellDim;
  win = inWin;
  shell = inShell;
  sDiv = insDiv;
  initialLocvarCalculation();
  for (int y = 0; y < yDim-1; y++)
  {
    for (int x = 1; x < xDim; x++)
    {
      sideshift(x,y);
    }
    upshift(0, y+1);
  }
  for (int x = 1; x < xDim; x++)
  {
    sideshift(x, yDim-1);
  }
}

void LocvarCalculatorPhysical::initialLocvarCalculation(){
  int xCoord, yCoord;
  xs = 0;
  x2s = 0;
  divisor = 0;
  for (int p = 0; p < winDim; p++)
  {
    xCoord = win[2*p];
    yCoord = win[2*p+1];

    bool isInWindow = (xCoord >= 0) && (yCoord >= 0) && (xCoord < xDim) && (yCoord < yDim);
    if (isInWindow)
    {
      xs += sDiv[xCoord+(xDim*yCoord)];
      x2s += sDiv[xCoord+(xDim*yCoord)]*sDiv[xCoord+(xDim*yCoord)];
      divisor += 1;
    }
  }
  xsLeft = xs;
  x2sLeft = x2s;
  divisorLeft = divisor;
  locvar.push_back((x2s - (xs*xs/divisor))/(divisor-1));
}

void LocvarCalculatorPhysical::sideshift(int x, int y){
  sideshift_sub(shell[0], shell[1], shell[3], shell[2*shellDim-1], x, y);
  for(int p = 1; p < shellDim-1; p++){
    sideshift_sub(shell[2*p], shell[2*p+1], shell[2*p+3], shell[2*p-1], x, y);
  }
  sideshift_sub(shell[2*shellDim-2], shell[2*shellDim-1], shell[1], shell[2*shellDim-3], x, y);
  locvar.push_back((x2s - (xs*xs/divisor))/(divisor-1));
}

void LocvarCalculatorPhysical::sideshift_sub(int sx_curr, int sy_curr, int sy_next, int sy_prev, int x, int y){
  int current_i = sx_curr+x;
  int current_j = sy_curr+y;
  int prev_i = sx_curr+x-1;

  bool nextInWindow, prevInWindow;
  nextInWindow = (0 <= sx_curr+x) && (sx_curr+x < xDim) && (0 <= sy_curr + y) && (sy_curr + y < yDim);
  prevInWindow = (0 < sx_curr+x) && (sx_curr+x <= xDim) && (0 <= sy_curr + y) && (sy_curr + y < yDim);

  bool addToXs, subtractFromXs;
  addToXs = (sx_curr >= 0) && (((sy_curr >= 0) && (sy_next < sy_curr)) || ((sy_curr < 0) && (sy_prev > sy_curr)));
  subtractFromXs = (sx_curr <= 0) && (((sy_curr <= 0) && (sy_curr < sy_next)) || ((sy_curr > 0) && (sy_prev < sy_curr)));

  if(addToXs && nextInWindow){
    xs += sDiv[current_i+xDim*current_j];
    x2s += sDiv[current_i+xDim*current_j]*sDiv[current_i+xDim*current_j];
    divisor += 1;
  }
  if(subtractFromXs && prevInWindow){
    xs -= sDiv[prev_i+xDim*current_j];
    x2s -= sDiv[prev_i+xDim*current_j]*sDiv[prev_i+xDim*current_j];
    divisor -= 1;
  }
}

void LocvarCalculatorPhysical::upshift(int x, int y){
  upshift_sub(shell[1], shell[0], shell[2], shell[2*shellDim-2], x, y);
  for(int p = 1; p < shellDim-1; p++){
    upshift_sub(shell[2*p+1], shell[2*p], shell[2*p+2], shell[2*p-2], x, y);
  }
  upshift_sub(shell[2*shellDim-1], shell[2*shellDim-2], shell[0], shell[2*shellDim-4], x, y);
  locvar.push_back((x2sLeft - (xsLeft*xsLeft/divisorLeft))/(divisorLeft-1));
  xs = xsLeft;
  x2s = x2sLeft;
  divisor = divisorLeft;
}

void LocvarCalculatorPhysical::upshift_sub(int sy_curr, int sx_curr, int sx_next, int sx_prev, int x, int y){
  int current_i = sx_curr+x;
  int current_j = sy_curr+y;
  int prev_j = sy_curr+y-1;

  bool nextInWindow, prevInWindow;
  nextInWindow = (0 <= sx_curr+x) && (sx_curr+x < xDim) && (sy_curr+y < yDim) && (0 <= sy_curr+y);
  prevInWindow = (0 <= sx_curr+x) && (sx_curr+x < xDim) && (sy_curr+y <= yDim) && (0 < sy_curr+y);

  bool addToXs, subtractFromXs;
  addToXs = (sy_curr >= 0) && (((sx_curr <= 0) && (sx_next > sx_curr)) || ((0 < sx_curr) && (sx_prev < sx_curr)));
  subtractFromXs = (sy_curr <= 0) && (((0 <= sx_curr) && (sx_next < sx_curr)) || ((sx_curr < 0) && (sx_curr < sx_prev)));

  if(addToXs && nextInWindow){
    xsLeft += sDiv[current_i+xDim*current_j];
    x2sLeft += sDiv[current_i+xDim*current_j]*sDiv[current_i+xDim*current_j];
    divisorLeft += 1;
  }
  if(subtractFromXs && prevInWindow){
    xsLeft -= sDiv[current_i+xDim*prev_j];
    x2sLeft -= sDiv[current_i+xDim*prev_j]*sDiv[current_i+xDim*prev_j];
    divisorLeft -= 1;
  }
}


Rcpp::NumericVector LocvarCalculatorPhysical::returnLocVar(){
  return Rcpp::wrap(locvar);
}

void LocvarCalculatorPhysical::rprintsdiv(){
  for (int y = 0; y < yDim; y++)
  {
    for (int x = 0; x < xDim; x++)
    {
      Rcpp::Rcout << " | " << sDiv[x+y*xDim];
    }
    Rcpp::Rcout << "\n";
  }
}

void LocvarCalculatorPhysical::rprintlcov(){
  for (int y = 0; y < yDim; y++)
  {
    for (int x = 0; x < xDim; x++)
    {
      Rcpp::Rcout << " | " << locvar[x+y*xDim];
    }
    Rcpp::Rcout << "\n";
  }
}

void LocvarCalculatorPhysical::rptintwin(){
  for (int p = 0; p < winDim; p++)
  {
    Rcpp::Rcout << win[2*p] << " | " << win[2*p+1] << "\n";
  }
  Rcpp::Rcout << "\n";
}

void LocvarCalculatorPhysical::rprintshell(){
  Rcpp::Rcout << "\n We show the values of the shell: \n\n";
  for (int p = 0; p < shellDim; p++)
  {
    Rcpp::Rcout << shell[2*p] << " | " << shell[2*p+1] << "\n";
  }
  Rcpp::Rcout << "\n";
}


//' Calculates the local variance with physical boundary conditions
//'
//' @param inWin Window
//' @param inShell Shell of window
//' @param insDiv field of which the local variance should be calculated
//' @param inDims Integer vector of dimensions for input, should be <xDim, yDim, winDim, shellDim>
//'
//' @return Vector with the local variance
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calculateLocVar_physical(Rcpp::NumericVector inWin, Rcpp::NumericVector inShell, Rcpp::NumericVector insDiv, Rcpp::IntegerVector inDims){
  // We extract the dimensions
  int xDim = inDims[0];
  int yDim = inDims[1];
  int winDim = inDims[2];
  int shellDim = inDims[3];

  // We extract the vectors
  std::vector<int> win = Rcpp::as< std::vector<int> >(inWin);
  std::vector<int> shell = Rcpp::as< std::vector<int> >(inShell);
  std::vector<double> sDiv = Rcpp::as< std::vector<double> >(insDiv);

  // We use the constructor
  LocvarCalculatorPhysical LocVarCalcInstance = LocvarCalculatorPhysical(xDim, yDim, winDim, shellDim, win, shell, sDiv);

  // Some tests for debugging
  //LocVarCalcInstance.rprintsdiv();
  //LocVarCalcInstance.rprintlcov();
  //LocVarCalcInstance.rptintwin();
  //LocVarCalcInstance.rprintshell();

  // We export the result
  // Rcpp::NumericVector res = LocVarCalcInstance.returnLocVar();
  // LocVarCalcInstance.~LocvarCalculatorPhysical();
  // return res;
return LocVarCalcInstance.returnLocVar();
}
