#include <Rcpp.h>
#include <vector>

// [[Rcpp::plugins("cpp11")]]

int myModulo(int a, int b){
  if ((a % b) >= 0)
  {
    return a % b;
  }
  else
  {
    return (a % b) + b;
  }

}

class LocvarCalculator {
private:
  int xDim, yDim, winDim, shellDim;
  double xs, x2s;
  std::vector<int> win;
  std::vector<int> shell;
  std::vector<double> sDiv;
  std::vector<double> locvar;
public:
  LocvarCalculator(int inxDim, int inyDim, int inwinDim, int inShellDim, std::vector<int> inWin, std::vector<int> inShell, std::vector<double> insDiv);
  //~LocvarCalculator();

  // The main shift functions
  void initialLocvarCalculation();
  void sideshift(int x, int y);
  void sideshift_hollow(int x, int y);
  void upshift(int x, int y);

  // Some sub-functions for shifting
  void sideshift_sub(int sx_curr, int sy_curr, int sy_next, int sy_prev, int x, int y);
  void upshift_sub(int sy_curr, int sx_curr, int sx_next, int sx_prev, int x, int y);

  // The export funciton
  Rcpp::NumericVector returnLocVar();

  // Some funcitons for debugging
  // void rprintsdiv();
  // void rprintlcov();
  // void rptintwin();
  // void rprintshell();
};


LocvarCalculator::LocvarCalculator(int inxDim, int inyDim, int inwinDim, int inShellDim, std::vector<int> inWin, std::vector<int> inShell, std::vector<double> insDiv){
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
    sideshift_hollow(0, y);
    upshift(0, y+1);
  }
  for (int x = 1; x < xDim; x++)
  {
    sideshift(x, yDim-1);
  }
}

/*
LocvarCalculator::~LocvarCalculator(){
  win.~vector();
  shell.~vector();
  sDiv.~vector();
  locvar.~vector();
}
*/

void LocvarCalculator::initialLocvarCalculation(){
  int xCoord, yCoord;
  xs = 0;
  x2s = 0;
  for (int p = 0; p < winDim; p++)
  {
    xCoord = myModulo(win[2*p],xDim);
    yCoord = myModulo(win[2*p+1], yDim);
    xs += sDiv[xCoord+(xDim*yCoord)];
    x2s += sDiv[xCoord+(xDim*yCoord)]*sDiv[xCoord+(xDim*yCoord)];
  }
  locvar.push_back((x2s - (xs*xs/winDim))/(winDim-1));
}

void LocvarCalculator::sideshift(int x, int y){
  sideshift_sub(shell[0], shell[1], shell[3], shell[2*shellDim-1], x, y);
  for(int p = 1; p < shellDim-1; p++){
    sideshift_sub(shell[2*p], shell[2*p+1], shell[2*p+3], shell[2*p-1], x, y);
  }
  sideshift_sub(shell[2*shellDim-2], shell[2*shellDim-1], shell[1], shell[2*shellDim-3], x, y);
  locvar.push_back((x2s - (xs*xs/winDim))/(winDim-1));
}

void LocvarCalculator::sideshift_hollow(int x, int y){
  sideshift_sub(shell[0], shell[1], shell[3], shell[2*shellDim-1], x, y);
  for(int p = 1; p < shellDim-1; p++){
    sideshift_sub(shell[2*p], shell[2*p+1], shell[2*p+3], shell[2*p-1], x, y);
  }
  sideshift_sub(shell[2*shellDim-2], shell[2*shellDim-1], shell[1], shell[2*shellDim-3], x, y);
}

void LocvarCalculator::sideshift_sub(int sx_curr, int sy_curr, int sy_next, int sy_prev, int x, int y){
  int current_i = myModulo(sx_curr+x, xDim);
  int current_j = myModulo(sy_curr+y, yDim);
  int prev_i = myModulo(sx_curr+x-1, xDim);

  bool addToXs, subtractFromXs;
  addToXs = (sx_curr >= 0) && (((sy_curr >= 0) && (sy_next < sy_curr)) || ((sy_curr < 0) && (sy_prev > sy_curr)));
  subtractFromXs = (sx_curr <= 0) && (((sy_curr <= 0) && (sy_curr < sy_next)) || ((sy_curr > 0) && (sy_prev < sy_curr)));

  if(addToXs){
    xs += sDiv[current_i+xDim*current_j];
    x2s += sDiv[current_i+xDim*current_j]*sDiv[current_i+xDim*current_j];
  }
  if(subtractFromXs){
    xs -= sDiv[prev_i+xDim*current_j];
    x2s -= sDiv[prev_i+xDim*current_j]*sDiv[prev_i+xDim*current_j];
  }
}

void LocvarCalculator::upshift(int x, int y){
  upshift_sub(shell[1], shell[0], shell[2], shell[2*shellDim-2], x, y);
  for(int p = 1; p < shellDim-1; p++){
    upshift_sub(shell[2*p+1], shell[2*p], shell[2*p+2], shell[2*p-2], x, y);
  }
  upshift_sub(shell[2*shellDim-1], shell[2*shellDim-2], shell[0], shell[2*shellDim-4], x, y);
  locvar.push_back((x2s - (xs*xs/winDim))/(winDim-1));
}

void LocvarCalculator::upshift_sub(int sy_curr, int sx_curr, int sx_next, int sx_prev, int x, int y){
  int current_i = myModulo(sx_curr+x, xDim);
  int current_j = myModulo(sy_curr+y, yDim);
  int prev_j = myModulo(sy_curr+y-1, yDim);

  bool addToXs, subtractFromXs;
  addToXs = (sy_curr >= 0) && (((sx_curr <= 0) && (sx_next > sx_curr)) || ((0 < sx_curr) && (sx_prev < sx_curr)));
  subtractFromXs = (sy_curr <= 0) && (((0 <= sx_curr) && (sx_next < sx_curr)) || ((sx_curr < 0) && (sx_curr < sx_prev)));

  if(addToXs){
    xs += sDiv[current_i+xDim*current_j];
    x2s += sDiv[current_i+xDim*current_j]*sDiv[current_i+xDim*current_j];
  }
  if(subtractFromXs){
    xs -= sDiv[current_i+xDim*prev_j];
    x2s -= sDiv[current_i+xDim*prev_j]*sDiv[current_i+xDim*prev_j];
  }
}

Rcpp::NumericVector LocvarCalculator::returnLocVar(){
  return Rcpp::wrap(locvar);
}

/*

void LocvarCalculator::rprintsdiv(){
  for (int y = 0; y < yDim; y++)
  {
    for (int x = 0; x < xDim; x++)
    {
      Rcpp::Rcout << " | " << sDiv[x+y*xDim];
    }
    Rcpp::Rcout << "\n";
  }
}

void LocvarCalculator::rprintlcov(){
  for (int y = 0; y < yDim; y++)
  {
    for (int x = 0; x < xDim; x++)
    {
      Rcpp::Rcout << " | " << locvar[x+y*xDim];
    }
    Rcpp::Rcout << "\n";
  }
}

void LocvarCalculator::rptintwin(){
  for (int p = 0; p < winDim; p++)
  {
    Rcpp::Rcout << win[2*p] << " | " << win[2*p+1] << "\n";
  }
  Rcpp::Rcout << "\n";
}

void LocvarCalculator::rprintshell(){
  Rcpp::Rcout << "\n We show the values of the shell: \n\n";
  for (int p = 0; p < shellDim; p++)
  {
    Rcpp::Rcout << shell[2*p] << " | " << shell[2*p+1] << "\n";
  }
  Rcpp::Rcout << "\n";
}
*/

//' Calculates the local variance with periodic boundary conditions
//'
//' @param inWin Window
//' @param inShell Shell of window
//' @param insDiv field of which the local variance should be calculated
//' @param inDims Integer vector of dimensions for input, should be <xDim, yDim, winDim, shellDim>
//'
//' @return Vector with the local variance
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calculateLocVar_periodic(Rcpp::NumericVector inWin, Rcpp::NumericVector inShell, Rcpp::NumericVector insDiv, Rcpp::IntegerVector inDims){
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
  LocvarCalculator LocVarCalcInstance = LocvarCalculator(xDim, yDim, winDim, shellDim, win, shell, sDiv);

  // Some tests for debugging
  //LocVarCalcInstance.rprintsdiv();
  //LocVarCalcInstance.rprintlcov();
  //LocVarCalcInstance.rptintwin();
  //LocVarCalcInstance.rprintshell();

  // We export the result
  // Rcpp::NumericVector res = LocVarCalcInstance.returnLocVar();
  // LocVarCalcInstance.~LocvarCalculator();
  // return res;
  return LocVarCalcInstance.returnLocVar();
}
