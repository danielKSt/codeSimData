#ifndef CODESIMDATA_LOCALSTATISTICS_H
#define CODESIMDATA_LOCALSTATISTICS_H

// The statistics we can calculate over the sliding window. The window keeps running sums of the
// values (xs), of the squared values (x2s) and of the absolute values (xabs), together with the
// pixel count, and every statistic below is calculated from these four.
enum LocalStatistic {
  LOCAL_VARIANCE = 0,
  LOCAL_MEAN_SQUARE = 1,
  LOCAL_MEAN = 2,
  LOCAL_MEAN_ABS = 3
};

inline double windowStatistic(int statistic, double xs, double x2s, double xabs, int divisor){
  switch (statistic)
  {
  case LOCAL_MEAN_SQUARE:
    return x2s/divisor;
  case LOCAL_MEAN:
    return xs/divisor;
  case LOCAL_MEAN_ABS:
    return xabs/divisor;
  default:
    return (x2s - (xs*xs/divisor))/(divisor-1);
  }
}

#endif
