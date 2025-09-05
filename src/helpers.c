#include "../include/helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function to read test signal from MATLAB */
void readTestSig(float complex *x, int length) {
  /* Open file for reading */
  FILE *fp;
  fp = fopen("../tests/testSignal.txt", "r");
  if (fp == NULL) {
    fprintf(stderr, "The file containing test signal does not exist!\n");
    exit(EXIT_FAILURE);
  }

  /* Read input file line by line and store in vector */
  float val[2];
  int ret;
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < 2; j++)
      ret = fscanf(fp, "%f\n", &val[j]);
    if (ret == 0) {
      fprintf(stderr, "File reading failed\n");
    }
    x[i] = val[0] + I * val[1];
  }

  /* Close file after use */
  fclose(fp);
}

/* Function to write signal to file */
void sig2file(void *x, int type, int length, char *path) {

  /* Open file for reading */
  FILE *fp;
  fp = fopen(path, "w");
  if (fp == NULL) {
    fprintf(stderr, "Unable to open file for writing!\n");
    exit(EXIT_FAILURE);
  }

  /* Write to file for all supported signal types */
  if (type == 0) { // Int
    int *y = (int *)x;
    for (int i = 0; i < length; i++) {
      fprintf(fp, "%d\n", y[i]);
    }
  } else if (type == 1) { // Unsigned char
    unsigned char *y = (unsigned char *)x;
    for (int i = 0; i < length; i++) {
      fprintf(fp, "%d\n", y[i]);
    }
  } else if (type == 2) { // float
    float *y = (float *)x;
    for (int i = 0; i < length; i++) {
      fprintf(fp, "%f\n", y[i]);
    }
  } else if (type == 3) { // Complex float
    float complex *y = (float complex *)x;
    float val[2];
    for (int i = 0; i < length; i++) {
      val[0] = crealf(y[i]);
      val[1] = cimagf(y[i]);
      for (int j = 0; j < 2; j++) {
        fprintf(fp, "%f\n", val[j]);
      }
    }
  }

  fclose(fp);
}
/* fliplr function for complex vectors */
void compFliplr(float complex *d_orig, float complex *flip, int length) {
  int l = length - 1;
  for (int i = 0; i < length; i++) {
    flip[l] = d_orig[i];
    l--;
  }
}

/* fliplr function for real vectors */
void reFliplr(float *d_orig, float *d_flip, int length) {
  int l = length - 1;
  for (int i = 0; i < length; i++) {
    d_flip[l] = d_orig[i];
    l--;
  }
}

/* fliplr function for int vectors */
void intFliplr(int *d_orig, int *d_flip, int length) {
  int l = length - 1;
  for (int i = 0; i < length; i++) {
    d_flip[l] = d_orig[i];
    l--;
  }
}

/* Find maximum value in input array and its index */
void max(float *arr, int length, float *mx, int *indx) {
  *mx = 0.0; /* Initial max */
  *indx = 0; /* Initial Index */
  for (int i = 0; i < length; i++) {
    if (arr[i] > *mx) {
      *mx = arr[i];
      *indx = i;
    }
  }
}

/* Function to find first or last K indices of non zero elements of vector
 */
void find(int *arr, int *indices, int length, int K, int dir) {
  /* If dir = 1, find 'first' K non zero elements. If dir = 0,
     find 'last' K non zero elements */
  if (dir == 1) {
    int k = 0;
    for (int i = 0; i < length; i++) {
      if (arr[i] != 0) {
        indices[k] = i;
        k++;
      }
      if (k == K) {
        break;
      }
    }
  } else {
    int k = 0;
    for (int i = length - 1; i > 0; i--) {
      if (arr[i] != 0) {
        indices[k] = i;
        k++;
      }
      if (k == K) {
        break;
      }
    }
    int indices_temp[K];
    memcpy(indices_temp, indices, sizeof(int) * K);
    intFliplr(indices_temp, indices, K); /* Sort into ascending order */
  }
}

/* Function to find max value in array and its index that also exceeds a
 * threshold */
void maxThresh(float *arr, float *thresh, int length, float *mx, int *indx) {
  /* The signal array and threshold array must be the same length */
  /* Array to hold values + indices that exceed threshold */
  float *vals = (float *)calloc(length, sizeof(float));
  int *indxs = (int *)calloc(length, sizeof(int));
  int j = 0;
  for (int i = 0; i < length; i++) {
    if (arr[i] > thresh[i]) {
      vals[j] = arr[i];
      indxs[j] = i;
      j++;
    }
  }

  if (j == 0) {
    /* If signal does not exceed thresholds */
    *mx = -1.0;
    *indx = -1;
  } else {
    /* Search through to find the max value and its index */
    int index;
    max(vals, length, mx, &index);
    *indx = indxs[index];
  }

  /* Clean up memory */
  free(vals);
  free(indxs);
}

/* Function to calculate magnitude squared of complex signal */
void mag2(float complex *x, float *y, int length) {
  for (int i = 0; i < length; i++)
    y[i] = crealf(x[i] * conjf(x[i]));
}

/* Function to calculate magnitude squared of complex signal */
float abs2(float complex x) { return crealf(x * conjf(x)); }

/* Function to estimate variance of complex signal */
float var(float complex *x, unsigned int len) {

  // Calculate mean
  float complex mu = compAvg(x, len);

  // Calculate un-biased variance
  float var = 0.0;
  for (int i = 0; i < len; i++) {
    var += abs2(x[i] - mu);
  }
  return var / (len - 1);
}

float complex compAvg(float complex *x, int length) {
  float avgRe = 0.0, avgIm = 0.0;
  float complex y;
  for (int i = 0; i < length; i++) {
    avgRe += crealf(x[i]);
    avgIm += cimagf(x[i]);
  }
  avgRe = avgRe / length;
  avgIm = avgIm / length;

  y = avgRe + avgIm * I;
  return y;
}

float reAvg(float *x, int length) {
  float avgRe = 0.0;
  for (int i = 0; i < length; i++) {
    avgRe += x[i];
  }
  float y = avgRe / (float)length;

  return y;
}

void getReal(float *y, float complex *x, int length) {
  for (int i = 0; i < length; i++)
    y[i] = crealf(x[i]);
}

void getImag(float *y, float complex *x, int length) {
  for (int i = 0; i < length; i++)
    y[i] = cimagf(x[i]);
}

/* Function to get abs value of real number */
float absre(float x) {

  float y = (x < 0) ? x * -1 : x;

  return y;
}

/* Sign function */
float sign(float x) {

  float y = (x < 0) ? -1.0 : 1.0;

  return y;
}

/* Find min of two positive numbers */
float min(float x, float y) {

  float m = (x < y) ? x : y;
  return m;
}

/* Function to convert decimal value to n-bit binary string */
void dec2bin(unsigned int x, unsigned int n, unsigned char *bin) {
  for (int i = n - 1; i >= 0; i--) {
    bin[i] = x & 1;
    x >>= 1;
  }
}
