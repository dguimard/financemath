package com.guimard.finance.math;

public abstract class RandomVector {

    int dim;

    public RandomVector(int d) {
        dim = d;
    }


    public int dim() {
        return dim;
    }

    public abstract double[] getValue(int t);


    public double[]
    conditionalExpectation(int t, int N) {
        double[] X;      //sample of X
        double[] sum = new double[dim];  //running sum of samples

        for (int n = 0; n < N; n++) {
            X = getValue(t);
            Vector.add(sum, X);         //add X to sum
        }

        return Vector.scalarMult(1.0 / N, sum);
    }


    public double[] expectation(int N) {
        return conditionalExpectation(0, N);
    }


    public double[][]
    conditionalMeanAndStandardDeviation(int t, int N, int sampleGroupSize) {
        int nSampleGroups = N / sampleGroupSize;

        double[] X,                           // sample vector
                group_sum = new double[dim],   // sum of samples in current group
                group_mean,
                sum = new double[dim],         // sum of sample group means
                sumSquares = new double[dim],  // sum of squares of group means
                mean,                        // mean of group means
                meanSquares,                 // mean of group mean squares
                variance,                    // variance of sample group means
                stdv;                        // standard deviation of same


        for (int n = 0; n < nSampleGroups; n++) {
            Vector.setZero(group_sum);

            // sum of samples in the current sample group
            for (int k = 0; k < sampleGroupSize; k++) {
                X = getValue(t);   // next sample
                Vector.add(group_sum, getValue(t));
            }
            group_mean = Vector.scalarMult(1.0 / sampleGroupSize, group_sum);

            //update sum of group means, group mean squares
            Vector.add(sum, group_mean);
            Vector.add(sumSquares, Vector.mult(group_mean, group_mean));
        }

        double[][] results = new double[2][];

        mean = Vector.scalarMult(1.0 / nSampleGroups, sum);
        // copy this into results[0] at once since the call
        // Vector.mult(mean,mean) below destroys the vector mean
        results[0] = new double[dim];
        for (int j = 0; j < dim; j++) results[0][j] = mean[j];

        meanSquares = Vector.scalarMult(1.0 / nSampleGroups, sumSquares);
        variance = Vector.subtract(meanSquares, Vector.mult(mean, mean));
        stdv = Vector.sqrt(variance);
        results[1] = stdv;

        return results;

    }


    public double[][] doMeanAndStandardDeviation(int N, int sampleGroupSize) {
        return conditionalMeanAndStandardDeviation(0, N, sampleGroupSize);
    }


    public double[][]
    meanAndStandardDeviation(int N, int m) {
        return conditionalMeanAndStandardDeviation(0, N, m);
    }


    public double[][]
    conditionalMeanAndStandardDeviation
            (int t, int N, int m, int sampleGroupSize) {
        int nSampleGroups = N / sampleGroupSize;  //note integer division

        double[] X,                           // sample vector
                group_sum = new double[dim],   // sum of samples in group
                group_mean,
                sum = new double[dim],         // sum of group means
                sumSquares = new double[dim],  // sum of group mean squares
                mean,                        // mean of group means
                meanSquares,                 // mean of group mean squares
                variance,                    // variance of group means
                stdv;


        long before = System.currentTimeMillis();


        for (int n = 0; n < nSampleGroups; n++) {


            Vector.setZero(group_sum);
            for (int k = 0; k < sampleGroupSize; k++) {
                X = getValue(t);
                Vector.add(group_sum, X);
            }
            group_mean = Vector.scalarMult(1.0 / sampleGroupSize, group_sum);

            //update sum of group means, group mean squares
            Vector.add(sum, group_mean);
            Vector.add(sumSquares, Vector.mult(group_mean, group_mean));
        }


        double[][] results = new double[2][];

        mean = Vector.scalarMult(1.0 / nSampleGroups, sum);
        // copy this into results[0] at once since the call
        // Vector.mult(mean,mean) below destroys the vector mean
        results[0] = new double[dim];
        for (int j = 0; j < dim; j++) results[0][j] = mean[j];

        meanSquares = Vector.scalarMult(1.0 / nSampleGroups, sumSquares);
        variance = Vector.subtract(meanSquares, Vector.mult(mean, mean));
        stdv = Vector.sqrt(variance);
        results[1] = stdv;

        return results;

    } //end meanAndStandardDeviation 


    public double[][] meanAndStandardDeviation
            (int N, int m, int sampleGroupSize) {
        return conditionalMeanAndStandardDeviation
                (0, N, m, sampleGroupSize);
    }


    public double conditionalCovariance(int i, int j, int t, int N) {
        double sum_Xi = 0,
                sum_Xj = 0,
                sum_XiXj = 0,
                mean_Xi,
                mean_Xj,
                mean_XiXj;

        for (int n = 0; n < N; n++) {
            double[] x = getValue(t);

            sum_Xi += x[i];
            sum_Xj += x[j];
            sum_XiXj += x[i] * x[j];
        }

        mean_Xi = sum_Xi / N;
        mean_Xj = sum_Xj / N;
        mean_XiXj = sum_XiXj / N;

        return mean_XiXj - mean_Xi * mean_Xj;

    }

    public double covariance(int i, int j, int N) {
        return conditionalCovariance(i, j, 0, N);
    }


    public double[][] conditionalCovarianceMatrix(int t, int N) {
        double[][] sum_XiXj = new double[dim][dim];
        double[] sum_Xi = new double[dim];


        long before = System.currentTimeMillis();
        int loops = 0,
                m = N / 200;  // number of samples completed for next report


        // sample means, covariances
        // note symmetry, operate on lower triangular half only
        for (int n = 0; n < N; n++) {
            double[] x = getValue(t);

            for (int i = 0; i < dim; i++) {
                sum_Xi[i] += x[i];
                for (int j = 0; j <= i; j++) sum_XiXj[i][j] += x[i] * x[j];
            }


        }

        for (int i = 0; i < dim; i++) sum_Xi[i] /= N;           // E_t(X_i)

        for (int i = 0; i < dim; i++)
            for (int j = 0; j <= i; j++) {

                sum_XiXj[i][j] /= N;                         // E_t(X_iX_j)
                sum_XiXj[i][j] -= sum_Xi[i] * sum_Xi[j];       // Cov_t(X_i,X_j)
                // fill out the upper triangular half
                sum_XiXj[j][i] = sum_XiXj[i][j];
            }

        return sum_XiXj;

    } // end


    public double[][] covarianceMatrix(int N) {
        return conditionalCovarianceMatrix(0, N);
    }

    public double conditionalCorrelation(int i, int j, int t, int N) {
        double sum_Xi = 0,
                sum_Xj = 0,
                sum_XiXi = 0,
                sum_XjXj = 0,
                sum_XiXj = 0;

        for (int n = 0; n < N; n++) {
            double[] x = getValue(t);

            sum_Xi += x[i];
            sum_Xj += x[j];
            sum_XiXi += x[i] * x[i];
            sum_XjXj += x[j] * x[j];
            sum_XiXj += x[i] * x[j];

        }

        //when divided by N^4: 
        double Var_XiVar_Xj =
                ((N * sum_XiXi - sum_Xi * sum_Xi) * (N * sum_XjXj - sum_Xj * sum_Xj));

        return (N * sum_XiXj - sum_Xi * sum_Xj) / Math.sqrt(Var_XiVar_Xj);

    }

    public double correlation(int i, int j, int N) {
        return conditionalCorrelation(i, j, 0, N);
    }


}