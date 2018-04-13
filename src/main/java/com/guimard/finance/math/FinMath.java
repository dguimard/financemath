package com.guimard.finance.math;import cern.colt.matrix.DoubleMatrix1D;import cern.colt.matrix.DoubleMatrix2D;import cern.colt.matrix.impl.DenseDoubleMatrix1D;import cern.colt.matrix.impl.DenseDoubleMatrix2D;import cern.colt.matrix.linalg.Algebra;import cern.colt.matrix.linalg.Blas;import cern.colt.matrix.linalg.SeqBlas;import static org.apache.commons.math3.special.Erf.erf;public class FinMath {    public static double round(double x, int n) {        int f = 1;        for (int j = 0; j < n; j++) f *= 10;        return (double) Math.round(f * x) / f;    }    public double normalDensity(double[] z) {        double f = 2.506628274631,    // sqrt(2pi)                q = 1;        for (int k = 0; k < z.length; k++) q *= (f * Math.exp(-z[k] * z[k] / 2));        return q;    }    public static double N(double x) {        return 0.5 * (1.0 + erf(x / Math.sqrt(2.0)));    }    public static double d_plus(double Q, double k, double Sigma) {        return (Math.log(Q / k) / Sigma) + Sigma / 2;    }    public static double d_minus(double Q, double k, double Sigma) {        return (Math.log(Q / k) / Sigma) - Sigma / 2;    }    public static double blackScholesFunction(double Q, double k, double Sigma) {        return Q * N(d_plus(Q, k, Sigma)) - k * N(d_minus(Q, k, Sigma));    }    public static double dBSF(double Q, double k, double Sigma) {        double Sqrt2Pi = 2.506628274631;        return 2 * k * Math.exp(-d_minus(Q, k, Sigma) * d_minus(Q, k, Sigma) / 2) / Sqrt2Pi;    }    public static double bsDiscountedCallPrice            (double S, double K, double tau, double sigma, double B) {        double Sigma = sigma * Math.sqrt(tau);        return blackScholesFunction(S, K, Sigma) / B;    }    public static double bsDiscountedPutPrice            (double S, double K, double tau, double sigma, double B) {        double Sigma = sigma * Math.sqrt(tau);        return (K + blackScholesFunction(S, K, Sigma) - S) / B;    }    public static double NewtonSolveBSF(double Q, double k, double y)            throws Exception {        double Sigma = 0.5,                              //initial value for solution                err;                                    //current error        try {            if (y >= Q) throw new Exception();            err = blackScholesFunction(Q, k, Sigma) - y;            while (Math.abs(err) > 0.001)              //Newton's algorithm            {                Sigma -= err / dBSF(Q, k, Sigma);                err = blackScholesFunction(Q, k, Sigma) - y;            }        } catch (Exception e) {            System.err.println("Impossible option price, returning zero.");            return 0;        }        return Sigma;    } //end BSF_Solve    public static double BisectionSolveBSF(double Q, double k, double y)            throws Exception {        double a = 0.00001, c, b = 1000.0, error;        int n = 0;        if (blackScholesFunction(Q, k, a) > y) return a;        try {            if (y >= Q) throw new Exception();            while (n < 100) {                c = 0.5 * (a + b);                error = blackScholesFunction(Q, k, c) - y;                if (Math.abs(error) < 0.00000001)                    return c;                else if (error > 0) b = c;                else a = c;                n++;            }        } catch (Exception e) {            System.err.println("Impossible option price, returning zero.");            return 0;        }        return a;    }    public static double dN(double x) {        double Sqrt2Pi = 2.506628274631;        return Math.exp(-x * x / 2) / Sqrt2Pi;    }    public static double N_Solve(double y) {        double x = 0.5,                              //initial value for solution                err;                                //current error        try {            if (y >= 1) throw new Exception();            err = 1;            while (Math.abs(err) > 0.000001)        //Newton's algorithm            {                x -= err / dN(x);                err = err = N(x) - y;            }        } catch (Exception e) {            System.err.println("Equation has no solution, returning zero,\n" +                    "Computation may be meaningless!");            return 0;        }        return x;    } //end BS_f_Solve    public static double N_Solve1(double y)            throws Exception {        double a = -100000, c, b = 100000, error;        int n = 0;        try {            if ((y < N(a)) || (y > N(b))) throw new Exception();            while (n < 50) {                c = 0.5 * (a + b);                error = N(c) - y;                if (Math.abs(error) < 0.000001)                    return c;                else if (error > 0) b = c;                else a = c;                n++;            }        } catch (Exception e) {            System.err.println                    ("Solution outside the range of this equation solver");            System.exit(0);        }        return a;    } //end N_Solve1    public static double N_Inverse(double x) {        final double SQRT_TWO_PI = 2.5066282746310;        // Coefficients for the rational approximation.        final double                e_1 = -3.969683028665376e+01,                e_2 = 2.209460984245205e+02,                e_3 = -2.759285104469687e+02,                e_4 = 1.383577518672690e+02,                e_5 = -3.066479806614716e+01,                e_6 = 2.506628277459239e+00;        final double                f_1 = -5.447609879822406e+01,                f_2 = 1.615858368580409e+02,                f_3 = -1.556989798598866e+02,                f_4 = 6.680131188771972e+01,                f_5 = -1.328068155288572e+01;        final double                g_1 = -7.784894002430293e-03,                g_2 = -3.223964580411365e-01,                g_3 = -2.400758277161838e+00,                g_4 = -2.549732539343734e+00,                g_5 = 4.374664141464968e+00,                g_6 = 2.938163982698783e+00;        final double                h_1 = 7.784695709041462e-03,                h_2 = 3.224671290700398e-01,                h_3 = 2.445134137142996e+00,                h_4 = 3.754408661907416e+00;        // Limits of the approximation region: ]-oo,x_l[, [x_l,x_u], ]x_u,+oo[        final double                x_l = 0.02425,                x_u = 1.0 - x_l;        double z, r;        // Lower region: 0 < x < x_l        if (x < x_l) {            z = Math.sqrt(-2.0 * Math.log(x));            z = (((((g_1 * z + g_2) * z + g_3) * z + g_4) * z + g_5) * z + g_6) /                    ((((h_1 * z + h_2) * z + h_3) * z + h_4) * z + 1.0);        }        // Central region: x_l <= x <= x_u        else if (x <= x_u) {            z = x - 0.5;            r = z * z;            z = (((((e_1 * r + e_2) * r + e_3) * r + e_4) * r + e_5) * r + e_6) * z /                    (((((f_1 * r + f_2) * r + f_3) * r + f_4) * r + f_5) * r + 1.0);        }        // Upper region. ( x_u < x < 1 )        else {            z = Math.sqrt(-2.0 * Math.log(1.0 - x));            z = -(((((g_1 * z + g_2) * z + g_3) * z + g_4) * z + g_5) * z + g_6) /                    ((((h_1 * z + h_2) * z + h_3) * z + h_4) * z + 1.0);        }        return z;    }    public static void CholeskyRoot(double[][] C, double[][] L) {        int n = C[0].length;    // dimension        double S;        for (int i = 0; i < n; i++)            for (int j = 0; j <= i; j++)             // lower triangular            {                S = 0;                for (int k = 0; k < j; k++) S += L[i][k] * L[j][k];                //now S=sum_{k=0}^{j-1}L_ik*L_jk=r_i(L)r_j(L)-L_ij*L_jj,                //and so S+L_ij*L_jj=r_i(L)r_j(L)=C[i][j]                if (j < i)   // L[j][j] already computed, compute L[i][j]                    L[i][j] = (C[i][j] - S) / L[j][j];                else      // j=i and so S+L^2_jj=C[j][j]                {                    if (C[j][j] - S <= 0) {                        System.err.println("Cholesky: matrix not positive definite.\n" +                                "Terminating...");                        System.exit(0);                    }                    L[j][j] = Math.sqrt(C[j][j] - S);                }            } //end for i    }    public static void PlusEquals            (DoubleMatrix1D w, double k, DoubleMatrix2D A, boolean transposeA,             DoubleMatrix1D z) {        // cern.colt.matrix.linalg.SeqBlas object to call        // basic linear algreba routines        Blas blas = SeqBlas.seqBlas;        // Sum+=h*gammaZ, set the upper triangular flag for g to true        blas.dgemv(transposeA, k, A, z, 1, w);    }    public static DoubleMatrix1D    linearSystemSolution(DoubleMatrix2D A, DoubleMatrix1D y) {        Algebra linAlg = Algebra.DEFAULT;        int m = A.rows(), n = A.columns();        DoubleMatrix2D B = A.like(m, 1);        for (int i = 0; i < m; i++) B.setQuick(i, 0, y.getQuick(i));        DoubleMatrix2D X = linAlg.solve(A, B);        return X.viewColumn(0);    }    public static void main(String[] args) {        DenseDoubleMatrix2D                A = new DenseDoubleMatrix2D(new double[][]{{1, 2, 1}, {2, 1, 1}, {1, 0, 1}});        System.err.println("A: " + A.toString());        DenseDoubleMatrix1D                y = new DenseDoubleMatrix1D(new double[]{2, 2, 2});        System.err.println("y: " + y.toString());        DoubleMatrix1D x = linearSystemSolution(A, y);        System.out.println("Solution\n" + x.toString());    }}