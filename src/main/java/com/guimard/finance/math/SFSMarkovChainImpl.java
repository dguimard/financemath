package com.guimard.finance.math;public class SFSMarkovChainImpl extends SFSMarkovChain {    double[][] Q;    public SFSMarkovChainImpl            (int T, int j_0, int N, int[] a, int[] b, double[][] Q) {        super(T, j_0);                // allocate MarkovChain        super.a = a;        super.b = b;        this.Q = Q;        super.N = N;                   // number of states        partition = new double[N][];        for (int i = 0; i < N; i++) {            // points in the the partition of [0,1) conditional on X(t)=i            // note how baseline a[i] is added to index j.            int ni = b[i] - a[i] + 2;         // number of indices a[i]<=j<=b[i]+1            partition[i] = new double[ni];            double sum = 0;            for (int j = 0; j < ni; j++) {                partition[i][j] = sum;                if (j < ni - 1) sum += Q[i][j];            }            if (sum < 0.9999999) {                String message;                message = "\n\nprobabilities don't sum to one, i=" + i +                        ", sum=" + sum +                        "\nthis can cause errors." +                        "\nprobabilities q(i,j):\n";                System.out.println(message);                for (int j = a[i]; j <= b[i]; j++)                    System.out.print(Q[i][j - a[i]] + ", ");            } // end if        } // end for i    } //end constructor    public int a(int i) {        return a[i];    }    public int b(int i) {        return b[i];    }    public double q(int i, int j) {        if ((0 <= i) && (i < N) && (a[i] <= j) && (j <= b[i])) return Q[i][j - a[i]];        else return 0;    }}