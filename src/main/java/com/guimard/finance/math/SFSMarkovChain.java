package com.guimard.finance.math;public abstract class SFSMarkovChain extends MarkovChain {    int N;    int[] a;    int[] b;    double[][] partition;    public SFSMarkovChain(int T, int j_0) {        super(T, j_0);    }    public SFSMarkovChain(int T, int j_0, int N) {        super(T, j_0);                // allocate MarkovChain        this.N = N;                    // dimension of square matrix P        partition = new double[N][];        a = new int[N];        b = new int[N];        for (int i = 0; i < N; i++) {            a[i] = a(i);            b[i] = b(i);            int ni = b[i] - a[i] + 2;         // number of indices a[i]<=j<=b[i]+1            partition[i] = new double[ni];            double sum = 0;            for (int j = 0; j < ni; j++) {                partition[i][j] = sum;                sum += q(i, j + a[i]);            }            if (sum < 0.9999999) {                String message;                message = "\n\nprobabilities don't sum to one, i=" + i +                        ", sum=" + sum +                        "\nthis can cause errors." +                        "\nprobabilities q(i,j):\n";                System.out.println(message);                for (int j = 0; j < N; j++)                    System.out.print(q(i, j) + ", ");            } // end if        } // end for i    } //end constructor    public SFSMarkovChain(int T, int j_0, int N, int[] a, int b[]) {        super(T, j_0);                // allocate MarkovChain        this.N = N;                    // dimension of square matrix P        partition = new double[N][];        this.a = a;        this.b = b;        for (int i = 0; i < N; i++) {            // points in the the partition of [0,1) conditional on X(t)=i            // note how baseline a[i] is added to index j.            int ni = b[i] - a[i] + 2;         // number of indices a[i]<=j<=b[i]+1            partition[i] = new double[ni];            double sum = 0;            for (int j = 0; j < ni; j++) {                partition[i][j] = sum;                sum += q(i, j + a[i]);            }            if (sum < 0.9999999) {                String message;                message = "\n\nprobabilities don't sum to one, i=" + i +                        ", sum=" + sum +                        "\nthis can cause errors." +                        "\nprobabilities q(i,j):\n";                System.out.println(message);                for (int j = 0; j < N; j++)                    System.out.print(q(i, j) + ", ");            } // end if        } // end for i    } //end constructor    public abstract double q(int i, int j);    public abstract int a(int i);    public abstract int b(int i);    public double q(int t, int i, int j) {        return 0;    }    public double I(int i, int j) {        return partition[i][j - a[i]];    }    private int j(int i, double u) {        if (i >= N) System.out.println("bad i=" + i + " u=" + u);        int j = a[i], k = b[i] + 1, m;        while (k - j > 1) {            m = (k + j) / 2;            if (I(i, m) > u) k = m;            else j = m;        }        return j;    } // end jj    public void timeStep(int t) {        int i = (int) path[t];         //current state        double u = RandomValue.U1();       //uniform draw from [0,1)        path[t + 1] = j(i, u);           //next state    }}