#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]);

int congruence(int a, int b, int c, bool* error);
int GCD(int i, int j);
int MAX(int i, int j);
int MIN(int i, int j);
int SIGN(int i);
int LCRG_EVALUATE(int a, int b, int c, int x);
void LCRG_ANBN(int a, int b, int c, int n, int* an, int* bn);
int POWER_MOD(int a, int n, int m);

int main(int argc, char* argv[]) {
    // ���� �������� ������, �� P-������� ��������� ���� �����������
    // ���� ���� ����������� �� � ���� ��������
    int a, an, b, bn, c, id, j, k, k_hi, p, u, v;
    int numP;
    double Start, Finish;
    double duration = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);//������� ���������
    MPI_Comm_rank(MPI_COMM_WORLD, &id);//���� ���������
    MPI_Status status;

    if (id == 0) {

        /*cout << "  The number of processors is P = " << p << "\n";
        cout << "\n";
        cout << "  This program shows how a stream of random numbers\n";
        cout << "  can be computed 'in parallel' in an MPI program.\n";
        cout << "\n";
        cout << "  We assume we are using a linear congruential\n";
        cout << "  random number generator or LCRG, which takes\n";
        cout << "  an integer input and returns a new integer output:\n";
        cout << "\n";*/
        cout << "    U = ( A * V + B ) mod C\n";//�������, �� ���� �������� ����������� ���������������� �����
        cout << "\n";
    }
    //��������� �������-�������
    a = 16807;
    b = 45789;
    c = 2147483637;

    if (id == 0) {
        cout << "  LCRG parameters:\n";
        cout << "\n";
        cout << "  A  = " << a << "\n";
        cout << "  B  = " << b << "\n";
        cout << "  C  = " << c << "\n";
    }

    //������� ������������ �����
    k_hi = 40000000;
    Start = MPI_Wtime();
    // �-�������� �������� ��������� ��� �-� ������� �����������
    //��������� �������� ��� ������� �����������
    v = 12345;
    LCRG_ANBN(a, b, c, p, &an, &bn);
    if (id == 0)
    {
        
        cout << "\n";
        cout << "  LCRG parameters for P processors:\n";
        cout << "\n";
        cout << "  AN = " << an << "\n";
        cout << "  BN = " << bn << "\n";
        cout << "  C  = " << c << "\n";
        cout << "\n";
        
        //���������� ������ �������� � ��������������� ���������������� �����
        for (j = 1; j <= k_hi/p; j++)
        {
            an = a;
            bn = b;
            for (int k = 1; k < p;k++) {
                LCRG_ANBN(an, bn, c, k, &an, &bn);
                u = v;
                v = LCRG_EVALUATE(an, bn, c, u);           
                /*cout << v<<"\n";*/
                MPI_Send(&v, 1, MPI_INT, k, 1, MPI_COMM_WORLD);                  
            }         
        }                
    }
    else {
        k = id;

        //���������� ��� ����� �������� �-� ������������� �����
        for (k = id + p; k <= k_hi; k = k + p)
        {           
            MPI_Recv(&v, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status); 
            u = v;
            v = LCRG_EVALUATE(an, bn, c, u);
        }
    }
    
    Finish = MPI_Wtime();
    duration = Finish - Start;
    MPI_Finalize();

    if (id == 0)
    {
        cout << "\n";
        cout <<"Time of parallel execution: "<< duration;//��� ��������� ������������ ����
    }
    return 0;
}

//����'������ ������� ( A * X ) mod B = C
int congruence(int a, int b, int c, bool* error) {
    #define N_MAX 100
    int a_copy;
    int a_mag;
    int a_sign;
    int b_copy;
    int b_mag;
    int c_copy;
    int g;
    int k;
    int n;
    int q[N_MAX];
    bool swap;
    int x, y, z;

    *error = false;
    x = 0;
    y = 0;

    if (a == 0 && b == 0 && c == 0) {
        x = 0;
        return x;
    }
    else if (a == 0 && b == 0 && c != 0) {
        *error = true;
        x = 0;
        return x;
    }
    else if (a == 0 && b != 0 && c == 0) {
        x = 0;
        return x;
    }
    else if (a == 0 && b != 0 && c != 0) {
        x = 0;
        if ((c % b) != 0) {
            *error = true;
        }
        return x;
    }
    else if (a != 0 && b == 0 && c == 0) {
        x = 0;
        return 0;
    }
    else if (a != 0 && b == 0 && c != 0) {
        x = c / a;
        if ((c % a) != 0) {
            *error = true;
        }
        return x;
    }
    else if (a != 0 && b != 0 && c == 0) {
        x = 0;
        return x;
    }
    //1. ���������� ��� ����� a,b, �� ����� ������� �� �
    g = GCD(a, b);

    if ((c % g) != 0) {
        *error = true;
        return x;
    }
    a_copy = a / g;
    b_copy = b / g;
    c_copy = c / g;
    //2. �������� ����� A,B �� �������� �� ������� � ����
    a_mag = abs(a_copy);
    a_sign = SIGN(a_copy);
    b_mag = abs(b_copy);
    //�������� �� a==1 ��� b==1
    if (a_mag == 1) {
        x = a_sign * c_copy;
        return x;
    }
    else if (b_mag == 1) {
        x = 0;
        return x;
    }
    //���������� ����������� �������� ������
    if (b_mag <= a_mag) {
        swap = false;
        q[0] = a_mag;
        q[1] = b_mag;
    }
    else {
        swap = true;
        q[0] = b_mag;
        q[1] = a_mag;
    }
    n = 3;

    for (;;) {
        q[n - 1] = (q[n - 3] % q[n - 2]);

        if (q[n - 1] == 1) {
            break;
        }

        n += 1;
        if (N_MAX < n) {
            *error = true;
            cout << "\n";
            cout << "CONGRUENCE - Fatal error!\n";
            cout << "  Exceeded number of iterations.\n";
            exit(1);
        }
    }

    //4. ����'������ ������� X*A_MAG + Y*B_MAG = 1
    y = 0;
    for (k = n; k >= 2;k--) {
        x = y;
        y = (1 - x * q[k - 2]) / q[k - 1];
    }

    //5. ���������� �� ����� ��������� �������� �� �����
    if (swap) {
        z = x;
        x = y;
        y = z;
    }
    //6. ��������� ����, ��� �� ���� ������� X * A + Y * B = 1
    x = x * a_sign;
    //7. ������� �� �, ��� �������� X * A + Y * B = C
    x = x * c_copy;
    //8. ������ �������� �� ������� �
    x = x % b;
    //9. ���� ����� ��'�����, ������ ���� �������
    if (x < 0) {
        x += b;
    }

    return x;
#undef N_MAX;
}

//��������� ��� ����� i , j
int GCD(int i, int j) {
    int ip;
    int iq;
    int ir;
    //���������� �� ����� �� � ������ 0
    if (i == 0) {
        return MAX(1, abs(j));
    }
    else if (j == 0) {
        return MAX(1, abs(i));
    }
    //������� ����� � ����� �����
    ip = MAX(abs(i), abs(j));
    iq = MIN(abs(i), abs(j));
    //��������� �������� ������
    for (;;) {
        ir = ip % iq;
        if (ir == 0) {
            break;
        }
        ip = iq;
        iq = ir;
    }

    return iq;
}
//���� ����������� �������� � ����
int MAX(int i, int j) {
    int value;
    if (j < i) {
        value = i;
    }
    else {
        value = j;
    }
    return value;
}
//���� �������� �������� � ����
int MIN(int i, int j) {
    int value;
    if (i < j) {
        value = i;
    }
    else {
        value = j;
    }
    return value;
}
//���� ���� �������� �����
int SIGN(int i) {
    int value;
    if (i < 0) {
        value = -1;
    }
    else if(i >0){
        value = 1;
    }
    else {
        value = 0;
    }
    return value;
}
//���������� ( A^N ) mod M 
//�� ��������� �������� ��������� �� ������� �� ������
int POWER_MOD(int a, int n, int m) {
    long long int aa; //square
    int d;
    long long int m2;
    int x;
    long long int x2;

    if (a < 0) {
        return -1;
    }

    if (m <= 0) {
        return -1;
    }
    if (n < 0) {
        return -1;
    }
    //�������� ������� ����� ����� ������ 
    aa = (long long int) a;
    m2 = (long long int) m;
    x2 = (long long int) 1;

    while (n > 0) {
        d = n % 2;
        if (d == 1) {
            x2 = (x2 * aa) % m2;
        }
        
        aa = (aa * aa) % m2;
        n = (n - d) / 2;
    }
    //������������, �� � ����� ���� 0
    if (x2 < 0) {
        x2 += m2;
    }

    x = (int)x2;

    return x;
}
//�������� "N-������" ������� ���������� ���������������� �����
//F(X_i) = (a * F(X_i-1) + b) mod c
//F(N) = a^N * F(0) + ( a^N - 1) / ( a - 1 ) * b
//     = AN * F(0) + BN
void LCRG_ANBN(int a, int b, int c, int n, int *an, int *bn) {
    int am1;
    int anm1tb;
    bool ierror;

    if (n < 0) {
        cerr << "\n";
        cerr << "LCRG_ANBN - Fatal error!\n";
        cerr << "Illegal input value of N = " << n << "\n";
        exit(1);
    }

    if (c <= 0) {
        cerr << "\n";
        cerr << "LCRG_ANBN - Fatal error!\n";
        cerr << "Illegal input value of C = " << c << "\n";
        exit(1);
    }

    if (n == 0) {
        *an = 1;
        *bn = 0;
    }
    else if (n == 1) {
        *an = a;
        *bn = b;
    }
    else
    {
        //��������� A^N
        *an = POWER_MOD(a, n, c);
        //����'����� (a - 1)* BN = (a ^ N - 1) mod B ������� for BN
        am1 = a - 1;
        anm1tb = (*an - 1) ;
        //����������� ����� ��� ����'������ ��������� ������
        *bn = congruence(am1, b, anm1tb, &ierror);

        /*if (ierror) 
        {
            cerr << "\n";
            cerr << "LCRG_ANBN - Fatal error!\n";
            cerr << "  An error occurred in the CONGRUENCE routine.\n";
            exit(1);
        }*/
    }
    return;
}
//�������� ���� ������� ���������� y = ( A * x + B ) mod C �� ���� ���������� 
int LCRG_EVALUATE(int a, int b, int c, int x) {
    long long int aLong;
    long long int bLong;
    long long int cLong;
    long long int xLong;
    int y;
    long long int yLong;

    aLong = (long long int) a;
    bLong = (long long int) b;
    cLong = (long long int) c;
    xLong = (long long int) x;

    yLong = (aLong * xLong + bLong) % cLong;
    y = (int)yLong;

    if (y < 0) {
        y += c;
    }

    return y;
}
