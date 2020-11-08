# Optimize Code


关于代码优化的一些启示

<!--more-->

今天看到一篇博文[Optimize Your Code: Matrix Multiplication](https://docs.microsoft.com/zh-cn/archive/blogs/xiangfan/optimize-your-code-matrix-multiplication)，里面提出了几个对于矩阵乘法的优化方法，于是试了一下。

## 优化方法

**1. 原始写法**

```c++
template<typename T>
void matrixMul1(int size, T** m1, T** m2, T** result) {
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            result[i][j] = 0;
            for (int k=0; k<size; k++) {
                result[i][j] += m1[i][k] *m2[k][j];
            }
        }
    }
}
```

**2. 采用临时变量**

```c++
template<typename T>
void matrixMul2(int size, T** m1, T** m2, T** result) {
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            T c = 0;
            for (int k=0; k<size; k++) {
                c += m1[i][k] *m2[k][j];
            }
            result[i][j] = c;
        }
    }
}
```

很简单，但确实很有效。采用临时变量节省了在`result`中寻址并写入值的时间

**3. 连续寻址**

```c++
template<typename T>
void transpose(int size, T** m) {
    for (int i = 0; i < size; i++) {
        for (int j = i+1; j < size; j++) {
            swap(m[i][j],m[j][i]);
        }
    }
}

template<typename T>
void matrixMul3(int size, T** m1, T** m2, T** result) {
    transpose(size, m2);
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            T c = 0;
            for (int k=0; k<size; k++) {
                c += m1[i][k] *m2[j][k];
            }
            result[i][j] = c;
        }
    }
    transpose(size, m2);
}
```

转置矩阵，使`m1[i]`和`m2[i]`都可以按顺序访问，这会极大提高性能。

## 开启优化对于性能的影响

采用下面的测试代码，比较上面的三种方案计算速度，并和Eigen库做对比。

```c++
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

clock_t start_t,end_t;

template<typename T>
T** randomMatrix(int size) {
    T** a;
    a = new T* [size];
    for (int i=0; i<size; i++) {
        a[i] = new T [size];
    }
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            a[i][j] = (T)(1.0 - rand()/(RAND_MAX/2.0));
        }
    }
    return a;
}

template<typename T>
void matrixMul1(int size, T** m1, T** m2, T** result) {
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            result[i][j] = 0;
            for (int k=0; k<size; k++) {
                result[i][j] += m1[i][k] *m2[k][j];
            }
        }
    }
}

template<typename T>
void matrixMul2(int size, T** m1, T** m2, T** result) {
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            T c = 0;
            for (int k=0; k<size; k++) {
                c += m1[i][k] *m2[k][j];
            }
            result[i][j] = c;
        }
    }
}


template<typename T>
void transpose(int size, T** m) {
    for (int i = 0; i < size; i++) {
        for (int j = i+1; j < size; j++) {
            swap(m[i][j],m[j][i]);
        }
    }
}

template<typename T>
void matrixMul3(int size, T** m1, T** m2, T** result) {
    transpose(size, m2);
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            T c = 0;
            for (int k=0; k<size; k++) {
                c += m1[i][k] *m2[j][k];
            }
            result[i][j] = c;
        }
    }
    transpose(size, m2);
}


int main(int argc, char *argv[]) {
    int n = atoi(argv[1]);

    double** m1 = randomMatrix<double>(n);
    double** m2 = randomMatrix<double>(n);
    double** m3 = randomMatrix<double>(n);

    double endTime = 0;
    for (int i=0; i<10; i++) {
        start_t = clock();
        matrixMul1(n,m1,m2,m3);
        end_t = clock();
        endTime += (double)(end_t-start_t)/CLOCKS_PER_SEC;
    }
    cout << endTime*100 << "ms" << endl;

    endTime = 0;
    for (int i=0; i<10; i++) {
        start_t = clock();
        matrixMul2(n,m1,m2,m3);
        end_t = clock();
        endTime += (double)(end_t-start_t)/CLOCKS_PER_SEC;
    }
    cout << endTime*100 << "ms" << endl;

    endTime = 0;
    for (int i=0; i<10; i++) {
        start_t = clock();
        matrixMul3(n,m1,m2,m3);
        end_t = clock();
        endTime += (double)(end_t-start_t)/CLOCKS_PER_SEC;
    }
    cout << endTime*100 << "ms" << endl;

    MatrixXd em1 = MatrixXd::Random(n,n);
    MatrixXd em2 = MatrixXd::Random(n,n);
    MatrixXd em3 = MatrixXd::Random(n,n);

    endTime = 0;
    for (int i=0; i<10; i++) {
        start_t = clock();
        em3 = em1*em2;
        end_t = clock();
        endTime += (double)(end_t-start_t)/CLOCKS_PER_SEC;
    }
    cout << endTime*100 << "ms" << endl;

    return 0;
}
```

对两个nxn的随机方阵做相乘运算，随机数范围-1~1，双精度，取十次计算的平均耗时。得到如下结果

| size = 500x500 | 无编译优化   | -O1         | -O2          | -O3         |
| -------------- | ------------ | ----------- | ------------ | ----------- |
| matrixMul1     | 729.774 ms   | 366.574 ms  | 159.6     ms | 159.71   ms |
| matrixMul2     | 532.769 ms   | 199.984 ms  | 172.28   ms  | 166.058 ms  |
| matrixMul3     | 375.81   ms  | 151.089 ms  | 149.625 ms   | 149.739 ms  |
| Eigen          | 1128      ms | 31.415   ms | 30.5796 ms   | 30.8327 ms  |

在没有编译优化的情况下，Eigen的速度甚至比最朴素的写法还要慢得多，而只要开启O1优化，就可以吊打自己写的代码了，一定是用了很多的黑科技比如指令集优化或者并行化？看来还是不要轻易自己造轮子，用别人轮子的时候也要小心啊。


