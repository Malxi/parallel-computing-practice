#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pthread.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <math.h>

#define DEBUG 0

int setup(int argc, char** argv, char* &fstreama, char* &fstreamb, int* dim){
    /* 错误编号 0表示打开文件成功, 1表示打开失败 */
    int error = 0;
    /* 文件类型 */
    FILE *fha, *fhb;
    /* 文件大小 */
    int fsizea, fsizeb;

    /* 参数检查 */
    if (argc < 4) {
        printf("Invalid arguments!\n"); 
        printf("Usage: ./serial filea fileb filec\n");
        printf("filea, fileb and filec are file names for matrix A, B and C\n");
        return 1;
    }

    /* 打开矩阵a文件 */
    if (!(fha = fopen(argv[1], "r"))) {
        printf("Can't open matrix file %s, Errno=%d\n", argv[1], errno);
        return 1;
    }
    /* 打开矩阵b文件 */
    if (!(fhb = fopen(argv[2], "r"))) {
        printf("Can't open matrix file %s, Errno=%d\n", argv[2], errno);
        return 1;
    }

    /* 获得文件描述信息 */
    struct stat fstata, fstatb;
    stat(argv[1], &fstata);
    stat(argv[2], &fstatb);
    fsizea = fstata.st_size;
    fsizeb = fstatb.st_size;
    /* 分配文件流缓存 */
    fstreama = (char *)malloc(fsizea);
    fstreamb = (char *)malloc(fsizeb);
    /* 读取矩阵A, B */
    fread(fstreama, sizeof(char), fsizea, fha);
    fread(fstreamb, sizeof(char), fsizeb, fhb);

    /* 读取矩阵的维度 */
    int n1, n2, n3, n4;
    n1 = ((int*)fstreama)[0]; 
    n2 = ((int*)fstreama)[1]; 
    n3 = ((int*)fstreamb)[0]; 
    n4 = ((int*)fstreamb)[1]; 
    
    /* 检查矩阵数值 */
    if (n1 <=0 || n2 <=0 || n3 <= 0 || n4 <=0 || n2 != n3) {
        printf("Matrix size error, %dx%d with %dx%d\n", n1, n2, n3, n4);
        return 1;
    } 

    /* 检查矩阵A的大小是否合法 */
    if (fsizea < (sizeof(int)*2 + sizeof(double)*n1*n2)) {
        printf("Actual size of A mismatches with stated size\n");
        return 1;
    }
    /* 检查矩阵B的大小是否合法 */
    if (fsizeb < (sizeof(int)*2 + sizeof(double)*n3*n4)) {
        printf("Actual size of B mismatches with stated size\n");
        return 1;
    }

    /* 返回维度信息 */
    dim[0] = n1;
    dim[1] = n2;
    dim[2] = n4;

    fclose(fha); 
    fclose(fhb); 
    return 0;
}

void scatter_matrix(double *fstream, int m, int n, double *matrix, int rootp){
    /* 用于区分进行的变量 */
    int myid, numprocs;
    int i, j, k;
    int pos;
    MPI_Status status;
    int maxrows = (m + rootp - 1)/rootp;
    int maxcols = (n + rootp - 1)/rootp;
    int i_min, i_max, j_min, j_max;
    int rank;
    double *matrix_v = (double *)((int *)matrix+2);

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    /* 0号进程将 fstream 中的数据分发给其余进程 */
    if(myid == 0){
        /* 依次处理每一个进程, 从1号进程开始 */
        for(rank = 1; rank < numprocs; rank++){
            /* 计算进程所分得的矩阵块在总矩阵中的坐标范围 */
            i_min = (rank / rootp) * maxrows;
            j_min = (rank % rootp) * maxcols;
            /* 需要处理最后一块的矩阵大小不平衡 */
            i_max = (rank / rootp + 1) * maxrows > m ? m : (rank / rootp  + 1) * maxrows;
            j_max = (rank % rootp + 1) * maxcols > n ? n : (rank % rootp + 1) * maxcols;
            /* 0号进程将需要发送的数据块拷贝到缓冲区 */
            k = 0;
            for(i = i_min; i < i_max; i++){
                for(j = j_min; j < j_max; j++){
                    pos = i * n + j;
                    matrix_v[k++] = fstream[pos];
                }
            }
            /* 将相应的数据块发送给相应的进程 */
            if(DEBUG){
                printf("%d send matrix with dim %d x %d!\n", myid, i_max - i_min, j_max - j_min);
            }
            MPI_Send(matrix_v, maxrows * maxcols, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD);
        }
        /* 0号进程读取自己的数据块 */
        ((int*)matrix)[0] = maxrows;
        ((int*)matrix)[1] = maxcols;
        k = 0;
        for(i = 0; i < maxrows; i++){
            for(j = 0; j < maxcols; j++){
                pos = i * n + j;
                matrix_v[k++] = fstream[pos];
            }
        }
    }else{
        /* 其余进程接收数据 */
        /* 计算进程所分得的矩阵块在总矩阵中的坐标范围 */
        i_min = (myid / rootp) * maxrows;
        j_min = (myid % rootp) * maxcols;
        /* 需要处理最后一块的矩阵大小不平衡 */
        i_max = (myid / rootp + 1) * maxrows > m ? m : (myid / rootp + 1) * maxrows;
        j_max = (myid % rootp + 1) * maxcols > n ? n : (myid % rootp + 1) * maxcols;
        /* 将矩阵的尺寸写入 */
        ((int*)matrix)[0] = i_max - i_min;
        ((int*)matrix)[1] = j_max - j_min;
        MPI_Recv(matrix_v, maxrows * maxcols, MPI_DOUBLE, 0, myid, MPI_COMM_WORLD, &status);
        if(DEBUG){
            printf("%d receive matrix with dim %d x %d!\n", myid, ((int*)matrix)[0], ((int*)matrix)[1]);
            for(i = 0; i < i_max - i_min; i++){
                for(j = 0; j < j_max - j_min; j++){
                    printf("%.4f  ", matrix_v[i * ((int*)matrix)[1] + j]);
                }
                printf("\n");
            }
        }
    }
}

void cannon(double *A, double *bufA, int bufA_size, double *B, double *bufB, int bufB_size, \
double *C, int bufC_size, int n1, int n3, int rootp, int myid){
    int cur_row, cur_col;
    int desproc;
    int srcproc;
    int i, j, k, p;
    int d1, d2, d3;
    double *A_v = (double *)((int *)A+2);
    double *B_v = (double *)((int *)B+2);
    double *C_v = (double *)((int *)C+2);
    MPI_Status status;
    MPI_Request request;
    
    /* 计算每个进程的位置 */
    cur_row = myid / rootp;
    cur_col = myid % rootp;

    /* 首先将矩阵块对齐 */
    /* 将块A中坐标为(i,j)的分块A(i,j)向左循环移动i步 */
    desproc = cur_row * rootp + (cur_col - cur_row + rootp) % rootp;
    if(desproc != myid){
        MPI_Sendrecv(A, bufA_size, MPI_BYTE, desproc, 1, bufA, bufA_size, MPI_BYTE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
        memcpy(A, bufA, bufA_size);
    }
    /* 将块B中坐标为(i,j)的分块B(i,j)向上循环移动j步 */
    desproc = ((cur_row - cur_col + rootp) % rootp) * rootp + cur_col;
    if(desproc != myid){
        MPI_Sendrecv(B, bufB_size, MPI_BYTE, desproc, 2, bufB, bufB_size, MPI_BYTE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
        memcpy(B, bufB, bufB_size);
    }
    
    /* 写入块C的维度 */
    ((int*)C)[0] = ((int*)A)[0];
    ((int*)C)[1] = ((int*)B)[1];

    /* 对矩阵C进行初始化 */
    d1 = ((int*)C)[0];
    d3 = ((int*)C)[1];
    for(i = 0; i < d1; i++){
        for(j = 0; j < d3; j++){
            C_v[i * d3 + j] = 0.0;
        }
    }
    
    if(DEBUG && myid == 0){
        printf("matrix alignment success!\n");
        printf("dim 1: %d, dim 2: %d, dim 3: %d\n", ((int*)A)[0], ((int*)B)[0], ((int*)B)[1]);
    }

    /* 进行乘法rootp操作后, 每次单位移动矩阵块 */
    for(p = 0; p < rootp; p++){
        /* 矩阵块相乘 C = A * B 需要考虑矩阵的尺寸 */
        if(d1 != ((int*)A)[0]){
            printf("block dim 1 error!\n");
        }
        d2 = ((int*)A)[1];
        /* 实际情况是每个进程的尺寸恰好满足乘法要求 */
        if(d2 != ((int*)B)[0]){
            printf("block dim 2 error!\n");
        }
        if(d3 != ((int*)B)[1]){
            printf("block dim 3 error!\n");
        }
        /* 矩阵相乘，这部分还可以进行多线程处理 */
        for(i = 0; i < d1; i++){
            for(j = 0; j < d3; j++){
                for(k = 0; k < d2; k++){
                    C_v[i * d3 + j] += A_v[i*d2 + k] * B_v[k*d3 + j];
                }
            }
        }
        if(DEBUG && myid == 0) printf("loop %d end\n", p);
        /* 将分块A左移1块, 避免死锁 */
        desproc = cur_row * rootp + (cur_col - 1 + rootp) % rootp;
        srcproc = cur_row * rootp + (cur_col + 1) % rootp;
        MPI_Issend(A, bufA_size, MPI_BYTE, desproc, 1, MPI_COMM_WORLD, &request);
        MPI_Recv(bufA, bufA_size, MPI_BYTE, srcproc, 1, MPI_COMM_WORLD, &status);
        MPI_Wait(&request, &status);
        memcpy(A, bufA, bufA_size);
        /* 将分块B上移1块, 避免死锁 */
        desproc = ((cur_row - 1 + rootp) % rootp) * rootp + cur_col;
        srcproc = ((cur_row + 1) % rootp) * rootp + cur_col;
        MPI_Issend(B, bufB_size, MPI_BYTE, desproc, 1, MPI_COMM_WORLD, &request);
        MPI_Recv(bufB, bufB_size, MPI_BYTE, srcproc, 1, MPI_COMM_WORLD, &status);
        MPI_Wait(&request, &status);
        memcpy(B, bufB, bufB_size);
    }
}

void gather_matrix(double *fstreamc, int n1, int n3, double *C, int rootp){
    /* 用于区分进行的变量 */
    int myid, numprocs;
    int i, j, k;
    int pos;
    MPI_Status status;
    int maxrows = (n1 + rootp - 1)/rootp;
    int maxcols = (n3 + rootp - 1)/rootp;
    int i_min, i_max, j_min, j_max;
    int rank;
    double *C_v = (double *)((int *)C+2);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    /* 0号进程接收来自其余进程的数据 */
    if(myid == 0){
        /* 0号进程首先将自己的部分存入 fstreamc */
        for(i = 0; i < maxrows; i++){
            for(j = 0; j < maxcols; j++){
                fstreamc[i * n3 + j] = C_v[i * maxcols + j];
            }
        }
        /* 依次处理每一个进程, 从1号进程开始 */
        for(rank = 1; rank < numprocs; rank++){
            /* 接收来自rank进程的数据 */
            MPI_Recv(C_v, maxrows * maxcols, MPI_DOUBLE, rank, rank, MPI_COMM_WORLD, &status);
            /* 计算进程rank C矩阵块在总矩阵中的坐标范围 */
            i_min = (rank / rootp) * maxrows;
            j_min = (rank % rootp) * maxcols;
            /* 需要处理最后一块的矩阵大小不平衡 */
            i_max = (rank / rootp + 1) * maxrows > n1 ? n1 : (rank / rootp + 1) * maxrows;
            j_max = (rank % rootp + 1) * maxcols > n3 ? n3 : (rank % rootp + 1) * maxcols;
            /* 0号进程将数据拷贝到 fstreamc */
            for(i = i_min; i < i_max; i++){
                for(j = j_min; j < j_max; j++){
                    pos = (i - i_min) * (j_max - j_min) + j - j_min;
                    fstreamc[i * n3 + j] = C_v[pos];
                }
            }
            if(DEBUG){
                printf("%d receive C matrix with dim %d x %d!\n", myid, i_max - i_min, j_max - j_min);
                for(i = 0; i < i_max - i_min; i++){
                    for(j = 0; j < j_max - j_min; j++){
                        printf("%.4f  ", C_v[i * (j_max - j_min) + j]);
                    }
                    printf("\n");
                }
            }
        }
    }else{
        /* 其余进程发送数据 */
        MPI_Send(C_v, maxrows * maxcols, MPI_DOUBLE, 0 , myid, MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv){
    /* 用于区分进行的变量 */
    int myid, numprocs;
    int rootp;
    /* 矩阵维数 A: n1 x n2, B: n2 x n3, C: n1 x n3 */
    int n1, n2, n3;
    /* 为矩阵 A, B, C 分配缓存, 由于A, B需移动故有两份缓存 */
    double *A, *B, *C, *bufA, *bufB;
    char *buf;
    int bufA_size, bufB_size, bufC_size;
    /* 在进程0上缓存矩阵A, B, C */
    char *fstreama, *fstreamb, *fstreamc;
    double elapsed_time;
    double read_time;
    double save_time;
    /* 存储矩阵维数 */
    int dim[3];
    /* 分块大小信息 */
    int maxrows_a, maxcols_a, maxrows_b, maxcols_b;
    /* 文件存取变量 */
    FILE* fhc;
    int fsizec;

    /* MPI初始化 */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    /* 算法限制需检查启动进程数是否为平方数 */
    rootp = sqrt(numprocs);
    if (numprocs != rootp*rootp) {
        printf("Processor number must be a square!\n");
        MPI_Finalize(); // Something error during preprocessing
        exit(-1);
    }
    read_time = MPI_Wtime();
    /* 0号进程预处理命令行, 从文件读取矩阵A, B, 并且将其维数存到 dim */
    if (myid == 0) {
        if (setup(argc, argv, fstreama, fstreamb, dim)) {
            MPI_Finalize(); // Something error during preprocessing
            exit(-1);
        }
    }
    /* 广播矩阵A, B的维度信息 */
    MPI_Bcast(dim, 3, MPI_INT, 0, MPI_COMM_WORLD);
    /* 每个进程都获得矩阵A, B的维度信息 */
    n1 = dim[0];
    n2 = dim[1];
    n3 = dim[2];

    if(DEBUG && myid == 0){
        printf("Read file success!\n");
    }

    /* 为缓存A, B, C, bufA 和 bufB分配空间 */
    /* 一个 m x n 的矩阵 2D 分别在 rootp x rootp 的处理器网格上 */
    /* 若 rootp 不能被 m 或 n 整除, 则子矩阵的尺寸不相同, 因为要移动 A , B 所以，根据分配为 A, B 行和列的最大值 */ 
    /* 取上整除 */
    maxrows_a = (n1 + rootp - 1)/rootp;
    maxcols_a = (n2 + rootp - 1)/rootp;
    /* 矩阵乘要求 B的rows与A的cols值相同 */
    maxrows_b = maxcols_a;
    maxcols_b = (n3 + rootp - 1)/rootp;

    /* 设置子矩阵内存大小, 多两个整型存储矩阵编号 */
    bufA_size = sizeof(int)*2 + sizeof(double)*maxrows_a*maxcols_a;
    bufB_size = sizeof(int)*2 + sizeof(double)*maxrows_b*maxcols_b;
    bufC_size = sizeof(int)*2 + sizeof(double)*maxrows_a*maxcols_b;
    
    /* 每个进程都分配存储空间 */
    if(!(buf = (char *)malloc(bufA_size*2 + bufB_size*2 + bufC_size))) {
        printf("Memory allocation failed\n");
        MPI_Finalize(); // Something error during preprocessing
        exit(-1);
    }
    /* 将相应变量指向内存对应位置 */
    A = (double*)buf;
    bufA = (double*) (buf + bufA_size);
    B = (double*) (buf + bufA_size*2);
    bufB = (double*) (buf + bufA_size*2 + bufB_size);
    C = (double*) (buf + bufA_size*2 + bufB_size*2);

    /* 0号进程将A, B 二维块分发给其他进程 */
    scatter_matrix((double*)(fstreama + sizeof(int)*2), n1, n2, A, rootp);
    scatter_matrix((double*)(fstreamb + sizeof(int)*2), n2, n3, B, rootp);
    
    read_time = MPI_Wtime() - read_time;

    /* 同步所有进程 */
    MPI_Barrier(MPI_COMM_WORLD);
    if(DEBUG && myid == 0){
        printf("scatter matrix success!\n");
    }
    /* 获取开始时间 */
    elapsed_time = MPI_Wtime();
    /* 调用cannon算法计算矩阵相乘 */
    cannon(A, bufA, bufA_size, B, bufB, bufB_size, C, bufC_size, n1, n3, rootp, myid);

    /* 同步所有进程 */
    MPI_Barrier(MPI_COMM_WORLD);

    if(DEBUG && myid == 0){
        printf("cannon success!\n");
    }
    
    /* 计算花费时间 */
    elapsed_time = MPI_Wtime() - elapsed_time;
    /* 0号进程收集计算结果并输出到文件 */
    save_time = MPI_Wtime();

    fsizec = sizeof(int)*2 + sizeof(double)*n1*n3;
    /* 打开保存文件描述符 */
    if(myid == 0) {
        if (!(fhc = fopen(argv[3], "w"))) {
            printf("Can't open file %s, Errno=%d\n", argv[3], errno);
            MPI_Finalize();
        }
        fstreamc = (char *)malloc(fsizec);
        ((int*)fstreamc)[0] = n1;
        ((int*)fstreamc)[1] = n3;
    }
    /* 收集计算结果 */
    gather_matrix((double*)(fstreamc + sizeof(int)*2), n1, n3, C, rootp);

    if(DEBUG && myid == 0){
        printf("gather matrix success!\n");
    }

    /* 确保所有内容都已收到 */
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* 0号进程打印运算结果, 并输出到文件 */
    if(myid == 0) {
        printf("Cannon algrithm: multiply a %dx%d with a %dx%d, use %.2f(s)\n", n1, n2, n2, n3, elapsed_time);
        fwrite(fstreamc, sizeof(char), fsizec, fhc);
        fclose(fhc);
        free(fstreama);
        free(fstreamb);
        free(fstreamc);
        save_time = MPI_Wtime() - save_time;
        printf("read_time = %.2f(s), save_time = %.2f(s), compute_time = %.2f(s)\n", read_time, save_time, elapsed_time);
    }
    free(buf);
    MPI_Finalize();
    return 0;
}
