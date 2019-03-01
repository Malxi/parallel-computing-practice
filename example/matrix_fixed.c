#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pthread.h>
#include <memory.h>

/* 数组大小宏定义 */
#define M 30
#define N 600
#define P 8

struct threadArg{
    int tid;
    float (*B)[P];
    float *A_row;
    float *C_row;
    int numthreads;
};

void *worker(void* arg){
    int i, j;
    struct threadArg *myarg = (struct threadArg *)arg;
    
    /* 平均分配B的所有列 */
    for(i = myarg->tid; i < P; i += myarg->numthreads){
        myarg->C_row[i] = 0.0;
        for(j = 0; j < N; j++){
            myarg->C_row[i] += myarg->A_row[j] * myarg->B[j][i];
        }
    }
    
    return NULL;
}

int main(int argc, char *argv[]){
    float A[M][N], B[N][P], C[M][P];
    float *A_row;
    float *C_row;
    pthread_t *tids;
    struct threadArg *targs;
    float buf[P];
    int taskNum;
    int taskMap[M + 5];
    int myid, numprocs;
    int numthreads;
    MPI_Status status;
    int numsend;
    int sender;
    int i, j;
    double starttime, endtime;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    
    /* 在0号进程内对矩阵A和B进程初始化 */
    if(!myid){
        for(i = 0; i < M; i++)
            for(j = 0; j < N; j++)
                    A[i][j] = i * j + 1;
        for(i = 0; i < N; i++)
            for(j = 0; j < P; j++)
                    B[i][j] = i * j + 1;
    }
    
    /* 0号进程将矩阵B广播给其他进程 */
    MPI_Bcast(B[0], N * P, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    if(!myid){
        starttime = MPI_Wtime();
        /* 0号进程进行任务分配(一共M个任务)，并回收计算结果 */
        /* 首次任务分配，若从进程个数小于任务数M，则取进程数，否则取M */
        j = (numprocs - 1) < M ? (numprocs - 1) : M;
        /* 将矩阵A的行数据发送给从进程，长度为N */
        for(i = 1; i <= j; i++){
            //MPI_Send(A[i-1], N, MPI_FLOAT, i, 99, MPI_COMM_WORLD);
            /* 将矩阵A的行数据发送给从进程，长度为N，tag设置为任务编号，并记录任务与进程映射 */
            MPI_Send(A[i-1], N, MPI_FLOAT, i, i, MPI_COMM_WORLD);
            taskMap[i] = i;
        }
        
        /* 已经分配的任务个数 */ 
        numsend = j;
        
        /* 此时已经发送j个任务，循环接收M个任务 */
        for(i = 1; i <= M; i++){
            /* 计算从进程发送的编号 */
            //sender = (i - 1) % (numprocs - 1) + 1;
            /* 接收第i个任务的数据 */
            //MPI_Recv(C[i-1], P, MPI_FLOAT, sender, 100, MPI_COMM_WORLD, &status);
            /* 接收数据到buf，根据tag来判断是那个任务完成，再根据任务映射来获得发送者编号处理数据 */
            MPI_Recv(buf, P, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            taskNum = status.MPI_TAG;
            sender = taskMap[taskNum];
            //printf("Master0 recieve task %d result from process %d\n", taskNum, sender);
            for(j = 0; j < P; j++){
                C[taskNum-1][j] = buf[j];
            }
            /* 任务未分发完毕，则发送给当前空出的进程, 否则发送结束信息 */
            if(numsend < M){
                numsend++;
                //MPI_Send(A[numsend-1], N, MPI_FLOAT, sender, 99, MPI_COMM_WORLD);
                MPI_Send(A[numsend-1], N, MPI_FLOAT, sender, numsend, MPI_COMM_WORLD);
                taskMap[numsend] = sender;
            }else{
                MPI_Send(&j, 0, MPI_INT, sender, 0, MPI_COMM_WORLD);
            }
           
        }
        
        /* 当进程数大于M时，需要对编号为M + 1和以后的进程发送结束信息 */
        for(i = M + 1; i < numprocs; i++){
            MPI_Send(&j, 0, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        
        endtime = MPI_Wtime();
        printf("Took %f secodes\n", endtime - starttime);
        /* 输出计算结果 */
        for(i = 0; i < M; i++)
            for(j = 0; j < P; j++)
                    printf(j == P - 1 ? "%10.0f\n":"%10.0f ", C[i][j]);
    
    }else{
        /* 从进程计算来自0号进程发送的数据 */
        /* 从进程获取节点CPU数 */
        numthreads = get_nprocs();
        /* 从进程分配pthread_t * tids;tids数组用来存放线程ID */
        tids = (pthread_t *)malloc(numthreads*sizeof(pthread_t));
        /* A_row存放主进程分配的A中一行 */
        A_row = (float *)malloc(N*sizeof(float));
        /* C_row存放计算结果C中的一行*/
        C_row = (float *)malloc(P*sizeof(float));
        /* 线程参数，由于传入参数比较多，采用结果传递 */
        targs = (struct threadArg *)malloc(numthreads*sizeof(struct threadArg));
        for(i = 0; i < numthreads; i++){
            targs[i].tid = i;
            targs[i].B = B;
            targs[i].A_row = A_row;
            targs[i].C_row = C_row;
            targs[i].numthreads = numthreads;
        }
        
        /* 从进程等待来自0号进程的任务 */
        while(1){
            /* 接收0进程发送来A的一行 */
            MPI_Recv(A_row, N, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            /* 若接收到标识为0的消息则退出执行 */
            if(status.MPI_TAG == 0)
                break;
            taskNum = status.MPI_TAG;
            //printf("Slave process %d recieve task %d result from master\n", myid, taskNum);

            /* 创建线程执行计算 */
            for(i = 0; i < numthreads; i++){
                pthread_create(&tids[i], NULL, worker, &targs[i]);
            }
            /* 等待线程内的计算完成 */
            for(i = 0; i < numthreads; i++){
                pthread_join(tids[i], NULL);
            }
            /* 将计算的结果发送给0号进程 */
            //MPI_Send(C_row, P, MPI_FLOAT, 0, 100, MPI_COMM_WORLD);
            MPI_Send(C_row, P, MPI_FLOAT, 0, taskNum, MPI_COMM_WORLD);
        }
    }
    
}
