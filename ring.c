#include <stdio.h>
#include <string.h>
#include "mpi.h"

#define DEBUG 0

int main(int argc, char *argv[]){
    int myid, numprocs, source;                          /* 进程相关变量 */
    MPI_Status status;
    char message[100];                                   /* 消息缓冲区 */
    double starttime, endtime;                           /* 计时起始时间和结束时间 */
    
    MPI_Init(&argc, &argv);                                   
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);                /* 获取当前进程id */
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);            /* 获取总进程数 */
    
    if(myid == 0){
        /* 0号进程行为 */
        strcpy(message, "Hello world");
        /* 起始时间 */
        starttime = MPI_Wtime();
        /* 0号进程发送message给下一个进程1号，tag设置为当前进程id+1 */
        MPI_Send(message, strlen(message) + 1, MPI_CHAR, myid + 1, myid + 1, MPI_COMM_WORLD);
        /* 0号进程接收来自最后一个进程id为numprocs - 1的消息，tag为最后一个进程id+1*/
        MPI_Recv(message, 100, MPI_CHAR, numprocs - 1, numprocs, MPI_COMM_WORLD, &status);
        /* 结束时间 */
        endtime = MPI_Wtime();
        /* 打印循环一圈所消耗的时间 */
        printf("A ring took %f secodes\n", endtime - starttime);
    }else{
        /* 非0号进程行为 */
        /* 接收来自前一个进程的消息，id为myid-1，tag为myid */
        MPI_Recv(message, 100, MPI_CHAR, myid - 1, myid, MPI_COMM_WORLD, &status);
        if(DEBUG) printf("Process %d of %d receive %s from %d\n", myid, numprocs, message, myid - 1);
        /* 将message发送给下一个进程，如果当前进程时最后一个，则发给0号进程 */
        MPI_Send(message, strlen(message) + 1, MPI_CHAR, (myid + 1) % numprocs , myid + 1, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    
    return 0;
}